"""生データフォルダの入力隔離・メソッド選択・解析専用配置（spec §4）。

このモジュールはpipelineの入力固定層。4つの公開関数を持つ。

- `inspect_inputs(source_root, request, *, exe_path)`: 元フォルダを読むだけで
  何一つ変更せず、採用する形式・メソッド・LBM・実行体・極性を決定し、
  「入力計画」(plan) を返す。
- `select_method(candidates)`: メソッド候補を内容（sha256 と、宣言された参照
  キーの解決先ハッシュ）でグルーピングし、一意なら代表を返す。異なる内容が
  複数あれば `METHOD_FILE_CHOICE_REQUIRED`。mtimeは一切使わない。
- `stage_inputs(plan, pipeline_root)`: planに基づき、raw・随伴・実効メソッドを
  `pipeline_root/input` へ実際に配置し、STAT fingerprint等を固定した
  「入力スナップショット」(snapshot) を返す。
- `verify_inputs(snapshot)`: 元rawと配置済みコピーのstat fingerprint・固定
  hashを再検査し、変化していれば例外にする（実行開始・上流終了時の再検査、
  および再開時の照合の両方から呼ばれる想定）。

副作用に注意: mcp_core・pipeline.store等の永続化層をimportしない。返す辞書は
JSON互換のプリミティブのみ（Pathは引数でだけ受け取り、戻り値はすべてstr）。

パスの絶対・相対（Task14への申し送り）:
    - `source_root` / `method.source_path` / `lbm.path` / `exe.path` /
      `pipeline_root` は**絶対文字列**（元フォルダ・元メソッド・元LBM・実行体・
      pipeline_rootそのものは、互いの外側にありうるため相対化できない）。
    - `raw_stat[].relative_path`（source_root基準）・`entries[].name`
      （同）・`companions`（同）・`staged_files[].relative_path`
      （pipeline_root基準）・`staged_files[].source_relative_path`
      （source_root基準）・`method.effective_relative_path`
      （pipeline_root基準）は**相対文字列**。
    pipeline-run.json等へ永続化する際の相対化（絶対パス側を含む）は
    Task14/17の責務であり、このモジュールは何もrelativizeしない。
"""
from __future__ import annotations

import copy
import hashlib
import os
import shutil
from pathlib import Path

from metabolomix.console import job_manager
from metabolomix.console import method_file as method_file_mod
from metabolomix.console import runner as console_runner
from metabolomix.console.input_prep import _companions_of
from metabolomix.console.output_collector import is_upstream_artifact
from metabolomix.core import app_control
from metabolomix.core import user_config
from metabolomix.core.atomic_io import DomainError, canonical_hash

__all__ = ["DEFAULT_MANIFEST_NAME", "inspect_inputs", "manifest_source_record",
           "resolve_manifest_path", "resolve_raw_inventory", "select_method",
           "stage_inputs", "verify_inputs"]

#: request.sample_manifest省略時に探す既定シート名（spec §7.1「元フォルダ直下の
#: `analysis-request.json`と`sample-manifest.tsv`を既定名として探索する」）。
#: 受付（`service.start_pipeline`の事前検査）・`resolve_metadata` handler・
#: 再開時のシート内容照合（`recovery.prepare_resume`）が同じ1つの定数を共有する。
DEFAULT_MANIFEST_NAME = "sample-manifest.tsv"

#: 初期版はlipidomics固定（Global Constraints）。pipeline-request.v1にomicsは無い。
_OMICS = "lipidomics"

_INPUT_SUBDIR = "input"
_INPUTS_META_SUBDIR = "inputs"
_EFFECTIVE_METHOD_NAME = "effective-method.txt"


# ---------- 小さなユーティリティ ----------

def _sha256_file(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def resolve_manifest_path(source_root: Path, request: dict) -> Path | None:
    """`request["sample_manifest"]`を実パスへ解決する。相対はsource_root基準（spec §7.1）。

    省略時は既定名（`DEFAULT_MANIFEST_NAME`）を探し、無ければ`None`
    （＝自動一覧生成へ回す）。明示された値は存在しなくてもそのまま返す
    ——「指定したシートが無い」ことは呼び出し側が
    `SAMPLE_MANIFEST_NOT_FOUND`として報告すべき事実で、ここで既定名へ
    黙って落とすと誤記入が別のシートで走ってしまう。
    """
    manifest_arg = request.get("sample_manifest")
    if manifest_arg is not None:
        candidate = Path(manifest_arg)
        return candidate if candidate.is_absolute() else Path(source_root) / candidate
    default_candidate = Path(source_root) / DEFAULT_MANIFEST_NAME
    return default_candidate if default_candidate.is_file() else None


def manifest_source_record(source_root: Path, request: dict) -> dict:
    """実験情報シートの「今の中身」を表す小さな記録を返す。

    `{"path": <str|None>, "sha256": <str|None>}`。パスが解決できない
    （自動一覧生成）ときは両方None、解決できても読めないときはsha256だけNone。

    **パス文字列だけでは訂正を検出できない**のがこの関数の存在理由。同じ
    `sample-manifest.tsv`を書き直して`pipeline_resume`しても、要求の
    `sample_manifest`は同じ文字列のままなので「変更なし」に見え、旧群割当の
    まま再開してしまう（`recovery._stages_to_reset`が読む側）。内容hashで
    比べれば、書き直しはそれだけで下流の差し戻しになる。
    """
    path = resolve_manifest_path(source_root, request)
    if path is None:
        return {"path": None, "sha256": None}
    try:
        digest = _sha256_file(path)
    except OSError:
        digest = None
    return {"path": str(path), "sha256": digest}


def _resolves_outside_root(path: Path, root: Path) -> bool:
    """path（symlink/reparse pointを含みうる）がroot外を指すかを判定する。

    ジャンクション・シンボリックリンクのどちらも`os.path.realpath`が解決するため、
    「symlinkかどうか」を個別判定する必要はない —— 通常のファイルはrealpathが
    自分自身のままなので常にroot配下、リンクだけが外を指しうる。
    """
    try:
        real = Path(os.path.realpath(path))
        real_root = Path(os.path.realpath(root))
    except OSError:
        return True
    try:
        real.relative_to(real_root)
    except ValueError:
        return True
    return False


def _iter_dir_raw_files(root_entry: Path, source_root: Path) -> list[Path]:
    """ディレクトリ形式rawの内部ファイルを列挙し、reparse point脱出も検査する。"""
    out: list[Path] = []
    for dirpath, _dirnames, filenames in os.walk(root_entry):
        dpath = Path(dirpath)
        if _resolves_outside_root(dpath, source_root):
            raise DomainError(
                "INPUT_ESCAPES_SOURCE_ROOT",
                f"生データ内のディレクトリが元フォルダ外を指しています: {dpath}",
                {"path": str(dpath)},
            )
        for name in filenames:
            fpath = dpath / name
            if _resolves_outside_root(fpath, source_root):
                raise DomainError(
                    "INPUT_ESCAPES_SOURCE_ROOT",
                    f"生データ内のファイルが元フォルダ外を指しています: {fpath}",
                    {"path": str(fpath)},
                )
            out.append(fpath)
    return sorted(out, key=lambda p: str(p).lower())


def _stat_entry(path: Path, relative_path: str, role: str) -> dict:
    st = path.stat()
    return {"relative_path": relative_path, "size": st.st_size,
            "mtime_ns": st.st_mtime_ns, "role": role}


# ---------- 形式解決 ----------

#: プロジェクト保存（`-p`）のときだけ読まれるアセンブリ。Application Control が
#: これを弾くと、MS-DIAL は**全検体の解析を終えた後**に落ちる。
_PROJECT_SAVE_ASSEMBLY = "MsdialLcImMsApi.dll"


def _assert_project_save_possible(request: dict, exe_path) -> None:
    """プロジェクト保存が Application Control に塞がれるなら計画時に止める。

    Console を起動してからでは、実測で 13 分半を費やした後に落ちる。判定は
    「ポリシーが Enforce」かつ「当該アセンブリが未署名」の両方が揃ったときだけ
    ——ポリシー単独で塞ぐと、署名済みの公式配布版を使う正当な構成まで止まる。
    """
    if not request.get("save_project"):
        return
    assembly = Path(exe_path).parent / _PROJECT_SAVE_ASSEMBLY
    if not app_control.project_save_blocked(assembly):
        return
    raise DomainError(
        "PROJECT_SAVE_BLOCKED",
        "このマシンの Application Control（Smart App Control）が、プロジェクト保存に"
        f"要る未署名のアセンブリ {_PROJECT_SAVE_ASSEMBLY} の読み込みを拒否します。"
        "このまま実行すると、MS-DIAL は全検体の解析を終えた**後**に失敗します。"
        "要求に save_project=false を指定してください（.mdproject は解析成果物では"
        "なく GUI で開くための便宜で、必須成果物には含まれません）。",
        {"assembly": str(assembly), "remedy": {"save_project": False}})


def _resolve_raw_format(source_root: Path, requested_extension: str | None) -> tuple[str, dict]:
    """採用する計測拡張子を決める（spec §4.1「1フォルダ直下の1形式」）。

    明示があればそれを使う（フォルダに実在しなければエラー）。省略時は、
    拡張子が1種類ならそれを採用。`.wiff`/`.wiff2`だけが混在し、同一stem集合が
    1対1で一致するときだけ`.wiff`を既定にする（spec §4.2）。それ以外の混在は
    `MIXED_RAW_FORMATS`で形式指定を要求する。rawが無く子フォルダだけがあれば
    `DATASET_SELECTION_REQUIRED`。
    """
    formats = job_manager.raw_input_summary(source_root)
    if requested_extension:
        ext = requested_extension.lower().lstrip(".")
        if ext not in formats:
            raise DomainError(
                "MIXED_RAW_FORMATS",
                f"要求された形式 '{ext}' の計測ファイルがありません: {source_root}",
                {"formats": formats, "requested": ext},
            )
        return ext, formats
    if not formats:
        candidates = [p.name for p in sorted(source_root.iterdir()) if p.is_dir()]
        if candidates:
            raise DomainError(
                "DATASET_SELECTION_REQUIRED",
                f"{source_root} 直下に計測ファイルがなく、子フォルダのみがあります。"
                "解析対象フォルダを明示してください。",
                {"candidates": candidates},
            )
        raise DomainError(
            "MIXED_RAW_FORMATS",
            f"データフォルダに MS-DIAL が読める計測ファイルがありません: {source_root}",
            {"formats": formats},
        )
    if len(formats) == 1:
        return next(iter(formats)), formats
    if set(formats) == {"wiff", "wiff2"}:
        wiff_stems = {p.stem for p in source_root.iterdir() if p.suffix.lower() == ".wiff"}
        wiff2_stems = {p.stem for p in source_root.iterdir() if p.suffix.lower() == ".wiff2"}
        if wiff_stems and wiff_stems == wiff2_stems:
            return "wiff", formats
    raise DomainError(
        "MIXED_RAW_FORMATS",
        "データフォルダに MS-DIAL が対象とする拡張子が2種類以上あり、自動選択できません: "
        + ", ".join(f"{ext}x{n}" for ext, n in sorted(formats.items())),
        {"formats": formats},
    )


def _collect_primaries_and_companions(
    source_root: Path, ext: str,
) -> tuple[list[Path], dict[Path, list[Path]]]:
    """選択済み拡張子の主ファイル一覧と、随伴ファイル一覧を返す（escapeチェック込み）。"""
    entries = sorted(source_root.iterdir())
    file_entries = [p for p in entries if p.is_file()]
    suffix = "." + ext
    primaries = [p for p in entries
                 if p.name.lower().endswith(suffix)
                 and p.suffix.lower().lstrip(".") == ext]
    if not primaries:
        raise DomainError(
            "MIXED_RAW_FORMATS",
            f"入力フォルダに .{ext} がありません: {source_root}",
            {"extension": ext},
        )

    companions_map: dict[Path, list[Path]] = {}
    for primary in primaries:
        if _resolves_outside_root(primary, source_root):
            raise DomainError(
                "INPUT_ESCAPES_SOURCE_ROOT",
                f"生データが元フォルダ外を指しています: {primary}",
                {"path": str(primary)},
            )
        if primary.is_dir():
            _iter_dir_raw_files(primary, source_root)  # escapeチェックのためだけに歩く
            companions_map[primary] = []
            continue
        found = _companions_of(file_entries, primary, ext)
        for companion in found:
            if _resolves_outside_root(companion, source_root):
                raise DomainError(
                    "INPUT_ESCAPES_SOURCE_ROOT",
                    f"随伴ファイルが元フォルダ外を指しています: {companion}",
                    {"path": str(companion)},
                )
        companions_map[primary] = found
    return primaries, companions_map


def _build_entries_and_stat(
    source_root: Path, primaries: list[Path], companions_map: dict[Path, list[Path]],
) -> tuple[list[dict], list[dict]]:
    """物理配置に使う`entries`（トップレベル名の一覧）とSTAT fingerprint一覧を作る。"""
    entries: list[dict] = []
    raw_stat: list[dict] = []
    for primary in primaries:
        if primary.is_dir():
            entries.append({"name": primary.name, "kind": "dir", "role": "primary"})
            for inner in _iter_dir_raw_files(primary, source_root):
                rel = str(inner.relative_to(source_root)).replace(os.sep, "/")
                raw_stat.append(_stat_entry(inner, rel, "primary"))
        else:
            entries.append({"name": primary.name, "kind": "file", "role": "primary"})
            raw_stat.append(_stat_entry(primary, primary.name, "primary"))
        for companion in companions_map.get(primary, []):
            entries.append({"name": companion.name, "kind": "file", "role": "companion",
                            "primary": primary.name})
            raw_stat.append(_stat_entry(companion, companion.name, "companion"))
    return entries, raw_stat


def resolve_raw_inventory(
    source_root: Path, requested_extension: str | None = None,
) -> tuple[str, list[dict]]:
    """1フォルダ直下の計測ファイル（主+随伴）を列挙する（新規・純粋関数）。

    既存の形式選択（`_resolve_raw_format`）とsidecar列挙
    （`_collect_primaries_and_companions` / `_build_entries_and_stat`）を
    そのまま再利用する——lipidomics v1（`inspect_inputs`）と別の選択規則を
    metabolomics側だけに新設しない。`inspect_inputs`自身は呼ばない・呼ばれない
    （既存関数は一切変更しない、純粋な追加）。

    戻り値は`(選択した拡張子, raw_stat一覧)`。raw_stat各要素は
    `{"relative_path", "size", "mtime_ns", "role"}`——内容hashはここでは
    計算しない（呼び出し側の責務。`metabolomix.console.profiles.hash_files`が担う。
    rawファイルは巨大なことがあるため、全量hashを要求する側だけがそのコストを
    負う設計にする）。
    """
    source_root = Path(source_root)
    ext, _formats = _resolve_raw_format(source_root, requested_extension)
    primaries, companions_map = _collect_primaries_and_companions(source_root, ext)
    _entries, raw_stat = _build_entries_and_stat(source_root, primaries, companions_map)
    return ext, raw_stat


# ---------- メソッド選択 ----------

def _method_candidate_dict(path: Path, ion_mode: str | None, mtime: float) -> dict:
    """1件のメソッドファイルからselect_method用のdictを作る（sha256 + 参照hash）。"""
    keys = method_file_mod.read_method_keys(path)
    fingerprint = method_file_mod.method_reference_fingerprint(keys, path)
    reference_hash = canonical_hash({
        key: (info["sha256"] if info["resolved"] else f"UNRESOLVED:{info['declared']}")
        for key, info in fingerprint.items()
    })
    return {
        "path": str(path),
        "sha256": _sha256_file(path),
        "mtime": mtime,
        "polarity": ion_mode,
        "reference_hash": reference_hash,
    }


def select_method(candidates: list[dict]) -> dict:
    """メソッド候補を内容でグルーピングし、一意なら代表を返す（spec §4.1/4.3, D07）。

    グルーピング鍵は `(sha256, reference_hash)`。生バイトが同一でも、宣言された
    参照キー（`Lbm file path`等）の解決先ハッシュが候補ごとに異なれば別の実効
    メソッドとして扱う —— 「異なる相対パス基準で同じ文字列を含む候補」を
    内容ハッシュだけでは見分けられないため。mtimeはグルーピングに一切使わない。
    """
    if not candidates:
        raise DomainError("METHOD_FILE_NOT_GIVEN", "メソッド候補がありません。", {})
    groups: dict[tuple, list[dict]] = {}
    for candidate in candidates:
        key = (candidate["sha256"], candidate.get("reference_hash"))
        groups.setdefault(key, []).append(candidate)
    if len(groups) != 1:
        raise DomainError(
            "METHOD_FILE_CHOICE_REQUIRED", "解析条件を一つ選んでください",
            {"candidates": candidates})
    chosen = sorted(next(iter(groups.values())), key=lambda c: c["path"])[0]
    return chosen


def _discover_candidate_dicts(source_root: Path, polarity: str | None) -> list[dict]:
    method_candidates, _searched = method_file_mod.discover_method_candidates(
        source_root, polarity=polarity, omics=_OMICS)
    if polarity is not None:
        # 明示極性が要求極性と食い違う候補は選択対象から外す（spec §4.1-2）。
        method_candidates = [c for c in method_candidates if c.ion_mode == polarity]
    return [_method_candidate_dict(Path(c.path), c.ion_mode, c.mtime) for c in method_candidates]


def _resolve_method(source_root: Path, request: dict) -> dict:
    """method_fileが明示されていればそれを唯一の候補として使い、省略時は
    discover_method_candidates + select_method で決める（spec §4.1-1/2/3）。"""
    explicit = request.get("method_file")
    if explicit:
        path = Path(explicit).expanduser()
        if not path.is_absolute():
            path = source_root / path
        if not path.is_file():
            raise DomainError(
                "METHOD_FILE_NOT_FOUND",
                f"method_fileが指すファイルがありません: {explicit}",
                {"method_file": explicit})
        keys = method_file_mod.read_method_keys(path)
        ion_mode = (keys.get("ion mode") or "").strip().lower() or None
        return _method_candidate_dict(path, ion_mode, path.stat().st_mtime)

    candidates = _discover_candidate_dicts(source_root, request.get("polarity"))
    if not candidates:
        raise DomainError(
            "METHOD_FILE_NOT_GIVEN",
            f"使えるメソッドファイルが見つかりません: {source_root}",
            {"source_root": str(source_root)})
    return select_method(candidates)


def _resolve_polarity(request: dict, method_path: Path) -> dict:
    """spec §4.1-4: 極性省略時は一意に選べたメソッドのIon mode宣言を採る。

    これは生データからの極性検証ではないので、呼び出し側は`unverified`へ
    その旨を記録する。
    """
    explicit = request.get("polarity")
    if explicit is not None:
        return {"value": explicit, "source": "request_explicit"}
    keys = method_file_mod.read_method_keys(method_path)
    declared = (keys.get("ion mode") or "").strip().lower()
    if declared not in ("positive", "negative"):
        raise DomainError(
            "POLARITY_UNDETERMINED",
            f"極性が未指定で、メソッドファイルの Ion mode も有効な値ではありません: "
            f"{declared!r}（{method_path}）",
            {"method_file": str(method_path), "ion_mode": declared})
    return {"value": declared, "source": "method_declaration"}


def _resolve_lbm_pinned(method_keys: dict, method_path: Path, exe_path: str,
                        request: dict, source_root: Path) -> dict:
    override = request.get("lbm_file")
    if override and not Path(override).is_absolute():
        override = str(source_root / override)
    try:
        lbm_setting = user_config.get_setting("msdial.lbm")
    except user_config.ConfigInvalidError as exc:
        raise DomainError(exc.code, exc.message, exc.details()) from exc
    lbm = method_file_mod.resolve_lbm(
        method_keys, method_path, omics=_OMICS, exe_path=exe_path,
        lbm_setting=lbm_setting, override=override)
    if lbm.error_code:
        raise DomainError(lbm.error_code, lbm.message or "",
                          {"candidates": list(lbm.candidates)} if lbm.candidates else {})
    if lbm.path is None:
        return {"path": None, "sha256": None, "source": lbm.source}
    resolved = Path(lbm.path).resolve()
    return {"path": str(resolved), "sha256": _sha256_file(resolved), "source": lbm.source}


def _resolve_exe(exe_path: Path) -> dict:
    exe_path = Path(exe_path)
    if not exe_path.is_file():
        raise DomainError("MSDIAL_EXE_NOT_FOUND", f"実行体がありません: {exe_path}",
                          {"exe_path": str(exe_path)})
    if not console_runner.is_console_exe(str(exe_path)):
        raise DomainError(
            "MSDIAL_EXE_NOT_CONSOLE",
            f"Console の実行体（[msdial] exe / MSDIAL_EXE）が MS-DIAL Console ではありません: {exe_path}",
            {"exe_path": str(exe_path)})
    return {"path": str(exe_path.resolve()), "sha256": _sha256_file(exe_path), "version": None}


# ---------- 公開API: inspect_inputs ----------

def inspect_inputs(source_root: Path, request: dict, *, exe_path: Path) -> dict:
    """元フォルダを読むだけで、採用する形式・メソッド・LBM・実行体・極性を決める。

    何も書き込まない（隔離・固定は`stage_inputs`の責務）。戻り値の「入力計画」は
    そのまま`stage_inputs`へ渡す。
    """
    source_root = Path(source_root).expanduser()
    if not source_root.is_dir():
        raise DomainError("DATASET_ROOT_NOT_FOUND", f"source_rootが存在しません: {source_root}",
                          {"source_root": str(source_root)})

    ext, _formats = _resolve_raw_format(source_root, request.get("keep_extension"))
    primaries, companions_map = _collect_primaries_and_companions(source_root, ext)
    entries, raw_stat = _build_entries_and_stat(source_root, primaries, companions_map)
    companions = {p.name: [c.name for c in cs] for p, cs in companions_map.items() if cs}

    chosen_method = _resolve_method(source_root, request)
    method_path = Path(chosen_method["path"])
    method_keys = method_file_mod.read_method_keys(method_path)

    # 最終選択後にだけ厳格な参照解決を行う（select_methodのグルーピングは弱い版）。
    method_file_mod.resolve_method_references(method_keys, method_path)

    polarity = _resolve_polarity(request, method_path)
    exe_info = _resolve_exe(Path(exe_path))
    lbm_info = _resolve_lbm_pinned(method_keys, method_path, exe_info["path"], request, source_root)

    unverified: list[str] = []
    if polarity["source"] == "method_declaration":
        unverified.append("polarity_from_method_declaration_not_verified_from_raw")

    # 相対で宣言されたパスキー（LBM・MSP・Text DB 等）は原本基準の絶対パスへ
    # 書き換える。LC-MS の Console は宣言パスを解決せず、相対値を自分の cwd
    # 基準で読む（メソッド基準で解くのは GC-MS 経路だけ）ので、実効コピーへ
    # verbatim に写すと見つからずに黙って飛ばされる（spec §4.3「コピー後に
    # 相対参照の意味を変えない」）。実在は問わない——LBM 以外の古い宣言で
    # 既存 lipidomics 実行を止めない（REFERENCE_KEYS を広げない方針と同じ）。
    overrides = method_file_mod.relative_path_overrides(method_keys, method_path)
    if lbm_info["path"] and (lbm_info["source"] != "method_file"
                             or method_file_mod.LBM_KEY in overrides):
        # 相対宣言なら resolve_lbm の解決結果で置き換え、build_tree/env/exe_dir
        # へフォールバックした場合は（原本に行が無いので）新規追加する。
        # 絶対宣言はそのまま（Console がそのまま読める）。
        overrides[method_file_mod.LBM_KEY] = lbm_info["path"]

    # Console を起動する前に判定する（起動後だと全検体の解析を終えてから落ちる）。
    _assert_project_save_possible(request, exe_info["path"])

    return {
        "source_root": str(source_root.resolve()),
        "selected_format": ext,
        "entries": entries,
        "raw_stat": raw_stat,
        "companions": companions,
        "method": {
            "source_path": str(method_path.resolve()),
            "sha256": chosen_method["sha256"],
            "effective_relative_path": None,
            "effective_sha256": None,
            "overrides": overrides,
        },
        "lbm": {"path": lbm_info["path"], "sha256": lbm_info["sha256"]},
        "exe": exe_info,
        "polarity": polarity,
        "unverified": unverified,
    }


def plan_from_profile(source_root: Path, request: dict, profile: dict) -> dict:
    """profileだけを情報源に入力配置計画を作る（v2）。

    `inspect_inputs`と**同じ形**のplanを返すので、`stage_inputs`・
    `_plan_fingerprint`・`_handle_prepare_inputs_v2`はそのまま共有できる。
    違うのは決め方だけ——methodをフォルダから推定せず、LBMを必須にせず、
    環境設定の実行体へフォールバックしない（profileが唯一の情報源）。

    解決とhash照合はやり直さない。`resolve_profile_inputs`が既にmethod・依存・
    実行体・rawの実在とhash一致を検証しているので、ここはその結果を形へ移すだけ。
    """
    # profilesはmodule先頭でこのモジュールをimportしている。関数内で読む
    # （module先頭に書くと循環import）。
    from metabolomix.console import profiles as profiles_mod

    source_root = Path(source_root).expanduser()
    if not source_root.is_dir():
        raise DomainError("DATASET_ROOT_NOT_FOUND", f"source_rootが存在しません: {source_root}",
                          {"source_root": str(source_root)})

    profile_path = Path(request["profile_file"])
    resolved = profiles_mod.resolve_profile_inputs(
        profile, profile_path.parent, raw_root=source_root)

    # raw形式はprofileが宣言しない（データ由来であってmethod由来ではない）。
    ext, _formats = _resolve_raw_format(source_root, request.get("keep_extension"))
    primaries, companions_map = _collect_primaries_and_companions(source_root, ext)
    entries, raw_stat = _build_entries_and_stat(source_root, primaries, companions_map)
    companions = {p.name: [c.name for c in cs]
                  for p, cs in companions_map.items() if cs}

    environment = resolved["execution_environment"]
    exe_path = Path(environment["executable_path"])
    if not console_runner.is_console_exe(str(exe_path)):
        # profileがGUIのMSDIAL.exeを固定していても、hashは一致してしまう。
        raise DomainError(
            "MSDIAL_EXE_NOT_CONSOLE",
            f"profileが宣言した実行体がMS-DIAL Consoleではありません: {exe_path}",
            {"exe_path": str(exe_path)})

    # v1がLBM1件に使っていた任意キーdictを、全依存へそのまま一般化し、
    # profile外で原本が相対宣言したパスキーも絶対化する（記録用の
    # `snapshot_profile`と同じ関数を通し、実行と記録をずらさない）。
    overrides = profiles_mod.effective_method_overrides(resolved)
    lbm = next((dep for dep in resolved["dependencies"]
                if dep["method_key"] == method_file_mod.LBM_KEY), None)

    _assert_project_save_possible(request, exe_path)

    return {
        "source_root": str(source_root.resolve()),
        "selected_format": ext,
        "entries": entries,
        "raw_stat": raw_stat,
        "companions": companions,
        "method": {
            "source_path": resolved["method"]["source_path"],
            "sha256": resolved["method"]["sha256"],
            "effective_relative_path": None,
            "effective_sha256": None,
            "overrides": overrides,
        },
        "lbm": ({"path": lbm["source_path"], "sha256": lbm["sha256"]} if lbm
                else {"path": None, "sha256": None}),
        "exe": {"path": str(exe_path), "sha256": environment["executable_sha256"],
                "version": environment["msdial_version"]},
        "polarity": {"value": resolved["polarity"], "source": "profile"},
        # 極性をrawから検証していない点はv1と同じ。黙って確定扱いにしない。
        "unverified": ["polarity_from_profile_not_verified_from_raw"],
    }


# ---------- 公開API: stage_inputs ----------

def _copy_or_link(src: Path, dst: Path, link_fn) -> None:
    try:
        link_fn(src, dst)
    except OSError:
        # ボリュームをまたぐ／リンク非対応のファイルシステム。リンクは書込隔離では
        # ないので、コピーでもハードリンクでも「元rawは読むだけ」という前提は
        # 変わらない（spec §4.2）。
        shutil.copy2(src, dst)


def _reconcile_staged_file(src: Path, dst: Path, name: str) -> None:
    src_stat = src.stat()
    dst_stat = dst.stat()
    if dst_stat.st_size != src_stat.st_size or dst_stat.st_mtime_ns != src_stat.st_mtime_ns:
        raise DomainError(
            "STAGED_INPUT_MISMATCH",
            f"配置済みのファイルが元入力と一致しません（再開時の照合）: {name}",
            {"name": name})


def _reconcile_staged_dir(dst: Path, expected_stats: list[dict], dirname: str) -> None:
    prefix = dirname + "/"
    for stat in expected_stats:
        rel = stat["relative_path"][len(prefix):]
        dst_file = dst / rel
        if not dst_file.is_file():
            raise DomainError(
                "STAGED_INPUT_MISMATCH",
                f"配置済みディレクトリに欠落ファイルがあります: {dirname}/{rel}",
                {"name": f"{dirname}/{rel}"})
        st = dst_file.stat()
        if st.st_size != stat["size"] or st.st_mtime_ns != stat["mtime_ns"]:
            raise DomainError(
                "STAGED_INPUT_MISMATCH",
                f"配置済みディレクトリのファイルが元入力と一致しません: {dirname}/{rel}",
                {"name": f"{dirname}/{rel}"})


def stage_inputs(plan: dict, pipeline_root: Path, *,
                  disk_usage=shutil.disk_usage, link_fn=os.link) -> dict:
    """planに基づき、raw・随伴・実効メソッドを`pipeline_root/input`へ実際に配置する。

    1. 合計サイズと空き容量を確認する（起動前に停止、spec §4.2）。
    2. 新規（空の配置先）ならos.link優先→copy2フォールバックで全ファイルを配置。
    3. 既存の配置先（resume）なら、期待される各ファイルの対応・fingerprintを
       照合し、不一致を単にskipしない（spec §4.2「同名ファイルは...検証し、
       不一致を単にskipしない」）。欠けているものだけ補って配置する。
    4. 実効メソッドファイルを`inputs/effective-method.txt`へ常に書く。
    """
    source_root = Path(plan["source_root"])
    pipeline_root = Path(pipeline_root)
    input_dir = pipeline_root / _INPUT_SUBDIR
    inputs_meta_dir = pipeline_root / _INPUTS_META_SUBDIR

    total_size = sum(entry["size"] for entry in plan["raw_stat"])
    usage_root = pipeline_root
    while not usage_root.exists():
        usage_root = usage_root.parent
    usage = disk_usage(usage_root)
    if usage.free < total_size:
        raise DomainError(
            "INSUFFICIENT_FREE_SPACE",
            f"配置先の空き容量が不足しています（必要 {total_size} バイト、"
            f"空き {usage.free} バイト）: {pipeline_root}",
            {"required": total_size, "available": usage.free})

    input_dir.mkdir(parents=True, exist_ok=True)
    entries = plan["entries"]
    expected_names = {e["name"] for e in entries}
    existing_names = {p.name for p in input_dir.iterdir()}
    # 上流が入力フォルダ側へ書いた生成物は「想定外」ではない。MS-DIAL は `-o` だけで
    # なく `-i` 側にも .arf / .arf2 / .pai2 / .dcl / .EIC.aef / _tags.xml を書き、
    # pipeline では `-i` が staged input そのもの（`service._handle_upstream` が
    # `pipeline_root/input` を job.dataset_root として渡す）。`console/execution.py`
    # が両ルートを snapshot して収集するのと同じ前提をここでも採る——採らないと、
    # Console が一度でも走った run では `rerun_upstream` が必ずここで止まる。
    # 判定規則は `output_collector` の role 表が唯一の出所（二重定義にしない）。
    unexpected = {name for name in (existing_names - expected_names)
                  if not is_upstream_artifact(name)}
    if unexpected:
        raise DomainError(
            "STAGED_INPUT_MISMATCH",
            f"配置済みの入力に想定外のファイルがあります: {sorted(unexpected)}",
            {"unexpected": sorted(unexpected)})

    raw_stat_by_prefix: dict[str, list[dict]] = {}
    for stat in plan["raw_stat"]:
        prefix = stat["relative_path"].split("/", 1)[0]
        raw_stat_by_prefix.setdefault(prefix, []).append(stat)

    staged_files: list[dict] = []
    for entry in entries:
        name = entry["name"]
        src = source_root / name
        dst = input_dir / name
        if entry["kind"] == "dir":
            if dst.exists():
                _reconcile_staged_dir(dst, raw_stat_by_prefix.get(name, []), name)
            else:
                shutil.copytree(src, dst)
            for stat in raw_stat_by_prefix.get(name, []):
                staged_files.append({
                    "relative_path": f"{_INPUT_SUBDIR}/{stat['relative_path']}",
                    "source_relative_path": stat["relative_path"]})
        else:
            if dst.exists():
                _reconcile_staged_file(src, dst, name)
            else:
                _copy_or_link(src, dst, link_fn)
            staged_files.append({
                "relative_path": f"{_INPUT_SUBDIR}/{name}",
                "source_relative_path": name})

    overrides = dict(plan["method"]["overrides"])
    method_path = Path(plan["method"]["source_path"])
    effective_path = inputs_meta_dir / _EFFECTIVE_METHOD_NAME
    try:
        method_file_mod.write_effective_method_file(method_path, effective_path, overrides)
    except UnicodeEncodeError as exc:
        raise DomainError(
            "METHOD_ENCODING_UNSUPPORTED",
            f"実効メソッドをASCIIで書き出せません（非ASCII文字を含みます）: {exc}",
            {"overrides": overrides}) from exc

    snapshot = copy.deepcopy(plan)
    snapshot["pipeline_root"] = str(pipeline_root.resolve())
    snapshot["staged_files"] = staged_files
    snapshot["method"]["effective_relative_path"] = f"{_INPUTS_META_SUBDIR}/{_EFFECTIVE_METHOD_NAME}"
    snapshot["method"]["effective_sha256"] = _sha256_file(effective_path)
    return snapshot


# ---------- 公開API: verify_inputs ----------

def verify_inputs(snapshot: dict) -> None:
    """元rawと配置済みコピーのstat fingerprint・固定hashを再検査する（spec D05）。

    変化していれば`INPUT_CHANGED`（元raw・method/LBM/exeの内容）または
    `STAGED_INPUT_MISMATCH`（配置済みコピー）にする。rawの全量内容ハッシュは
    検証しない（spec: 「rawの全量ハッシュは初期版の必須条件にしない」）。
    """
    source_root = Path(snapshot["source_root"])
    changed: list[str] = []
    for stat in snapshot["raw_stat"]:
        path = source_root / stat["relative_path"]
        try:
            st = path.stat()
        except OSError:
            changed.append(stat["relative_path"])
            continue
        if st.st_size != stat["size"] or st.st_mtime_ns != stat["mtime_ns"]:
            changed.append(stat["relative_path"])
    if changed:
        raise DomainError(
            "INPUT_CHANGED",
            f"元の生データが実行開始後に変化しています: {changed}",
            {"changed": changed})

    pipeline_root = snapshot.get("pipeline_root")
    if pipeline_root:
        pipeline_root = Path(pipeline_root)
        raw_stat_by_rel = {s["relative_path"]: s for s in snapshot["raw_stat"]}
        mismatched: list[str] = []
        for staged in snapshot.get("staged_files", []):
            source_stat_entry = raw_stat_by_rel.get(staged["source_relative_path"])
            if source_stat_entry is None:
                continue
            path = pipeline_root / staged["relative_path"]
            try:
                st = path.stat()
            except OSError:
                mismatched.append(staged["relative_path"])
                continue
            if (st.st_size != source_stat_entry["size"]
                    or st.st_mtime_ns != source_stat_entry["mtime_ns"]):
                mismatched.append(staged["relative_path"])
        if mismatched:
            raise DomainError(
                "STAGED_INPUT_MISMATCH",
                f"配置済みの入力が元データと一致しません: {mismatched}",
                {"mismatched": mismatched})

    for label, info, path_key in (
        ("method", snapshot.get("method"), "source_path"),
        ("lbm", snapshot.get("lbm"), "path"),
        ("exe", snapshot.get("exe"), "path"),
    ):
        if not info or not info.get(path_key) or not info.get("sha256"):
            continue
        path = Path(info[path_key])
        try:
            digest = _sha256_file(path)
        except OSError as exc:
            raise DomainError("INPUT_CHANGED", f"{label}のファイルが読めません: {path}",
                              {"which": label, "path": str(path)}) from exc
        if digest != info["sha256"]:
            raise DomainError("INPUT_CHANGED", f"{label}の内容が変化しています: {path}",
                              {"which": label, "path": str(path)})
