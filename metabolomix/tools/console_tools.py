"""Console 実行層の MCP ツール: メソッドファイルの候補列挙と雛形作成。

Console の計画・実行・状態は pipeline（`pipeline_plan` / `pipeline_run` / `pipeline_status`）が
担う。ここに残すのは、メソッドファイルを用意するための 2 ツールだけ。
"""
from __future__ import annotations

from pathlib import Path

from mcp.types import ToolAnnotations

from metabolomix.core.mcp_core import mcp
from metabolomix.core.mcp_errors import console_error
from metabolomix.core.serialization import json_payload

__all__ = ["console_method_template", "console_method_candidates"]


def _non_ascii_override_error(overrides: dict[str, str]) -> str | None:
    """実効メソッドへ書く値に ASCII で表せない文字があれば封筒を返す。

    Console はメソッドを `Encoding.ASCII` で読む（`ConfigParser`）。非 ASCII の
    フォルダ名は `?` に化けてライブラリが見つからず、黙って同定 0 件になる。
    pipeline の `METHOD_ENCODING_UNSUPPORTED` と同じコードで止める。
    """
    bad = {key: value for key, value in overrides.items() if not value.isascii()}
    if not bad:
        return None
    return console_error(
        "METHOD_ENCODING_UNSUPPORTED",
        "実効メソッドへ書くパスに ASCII 以外の文字があります。MS-DIAL Console は"
        "メソッドファイルを ASCII で読むため、このパスのライブラリは見つからず"
        "黙って同定 0 件になります。ASCII だけのフォルダへ置き直してください: "
        + ", ".join(f"{key}={value}" for key, value in sorted(bad.items())),
        {"keys": sorted(bad)})


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=False, destructiveHint=False, idempotentHint=True),
          structured_output=False)
def console_method_template(
    out_path: str,
    polarity: str,
    based_on: str | None = None,
    dataset_root: str | None = None,
    omics: str = "lipidomics",
    library_mode: str = "auto",
    msp_file: str | None = None,
) -> str:
    """既存のパラメータから、別極性用のメソッドファイルを作ります（最後の手段）。

    **まず pipeline_plan の request で method_file を省略して試してください。** MS-DIAL GUI は解析の
    たびに `<project>_param_<終了時刻>.txt` を自動保存するので、その極性で一度でも
    GUI 実行があれば、作る必要はありません。

    このツールが要るのは「その極性で一度も処理していない」場合だけです
    （別極性は処理済み、というのが典型）。やることは 2 行の差し替えです:
    `Ion mode` と `Searched adduct ions` をその極性の標準セットにする。
    **検出・アライメント条件は元のまま引き継ぎます**（勝手に変えると別の解析になる）。
    加えて、GUI 由来では必ず空の `Lbm file path` を解決して埋めます。

    out_path: 書き出し先。
    polarity: "positive" / "negative"（作りたい側）。
    based_on: 元にするパラメータファイル。省略時は dataset_root から探します
        （**極性は問いません** — 別極性から作るのがこのツールの用途なので）。
    dataset_root: based_on 省略時の探索先。dataset_root 直下・兄弟フォルダ・
        過去 run まで探します。候補が複数あれば based_on で選ぶよう
        METHOD_FILE_CHOICE_REQUIRED で止まります（土台の選択は解析条件そのもの
        なので、最新を黙って採りません）。
    library_mode: "auto"（既定）か "msp_only"。"msp_only" は pipeline_plan と同じく
        `Lbm file path` を空にし、`Msp file path` に作りたい極性の研究室 MSP を書き、
        Text DB と注釈器設定表の宣言を空にする（解決順・止まる条件も pipeline_plan と同じ）。
    msp_file: library_mode="msp_only" で使う MSP の明示パス。
    """
    if polarity not in ("positive", "negative"):
        return console_error("JOB_NOT_PLANNED",
                             f"polarity は 'positive' または 'negative' です: {polarity!r}")
    library_error = _library_mode_error(library_mode, None, msp_file)
    if library_error:
        return library_error

    from metabolomix.console import method_file as method_file_mod

    if based_on:
        src = Path(based_on).expanduser()
    elif dataset_root:
        # 極性で絞らない。「別極性から作る」のがこのツールの用途なので、
        # polarity=None のまま候補を集める。
        candidates, searched = method_file_mod.discover_method_candidates(
            Path(dataset_root).expanduser(), polarity=None, omics=omics)
        if not candidates:
            return console_error(
                "METHOD_FILE_NOT_GIVEN",
                "元にできるパラメータファイル（`*_param_<ts>.txt`）が見つかりません: "
                f"{dataset_root}  MS-DIAL GUI で一度も解析していないフォルダには"
                "存在しません（兄弟フォルダと過去 run も探しました）。"
                "他のデータセットのパラメータを based_on で明示してください。",
                {"dataset_root": dataset_root, "searched": searched})
        if len(candidates) > 1:
            return console_error(
                "METHOD_FILE_CHOICE_REQUIRED",
                f"元にできる候補が {len(candidates)} 件あります。based_on で 1 つ"
                "選んでください。検出・アライメント条件は選んだファイルのものが"
                "そのまま引き継がれるため、どれを土台にするかは解析条件の選択です。",
                {"dataset_root": dataset_root, "searched": searched,
                 **_capped_candidates(candidates)},
                required_tools=["console_method_template"])
        src = Path(candidates[0].path)
    else:
        return console_error("METHOD_FILE_NOT_GIVEN",
                             "based_on か dataset_root のどちらかを指定してください。")

    if not src.is_file():
        return console_error("METHOD_FILE_NOT_FOUND", f"元ファイルが見つかりません: {src}")
    if not _looks_like_method_text(src):
        return console_error(
            "METHOD_FILE_NOT_TEXT",
            f"元ファイルが MS-DIAL Console の読める ASCII テキストではありません: {src}  "
            ".mdproject / .mddata は ZIP で、解析パラメータを含んでいません"
            "（中身は .mddata へのポインタだけです）。",
            {"based_on": str(src)})

    src_keys = method_file_mod.read_method_keys(src)

    from metabolomix.console import runner as console_runner
    try:
        exe = console_runner.get_exe_path()
    except EnvironmentError as exc:
        # exe はビルド生成物・exe フォルダからの LBM 推定にだけ使う。未設定でも
        # メソッドの宣言や [library] lbm で解決できれば足りるので、ここでは止めない。
        if getattr(exc, "code", "") == "CONFIG_INVALID":
            return _exe_error(exc)
        exe = None
    # 別フォルダへ書き出すので、相対宣言は原本基準の絶対パスに固定する
    # （_resolve_libraries が作る上書きに含まれる）。
    resolved_libraries = _resolve_libraries(
        src_keys, src, omics=omics, exe=exe, polarity=polarity, library_mode=library_mode,
        msp_file=msp_file)
    if isinstance(resolved_libraries, str):
        return resolved_libraries
    library_overrides, library_payload = resolved_libraries
    overrides = {
        **library_overrides,
        method_file_mod.ION_MODE_KEY: polarity.capitalize(),
        method_file_mod.ADDUCT_KEY: method_file_mod.STANDARD_ADDUCTS[polarity],
    }
    encoding_error = _non_ascii_override_error(overrides)
    if encoding_error:
        return encoding_error

    dest = Path(out_path).expanduser()
    try:
        method_file_mod.write_effective_method_file(src, dest, overrides)
    except OSError as exc:
        return console_error("METHOD_FILE_NOT_FOUND", f"書き出せませんでした: {exc}")

    return json_payload({
        "status": "written",
        "out_path": str(dest),
        "based_on": str(src),
        "polarity": polarity,
        **library_payload,
        "changed_keys": sorted(overrides),
        "caveat": "検出・アライメント条件は元ファイルのまま引き継いでいます。"
                  "その極性に妥当かは実行前に確認してください。",
        # msp_only の雛形は `Lbm file path` が空なので、auto で計画すると LBM が自動補完され
        # 契約（LBM 不使用）が黙って破れる。次の一手にもモードを引き継ぐ。
        "next": (f"pipeline_plan(dataset_root=..., request={{'method_file': {str(dest)!r}, "
                 f"'polarity': {polarity!r}"
                 + (", 'library_mode': 'msp_only'}" if library_mode == "msp_only" else "}")
                 + ")"),
    })


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=True, destructiveHint=False, idempotentHint=True),
          structured_output=False)
def console_method_candidates(
    dataset_root: str,
    polarity: str | None = None,
    omics: str | None = "lipidomics",
    search_dirs: list[str] | None = None,
) -> str:
    """使えるメソッドファイルの候補を列挙します（pipeline_plan が失敗する前に呼べます）。

    MS-DIAL GUI は解析のたびに `<project>_param_<終了時刻>.txt` をプロジェクト
    フォルダへ自動保存します。その極性で一度も GUI 実行が無いフォルダには存在
    しないため、**兄弟フォルダ**（POS の隣の NEG 等）と**過去 run** まで探します。

    dataset_root: 生データフォルダ。
    polarity: "positive" / "negative"。省略すると極性で区別せず全件返します。
        指定すると各候補に `usable` が付き、一致しないものは
        `needs_polarity_conversion`（console_method_template を経由させる）。
    omics: 既定 "lipidomics"。None で絞り込みません。
    search_dirs: 追加で探すフォルダ。再帰はしません。

    候補には比較用の `key_params`（検出・アライメント条件のうち結果を変える少数）が
    付きます。候補が 10 件を超えるときは付きません（戻り値が肥大するため）。
    """
    root = Path(dataset_root).expanduser()
    if not root.is_dir():
        return console_error("DATASET_ROOT_NOT_FOUND",
                             f"データフォルダが見つかりません: {dataset_root}",
                             {"dataset_root": str(dataset_root)})

    from metabolomix.console import method_file as method_file_mod

    candidates, searched = method_file_mod.discover_method_candidates(
        root, polarity=polarity, omics=omics, search_dirs=search_dirs)
    return json_payload({
        "dataset_root": str(root),
        "polarity": polarity,
        "omics": omics,
        "searched": searched,
        "n_candidates": len(candidates),
        **_capped_candidates(candidates),
        "next": _candidates_next_hint(polarity),
    })


def _looks_like_method_text(path: Path) -> bool:
    """MS-DIAL Console の ConfigParser が読める形かを判定する。

    ConfigParser は ASCII のテキストを 1 行ずつ読み、`#` 始まりを飛ばして
    最初の `:` または `=` で key/value に割る。ここでは先頭 8192 バイトだけを
    調べ、テキストであることと key/value 行が 1 つ以上あることを確認する。

    **先頭だけを読む**。この判定が弾く相手（.mdproject / .mddata）は 60 サンプル
    規模で GB 級になり、全体を読むとエラー封筒を返す前にメモリを潰す。
    """
    try:
        with path.open("rb") as fh:
            head = fh.read(8192)
    except OSError:
        return False
    if b"\x00" in head:
        return False
    text = head.decode("ascii", errors="replace")
    for line in text.splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        positions = [i for i in (stripped.find(":"), stripped.find("=")) if i > 0]
        if positions:
            return True
    return False


def _candidate_payload(candidate) -> dict:
    """MethodCandidate を UI が読める辞書にする。mtime は ISO 文字列で返す。

    `usable` は polarity を指定せずに探索したとき None になる（「一致するか」を
    判定しようがないため）。"direct" と取り違えられないよう、`key_params` と
    同じくキーごと省く。
    """
    from datetime import datetime
    payload = {
        "path": candidate.path,
        "origin": candidate.origin,
        "ion_mode": candidate.ion_mode,
        "omics": candidate.omics,
        "has_lbm": candidate.has_lbm,
        "mtime": datetime.fromtimestamp(candidate.mtime).isoformat(timespec="seconds"),
    }
    if candidate.usable is not None:
        payload["usable"] = candidate.usable
    if candidate.key_params is not None:
        payload["key_params"] = candidate.key_params
    return payload


def _candidates_next_hint(polarity: str | None) -> str:
    """console_method_candidates の `next` 文面。polarity 省略時は usable が
    候補に付かない（finding 2）ので、判断材料を ion_mode に差し替える。
    """
    if polarity is None:
        return ("polarity を省略したため候補に usable は付きません。"
                "ion_mode を見て、目的の極性と一致するものを "
                "pipeline_plan の request.method_file に、違うものは "
                "console_method_template(based_on=..., polarity=...) を"
                "通してから pipeline_plan の request.method_file に渡してください。")
    return ("usable=direct なら pipeline_plan の request.method_file に、"
            "usable=needs_polarity_conversion なら console_method_template("
            "based_on=..., polarity=...) を通してから pipeline_plan の request.method_file に渡す")


def _capped_candidates(candidates) -> dict:
    """候補一覧を `MAX_REPORTED_CANDIDATES` で切る。

    候補は `discover_method_candidates` が direct 優先・新しい順にソート済みなので、
    先頭から切れば最も有用な候補が残る。切ったときだけ `truncated: true` を立て、
    切っていないときはキー自体を出さない。呼び出し側が別途持つ「総数」
    （`n_candidates` や封筒メッセージの件数）はここでは変えない —
    切った件数と混同させないため。
    """
    from metabolomix.console import method_file as method_file_mod

    emitted = candidates[:method_file_mod.MAX_REPORTED_CANDIDATES]
    section: dict = {"candidates": [_candidate_payload(c) for c in emitted]}
    if len(candidates) > method_file_mod.MAX_REPORTED_CANDIDATES:
        section["truncated"] = True
    return section


def _exe_error(exc: EnvironmentError) -> str:
    from metabolomix.console.runner import msdial_exe_setup_help
    """`get_exe_path` の失敗を封筒にする。

    未設定（`MSDIAL_EXE_NOT_FOUND`）には人間がすべき手順を、設定ファイルが読めない
    （`CONFIG_INVALID`）ときはその行・列を載せる。
    """
    code = getattr(exc, "code", "MSDIAL_EXE_NOT_FOUND")
    details = dict(getattr(exc, "details", {}) or {})
    if code == "MSDIAL_EXE_NOT_FOUND":
        details = {**msdial_exe_setup_help(), **details}
    return console_error(code, str(exc), details or None)


def _lbm_setting():
    """`[library] lbm` の設定。設定ファイルが読めなければ封筒（str）を返す。"""
    from metabolomix.core import user_config
    try:
        return user_config.get_setting(user_config.LIBRARY_LBM_SETTING_KEY)
    except user_config.ConfigInvalidError as exc:
        return console_error(exc.code, exc.message, exc.details())


_LIBRARY_MODES = ("auto", "msp_only")


def _library_mode_error(library_mode: str, lbm_file: str | None, msp_file: str | None) -> str | None:
    """library_mode の値と、lbm_file / msp_file との組み合わせを検査する。"""
    if library_mode not in _LIBRARY_MODES:
        return console_error("JOB_NOT_PLANNED",
                             f"library_mode は 'auto' または 'msp_only' です: {library_mode!r}")
    if library_mode == "msp_only" and lbm_file:
        return console_error(
            "LIBRARY_MODE_CONFLICT",
            "library_mode='msp_only' は LBM を使わないので、lbm_file と同時には指定できません。",
            {"library_mode": library_mode})
    if library_mode == "auto" and msp_file:
        return console_error(
            "LIBRARY_MODE_CONFLICT",
            "msp_file は library_mode='msp_only' のときだけ使えます。library_mode='auto' では"
            "MSP を解決しないので、渡しても使われません。",
            {"library_mode": library_mode})
    return None


def _msp_setting(polarity: str):
    """極性の研究室 MSP の設定。設定ファイルが読めなければ封筒（str）を返す。"""
    from metabolomix.core import user_config
    try:
        return user_config.get_setting(user_config.MSP_SETTING_KEYS[polarity])
    except user_config.ConfigInvalidError as exc:
        return console_error(exc.code, exc.message, exc.details())


def _resolve_libraries(method_keys: dict, method_path: Path, *, omics: str, exe: str | None,
                       polarity: str, library_mode: str, lbm_file: str | None = None,
                       msp_file: str | None = None):
    """実効メソッドの上書き（相対宣言の絶対化＋ライブラリ）と戻り値の断片を返す。

    失敗したら封筒（str）。研究室 MSP の置き場所は戻り値・エラー文に載せない
    （ファイル名だけ。core/path_resolvers.LibraryPathError と同じ規約）。
    """
    from metabolomix.console import method_file as method_file_mod
    from metabolomix.core import user_config

    overrides = method_file_mod.relative_path_overrides(method_keys, method_path)
    if library_mode == "msp_only":
        msp_setting = _msp_setting(polarity)
        if isinstance(msp_setting, str):
            return msp_setting
        msp = method_file_mod.resolve_msp(msp_file, polarity, msp_setting)
        if msp.error_code:
            details = ({"msp_file": Path(msp_file).name} if msp.source == "argument"
                       else user_config.describe_missing(
                           user_config.MSP_SETTING_KEYS[polarity], msp_setting))
            return console_error(msp.error_code, msp.message or "", details)
        if not msp.path.isascii():
            return console_error(
                "METHOD_ENCODING_UNSUPPORTED",
                f"研究室 MSP（{Path(msp.path).name}）のパスに ASCII 以外の文字があります。"
                "MS-DIAL Console はメソッドファイルを ASCII で読むため、このライブラリは"
                "見つからず黙って同定 0 件になります。ASCII だけのフォルダへ置き直してください。",
                {"keys": [method_file_mod.MSP_KEY]})
        extra, removed = method_file_mod.msp_only_overrides(method_keys, msp.path)
        overrides.update(extra)
        return overrides, {
            "library_mode": "msp_only",
            "lbm": {"path": None, "source": "disabled"},
            "msp": {"file": Path(msp.path).name, "source": msp.source},
            "removed_declarations": removed,
        }

    lbm_setting = _lbm_setting()
    if isinstance(lbm_setting, str):
        return lbm_setting
    lbm = method_file_mod.resolve_lbm(
        method_keys, method_path, omics=omics, exe_path=exe, lbm_setting=lbm_setting,
        override=lbm_file)
    if lbm.error_code:
        return console_error(lbm.error_code, lbm.message or "",
                             {"candidates": list(lbm.candidates)} if lbm.candidates else None)
    # `[library] lbm` から解いた LBM は研究室の資産なので、戻り値・エラー文には
    # ファイル名だけを出す（どの経路でも。spec 2026-10-07 §3.3）。実効メソッドへは
    # 絶対パスを書く（Console が読むため）。ほかの出所（宣言・ビルド生成物・exe フォルダ）は
    # これまでどおり絶対パスを示す。
    lab_lbm = lbm.source == method_file_mod.SETTING_SOURCE and bool(lbm.path)
    if lab_lbm and not lbm.path.isascii():
        return console_error(
            "METHOD_ENCODING_UNSUPPORTED",
            f"研究室 LBM（{Path(lbm.path).name}）のパスに ASCII 以外の文字があります。"
            "MS-DIAL Console はメソッドファイルを ASCII で読むため、このライブラリは"
            "見つからず黙って同定 0 件になります。ASCII だけのフォルダへ置き直してください。",
            {"keys": [method_file_mod.LBM_KEY]})
    if lbm.path:
        overrides[method_file_mod.LBM_KEY] = lbm.path
    return overrides, {
        "library_mode": "auto",
        "lbm": {"path": Path(lbm.path).name if lab_lbm else lbm.path, "source": lbm.source},
        "msp": {"file": None, "source": "not_used"},
        "removed_declarations": [],
    }
