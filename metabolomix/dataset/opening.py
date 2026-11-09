"""データセットを開く 1 つの入口（解析フォルダ・mzTab-M・Console の job・pipeline の run）。

`load_dataset`（対話も再解析も）はここを呼ぶ。入口は 4 つのうちちょうど 1 つで、読むアダプタは `source` で決まる:

| 入口 | `auto`（既定） | `folder` | `mztab` |
|---|---|---|---|
| `directory`（解析フォルダ、またはその中の `.arf` / `.arf2`） | 解析フォルダ | 解析フォルダ | 拒否 |
| `mztab_path` | mzTab-M | 拒否 | mzTab-M |
| `job_path` / `pipeline_path` | 解析フォルダ | 解析フォルダ | mzTab-M（今までの経路そのまま） |

job / pipeline を解析フォルダで読むときも、ジョブの門（`mztab.loading.job_gate`）は
mzTab-M で読むときと同じものを通る——未完了のジョブを既定で拒否し、出所の裏取りと
探索専用の印を同じ規則で書く。読む `.arf` は、正準 mzTab-M の隣（pipeline は run の
`input/` も）の候補のうち、mzTab-M のアダプタが検出状態を取り込むときと同じ検証
（`mztab.evidence.matching_arfs`。特徴数・m/z・試料名）を通ったちょうど 1 本。通るものが
無ければ `ALIGNMENT_NOT_FOUND`、2 本以上なら `ALIGNMENT_AMBIGUOUS` で止め、mzTab-M へ黙って
切り替えない（どちらで読んだかで数字の出どころが変わるため）。

読み込んだ基はハンドルの登録簿（`core.handles.registry`）に覚え、同じ入力を開き直したときは
読み直さない。「同じ」は、読み込みが読んだファイルと探したファイル（後から現れうる `.arf`
の候補など）の size と更新時刻がすべて同じこと（`_signature`）。返す `dataset` は毎回その
複製で、呼び出し側どうしで状態は共有しない。
"""
from __future__ import annotations

import copy
import dataclasses
import threading
import time
import uuid
from dataclasses import dataclass, field
from pathlib import Path

from metabolomix.core import handles
from metabolomix.core.atomic_io import DomainError
from metabolomix.dataset.alignment_files import AlignmentSelection, select_alignment
from metabolomix.mztab.dataset_state import DatasetState

__all__ = ["KINDS", "SOURCES", "OpenedDataset", "open_dataset", "pipeline_arf2",
           "pipeline_arf2_candidates"]

#: 入口の引数名（`kind`）。
KINDS = ("directory", "mztab_path", "job_path", "pipeline_path")
#: `source` に渡せる値。
SOURCES = ("auto", "folder", "mztab")
#: mzTab-M のアダプタが受ける `measure`（mzTab-M の定量種別はファイルとジョブの宣言で決まる）。
_MZTAB_MEASURE = "peak_height"
#: 解析フォルダが見つからないときの案内。
_MZTAB_HINT = ("mzTab-M のアダプタで読むなら source=\"mztab\" を指定してください"
               "（解析フォルダが決まらないときに mzTab-M へ黙って切り替えることはしない）。")
#: ジョブの正準 mzTab-M を読めず、解析フォルダの経路がバッチを選べないときの案内。
_BROKEN_MZTAB_HINT = ("解析フォルダの経路は、どのバッチの .arf を読むかを正しい mzTab-M との照合で"
                      "決める。ジョブの検査（出所の裏取り・探索専用の印）なしで解析フォルダを"
                      "開くなら directory=<run の msdial フォルダ> を渡してください。")


@dataclass
class OpenedDataset:
    """開いた結果。`dataset` は呼び出し側のもの（基の複製）。

    kind: `directory` / `mztab_path` / `job_path` / `pipeline_path`
    path: 入口に渡されたパスの絶対パス
    adapter: `msdial_folder`（`.arf` / `.arf2` を直接）/ `mztab`（mzTab-M）
    selection: 解析フォルダのアダプタが選んだ `.arf` 一式（mzTab-M のアダプタは None）
    notes: どう開いたかの注記（読み込み自体の注記は `dataset.load_notes`）
    reused_base: 読み込み済みの基を使い回したか
    mztab_path: `.arf` と照合した正準 mzTab-M（mzTab-M のアダプタと、job / pipeline を
        解析フォルダで読んだとき。基の指紋が今の `.arf` の候補を数えるのに使う）
    """

    dataset: DatasetState
    kind: str
    path: str
    adapter: str
    selection: AlignmentSelection | None
    notes: list[str] = field(default_factory=list)
    reused_base: bool = False
    mztab_path: str | None = None


def open_dataset(*, directory=None, mztab_path=None, job_path=None, pipeline_path=None,
                 source: str = "auto", measure: str = "peak_height",
                 allow_incomplete: bool = False) -> OpenedDataset:
    """4 つの入口のちょうど 1 つを、`source` で決まるアダプタで開く（モジュールの説明の表）。

    measure: 解析フォルダのアダプタが値に使う `.arf` の欄（`peak_height` / `peak_area_above_zero`）。
        mzTab-M のアダプタは既定の `peak_height` だけを受ける（定量種別はファイルで決まる）。
    allow_incomplete: 完了していないジョブ・run の出力を探索専用として読む（どちらのアダプタでも）。

    エラーは `DomainError`（`DATASET_BAD_REQUEST` / `ALIGNMENT_NOT_FOUND` / `ALIGNMENT_INCOMPLETE` /
    `INCOMPLETE_ANALYSIS_JOB` / `MZTAB_NOT_FOUND` など）。
    """
    from metabolomix.dataset.msdial_adapter import MEASURES

    if type(allow_incomplete) is not bool:
        # 真偽値以外を真偽として読むと、"false" という文字列が True になる。
        raise DomainError("DATASET_BAD_REQUEST",
                          f"allow_incomplete は bool です: {allow_incomplete!r}",
                          {"allow_incomplete": repr(allow_incomplete)})
    if source not in SOURCES:
        raise DomainError("DATASET_BAD_REQUEST",
                          f"source は {' / '.join(SOURCES)} のどれかです: {source!r}",
                          {"field": "source"})
    given = {kind: value for kind, value in zip(
        KINDS, (directory, mztab_path, job_path, pipeline_path)) if value}
    if len(given) != 1:
        raise DomainError(
            "DATASET_BAD_REQUEST",
            "directory / mztab_path / job_path / pipeline_path のどれか 1 つを指定してください"
            "（directory: MS-DIAL の解析フォルダ。mztab_path: .mzTab への直接パス。"
            "job_path: analysis-job.json。pipeline_path: pipeline-run.json かその親）。",
            {"given": sorted(given)})
    kind, raw = next(iter(given.items()))
    adapter = _adapter(kind, source)
    if measure not in MEASURES:
        raise DomainError("DATASET_BAD_REQUEST",
                          f"measure は {' / '.join(MEASURES)} のどちらかです: {measure!r}",
                          {"field": "measure"})
    if adapter == "mztab" and measure != _MZTAB_MEASURE:
        raise DomainError(
            "DATASET_BAD_REQUEST",
            f"measure={measure} は解析フォルダのアダプタだけが使います。mzTab-M の定量種別は"
            "ファイル（ジョブの宣言）で決まるので、別の値を読むなら解析フォルダで読んでください。",
            {"field": "measure", "adapter": adapter})

    input_path = Path(str(raw)).expanduser().resolve()
    key = f"{kind}:{adapter}:{measure}:{allow_incomplete}:{input_path}"
    cached = handles.registry.cached_base(key)
    if cached is not None and _signature(cached[0]) == cached[1]:
        return _copy(cached[0], reused=True)
    # 同じ入力の同時の読み込みは 1 回にする: キーのロックを取ってから基のキャッシュを見直し、
    # 無ければ読む。待っていた側は先に読んだ側の基を使う。別の入力はロックが別なので並ぶ。
    with _inflight_lock(key):
        cached = handles.registry.cached_base(key)
        if cached is not None and _signature(cached[0]) == cached[1]:
            return _copy(cached[0], reused=True)
        # 読み込みの間に更新された依存があれば、基は古い中身かもしれないので覚えない
        # （結果は返す）。読む前の時刻と、読んだ後の依存の更新時刻だけを比べる——どの経路・
        # どのパスの綴りでも効く。読む直前に書かれたファイルは覚えられないだけで、安全側。
        started = time.time_ns()
        base = _open(kind, str(raw), str(input_path), adapter, measure, allow_incomplete)
        signature = _signature(base)
        if _changed_during_load(started, signature):
            base.notes.append("読み込みの間に入力が変わった（読み込みの開始後に更新された"
                              "ファイルがある）ので、この基は覚えない。次の読み込みは読み直す。")
        else:
            handles.registry.put_base(key, base, signature)
    return _copy(base, reused=False)


def _changed_during_load(started_ns: int, signature: tuple) -> bool:
    """読み込みの開始後に更新された依存があるか（更新時刻が開始時刻以降）。"""
    return any(mtime is not None and mtime >= started_ns for _path, _size, mtime in signature)


#: 基のキャッシュのキーごとの「読み込み中」ロック。辞書は `_INFLIGHT_GUARD` で守る。
#: キーは入力のパスと読み方の組なので、増える数は開いた入力の数まで（ロックは軽い）。
_INFLIGHT: dict[str, threading.Lock] = {}
_INFLIGHT_GUARD = threading.Lock()


def _inflight_lock(key: str) -> threading.Lock:
    """基のキャッシュのキーごとのロック（同じキーには同じロックを返す）。"""
    with _INFLIGHT_GUARD:
        lock = _INFLIGHT.get(key)
        if lock is None:
            lock = _INFLIGHT[key] = threading.Lock()
        return lock


def _adapter(kind: str, source: str) -> str:
    if kind == "directory":
        if source == "mztab":
            raise DomainError("DATASET_BAD_REQUEST",
                              "directory は解析フォルダのアダプタでしか読めません"
                              "（source=\"mztab\" を使うなら mztab_path / job_path / pipeline_path）。",
                              {"field": "source", "kind": kind})
        return "msdial_folder"
    if kind == "mztab_path":
        if source == "folder":
            raise DomainError("DATASET_BAD_REQUEST",
                              "mztab_path は mzTab-M のアダプタでしか読めません"
                              "（解析フォルダで読むなら directory）。",
                              {"field": "source", "kind": kind})
        return "mztab"
    return "mztab" if source == "mztab" else "msdial_folder"


def _copy(base: OpenedDataset, *, reused: bool) -> OpenedDataset:
    """基の複製を返す。`dataset_id` は DatasetState の実体ごとの ID なので振り直す。"""
    ds = copy.deepcopy(base.dataset)
    ds.dataset_id = f"ds_{uuid.uuid4().hex}"
    return dataclasses.replace(base, dataset=ds, notes=list(base.notes), reused_base=reused)


def _signature(base: OpenedDataset) -> tuple:
    """基の指紋。読み込みが読んだファイルと探したファイルの (パス, size, mtime)。

    中身が変われば基を使い回さない。無いファイルは「無い」と数えるので、後から現れたら
    読み直す（mzTab-M の隣に後から置かれた `.arf` は検出状態を変え、解析フォルダに後から
    置かれた `.dcl` やライブラリはハンドルの照合を変える）。探したファイルの一覧（`.arf` の
    候補など）は呼ばれるたびに今のフォルダから数え直す。

    フォルダ（解析フォルダ・pipeline の run）は在るか無いかだけを数え、更新時刻は入れない。中で効く
    ファイルは上の一覧と付随ファイルが個別に数えており、フォルダの更新時刻は記録（`reports/`）が
    できただけでも変わる（読み込みの最初の概観がそこへ書く）ので、基を無駄に捨てる。
    """
    out = []
    for path in sorted(_dependencies(base)):
        try:
            stat = Path(path).stat()
            if Path(path).is_dir():
                out.append((path, 0, 0))
                continue
            out.append((path, stat.st_size, stat.st_mtime_ns))
        except OSError:
            out.append((path, None, None))
    return tuple(out)


def _dependencies(base: OpenedDataset) -> set[str]:
    """読み込みが読んだ・探したファイル（`_signature` の材料）。"""
    ds = base.dataset
    paths = {base.path, *(getattr(ds, "source_files", None) or {})}
    record_path = (getattr(ds, "library_record", None) or {}).get("path")
    if record_path:
        paths.add(str(record_path))
    if getattr(ds, "source_format", None) == "msdial_folder":
        paths.update(_folder_sidecars(ds))
    join = ((getattr(ds, "feature_qc", None) or {}).get("join") or {}).get("path")
    if join:
        paths.add(str(join))
    if base.mztab_path:
        from metabolomix.mztab import evidence

        # 正準 mzTab-M そのもの（解析フォルダの経路も `.arf` の照合と出所の裏取りに読む）。
        paths.add(str(base.mztab_path))
        paths.update(str(p) for p in evidence.arf_candidates(base.mztab_path,
                                                             ds.artifact_paths))
    if base.selection is not None:
        paths.update(_folder_listing(base.selection.arf.parent))
    if getattr(ds, "job_path", None):
        paths.update(_job_files(Path(ds.job_path)))
    if base.kind == "pipeline_path":
        from metabolomix.pipeline.recovery import normalize_pipeline_root

        root = normalize_pipeline_root(Path(base.path))
        paths.add(str(root / "pipeline-run.json"))
        paths.update(str(p) for p in (root / "input").glob("AlignResult-*.arf2"))
    return paths


def _folder_sidecars(ds) -> list[str]:
    """`.arf2` ごとに、解析フォルダのアダプタが読む・探す同じバッチの付随ファイル。

    タグ（`_tags.xml`）・キュレーションのフラグ（`curation/flags.jsonl`）・`.dcl`・
    `.EIC.aef`。無くても数える（後から現れたら読み直す）。
    """
    from metabolomix.curation.flags import FLAGS_FILENAME, curation_dir
    from metabolomix.msdial.alignment_spot import alignment_stem, find_by_name
    from metabolomix.msdial.tags import alignment_tag_path

    out = []
    for path in getattr(ds, "source_files", None) or {}:
        if not path.lower().endswith(".arf2"):
            continue
        out += [str(alignment_tag_path(path)), str(curation_dir(path) / FLAGS_FILENAME)]
        stem = alignment_stem(path)
        if stem is not None:
            for suffix in (".dcl", ".EIC.aef"):
                # 大文字小文字違いの綴りで実在するものはその名前で数える（無ければ標準の綴り）。
                found = find_by_name(stem.parent, stem.name + suffix)
                out.append(str(found or stem.with_name(stem.name + suffix)))
    return out


def _folder_listing(folder: Path) -> list[str]:
    """`select_alignment` が探すもの: フォルダの `.arf` と `*_Loaded.*.dbs` と `.mddata` / `.mdproject`（今の一覧）。

    `.mddata` は後から現れうる（試料の役割・群・バッチの宣言元で、2 本目が現れて選択が変わることも
    ある）ので、無かった基を使い回さない。`.mdproject` は `.mddata` が無いときの宣言元（`msdial/classes.py`）。
    フォルダの更新時刻は指紋に入れないので、ここに載らない種類は後から現れても気づかない。"""
    try:
        return [str(p) for p in folder.iterdir()
                if p.is_file() and (p.name.lower().endswith((".arf", ".mddata", ".mdproject"))
                                    or p.match("*_Loaded.*.dbs"))]
    except OSError:
        return []


def _job_files(job_path: Path) -> list[str]:
    """ジョブ経由の読み込みが読むファイル: ジョブの記録・実効メソッド・それが宣言する MSP・
    終了証跡（呼ばれるたびに今の記録から数え直す。記録を読めなければ記録だけ）。"""
    from metabolomix.console.execution import receipt_path
    from metabolomix.handoff.schema import AnalysisJob
    from metabolomix.mztab import loading

    out = [str(job_path)]
    try:
        job = AnalysisJob.load(job_path)
    except (OSError, ValueError, KeyError, TypeError):
        return out
    method = (getattr(job, "method_file", "") or "").strip()
    if method:
        out.append(method)
    msp = loading.effective_method_msp(job)
    if msp is not None:
        out.append(str(msp))
    out.append(str(receipt_path(Path(job.run_dir))))
    return out


# ---------- 読み込み ----------

def _open(kind: str, raw: str, path: str, adapter: str, measure: str,
          allow_incomplete: bool) -> OpenedDataset:
    from metabolomix.dataset import library_memory

    selection = None
    notes: list[str] = []
    mztab_path = None
    if adapter == "mztab":
        ds = _load(kind, raw, allow_incomplete)
        # 直接指定なら入力そのもの（mzTab-M の読み取りは拡張子を問わない）。job / pipeline は
        # 読み込んだ mzTab-M（この経路の source_files は mzTab-M だけ）。
        mztab_path = path if kind == "mztab_path" else next(
            (p for p in ds.source_files if p.lower().endswith(".mztab")),
            next(iter(ds.source_files), None))
    elif kind == "directory":
        from metabolomix.dataset.msdial_adapter import load_msdial_folder

        selection = select_alignment(Path(raw).expanduser())
        ds = load_msdial_folder(selection, measure=measure)
    elif kind == "job_path":
        ds, selection, gate = _folder_via_job(Path(raw).expanduser(), measure=measure,
                                              allow_incomplete=allow_incomplete)
        notes.append(_folder_note(gate, "ジョブ"))
        mztab_path = str(gate.mztab_path)
    else:
        ds, selection, gate, notes = _folder_via_run(raw, measure=measure,
                                                     allow_incomplete=allow_incomplete)
        mztab_path = str(gate.mztab_path)

    if kind == "directory":
        record = library_memory.remember_from_folder(selection)
    elif kind == "job_path":
        record = library_memory.remember_from_job(Path(raw).expanduser())
    elif kind == "pipeline_path":
        from metabolomix.pipeline.recovery import normalize_pipeline_root

        record = library_memory.remember_from_run(normalize_pipeline_root(Path(raw).expanduser()))
    else:
        record = library_memory.unknown_record()
    ds.library_record = record
    # 記録の置き場所の基準: 解析フォルダは選んだ `.arf` のフォルダ（job / pipeline もここ）、
    # mzTab-M は mzTab-M のフォルダ。
    if selection is not None:
        ds.analysis_dir = str(Path(selection.arf).parent)
    elif mztab_path:
        ds.analysis_dir = str(Path(mztab_path).parent)
    return OpenedDataset(dataset=ds, kind=kind, path=path, adapter=adapter,
                         selection=selection, notes=notes, mztab_path=mztab_path)


def _load(kind: str, path: str, allow_incomplete: bool):
    """mzTab-M のアダプタ（今までの読み込みそのまま）。"""
    if kind == "pipeline_path":
        from metabolomix.pipeline import session_handoff
        return session_handoff.load_run_dataset(path, allow_incomplete=allow_incomplete)
    from metabolomix.mztab.loading import load_dataset_state
    return load_dataset_state(
        mztab_path=Path(path) if kind == "mztab_path" else None,
        job_path=Path(path) if kind == "job_path" else None,
        allow_incomplete=allow_incomplete)


def _folder_via_job(job_path: Path, *, measure: str, allow_incomplete: bool,
                    pipeline_root: Path | None = None):
    """ジョブの門を通し、正準 mzTab-M と同じアラインメントの解析フォルダを読む。(ds, selection, gate)。

    読む `.arf` は `_matching_arf` が mzTab-M との照合で選ぶ。ジョブが記録した付随
    アーティファクトは mzTab-M のアダプタと同じく `ds.artifact_paths` に載せる。
    """
    from metabolomix.analysis import sample_manifest
    from metabolomix.dataset.msdial_adapter import load_msdial_folder
    from metabolomix.mztab import loading

    gate = loading.job_gate(job_path, allow_incomplete=allow_incomplete)
    artifact_paths = loading.job_artifact_paths(gate.job)
    target, raw_by_name, meta_by_name, scan = _matching_arf(gate, artifact_paths, pipeline_root)
    try:
        selection = select_alignment(target)
    except DomainError as exc:
        raise _with_hint(exc) from exc

    # 照合に読んだ `.arf` をそのまま使う（同じファイルを 2 度なめない）。
    ds = load_msdial_folder(selection, measure=measure, scan=scan)
    # 試料ごとの raw のパスは正準 mzTab-M（`ms_run[N]-location`）から。`.mddata` が無い
    # Console の出力でも、run の試料対応表と結べる。
    raws = [raw_by_name.get((name or "").strip()) for name in ds.sample_names]
    ds.sample_raw_paths = raws if any(raws) else None
    # バッチと注入順も同じ mzTab-M の MTD から（`.mddata` が宣言した値は勝つ）。`.mddata` が無い
    # Console の出力では、ここで引き継がないと注入順が無いことになりドリフト補正が使えない。
    ds.assay_metadata = {
        assay_id: dict(meta_by_name[name.strip()])
        for assay_id, name in zip(ds.sample_assay_ids, ds.sample_names)
        if name.strip() in meta_by_name}
    # raw のパスもバッチも既定の行に入るので、読み込み時に作った行を作り直す（読み込みの
    # 一部なので metadata_revision は進めない）。
    revision = ds.metadata_revision
    sample_manifest.apply_defaults_and_inputs(ds, sample_manifest.resolve_metadata(ds, None))
    ds.metadata_revision = revision
    for role, paths in artifact_paths.items():
        known = ds.artifact_paths.setdefault(role, [])
        known.extend(p for p in paths if p not in known)
    ds.job_path = str(job_path)
    ds.exploratory_only = gate.exploratory_only
    ds.source_verification = "verified" if gate.verified else "legacy_unverified"
    if gate.incomplete_warning:
        ds.load_notes.insert(0, gate.incomplete_warning)
    if gate.unverified_warning:
        ds.load_notes.append(gate.unverified_warning)
    return ds, selection, gate


def _matching_arf(gate, artifact_paths: dict[str, list[str]],
                  pipeline_root: Path | None):
    """正準 mzTab-M と同じアラインメントの `.arf` をちょうど 1 本選ぶ。

    戻り値は (その `.arf`, 試料名（strip 済み） -> raw のパス, 試料名 -> {batch, run_order},
    その `.arf` を照合に読んだ `ArfScan`（解析フォルダのアダプタが読み直さずに使う。
    定量値は `MEASURES` の両方を持つ）)。
    raw のパスは照合に読んだ mzTab-M の `ms_run[N]-location`、バッチと注入順は同じ mzTab-M の
    MTD（`MS:4000088` / `MS:4000089`。どちらも読み直さずに使う）。

    候補は mzTab-M のアダプタの evidence と同じ順（`arf_candidates`。handoff の
    `peak_matrix_source` → 隣の `*PeakProperties.arf` → 隣のその他）。どれも検証
    （`matching_arfs`。特徴数・m/z・試料名）を通らず pipeline の run なら、`input/` の
    `AlignResult-*.arf2`（`pipeline_arf2`）と同じバッチの `.arf` も同じ検証にかける。
    名前の順で先頭を取らない——同じフォルダに別のバッチが居るのは普通で、名前の古い
    バッチを「同じアラインメント」と呼んでしまう。
    """
    from metabolomix.dataset.msdial_adapter import MEASURES
    from metabolomix.mztab import evidence, loading

    axes: list = []
    raw_by_name: dict[str, str] = {}
    meta_by_name: dict[str, dict] = {}

    def verify(paths: list[Path]):
        if not axes:
            try:
                state, parsed = loading.read_mztab_parsed(gate.mztab_path)
            except DomainError as exc:
                message = exc.message if exc.message.endswith("。") else f"{exc.message}。"
                reason = (f"読めなかった理由: {exc.code}。この mzTab-M 自体が読めないので、"
                          "source=\"mztab\" でも読めない。")
                raise DomainError(exc.code, f"{message}{reason}{_BROKEN_MZTAB_HINT}",
                                  exc.details) from exc
            axes.extend(evidence.join_axes(state))
            raw_by_name.update(
                (name.strip(), raw) for name, raw in zip(
                    state.sample_names, loading.sample_raw_paths(gate.job, state, parsed))
                if raw)
            for name, assay_id in zip(state.sample_names, state.sample_assay_ids):
                meta = (state.assay_metadata or {}).get(assay_id) or {}
                if meta.get("batch") is not None or meta.get("run_order") is not None:
                    meta_by_name[name.strip()] = {"batch": meta.get("batch"),
                                                  "run_order": meta.get("run_order")}
        return evidence.matching_arfs(paths, feature_mz=axes[0], sample_names=axes[1],
                                      measures=tuple(MEASURES.values()))

    candidates = evidence.arf_candidates(gate.mztab_path, artifact_paths)
    passed, rejected, scan = verify(candidates) if candidates else ([], {}, None)
    if not passed and pipeline_root is not None:
        # run の input/ のバッチを新しい順に、合うものが見つかるまで試す。
        # すでに候補に入っている `.arf`（解決後のパスで比べる）は二重に検証・列挙しない。
        for arf2 in pipeline_arf2_candidates(pipeline_root):
            try:
                fallback = select_alignment(arf2).arf
            except DomainError as exc:
                raise _with_hint(exc) from exc
            if fallback.resolve() in {p.resolve() for p in candidates}:
                continue
            candidates.append(fallback)
            more, more_rejected, more_scan = verify([fallback])
            passed += more
            rejected.update(more_rejected)
            scan = scan or more_scan
            if passed:
                break

    details = {"mztab": gate.mztab_path.name, "candidates": [p.name for p in candidates],
               "source_hint": "mztab"}
    if len(passed) == 1:
        return passed[0], raw_by_name, meta_by_name, scan
    if passed:
        # directory で開くとジョブの門を通らない。未完了のジョブでそれを勧めると、探索専用の
        # 印が外れたまま比較まで進めてしまうので、そのときは mzTab-M だけを案内する。
        way_out = _MZTAB_HINT if gate.exploratory_only else (
            "どれを読むかは directory に .arf を直接渡して決められるが、その読み方はジョブの検査"
            "（出所の裏取り・探索専用の印）を通らない。" + _MZTAB_HINT)
        raise DomainError(
            "ALIGNMENT_AMBIGUOUS",
            f"ジョブの mzTab-M（{gate.mztab_path.name}）と同じアラインメントだと検証できた .arf が"
            f"複数あります（{', '.join(p.name for p in passed)}）。{way_out}",
            {**details, "matching": [p.name for p in passed]})
    if not candidates:
        where = (f"ジョブの mzTab-M（{gate.mztab_path.name}）の隣"
                 + ("にも run の input/ にも" if pipeline_root is not None else "に"))
        raise DomainError(
            "ALIGNMENT_NOT_FOUND",
            f"{where} MS-DIAL のアラインメント（.arf）が見つかりません。{_MZTAB_HINT}", details)
    listed = ", ".join(f"{name}（{reason}）" for name, reason in rejected.items())
    raise DomainError(
        "ALIGNMENT_NOT_FOUND",
        f"ジョブの mzTab-M（{gate.mztab_path.name}）と同じアラインメントの .arf がありません"
        f"（特徴数・m/z・試料名で照合して、どの候補も合わなかった: {listed}）。{_MZTAB_HINT}",
        {**details, "rejected": rejected})


def _with_hint(exc: DomainError) -> DomainError:
    """`select_alignment` のエラーに `source="mztab"` の案内を足す。"""
    message = exc.message if exc.message.endswith("。") else f"{exc.message}。"
    return DomainError(exc.code, f"{message}{_MZTAB_HINT}",
                       {**(exc.details or {}), "source_hint": "mztab"})


def _folder_via_run(pipeline_path: str, *, measure: str, allow_incomplete: bool):
    """pipeline の run を解析フォルダで読む。run から引き継ぐのは解決済みの試料対応表だけ。"""
    from metabolomix.analysis import sample_manifest
    from metabolomix.pipeline import session_handoff

    root, record, job_path = session_handoff.run_console_job(pipeline_path)
    ds, selection, gate = _folder_via_job(job_path, measure=measure,
                                          allow_incomplete=allow_incomplete, pipeline_root=root)
    rows = session_handoff.run_manifest_rows(root, record)
    if rows is not None:
        # roles / sample_meta も組み直す（run の試料対応表の役割を、前処理の前から見せる）。
        sample_manifest.apply_defaults_and_inputs(ds, sample_manifest.resolve_metadata(ds, rows))
    notes = [_folder_note(gate, "run"),
             "run の解析行列・feature binding・差次的結果は読んでいない（解析フォルダの"
             "アダプタが run から引き継ぐのは解決済みの試料対応表だけ）。それらも載せるなら"
             " source=\"mztab\" で読む。"]
    return ds, selection, gate, notes


def _folder_note(gate, where: str) -> str:
    return (f"{where}の mzTab-M（{gate.mztab_path.name}）と同じアラインメントの解析フォルダを"
            "読んだ。mzTab-M を読むなら source=\"mztab\"。")


def pipeline_arf2_candidates(pipeline_path) -> list[Path]:
    """pipeline の run の `input/AlignResult-*.arf2` を、名前の新しい順に並べる。

    run では mzTab が `console/attempt-*/msdial/` の下、`.arf`・`.arf2`・`.dcl` が
    `input/` にあるので、mzTab の隣を探す既定の経路では見つからない。バッチが複数あるときは
    名前の降順（決定的にするため。`AlignResult-<日時など>` は名前の順が新しさの順になる）。
    `pipeline_path` は `pipeline-run.json` でもその親でもよい。
    """
    from metabolomix.pipeline.recovery import normalize_pipeline_root

    root = normalize_pipeline_root(Path(pipeline_path))
    return [p for p in sorted((root / "input").glob("AlignResult-*.arf2"),
                              key=lambda p: p.name, reverse=True) if p.is_file()]


def pipeline_arf2(pipeline_path) -> Path | None:
    """pipeline の run の `input/AlignResult-*.arf2` のうち名前が最も新しいもの。無ければ None。"""
    return next(iter(pipeline_arf2_candidates(pipeline_path)), None)
