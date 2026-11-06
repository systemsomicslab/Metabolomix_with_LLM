"""データセットを開く 1 つの入口（解析フォルダ・mzTab-M・Console の job・pipeline の run）。

`dataset_load`（対話のセッション）と `dataset_handle_open`（再解析のハンドル）はどちらも
ここを呼ぶ。入口は 4 つのうちちょうど 1 つで、読むアダプタは `source` で決まる:

| 入口 | `auto`（既定） | `folder` | `mztab` |
|---|---|---|---|
| `directory`（解析フォルダ、またはその中の `.arf` / `.arf2`） | 解析フォルダ | 解析フォルダ | 拒否 |
| `mztab_path` | mzTab-M | 拒否 | mzTab-M |
| `job_path` / `pipeline_path` | 解析フォルダ | 解析フォルダ | mzTab-M（今までの経路そのまま） |

job / pipeline を解析フォルダで読むときも、ジョブの門（`mztab.loading.job_gate`）は
mzTab-M で読むときと同じものを通る——未完了のジョブを既定で拒否し、出所の裏取りと
探索専用の印を同じ規則で書く。正準 mzTab-M の隣に `.arf` が無ければ `ALIGNMENT_NOT_FOUND`
で止め、mzTab-M へ黙って切り替えない（どちらで読んだかで数字の出どころが変わるため）。

読み込んだ基はハンドルの登録簿（`core.handles.registry`）に覚え、同じ入力（読んだファイルの
size と更新時刻が同じ）を開き直したときは読み直さない。返す `dataset` は毎回その複製で、
呼び出し側どうしで状態は共有しない。
"""
from __future__ import annotations

import copy
import dataclasses
import uuid
from dataclasses import dataclass, field
from pathlib import Path

from metabolomix.core import handles
from metabolomix.core.atomic_io import DomainError
from metabolomix.dataset.alignment_files import AlignmentSelection, select_alignment
from metabolomix.mztab.dataset_state import DatasetState

__all__ = ["KINDS", "SOURCES", "OpenedDataset", "open_dataset"]

#: 入口の引数名（`kind`）。
KINDS = ("directory", "mztab_path", "job_path", "pipeline_path")
#: `source` に渡せる値。
SOURCES = ("auto", "folder", "mztab")
#: mzTab-M のアダプタが受ける `measure`（mzTab-M の定量種別はファイルとジョブの宣言で決まる）。
_MZTAB_MEASURE = "peak_height"
#: 解析フォルダが見つからないときの案内。
_MZTAB_HINT = ("mzTab-M のアダプタで読むなら source=\"mztab\" を指定してください"
               "（解析フォルダが無いときに mzTab-M へ黙って切り替えることはしない）。")


@dataclass
class OpenedDataset:
    """開いた結果。`dataset` は呼び出し側のもの（基の複製）。

    kind: `directory` / `mztab_path` / `job_path` / `pipeline_path`
    path: 入口に渡されたパスの絶対パス
    adapter: `msdial_folder`（`.arf` / `.arf2` を直接）/ `mztab`（mzTab-M）
    selection: 解析フォルダのアダプタが選んだ `.arf` 一式（mzTab-M のアダプタは None）
    notes: どう開いたかの注記（読み込み自体の注記は `dataset.load_notes`）
    reused_base: 読み込み済みの基を使い回したか
    """

    dataset: DatasetState
    kind: str
    path: str
    adapter: str
    selection: AlignmentSelection | None
    notes: list[str] = field(default_factory=list)
    reused_base: bool = False


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
    from metabolomix.dataset.msdial_adapter import _MEASURES

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
    if measure not in _MEASURES:
        raise DomainError("DATASET_BAD_REQUEST",
                          f"measure は {' / '.join(_MEASURES)} のどちらかです: {measure!r}",
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
    if cached is not None and _signature(cached[0].dataset, input_path) == cached[1]:
        return _copy(cached[0], reused=True)
    base = _open(kind, str(raw), str(input_path), adapter, measure, allow_incomplete)
    handles.registry.put_base(key, base, _signature(base.dataset, input_path))
    return _copy(base, reused=False)


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


def _signature(ds, input_path: Path) -> tuple:
    """入力の指紋（読んだファイルの size と mtime）。中身が変われば基を使い回さない。

    解析フォルダのアダプタは値のほかに `.arf2` の隣のタグ（`_tags.xml`）とキュレーションの
    フラグ（`curation/flags.jsonl`）も読むので、それも数える（無いファイルは「無い」と数え、
    後から現れたら読み直す）。
    """
    paths = sorted({str(input_path), *(getattr(ds, "source_files", None) or {}),
                    *_folder_sidecars(ds)})
    out = []
    for path in paths:
        try:
            stat = Path(path).stat()
            out.append((path, stat.st_size, stat.st_mtime_ns))
        except OSError:
            out.append((path, None, None))
    return tuple(out)


def _folder_sidecars(ds) -> list[str]:
    if getattr(ds, "source_format", None) != "msdial_folder":
        return []
    from metabolomix.curation.flags import FLAGS_FILENAME, curation_dir
    from metabolomix.msdial.tags import alignment_tag_path

    out = []
    for path in getattr(ds, "source_files", None) or {}:
        if path.lower().endswith(".arf2"):
            out += [str(alignment_tag_path(path)), str(curation_dir(path) / FLAGS_FILENAME)]
    return out


# ---------- 読み込み ----------

def _open(kind: str, raw: str, path: str, adapter: str, measure: str,
          allow_incomplete: bool) -> OpenedDataset:
    from metabolomix.dataset import library_memory

    selection = None
    notes: list[str] = []
    if adapter == "mztab":
        ds = _load(kind, raw, allow_incomplete)
    elif kind == "directory":
        from metabolomix.dataset.msdial_adapter import load_msdial_folder

        selection = select_alignment(Path(raw).expanduser())
        ds = load_msdial_folder(selection, measure=measure)
    elif kind == "job_path":
        ds, selection, gate = _folder_via_job(Path(raw).expanduser(), measure=measure,
                                              allow_incomplete=allow_incomplete)
        notes.append(_folder_note(gate, "ジョブ"))
    else:
        ds, selection, notes = _folder_via_run(raw, measure=measure,
                                               allow_incomplete=allow_incomplete)

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
    return OpenedDataset(dataset=ds, kind=kind, path=path, adapter=adapter,
                         selection=selection, notes=notes)


def _load(kind: str, path: str, allow_incomplete: bool):
    """mzTab-M のアダプタ（今までの `dataset_load` / `dataset_handle_open` の読み込みそのまま）。"""
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
    """ジョブの門を通し、正準 mzTab-M の隣の解析フォルダを読む。(ds, selection, gate)。

    `.arf` は mzTab-M の evidence と同じ候補順（`arf_candidates`。handoff の
    `peak_matrix_source` → 隣の `*PeakProperties.arf` → 隣のその他）の先頭を使う。
    pipeline の run では `.arf` 一式が `input/` にあるので、隣に無ければそこを探す。
    """
    from metabolomix.dataset.msdial_adapter import load_msdial_folder
    from metabolomix.mztab import evidence, loading

    gate = loading.job_gate(job_path, allow_incomplete=allow_incomplete)
    job = gate.job
    artifact_paths: dict[str, list[str]] = {}
    for art in job.artifacts:
        artifact_paths.setdefault(art.role, []).append(
            str(loading.artifact_abs_path(job, getattr(art, "root", "run_dir"), art.path)))
    candidates = evidence.arf_candidates(gate.mztab_path, artifact_paths)
    target = candidates[0] if candidates else None
    if target is None and pipeline_root is not None:
        target = _pipeline_arf2(pipeline_root)
    if target is None:
        where = (f"ジョブの mzTab-M（{gate.mztab_path.name}）の隣"
                 + ("にも run の input/ にも" if pipeline_root is not None else "に"))
        raise DomainError(
            "ALIGNMENT_NOT_FOUND",
            f"{where} MS-DIAL のアラインメント（.arf）が見つかりません。{_MZTAB_HINT}",
            {"mztab": gate.mztab_path.name, "source_hint": "mztab"})
    try:
        selection = select_alignment(target)
    except DomainError as exc:
        message = exc.message if exc.message.endswith("。") else f"{exc.message}。"
        raise DomainError(exc.code, f"{message}{_MZTAB_HINT}",
                          {**(exc.details or {}), "source_hint": "mztab"}) from exc

    ds = load_msdial_folder(selection, measure=measure)
    ds.job_path = str(job_path)
    ds.exploratory_only = gate.exploratory_only
    ds.source_verification = "verified" if gate.verified else "legacy_unverified"
    if gate.incomplete_warning:
        ds.load_notes.insert(0, gate.incomplete_warning)
    if gate.unverified_warning:
        ds.load_notes.append(gate.unverified_warning)
    return ds, selection, gate


def _folder_via_run(pipeline_path: str, *, measure: str, allow_incomplete: bool):
    """pipeline の run を解析フォルダで読む。run から引き継ぐのは解決済みの試料対応表だけ。"""
    from metabolomix.analysis import sample_manifest
    from metabolomix.pipeline import session_handoff

    root, record, job_path = session_handoff.run_console_job(pipeline_path)
    ds, selection, gate = _folder_via_job(job_path, measure=measure,
                                          allow_incomplete=allow_incomplete, pipeline_root=root)
    rows = session_handoff._run_manifest_rows(root, record)
    if rows is not None:
        sample_manifest.apply_metadata(ds, sample_manifest.resolve_metadata(ds, rows))
    notes = [_folder_note(gate, "run"),
             "run の解析行列・feature binding・差次的結果は読んでいない（解析フォルダの"
             "アダプタが run から引き継ぐのは解決済みの試料対応表だけ）。それらも載せるなら"
             " source=\"mztab\" で読む。"]
    return ds, selection, notes


def _folder_note(gate, where: str) -> str:
    return (f"{where}の mzTab-M（{gate.mztab_path.name}）と同じアラインメントの解析フォルダを"
            "読んだ。mzTab-M を読むなら source=\"mztab\"。")


def _pipeline_arf2(pipeline_path) -> Path | None:
    """pipeline の run の `input/AlignResult-*.arf2`。無ければ None。

    run では mzTab が `console/attempt-*/msdial/` の下、`.arf`・`.arf2`・`.dcl` が
    `input/` にあるので、mzTab の隣を探す既定の経路では見つからない。`AlignResult-*.arf2`
    が複数あるときは名前が最も新しい（降順で先頭の）ものを選ぶ（決定的にするため。
    `AlignResult-<日時など>` は名前の順が新しさの順になる）。`pipeline_path` は
    `pipeline-run.json` でもその親でもよい。
    """
    from metabolomix.pipeline.recovery import normalize_pipeline_root

    root = normalize_pipeline_root(Path(pipeline_path))
    candidates = sorted((root / "input").glob("AlignResult-*.arf2"),
                        key=lambda p: p.name, reverse=True)
    return next((p for p in candidates if p.is_file()), None)
