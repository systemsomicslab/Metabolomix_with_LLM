"""pipeline関連テストの合成fixture置き場。

実rawや既存成果物をfixtureにしない（CLAUDE.md方針）。テストごとに独自の形を作らせない
ための唯一の正準とする。Task 3・6・9・13・19がここへhelperを足していく前提なので、
既存helperの必須フィールドは減らさない。

中身は3層。

1. 合成入力ビルダー（`make_source` / `write_mztab` / `make_dataset` /
   `metadata_rows` / `write_sample_manifest`）——純粋なファイル・dictの生成。
2. 偽Console起動の差し替え口（`fake_console_command` / `use_fake_console` /
   `mztab_text`）——「何を起動するか」だけを替える。
3. Task 19の受入harness（`PipelineHarness` / `read_contract_tsv`）——受付から
   実子プロセスworkerまでを本物どおり回し、状態と実成果物だけを観測する。
"""
from __future__ import annotations

import csv
import io
import json
import os
import signal
import subprocess
import sys
import time
from pathlib import Path

from lipidmix.core.process_control import resolve_python_launch

#: このcheckoutのルート（`tests/` の1階層上）。子プロセスのcwdをここへ固定する。
_REPO_ROOT = Path(__file__).resolve().parents[1]


def make_source(root: Path, *, extension: str = ".wiff", kind: str = "file") -> dict:
    """Task 13 (入力隔離・メソッド選択) 用の合成生データフォルダを作る。

    root直下にS0..S7.wiff・GUI自動保存形式のlab_param_202609050001.txt
    （Ion mode/Target omics/Lbm file path）・fake.lbm2・fake.exeを作り、
    それぞれのPathを返す。fake.exeは実行しない（中身は placeholder テキスト）
    ——単体テストではlipidmix.console.runner.is_console_exeを差し替えて使う。

    Parameters
    ----------
    extension, kind:
        Task 19 (D08「ファイルraw、ディレクトリraw」) 用。`kind="dir"`にすると
        1測定単位が**フォルダ**（Agilent `.d` 等）になり、各フォルダの中に
        内部ファイルを1件置く。既定は従来どおりファイル形式の`.wiff`で、
        既存の呼び出し側の意味は変えない。本数は常に8件（8検体×2群という
        下流の前提に合わせてある）。
    """
    root = Path(root)
    root.mkdir(parents=True, exist_ok=True)
    if kind not in ("file", "dir"):
        raise ValueError(f"kindは'file'か'dir'です: {kind!r}")
    for i in range(8):
        entry = root / f"S{i}{extension}"
        if kind == "dir":
            entry.mkdir(exist_ok=True)
            (entry / "AcqData.bin").write_text(f"raw-{i}", encoding="ascii")
        else:
            entry.write_text(f"raw-{i}", encoding="ascii")

    lbm = root / "fake.lbm2"
    lbm.write_bytes(b"fake-lbm-library")

    method = root / "lab_param_202609050001.txt"
    method.write_text(
        "Ion mode: negative\n"
        "Target omics: Lipidomics\n"
        "Lbm file path: fake.lbm2\n",
        encoding="ascii", newline="\n")

    exe = root / "fake.exe"
    exe.write_text("fake console executable placeholder", encoding="ascii")

    return {"root": root, "method": method, "lbm": lbm, "exe": exe}


def execution_record(**overrides):
    """console-execution.v1 の妥当な最小レコードを返す。

    overridesで一部フィールドを壊し、validate_execution_recordの拒否条件を
    テストするために使う。
    """
    record = {
        "schema": "console-execution.v1", "execution_id": "exec-1", "job_id": "job-1",
        "started_at": "2026-09-05T00:00:00Z", "ended_at": "2026-09-05T00:00:01Z",
        "pid": 1234, "process_identity": {"pid": 1234, "creation_time": 100},
        "command_sha256": "a" * 64, "method_sha256": "b" * 64, "exe_sha256": "c" * 64,
        "exit_code": 0, "termination": "exited", "timeout_s": 21600,
        "collection": {"status": "succeeded"}, "validation": {"status": "pending"},
    }
    record.update(overrides)
    return record


# ---------- write_mztab ----------

# 少なくとも3特徴、識別済み(with_inchikey時)・InChIKeyを想定した合成値。
# 実データの脂質名文法に寄せる必要はない（構造検証のための合成値のため）。
_FEATURES = (
    ("PC 34:1", "IPCSVZSSVZVIGE-UHFFFAOYSA-N", 760.5851, 620.0),
    ("PE 36:2", "XKMRRTOUMJRJIA-UHFFFAOYSA-N", 742.5395, 590.0),
    ("TG 52:3", "DGGXCMYPQAOAJC-UHFFFAOYSA-N", 878.7789, 910.0),
)


def write_mztab(path: Path, sources: list[Path], *, with_inchikey: bool = True) -> Path:
    """実際のparse_mztab/validate_mztabに通る合成mzTab-M 2.0ファイルを書く。

    sources の各要素を ms_run[N] として file URI 化し、assay[N]-ms_run_ref で
    結ぶ（1 source = 1 assay）。tests/test_mztab_tools.py の SFH 構文（SFH ヘッダ
    行 + SMF データ行）を再利用し、少なくとも3特徴・SML・SMEを含む。文字列
    "x" のような偽物ではなく、実物の構造検証・定量行列抽出が通る内容にする。

    同定（InChIKey・化合物名）は **SME セクション**に書き、SMF からは
    `SME_ID_REFS` で参照する。mzTab-M 2.0.0-M では構造・名称が SME にしか無く、
    `lipidmix.mztab.dataset_state` も SME 側だけを読むため——SMF 列に書いても
    `feature_metadata["inchikey"]` は全件 None になる。

    Parameters
    ----------
    with_inchikey:
        False なら SME の database_identifier / chemical_name を全行 "null"
        （未識別）にする。spec E01「InChIKey 0 件」の入力。
    """
    path = Path(path)
    n = len(sources)
    if n < 1:
        raise ValueError("sources には少なくとも1件のPathが必要です")

    lines = [
        "MTD\tmzTab-version\t2.0.0-M",
        "MTD\tmzTab-mode\tComplete",
        "MTD\tmzTab-type\tQuantification",
        "MTD\tsoftware[1]\t[MS, MS:1003082, MS-DIAL, Msdial console 5.5.241113]",
        "MTD\tquantification_method\t[MS, MS:1001829, label-free raw feature quantitation, ]",
        "MTD\tsmall_molecule-quantification_unit\t[PRIDE, PRIDE:0000429, Abundance, ]",
    ]
    for i, src in enumerate(sources, start=1):
        uri = Path(src).resolve().as_uri()
        lines.append(f"MTD\tms_run[{i}]-location\t{uri}")
        lines.append(f"MTD\tassay[{i}]-ms_run_ref\tms_run[{i}]")
        lines.append(f"MTD\tassay[{i}]\tS{i}")

    abundance_cols = "\t".join(f"abundance_assay[{i}]" for i in range(1, n + 1))
    # SMFは`SME_ID_REFS`で同定証拠（SME）を指す。**構造・名称はSMEにしか無い**のが
    # mzTab-M 2.0.0-Mの実際で、`lipidmix.mztab.dataset_state`もSME側だけを読む
    # （SMFのdatabase_identifier/chemical_nameは読まない）。ここでSMEへ書かないと
    # `with_inchikey=True`でもfeature_metadataのinchikeyが全件Noneになり、
    # 契約TSVのエクスポートがNO_ANNOTATED_FEATURESで落ちる。
    lines.append(
        "SFH\tSMF_ID\tSME_ID_REFS\tSML_ID_REFS\tdatabase_identifier\tchemical_name\t"
        f"smiles\tinchi\texp_mass_to_charge\tretention_time_in_seconds\t{abundance_cols}"
    )
    for idx, (name, inchikey, mz, rt) in enumerate(_FEATURES, start=1):
        identifier = inchikey if with_inchikey else "null"
        abundances = "\t".join(f"{1000.0 * idx + 10.0 * j:.1f}" for j in range(1, n + 1))
        lines.append(
            f"SMF\t{idx}\t{idx}\tSML:{idx}\t{identifier}\t{name}\tnull\tnull\t"
            f"{mz}\t{rt}\t{abundances}"
        )

    # アダクトは正極性1・負極性1・不明1件にして多数決を意図的に引き分け(None)に
    # する。極性を断定すると、write_mztab を負極性jobで使うテストと衝突する
    # （validate_outputsの極性多数決チェックが致命的不一致として弾いてしまう）。
    _sml_adducts = ("[M+H]1+", "[M-H]1-", "null")
    lines.append("SMH\tSML_ID\tSMF_ID_REFS\tdatabase_identifier\tadduct_ions")
    for idx, adduct in enumerate(_sml_adducts[: len(_FEATURES)], start=1):
        lines.append(f"SML\t{idx}\t{idx}\tnull\t{adduct}")

    lines.append("SEH\tSME_ID\tSMF_ID_REFS\tdatabase_identifier\tchemical_name\t"
                 "smiles\tinchi\trank")
    for idx, (name, inchikey, _mz, _rt) in enumerate(_FEATURES, start=1):
        identifier = inchikey if with_inchikey else "null"
        chemical_name = name if with_inchikey else "null"
        lines.append(f"SME\t{idx}\tSMF:{idx}\t{identifier}\t{chemical_name}\tnull\tnull\t1")

    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


# ---------- 実プロセスとして走る偽 Console ----------
# supervise は「本当に起動したプロセスが本当にどう終わったか」を証跡にする。
# subprocess をモックで置き換えると、監視・停止・収集のどれが壊れても緑のままに
# なるため、実行系のテストは実際に子プロセスを起こす。tests/fixtures/fake_console.py
# は決まったシナリオ用（import せずスクリプトパスとして渡す）で、こちらは
# 「任意のファイルを書いて任意の終了コードで終わる」だけの汎用版。

def fake_console_command(files: dict | None = None, *, exit_code: int = 0,
                         sleep_s: float = 0.0, binary: dict | None = None) -> list[str]:
    """指定のファイルを書き、必要なら待ってから、指定の終了コードで終わるコマンド。

    `binary` はバイト列をそのまま書くファイル（`.arf` のような MessagePack 系）。
    テキストとして書くと壊れ、読み側が「読めなかった」ことを警告に落として
    先へ進むため、**証拠が無いまま工程が成功する**——偽 Console の都合が
    検証結果に化けないよう、バイナリはバイナリのまま渡す。
    """
    payload = [(str(path), text) for path, text in
               sorted((files or {}).items(), key=lambda kv: str(kv[0]))]
    blobs = [(str(path), bytes(data)) for path, data in
             sorted((binary or {}).items(), key=lambda kv: str(kv[0]))]
    script = (
        "import pathlib, sys, time\n"
        f"for path, text in {payload!r}:\n"
        "    p = pathlib.Path(path)\n"
        "    p.parent.mkdir(parents=True, exist_ok=True)\n"
        "    p.write_text(text, encoding='utf-8')\n"
        f"for path, data in {blobs!r}:\n"
        "    p = pathlib.Path(path)\n"
        "    p.parent.mkdir(parents=True, exist_ok=True)\n"
        "    p.write_bytes(data)\n"
        f"time.sleep({sleep_s!r})\n"
        f"sys.exit({exit_code!r})\n"
    )
    return [sys.executable, "-c", script]


def use_fake_console(monkeypatch, command: list[str]) -> None:
    """Console のコマンドライン組み立てだけを偽 Console へ差し替える。

    差し替えるのは「何を起動するか」だけで、起動・監視・収集・完了判定は本物を
    通る。`supervise` の `command` 注入は MCP の公開引数にできない（任意コマンドの
    実行口になる）ので、テストからは唯一の組み立て場所である build_msdial_cmd を
    差し替える。
    """
    monkeypatch.setattr("lipidmix.console.runner.build_msdial_cmd",
                        lambda *args, **kwargs: command)


def mztab_text(tmp_path: Path, sources) -> str:
    """構造検証と定量抽出を通る合成 mzTab-M の中身を返す。

    完了ゲートは「拡張子が .mzTab である」ことではなく、主 mzTab が一意に選べ・
    構造が妥当で・定量行列に有限値があり・予定した入力が全て assay に対応して
    いることを要求する。偽 Console にはこの中身を書かせる。
    """
    scratch = Path(tmp_path) / "_mztab_template"
    text = write_mztab(scratch, list(sources)).read_text(encoding="utf-8")
    scratch.unlink()
    return text


# ---------- DatasetState ----------

def make_dataset():
    """解析サービスのテスト用に、8 検体 × 6 特徴の DatasetState を作る。

    実 mzTab を読まずに済ませる（実データ・既存成果物を fixture にしない方針）。
    値は固定 seed の一様乱数で、前処理・PCA・2 群比較が実際に走る規模にしてある。
    """
    import numpy as np

    from lipidmix.mztab.dataset_state import DatasetState
    ds = DatasetState()
    ds.feature_matrix = np.random.default_rng(7).uniform(100, 200, (6, 8))
    ds.sample_names = [f"S{i}" for i in range(8)]
    ds.sample_assay_ids = [f"assay[{i + 1}]" for i in range(8)]
    ds.feature_ids = [f"F{i}" for i in range(6)]
    ds.quantification_measure = "peak_height"
    ds.validation_result = {"ok": True, "errors": [], "warnings": []}
    ds.feature_metadata = {f"F{i}": {"name": f"Lipid {i}", "mz": 500 + i,
                                     "rt": 2.0,
                                     "inchikey": "IPCSVZSSVZVIGE-UHFFFAOYSA-N",
                                     "inchikey_source": "database_identifier"}
                           for i in range(6)}
    return ds


# ---------- sample-manifest.v1 の合成行（Task 9 apply_metadata の契約） ----------

def metadata_rows(ds, *, n_qc=4, confirmed=True):
    """make_dataset()の8検体へ、resolve_metadata出力と同じ形の行を合成する。

    n_qcで健全なQC数（0〜4件）を振り、conservative-v1のQC閾値（Task 11）を
    跨いだ挙動をテストできるようにする。confirmed=Falseはbatch/injection_order/
    qc_poolだけをmzTab由来unverifiedへ落とし、他タスクが「出所で挙動が変わる」
    経路を確認するための契約。
    """
    assert len(ds.sample_names) == 8 and 0 <= n_qc <= 4
    qc_orders = [1, 3, 6, 8][:n_qc]
    orders = qc_orders + [i for i in range(1, 9) if i not in qc_orders]
    rows = []
    for i, name in enumerate(ds.sample_names):
        row = dict(sample_id=name, source_file=f"{name}.wiff",
                   role="qc" if i < n_qc else "sample", group=None,
                   batch="B1", injection_order=orders[i],
                   qc_pool="pool1" if i < n_qc else None, include=True)
        row["provenance"] = {
            key: {"value": value, "source": "user_manifest", "confidence": "confirmed"}
            for key, value in row.items()
        }
        if not confirmed:
            for key in ("batch", "injection_order", "qc_pool"):
                row["provenance"][key].update(source="mztab", confidence="unverified")
        row["conflicts"] = []
        rows.append(row)
    return rows


# ---------- Task 19: 実子プロセスworkerで一気通貫を回すharness ----------
#
# 偽物にするのは**workerが起動するConsoleのコマンドだけ**。受付（plan/start/
# resume）・engine・loading・metadata・preprocessing・differential・export・
# renderはすべて本物を通す。差し替え点は2つに限定してある。
#
# 1. `service.launch_pipeline_worker` — 起動する実行ファイルだけを
#    `tests/pipeline_worker_harness.py` へ替える。`start_pipeline`/`resume_pipeline`
#    自身の受付・再利用判定・handshakeは本物のまま通る（本番workerは試験用の
#    環境変数・任意import名を一切受け付けないため、差し替えるならここしかない）。
#    起こしたPopenは全部harnessが持ち、`close()`で確実に始末する。
# 2. harness worker内での `lipidmix.console.runner.build_msdial_cmd` —
#    偽Consoleスクリプト（`tests/fixtures/fake_console.py`）を起こすコマンドへ替える。
#    起動・監視・停止・収集・完了判定はすべて本物の`supervise`を通る。

_HARNESS_WORKER = Path(__file__).resolve().parent / "pipeline_worker_harness.py"

#: `wait`の観測上限（秒。brief拘束）。超えたら所有workerだけを止めて診断を出す。
WAIT_TIMEOUT_S = 20.0
_WAIT_POLL_S = 0.05

#: sample-manifest.v1 の8列（`lipidmix.analysis.sample_manifest.FIELDS`と同じ順）。
_MANIFEST_FIELDS = ("sample_id", "source_file", "role", "group", "batch",
                    "injection_order", "qc_pool", "include")

#: 8検体を4:4へ割る既定の群割当（S0..S3=control, S4..S7=treated）。
#: `write_mztab`の合成強度は assay 番号とともに単調増加するので、treated
#: （後半4検体）が必ず高い＝log2FCが正になる——「向き」を実数値でassertできる
#: ようにするための意図的な設計。
DEFAULT_GROUPS = ("control",) * 4 + ("treated",) * 4

#: 既定の比較（reference=control, test=treated）。正のlog2FCはtest群が高い方向
#: （Global Constraints / `export_contract.LOG2FC_SIGN`）。
DEFAULT_COMPARISON = {
    "comparison_id": "treated_vs_control",
    "reference_group": "control",
    "test_group": "treated",
}


def write_sample_manifest(path: Path, *, source_names, groups=DEFAULT_GROUPS,
                          roles=None) -> Path:
    """sample-manifest.v1 の完全なシート（全予定rawに1行ずつ）を書く。

    `source_names`は`source_root`基準の相対名（`S0.wiff` 等）。`parse_manifest`は
    列の過不足・行欠落・未知source_fileを一切許さないので、ここでは常に全件・
    8列・正しい見出しで書く（壊れたシートを使いたいテストは自分で書く）。
    """
    path = Path(path)
    names = list(source_names)
    roles = list(roles) if roles is not None else ["sample"] * len(names)
    groups = list(groups)
    assert len(groups) == len(names) and len(roles) == len(names)
    lines = ["# schema = sample-manifest.v1", "\t".join(_MANIFEST_FIELDS)]
    for index, name in enumerate(names):
        lines.append("\t".join([
            f"S{index}", name, roles[index], groups[index], "B1",
            str(index + 1), "", "true",
        ]))
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
    return path


def read_contract_tsv(path):
    """契約TSVを実際に読み戻し、(メタ行dict, DictReaderの行list, 列名list)を返す。

    「拡張子がTSVで中身は見ない」を避けるための唯一の読み手。メタ行（`#`始まり）と
    データ行を分け、データ行は`csv.DictReader`へ通す——列名・列順・空欄の保持を
    自分でパースし直して確かめる（`export_contract.EXPORT_COLUMNS`は別リポジトリ
    massbank-context との契約）。
    """
    text = Path(path).read_text(encoding="utf-8")
    meta: dict = {}
    for line in text.splitlines():
        if not line.startswith("#"):
            continue
        for field in line.lstrip("# ").split("\t"):
            if " = " in field:
                key, _, value = field.partition(" = ")
                meta[key.strip()] = value.strip()
    data_lines = [ln for ln in text.splitlines() if ln and not ln.startswith("#")]
    reader = csv.DictReader(io.StringIO("\n".join(data_lines)), delimiter="\t")
    rows = list(reader)
    return meta, rows, reader.fieldnames


class PipelineHarness:
    """`plan/start/resume` → 実子プロセスworker → 永続状態、を回す受入用harness。

    `wait`は**観測だけ**を行う（`pipeline_status`も`read_status`も呼ばず、
    `pipeline-run.json`を直接読む——「statusを一度も呼ばずに完走する」ことを
    確かめる経路を、harness自身が壊さないため）。

    **`wait`は「所有workerの終了」を待ってから1回だけ読む**（状態ファイルを
    ポーリングしない）。Windowsでは`pipeline-run.json`を読み込みで開いている
    だけで、workerの`atomic_write_json`の`os.replace`が
    `PermissionError(WinError 5)`で失敗する（実測。`FILE_SHARE_DELETE`を付けて
    開いても同じ）。つまり**観測者がworkerを落としうる**——本番の
    `store.load_run`/`recovery.read_status`も同じ開き方なので、これはharnessの
    都合ではなく製品側の危険であり、Task19のreportへ「発見」として上げてある。
    harness側は、状態ファイルへ触るのを「workerが確実に居なくなった後の1回」に
    限ることで、この危険を自分では踏まない。
    """

    def __init__(self, tmp_path, monkeypatch):
        self.tmp_path = Path(tmp_path)
        self.monkeypatch = monkeypatch
        self._sources: dict = {}
        self._workers: list = []          # [(pipeline名, Popen, ログファイルobject)]
        self._detached: list = []         # [(pipeline名, launch_detachedの戻り値)]
        self._launch_options: dict = {}   # pipeline_root名 → 起動引数
        self._pending_options: dict = {"scenario": "success"}
        self.counters_dir = self.tmp_path / "_fake_console_counters"
        self.logs_dir = self.tmp_path / "_worker_logs"
        self.counters_dir.mkdir(parents=True, exist_ok=True)
        self.logs_dir.mkdir(parents=True, exist_ok=True)
        self._configured = False

    # ---------- 環境 ----------

    def _configure(self) -> None:
        """受付層が要る環境（MSDIAL_EXE・console判定・受付索引）を一度だけ整える。"""
        if self._configured:
            return
        from lipidmix.pipeline import service

        exe = self.tmp_path / "fake-msdialcui.exe"
        exe.write_text("fake console executable placeholder", encoding="ascii")
        self.monkeypatch.setenv("MSDIAL_EXE", str(exe))
        self.monkeypatch.setattr("lipidmix.console.runner.is_console_exe",
                                 lambda *a, **k: True)
        # 受付索引はユーザーのLOCALAPPDATAではなくtmpへ。テスト間で共有させない。
        self.monkeypatch.setenv("LIPIDMIX_PIPELINE_INDEX_DIR",
                                str(self.tmp_path / "_pipeline_index"))
        # 起動するworkerだけを差し替える（受付・handshake・再利用判定は本物）。
        self.monkeypatch.setattr(service, "launch_pipeline_worker", self._launch_worker)
        # handshakeの意味論は`tests/test_pipeline_service.py`の単体testの担当。
        # ここで実workerのimport（numpy/scipy/matplotlib）を待つ意味は無い。
        self.monkeypatch.setattr(service, "_HANDSHAKE_TIMEOUT_S", 0.2)
        self.monkeypatch.setattr(service, "_HANDSHAKE_POLL_S", 0.02)
        self._configured = True

    def source(self, name: str = "source", **kwargs) -> dict:
        """合成生データフォルダを（無ければ作って）返す。名前ごとに1つ。"""
        if name not in self._sources:
            self._sources[name] = make_source(self.tmp_path / name, **kwargs)
        return self._sources[name]

    # ---------- パス補助 ----------

    @staticmethod
    def pipeline_root(run) -> Path:
        """`pipeline-run.json`でもpipeline_rootでも受け取ってrootを返す。"""
        path = Path(run)
        return path.parent if path.name == "pipeline-run.json" else path

    @staticmethod
    def run_json(pipeline_root) -> Path:
        return Path(pipeline_root) / "pipeline-run.json"

    def counter_path(self, run) -> Path:
        return self.counters_dir / f"{self.pipeline_root(run).name}.tsv"

    def record(self, run) -> dict:
        """`pipeline-run.json`をそのまま読む（読取専用・statusツールを通さない）。"""
        return json.loads(
            self.run_json(self.pipeline_root(run)).read_text(encoding="utf-8"))

    # ---------- 起動 ----------

    def _options_for(self, pipeline_root: Path, overrides: dict | None = None) -> dict:
        merged = dict(self._launch_options.get(pipeline_root.name)
                      or self._pending_options)
        merged.update(overrides or {})
        self._launch_options[pipeline_root.name] = merged
        return merged

    def _worker_command(self, pipeline_root: Path, options: dict) -> list:
        command = [sys.executable, str(_HARNESS_WORKER),
                   "--pipeline", str(pipeline_root),
                   "--scenario", options.get("scenario", "success"),
                   "--counter", str(self.counter_path(pipeline_root))]
        if options.get("no_inchikey"):
            command.append("--no-inchikey")
        if options.get("sleep_stage"):
            command += ["--sleep-stage", str(options["sleep_stage"]),
                        "--sleep-seconds", str(options.get("sleep_seconds", 0.0))]
        return command

    def _launch_worker(self, path) -> dict:
        """`service.launch_pipeline_worker`の差し替え先。実子プロセスを起こす。

        stdout/stderrは**実ファイルハンドル**へ結ぶ。パイプにすると、起動元が
        子の生存期間ぶん待たされる経路（Task2の実測不具合。`bInheritHandles=TRUE`
        で継承されたstdoutパイプ）を harness 自身が作ってしまう。
        """
        pipeline_root = self.pipeline_root(path)
        options = self._options_for(pipeline_root)
        command = self._worker_command(pipeline_root, options)
        log_path = self.logs_dir / f"{pipeline_root.name}-{len(self._workers)}.log"
        log = open(log_path, "ab")
        # venvのリダイレクタを挟まない。挟むと`proc`はリダイレクタになり、
        # `proc.wait()`が返った時点で本物のworker（engineがidentityを刻んだ
        # プロセス）はまだ終了処理中でありうる（GitHub #3）。
        command, env = resolve_python_launch(command)
        proc = subprocess.Popen(
            command, env=env, cwd=str(_REPO_ROOT), stdin=subprocess.DEVNULL,
            stdout=log, stderr=subprocess.STDOUT)
        self._workers.append((pipeline_root.name, proc, log))
        return {"launched": True, "pid": proc.pid, "log_path": str(log_path)}

    def set_launch_options(self, run, **options) -> dict:
        """このpipelineの次回worker起動引数を上書きする（scenario・sleep等）。"""
        return self._options_for(self.pipeline_root(run), options)

    def worker_processes(self, run) -> list:
        """このpipelineに対してharnessが起こしたworkerのPopen一覧（起動順）。"""
        name = self.pipeline_root(run).name
        return [proc for key, proc, _log in self._workers if key == name]

    def start(self, *, target="exploratory", comparisons=None, scenario="success",
              save_project=False, source_name="source", request_id=None,
              launch=True, extra_request=None, timeout_s=600, **launch_options):
        """`service.start_pipeline`（`launch=False`なら`plan_pipeline`）を本物どおり
        呼び、`pipeline-run.json`のPathを返す。

        `timeout_s`はConsoleの実行上限（要求フィールド）。`hang`シナリオの
        timeout検証で短くする。
        """
        self._configure()
        from lipidmix.pipeline import service

        source = self.source(source_name)
        request = {"target": target, "save_project": save_project,
                   "comparisons": list(comparisons or []), "timeout_s": timeout_s}
        request.update(extra_request or {})
        self._pending_options = {"scenario": scenario, **launch_options}
        entry = service.start_pipeline if launch else service.plan_pipeline
        receipt = entry(source["root"], request, request_id=request_id)
        pipeline_root = Path(receipt["pipeline_path"])
        self._options_for(pipeline_root)
        return self.run_json(pipeline_root)

    def start_scenario(self, scenario: str, **kwargs):
        """受入表のシナリオ1件を、既定の探索目標で起動する。"""
        return self.start(scenario=scenario, **kwargs)

    def relaunch(self, run, **options):
        """既存runに対してもう1つworkerを起こす（再開・並行・喪失の検証用）。"""
        pipeline_root = self.pipeline_root(run)
        self._options_for(pipeline_root, options)
        return self._launch_worker(pipeline_root)

    def launch_detached(self, run, **options) -> dict:
        """起動役プロセスを別に立て、そこから切り離してworkerを起こす。

        A05「MCP接続を切断、statusを一度も呼ばない」の再現。pytestプロセスから
        直接切り離すと「起動元は生きたまま」になってしまうので、**起動役そのものを
        短命な子プロセスにし、その終了を待ってから**workerの完走を観測する。
        """
        pipeline_root = self.pipeline_root(run)
        options = self._options_for(pipeline_root, options)
        info_path = self.logs_dir / f"{pipeline_root.name}-launch-info.json"
        command = self._worker_command(pipeline_root, options) + [
            "--launch-detached", "--launch-info", str(info_path)]
        completed = subprocess.run(
            command, cwd=str(_REPO_ROOT), stdin=subprocess.DEVNULL,
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=60)
        assert completed.returncode == 0, "切り離し起動役が異常終了しました"
        info = json.loads(info_path.read_text(encoding="utf-8"))
        self._detached.append((pipeline_root.name, info))
        return info

    # ---------- 観測 ----------

    def await_workers(self, run, *, timeout: float = WAIT_TIMEOUT_S) -> bool:
        """このpipelineへharnessが起こしたworkerが全部終わるまで待つ（観測のみ）。

        状態ファイルには一切触れない——待つ相手は自分の子プロセス（`Popen.wait`）と、
        切り離したworkerのprocess identity（`same_process`＝`OpenProcess`）だけ。
        """
        from lipidmix.core.process_control import same_process

        deadline = time.monotonic() + timeout
        name = self.pipeline_root(run).name
        for key, proc, _log in list(self._workers):
            if key != name:
                continue
            try:
                proc.wait(timeout=max(0.0, deadline - time.monotonic()))
            except subprocess.TimeoutExpired:
                return False
        for key, info in list(self._detached):
            if key != name:
                continue
            while same_process(info.get("identity") or {}):
                if time.monotonic() >= deadline:
                    return False
                time.sleep(_WAIT_POLL_S)
        return True

    def wait(self, run, *, expected: str, timeout: float = WAIT_TIMEOUT_S) -> dict:
        """所有workerの終了を待ってから`pipeline-run.json`を1回だけ読む。

        上限で打ち切ったら診断を出し、**harnessが起こしたworkerだけ**を止める。
        他プロセスへは一切触れない。
        """
        pipeline_root = self.pipeline_root(run)
        finished = self.await_workers(pipeline_root, timeout=timeout)
        record = self.record(pipeline_root)
        if record.get("status") == expected:
            return record
        diagnostics = self.diagnostics(pipeline_root, record)
        self.stop_workers()
        reason = ("worker終了後のstatusが期待と違います" if finished
                  else f"{timeout}秒以内にworkerが終了しませんでした")
        raise AssertionError(f"{reason}（expected={expected!r}）。\n{diagnostics}")

    def stage_marker(self, run, stage_id: str) -> Path:
        """`--sleep-stage`で待たせているstageへworkerが入った目印のパス。"""
        from tests.pipeline_worker_harness import stage_marker_path
        return stage_marker_path(self.pipeline_root(run), stage_id)

    def wait_for_stage(self, run, stage_id: str, *, timeout: float = WAIT_TIMEOUT_S) -> None:
        """workerが指定stageへ入るまで待つ（状態ファイルには触れない）。"""
        marker = self.stage_marker(run, stage_id)
        self.wait_for(marker.exists, timeout=timeout,
                      what=f"workerがstage {stage_id!r} へ到達すること")

    def wait_for(self, predicate, *, timeout: float = WAIT_TIMEOUT_S,
                 what: str = "条件") -> None:
        """任意の観測条件が成立するまで待つ（取消・喪失など状態以外の待ち合わせ用）。"""
        deadline = time.monotonic() + timeout
        while time.monotonic() < deadline:
            if predicate():
                return
            time.sleep(_WAIT_POLL_S)
        self.stop_workers()
        raise AssertionError(f"{timeout}秒以内に{what}が成立しませんでした。")

    def diagnostics(self, run, record: dict | None = None) -> str:
        pipeline_root = self.pipeline_root(run)
        record = record if record is not None else self.record(pipeline_root)
        stages = {sid: s.get("status") for sid, s in (record.get("stages") or {}).items()}
        lines = [
            f"pipeline_root = {pipeline_root}",
            f"status        = {record.get('status')!r}",
            f"needs_input   = {record.get('needs_input')}",
            f"stages        = {stages}",
            f"warnings      = {record.get('warnings')}",
        ]
        for key, proc, _log in self._workers:
            lines.append(f"worker[{key}] pid={proc.pid} returncode={proc.poll()}")
        for log_path in sorted(self.logs_dir.glob("*.log")):
            text = log_path.read_text(encoding="utf-8", errors="replace").strip()
            if text:
                lines.append(f"--- {log_path.name} ---\n{text[-2000:]}")
        return "\n".join(lines)

    def console_start_count(self, run) -> int:
        """偽Consoleが**実際に起動した回数**（起動ごとに1行が追記される）。"""
        return len(self.console_launch_rows(run))

    def console_launch_rows(self, run) -> list:
        path = self.counter_path(run)
        if not path.exists():
            return []
        return [ln.split("\t") for ln in path.read_text(encoding="utf-8").splitlines()
                if ln.strip()]

    def console_run_dir(self, run, *, attempt: int | None = None) -> Path:
        """pipelineが所有するConsole job のrun_dir（`_CONSOLE_RUN_SUBDIR`規約）。

        Console実行はattemptごとに別ディレクトリ（`console/attempt-NNNN/`）を
        持つ——`rerun_upstream=True`が前回の終了証跡・ログ・復旧用スナップ
        ショットを上書きしないため。`attempt`を省略すると**最新のattempt**を
        返す（`console/`直下しか無い古い形も一応受ける）。
        """
        base = self.pipeline_root(run) / "console"
        if attempt is not None:
            return base / f"attempt-{attempt:04d}"
        attempts = sorted(p for p in base.glob("attempt-*") if p.is_dir())
        return attempts[-1] if attempts else base

    def console_receipt(self, run) -> dict:
        """ディスク上の終了証跡（console-execution.v1）をそのまま読む。"""
        path = self.console_run_dir(run) / "execution-result.json"
        return json.loads(path.read_text(encoding="utf-8"))

    def console_job(self, run):
        """ディスク上のanalysis-job（収集された成果物一覧を読むため）。"""
        from lipidmix.console.job_manager import JOB_FILENAME, load_job
        return load_job(self.console_run_dir(run) / JOB_FILENAME)

    def output_ref(self, record: dict, kind: str) -> dict:
        """`output_name == kind`の最新ref（`results`はappend-only）。"""
        matches = [r for r in record.get("results") or []
                   if isinstance(r, dict) and r.get("output_name") == kind]
        if not matches:
            names = sorted({r.get("output_name") for r in record.get("results") or []})
            raise AssertionError(f"output_name={kind!r}のrefがありません（あるのは {names}）")
        return matches[-1]

    def output_path(self, record: dict, kind: str) -> str:
        """`output_name == kind`の最新refが指す実ファイルの絶対パス。"""
        root = Path(record["identity"]["pipeline_root"])
        return str(root / self.output_ref(record, kind)["relative_path"])

    def result_data(self, record: dict, kind: str) -> dict:
        """data形式で永続化された結果（PCA・差次的解析）の中身を読み戻す。"""
        payload = json.loads(
            Path(self.output_path(record, kind)).read_text(encoding="utf-8"))
        return payload["data"]

    # ---------- 訂正・再開 ----------

    def write_manifest(self, *, source_name="source", groups=DEFAULT_GROUPS,
                       name="sample-manifest.tsv") -> Path:
        """元フォルダへ完全なシートを置く（既定名なので自動探索でも拾われる）。"""
        source = self.source(source_name)
        root = Path(source["root"])
        raw_names = sorted(p.name for p in root.iterdir()
                           if p.name[0] == "S" and p.name[1:2].isdigit())
        return write_sample_manifest(root / name, source_names=raw_names, groups=groups)

    def resume_with_groups(self, run, *, comparisons=None, source_name="source",
                           groups=DEFAULT_GROUPS, **launch_options):
        """8検体control4/treated4の完全なシートと向き付き比較を与えて再開する。

        `service.resume_pipeline`を本物どおり呼ぶ（起動されるworkerだけharness版）。
        """
        from lipidmix.pipeline import service

        pipeline_root = self.pipeline_root(run)
        manifest = self.write_manifest(source_name=source_name, groups=groups)
        self._options_for(pipeline_root, launch_options)
        updates = {
            "sample_manifest": manifest.name,
            "comparisons": list(comparisons if comparisons is not None
                                else [DEFAULT_COMPARISON]),
        }
        return service.resume_pipeline(pipeline_root, updates=updates)

    # ---------- 後始末 ----------

    def stop_workers(self) -> None:
        """harnessが起こしたworkerだけを止める（他プロセスへは一切触らない）。"""
        for _key, proc, _log in self._workers:
            if proc.poll() is None:
                proc.kill()
        for _key, proc, log in self._workers:
            try:
                proc.wait(timeout=30)
            except subprocess.TimeoutExpired:
                pass
            if not log.closed:
                log.close()

    def survivors(self) -> list:
        """harnessが起こしたプロセスのうち、今も生きているもの。

        「子プロセスを漏らしていない」を、無関係なpython.exeの総数ではなく
        **自分が起こしたものが全部消えたか**で判定するための材料
        （偽Console本体と、その孫プロセスまで数える）。
        """
        from lipidmix.core.process_control import process_identity, same_process

        alive = []
        for _key, proc, _log in self._workers:
            if proc.poll() is None and process_identity(proc.pid) is not None:
                alive.append({"pid": proc.pid, "kind": "worker"})
        for _key, info in self._detached:
            if same_process(info.get("identity") or {}):
                alive.append({"pid": info.get("pid"), "kind": "detached_worker"})
        for counter in sorted(self.counters_dir.glob("*.tsv")):
            for row in counter.read_text(encoding="utf-8").splitlines():
                parts = row.split("\t")
                if len(parts) < 3:
                    continue
                for index, kind in ((1, "fake_console"), (2, "fake_console_grandchild")):
                    if parts[index].isdigit() and process_identity(int(parts[index])):
                        alive.append({"pid": int(parts[index]), "kind": kind})
        return alive

    def close(self) -> None:
        self.stop_workers()
        self._kill_detached()

    def _kill_detached(self) -> None:
        from lipidmix.core.process_control import same_process

        for _key, info in self._detached:
            identity = info.get("identity") or {}
            if not same_process(identity):
                continue
            try:
                os.kill(int(identity["pid"]), signal.SIGTERM)
            except OSError:
                pass
