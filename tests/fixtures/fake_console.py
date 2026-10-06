"""試験用の偽 MS-DIAL Console（実プロセスとして起動されるスクリプト）。

`metabolomix.console.execution.supervise` は「本当に起動したプロセスが本当に
どう終わったか」を証跡にする。そこをモックで置き換えると、監視・停止・収集の
どれが壊れても緑のままになるため、このファイルは **import せずスクリプトパス
として子プロセスへ渡す**（テストからの `import` は禁止。`tests.pipeline_fixtures`
の write_mztab だけをこちらから使う）。

シナリオ（`--scenario`）:

``success``
    Task 3 の合成 mzTab-M（`-i` 直下の全 raw に対応する assay）を `-o` へ書いて 0 で終了。
``nonzero``
    中間ファイルだけを `-o` へ残して 1 で終了。**終了コードが非ゼロでも成果物は残る**
    という現実を再現するためのもので、収集を飛ばさないことの検証に使う。
``hang``
    中間ファイル（`intermediate.pai2`）を `-o` へ出してから孫プロセス（実 Python）を
    起こし、自分も待ち続ける。timeout・取消で Job Object が一族ごと停止すること、
    その時点までの成果物が収集されることを本物のプロセスで確かめる。
``invalid``
    parse は通るが構造検証に落ちる mzTab を書いて 0 で終了。
``missing_sample``
    raw 1 本分の assay を欠いた合成 mzTab を書いて 0 で終了。

`-i` / `-o` / `-m` は実 Console と同じ意味。省略時は `--job` の analysis-job.json
から解決する。`--counter` に渡したファイルへは 1 回の起動につき 1 行だけ
TSV（scenario / 自 pid / 孫 pid）を追記する。**起動回数はこの行数で数える**
（「起動したつもり」ではなく実際の起動を数えるため）。

`--no-inchikey` は `success` の mzTab を全行未識別（`database_identifier = null`）に
する。E01「InChIKey 0 件で必須 TSV を作れない」を、壊れたファイルではなく
**構造としては完全に妥当な mzTab** で再現するためのもの（scenario を増やさない
のは、受入表の scenario 列を spec/brief のままに保つため）。
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from metabolomix.console.job_manager import list_raw_inputs  # noqa: E402
from metabolomix.core.process_control import resolve_python_launch  # noqa: E402
from tests.pipeline_fixtures import write_mztab  # noqa: E402

#: `success` が書く主 mzTab のファイル名。MS-DIAL GUI の `Height_` prefix を
#: 持たせ、定量種別が実物と同じ経路（ファイル名信号）で決まるようにする。
PRIMARY_MZTAB_NAME = "Height_AlignmentResult_fake.mzTab"

#: 孫プロセスの生存時間（秒）。テストは必ず Job Object 側で殺すので、
#: 取りこぼしたときに気付けるよう十分長くしておく。
GRANDCHILD_LIFETIME = 300


def _parse_args(argv: list[str]) -> argparse.Namespace:
    """実 Console の引数列（`lcms` サブコマンドや `-p`）を読み飛ばして解釈する。"""
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--scenario", required=True)
    parser.add_argument("-i", dest="input_dir")
    parser.add_argument("-o", dest="output_dir")
    parser.add_argument("-m", dest="method_file")
    parser.add_argument("--job", dest="job_path")
    parser.add_argument("--counter", dest="counter")
    parser.add_argument("--no-inchikey", dest="no_inchikey", action="store_true")
    args, _unknown = parser.parse_known_args(argv)
    return args


def _resolve_dirs(args: argparse.Namespace) -> tuple[Path, Path]:
    """(-i, -o) を決める。省略時は analysis-job.json から解決する。"""
    input_dir = Path(args.input_dir) if args.input_dir else None
    output_dir = Path(args.output_dir) if args.output_dir else None
    if (input_dir is None or output_dir is None):
        if not args.job_path:
            raise SystemExit("fake_console: -i/-o か --job のいずれかが必要です")
        data = json.loads(Path(args.job_path).read_text(encoding="utf-8"))
        if input_dir is None:
            input_dir = Path(data["source"]["dataset_root"])
        if output_dir is None:
            output_dir = Path(data["run_dir"]) / "msdial"
    output_dir.mkdir(parents=True, exist_ok=True)
    return input_dir, output_dir


def _raw_sources(input_dir: Path) -> list[Path]:
    """`-i` 直下の合成 raw を安定順で返す。

    監視側が固定した input inventory（`write_supervision_inputs` の
    `raw_inventory`）は `job_manager.list_raw_inputs` で作られる。ここでも
    **同じ関数**を使う——別の規則で列挙すると、ms_run location と予定入力が
    「テスト側の列挙差」で食い違い、意図しない SAMPLE_MAPPING_MISSING を
    生む。ディレクトリ形式 raw（Agilent `.d` 等）が 1 測定単位として
    そのまま 1 assay になるのも、この共有によって自動的に揃う。
    """
    return [p.resolve() for p in list_raw_inputs(input_dir)]


def _write_invalid_mztab(path: Path, sources: list[Path]) -> None:
    """parse は通るが構造検証に落ちる mzTab を書く（SMF データ行が 0 件）。

    「拡張子が .mzTab である」だけでは完了にならないことを、壊れた JSON では
    なく**実物のパーサに通る壊れ方**で確かめるための入力。
    """
    lines = [
        "MTD\tmzTab-version\t2.0.0-M",
        "MTD\tmzTab-mode\tComplete",
        "MTD\tmzTab-type\tQuantification",
        "MTD\tsmall_molecule-quantification_unit\t[PRIDE, PRIDE:0000429, Abundance, ]",
    ]
    for i, src in enumerate(sources, start=1):
        lines.append(f"MTD\tms_run[{i}]-location\t{Path(src).resolve().as_uri()}")
        lines.append(f"MTD\tassay[{i}]-ms_run_ref\tms_run[{i}]")
        lines.append(f"MTD\tassay[{i}]\tS{i}")
    abundance_cols = "\t".join(f"abundance_assay[{i}]" for i in range(1, len(sources) + 1))
    lines.append(
        "SFH\tSMF_ID\tSML_ID_REFS\tdatabase_identifier\tchemical_name\tsmiles\tinchi\t"
        f"exp_mass_to_charge\tretention_time_in_seconds\t{abundance_cols}"
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _spawn_grandchild() -> int:
    """孫の実 Python を起こして pid を返す。

    stdout/stderr/stdin はすべて DEVNULL にする。親（この偽 Console）の stdout は
    監視側が開いたログファイルのハンドルであり、パイプを介在させると待ち合わせが
    絡む。孫は「Job Object が一族を停止できるか」を確かめるためだけの存在なので、
    出力経路を一切持たせない。

    venvのリダイレクタは挟まない。挟むと記録するpidがリダイレクタになり、
    「孫が消えたか」の検査が孫本体を見なくなる（GitHub #3）。
    """
    command, env = resolve_python_launch(
        [sys.executable, "-c", f"import time; time.sleep({GRANDCHILD_LIFETIME})"])
    grandchild = subprocess.Popen(
        command, env=env, stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
    )
    return int(grandchild.pid)


def _append_counter(counter: str | None, scenario: str, grandchild_pid: int | None) -> None:
    """1 回の起動につき 1 行だけ追記する（scenario / 自 pid / 孫 pid）。"""
    if not counter:
        return
    path = Path(counter)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as fh:
        fh.write(f"{scenario}\t{os.getpid()}\t{grandchild_pid if grandchild_pid else '-'}\n")


def main(argv: list[str]) -> int:
    args = _parse_args(argv)
    scenario = args.scenario
    input_dir, output_dir = _resolve_dirs(args)

    if scenario == "hang":
        # 中間ファイルを先に置いてから孫を起こす。timeout・取消のあとで
        # 「その時点までの成果物が収集されているか」を確かめられるようにする
        # ——「止めた＝何も残らない」ではないのが現実のConsoleの振る舞い。
        (output_dir / "intermediate.pai2").write_bytes(b"fake intermediate peak file")
        grandchild_pid = _spawn_grandchild()
        _append_counter(args.counter, scenario, grandchild_pid)
        print(f"fake_console: hanging with grandchild {grandchild_pid}", flush=True)
        time.sleep(GRANDCHILD_LIFETIME)
        return 0

    _append_counter(args.counter, scenario, None)

    sources = _raw_sources(input_dir)
    if scenario == "success":
        write_mztab(output_dir / PRIMARY_MZTAB_NAME, sources,
                    with_inchikey=not args.no_inchikey)
        return 0
    if scenario == "missing_sample":
        if len(sources) < 2:
            raise SystemExit("fake_console: missing_sample には raw が 2 本以上必要です")
        write_mztab(output_dir / PRIMARY_MZTAB_NAME, sources[:-1])
        return 0
    if scenario == "invalid":
        _write_invalid_mztab(output_dir / PRIMARY_MZTAB_NAME, sources)
        return 0
    if scenario == "nonzero":
        (output_dir / "intermediate.pai2").write_bytes(b"fake intermediate peak file")
        print("fake_console: failing after writing an intermediate file", flush=True)
        return 1
    raise SystemExit(f"fake_console: 未知の scenario: {scenario}")


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
