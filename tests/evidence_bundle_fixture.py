"""massbank-context と共有する根拠バンドル v3 の見本を決定的に作る。

`python -m tests.evidence_bundle_fixture <出力先>` で `tests/fixtures/evidence_bundle_v3` を作り直す。
合成の mzTab-M だけから作り、研究室のライブラリ・ローカルの絶対パスは入れない。
"""
from __future__ import annotations

import sys
import tempfile
from pathlib import Path

import numpy as np

from metabolomix.analysis.evidence_export import export_dataset_evidence
from metabolomix.mztab import loading
from tests.engine_chain import register_two_group
from tests.reanalysis_fixtures import POS_FEATURES, write_mztab

_SECOND_CANDIDATE = {"name": "Betaine", "inchikey": "KWIUHFFTVRNATP-UHFFFAOYSA-N",
                     "formula": "C5H11NO2", "adduct": "[M+H]1+", "theo": 118.0863}
# 有意な行（1）・有意でない行（2・3）に加え、SME 候補が 2 件の特徴（1）を持つ。
_FEATURES = (
    {**POS_FEATURES[0], "values": [10, 12, 11, 40, 44, 42, 15, 16],
     "sme": [*POS_FEATURES[0]["sme"], _SECOND_CANDIDATE]},
    *POS_FEATURES[1:],
)


def build_fixture_bundle(output_dir: Path) -> Path:
    """見本のバンドルを `output_dir` に書き、manifest のパスを返す。"""
    with tempfile.TemporaryDirectory() as tmp:
        ds = loading.load_dataset_state(mztab_path=write_mztab(Path(tmp) / "fixture.mzTab",
                                                               features=_FEATURES))
    # source_files は {パス: sha256}。一時ディレクトリの絶対パスを、見本の名前に置き換える。
    ds.source_files = {"fixture.mzTab": next(iter(ds.source_files.values()))}
    mask = np.ones((len(ds.feature_ids), len(ds.sample_names)), dtype=bool)
    mask[ds.feature_ids.index("1"), ds.sample_names.index("RA1")] = False
    # 付加体の無い行（SME 候補も無い特徴 3）を作る。付加体は SML の adduct_ions から読まれる。
    for row in ds.sml_rows:
        if row.get("SML_ID") == "3":
            row["adduct_ions"] = None
    ds.detected_mask = mask
    ds.feature_matrix[ds.feature_ids.index("2"), ds.sample_names.index("HC1")] = np.nan
    # 名前からの推定は `QC1` を qc にしないので、役割と群は実験情報シートで与える。
    # 注入順（injection_order）は全試料に整数で与える（samples.tsv の injection_order が埋まる）。
    ds.sample_metadata_rows = [{"sample_id": name, "injection_order": order,
                                "role": "qc" if name.startswith("QC") else "sample",
                                "group": None if name.startswith("QC") else name[:2]}
                               for order, name in enumerate(ds.sample_names, start=1)]
    # 統計エンジンの 2 群比較（HC を基準、RA を比較）。行列と比較の記録を ds に登録する
    # （`dataset_preprocess` → `dataset_compare` と同じ登録。根拠バンドルの試料ごとの値は行列から読む）。
    record = register_two_group(ds, "HC", "RA")
    export_dataset_evidence(ds, record, Path(output_dir), polarity="positive")
    return Path(output_dir) / "manifest.json"


if __name__ == "__main__":
    print(build_fixture_bundle(Path(sys.argv[1])))
