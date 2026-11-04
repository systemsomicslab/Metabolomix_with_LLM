"""massbank-context と共有する根拠バンドル v3 の見本を決定的に作る。

`python -m tests.evidence_bundle_fixture <出力先>` で `tests/fixtures/evidence_bundle_v3` を作り直す。
合成の mzTab-M だけから作り、研究室のライブラリ・ローカルの絶対パスは入れない。
"""
from __future__ import annotations

import sys
import tempfile
from pathlib import Path

import numpy as np

from metabolomix.analysis import dataset_service
from metabolomix.analysis.evidence_export import export_dataset_evidence
from metabolomix.mztab import loading
from tests.reanalysis_fixtures import POS_FEATURES, write_mztab

HC = ["HC1", "HC2", "HC3"]
RA = ["RA1", "RA2", "RA3"]
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
    ds.detected_mask = mask
    ds.feature_matrix[ds.feature_ids.index("2"), ds.sample_names.index("HC1")] = np.nan
    # 名前からの推定は `QC1` を qc にしないので、役割は実験情報シートで与える。
    ds.sample_metadata_rows = [{"sample_id": name, "role": "qc" if name.startswith("QC") else "sample"}
                               for name in ds.sample_names]
    dataset_service.preprocess_dataset(ds, {"normalize": "none"})
    result = dataset_service.compare_dataset(ds, HC, RA, group_a_label="HC", group_b_label="RA")
    export_dataset_evidence(ds, result, Path(output_dir), polarity="positive")
    return Path(output_dir) / "manifest.json"


if __name__ == "__main__":
    print(build_fixture_bundle(Path(sys.argv[1])))
