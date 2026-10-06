"""データ探索ディレクトリの一元設定。

既定では本リポジトリ直下の `data/` を見るが、環境変数 `LIPIDMIX_DATA_DIR` が
設定されていればそのパスを優先する。これにより、コードを変更せずに任意の
MS-DIAL出力フォルダ（例: 生データと同じフォルダ）を解析対象にできる。

使い方:
    from metabolomix.core.data_config import get_data_dir
    DATA_DIR = get_data_dir()

セッション内で対象を切り替える例 (PowerShell):
    $env:LIPIDMIX_DATA_DIR = "C:\\Users\\yuu18\\datasets\\2_lipidome_lcms\\NEG"
"""

import os
from pathlib import Path

ENV_VAR = "LIPIDMIX_DATA_DIR"

# リポジトリ直下の data/ を既定とする
# mcp_core は data_config を import する側なので、mcp_core.BASE_DIR は参照できない
# （循環する）。独立に同じ値を計算し、一致は tests/test_package_layout.py で縛る。
DEFAULT_DATA_DIR = Path(__file__).resolve().parents[2] / "data"


def get_data_dir() -> Path:
    """解析対象データを探すディレクトリを返す。

    環境変数 LIPIDMIX_DATA_DIR があればそれを、無ければ <project>/data を使う。
    """
    override = os.environ.get(ENV_VAR)
    if override:
        return Path(override).expanduser()
    return DEFAULT_DATA_DIR
