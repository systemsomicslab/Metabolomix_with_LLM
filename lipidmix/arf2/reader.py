# lipidmix/arf2/reader.py
import io
import msgpack
import lz4.block
import pandas as pd
import numpy as np
from enum import Enum
from typing import List, Optional, Dict
import os
import glob
import sys
import pprint

from lipidmix.core.data_config import get_data_dir

class IonMode(Enum):
    Positive = 0
    Negative = 1
    Both = 2
    Unknown = -1

    @classmethod
    def from_int(cls, val: int) -> str:
        try:
            return cls(val).name
        except ValueError:
            return cls.Unknown.name

def _decode(val) -> str:
    if isinstance(val, bytes):
        return val.decode('utf-8', errors='ignore')
    elif isinstance(val, str):
        return val
    return ""

def _extract_times(data: list) -> dict:
    result = {}
    if not isinstance(data, list):
        return result
        
    for t in data[:4]:
        if isinstance(t, list) and len(t) >= 2 and isinstance(t[1], (list, tuple)) and len(t[1]) > 0:
            if t[0] == 1: result["rt"] = float(t[1][0])
            elif t[0] == 2: result["ri"] = float(t[1][0])
            elif t[0] == 3: result["m_z"] = float(t[1][0])
            elif t[0] == 4: result["dt"] = float(t[1][0])
    return result

def deserialize_lz4_packed_msgpack(data: bytes) -> list:
    _header, compressed_data = msgpack.unpackb(data, raw=False)
    stream = io.BytesIO(compressed_data)
    unpacker = msgpack.Unpacker(stream, raw=False)
    size = next(unpacker)
    read_size = unpacker.tell()
    decompressed_data = lz4.block.decompress(compressed_data[read_size:], uncompressed_size=size)

    all_spots_raw = []
    inner_unpacker = msgpack.Unpacker(io.BytesIO(decompressed_data), raw=False, strict_map_key=False)
    for item in inner_unpacker:
        all_spots_raw.append(item)
    
    return all_spots_raw

def _to_float(val, default: float = 0.0) -> float:
    """数値以外/Noneを安全に float へ変換する。"""
    try:
        return float(val)
    except (TypeError, ValueError):
        return default


def _nested_at(val, idx: int, default=None):
    """list/tuple の idx 番目を安全に取り出す（MsgPackの入れ子オブジェクト用）。"""
    if isinstance(val, (list, tuple)) and len(val) > idx:
        return val[idx]
    return default


def extract_arf2_data(data: list) -> Optional[dict]:
    """1つのスポットデータ(.arf2形式 = AlignmentSpotProperty)からカタログ情報を抽出する。

    Key番号は docs/schema/AlignmentSpotProperty.md の Key 番号表に対応。
    既存6項目に加え、同定情報・品質メトリクスを抽出する（後方互換のため既存キーは維持）。
    """
    if not isinstance(data, list) or len(data) < 10:
        return None

    def get(index, default=None):
        return data[index] if len(data) > index else default

    times = _extract_times(get(4)) if isinstance(get(4), list) else {}

    # Key13 Formula = ['C5H10O2', 質量, ...] の先頭が組成式文字列
    formula = _decode(_nested_at(get(13), 0, "")) or ""
    # Key54 AdductType = [質量差, charge, '[M-H]-', ...] の index2 が付加体文字列
    adduct = _decode(_nested_at(get(54), 2, "")) or ""

    return {
        # --- 既存(後方互換: summarize_arf2_data / server.py が依存) ---
        "MasterAlignmentID": int(_to_float(get(0))),
        "AlignmentID": int(_to_float(get(1))),
        "RT": float(times.get("rt", 0.0)),
        "MassCenter": _to_float(get(5)),
        "IonMode": IonMode.from_int(get(11)),
        "Name": _decode(get(12)) or "Unknown",
        "HeightAverage": _to_float(get(31)),
        # --- 同定 / 化学情報 (Key13-16, 54) ---
        "Formula": formula,
        "Ontology": _decode(get(14)) or "",       # 脂質クラス (FA, PC, TG...)
        "SMILES": _decode(get(15)) or "",
        "InChIKey": _decode(get(16)) or "",
        "AdductType": adduct,                      # [M-H]- など
        # --- 強度 / 品質メトリクス (Key32-37, 43-44, 49, 51) ---
        "HeightMin": _to_float(get(32)),
        "HeightMax": _to_float(get(33)),
        "PeakWidthAverage": _to_float(get(34)),
        "SignalToNoiseAve": _to_float(get(35)),
        "SignalToNoiseMax": _to_float(get(36)),
        "SignalToNoiseMin": _to_float(get(37)),
        "MassMin": _to_float(get(43)),
        "MassMax": _to_float(get(44)),
        "FillPercentage": _to_float(get(49)),      # 検出されたサンプルの割合
        "MonoIsotopicPercentage": _to_float(get(51)),
    }

def iter_raw_spots(datas: list):
    """`deserialize_lz4_packed_msgpack` の戻りから AlignmentSpotProperty の生配列を順に返す。

    `.arf2` は「スポットを直接並べる」形と「リストの入れ子に包む」形の両方がある
    （`deserialize` と同じ規則）。Key 番号を読む層（`match_results` など）は
    ここから生配列を受け取り、展開規則を重複して持たない。
    """
    for d in datas:
        if isinstance(d, list) and len(d) < 14:
            for item in d:
                if isinstance(item, list) and len(item) > 0:
                    for spot_raw in item:
                        if isinstance(spot_raw, list) and len(spot_raw) >= 10:
                            yield spot_raw
            continue
        if isinstance(d, list) and len(d) >= 10:
            yield d


def load_raw_spots(file_path) -> list:
    """`.arf2` を読み、AlignmentSpotProperty の生配列の一覧を返す。"""
    with open(os.fspath(file_path), "rb") as handle:
        datas = deserialize_lz4_packed_msgpack(handle.read())
    return list(iter_raw_spots(datas))


def deserialize(file_like_object) -> List[dict]:
    """バイナリストリームから .arf2 データをパースして辞書のリストを返す"""
    datas = deserialize_lz4_packed_msgpack(file_like_object.read())
    results = []
    for spot_raw in iter_raw_spots(datas):
        formatted = extract_arf2_data(spot_raw)
        if formatted:
            results.append(formatted)
    return results

def summarize_arf2_data(deserialized_list: List[dict]) -> dict:
    """ARF2データの全体像（メタデータ）の統計的要約を生成"""
    if not deserialized_list:
        return {"error": "データがありません。"}

    df = pd.DataFrame(deserialized_list)
    
    summary = {
        "total_spots": len(df)
    }
    
    if 'HeightAverage' in df.columns:
        heights = df[df['HeightAverage'] > 0]['HeightAverage']
        summary["height_average_median"] = float(heights.median()) if not heights.empty else 0.0
        summary["height_average_max"] = float(heights.max()) if not heights.empty else 0.0
        
    if 'RT' in df.columns:
        summary["rt_range"] = (float(df['RT'].min()), float(df['RT'].max()))
        
    if 'MassCenter' in df.columns:
        summary["mass_range"] = (float(df['MassCenter'].min()), float(df['MassCenter'].max()))
        
    if 'IonMode' in df.columns:
        # IonModeの文字列表現をカウント
        summary["ion_modes"] = df['IonMode'].apply(lambda x: x.name if isinstance(x, IonMode) else str(x)).value_counts().to_dict()
        
    if 'Name' in df.columns:
        # アノテーション済み（Unknown または 空白 でない）の数をカウント
        summary["annotated_count"] = int(df['Name'].apply(lambda x: 1 if x and str(x).lower() != "unknown" and str(x).strip() != "" else 0).sum())
        summary["annotation_rate"] = float(summary["annotated_count"] / len(df) * 100) if len(df) > 0 else 0.0

    if 'Ontology' in df.columns:
        # 脂質クラス(オントロジー)分布の上位
        onto = df['Ontology'].apply(lambda x: str(x).strip()).replace("", np.nan).dropna()
        if not onto.empty:
            summary["ontology_top"] = onto.value_counts().head(10).to_dict()

    if 'SignalToNoiseAve' in df.columns:
        sn = df[df['SignalToNoiseAve'] > 0]['SignalToNoiseAve']
        if not sn.empty:
            summary["sn_median"] = float(sn.median())

    return summary

def generate_text_summary(deserialized_list: List[dict]) -> str:
    """LLMが読みやすい自然言語での要約を生成"""
    summary = summarize_arf2_data(deserialized_list)
    if "error" in summary: return summary["error"]
    
    text = (
        f"### ARF2 カタログ要約\n"
        f"- **総スポット（ピーク）数**: {summary.get('total_spots', 0):,} 個\n"
    )
    
    annotated = summary.get('annotated_count', 0)
    rate = summary.get('annotation_rate', 0.0)
    text += f"- **アノテーション済み**: {annotated:,} 個 ({rate:.1f}%)\n"
    
    if summary.get('rt_range'):
        text += f"- **保持時間(RT)範囲**: {summary['rt_range'][0]:.2f} 〜 {summary['rt_range'][1]:.2f} min\n"
    if summary.get('mass_range'):
        text += f"- **質量(m/z)範囲**: {summary['mass_range'][0]:.4f} 〜 {summary['mass_range'][1]:.4f}\n"
    if summary.get('height_average_median'):
        text += f"- **平均強度の代表値(中央値)**: {summary['height_average_median']:,.1f}\n"
    if summary.get('ion_modes'):
        modes = ", ".join([f"{k}: {v:,}" for k, v in summary['ion_modes'].items()])
        text += f"- **検出イオンモード**: {modes}\n"
    if summary.get('sn_median'):
        text += f"- **S/N(平均)の代表値(中央値)**: {summary['sn_median']:.1f}\n"
    if summary.get('ontology_top'):
        onto = ", ".join([f"{k}: {v:,}" for k, v in summary['ontology_top'].items()])
        text += f"- **脂質クラス上位(Ontology)**: {onto}\n"

    return text

# 同一ファイルの再デシリアライズを避けるキャッシュ。.arf2 は実データで 714〜19,388
# スポットあり、差次的解析の名前補完・EIC の同定照合・同定注釈の3経路がそれぞれ
# 独立に開いて読み直していた。ファイルが MS-DIAL で作り直されたら読み直すよう、
# パスだけでなく mtime とサイズもキーに含める。保持は1件（大きいので溜めない）。
_CATALOG_CACHE: dict[tuple, list] = {}


def load_catalog(file_path) -> list:
    """`.arf2` のスポット一覧を返す（同一ファイルなら再パースしない）。

    返すリストは**共有**されるため、呼び出し側は変更しないこと（読み取り専用）。
    """
    path = os.fspath(file_path)
    stat = os.stat(path)
    key = (os.path.abspath(path), stat.st_mtime_ns, stat.st_size)
    cached = _CATALOG_CACHE.get(key)
    if cached is None:
        with open(path, "rb") as handle:
            cached = deserialize(io.BytesIO(handle.read()))
        _CATALOG_CACHE.clear()
        _CATALOG_CACHE[key] = cached
    return cached


def format_spots_as_table(deserialized_list: List[dict], delimiter: str = "\t",
                          columns: Optional[List[str]] = None) -> str:
    """スポット辞書のリストを CSV/TSV テーブル文字列に変換する（LLMへ渡す軽量符号化）。

    整形JSON(indent/キー反復)を避け、列名を1回だけ出すことでトークン量を大幅削減する。
    delimiter="\\t" でTSV、"," でCSV。columns省略時は全キーを先頭行のキー順で出力。
    """
    if not deserialized_list:
        return ""

    if columns is None:
        columns = list(deserialized_list[0].keys())

    def cell(v) -> str:
        if v is None:
            return ""
        if isinstance(v, float):
            # 不要な桁を抑える（m/zは小数4桁、その他は適度に丸め）
            return f"{v:.4f}".rstrip("0").rstrip(".") if v else "0"
        s = str(v)
        # 区切り文字・改行をエスケープ（CSV最小限対応）
        if delimiter in s or "\n" in s or '"' in s:
            s = '"' + s.replace('"', '""') + '"'
        return s

    lines = [delimiter.join(columns)]
    for row in deserialized_list:
        lines.append(delimiter.join(cell(row.get(c)) for c in columns))
    return "\n".join(lines)


def find_input_file(file_path: str | None = None, index: int = 0) -> str | None:
    if file_path and os.path.exists(file_path):
        return file_path

    # .arf2 ファイルのみを検索（探索先は data_config 経由で環境変数上書き可）
    candidates = glob.glob(os.path.join(str(get_data_dir()), "*.arf2"))
    return candidates[index] if candidates else None

if __name__ == "__main__":
    file_path = find_input_file()

    if not file_path:
        print(".arf2 ファイルが見つかりません。")
        sys.exit(1)

    print(f"Loading {file_path}...")
    with open(file_path, 'rb') as f:
        deserialized_data = deserialize(f)

    if deserialized_data:
        # ===== デシリアライズ結果をCLI表示 =====
        df = pd.DataFrame(deserialized_data)
        print(f"\n===== デシリアライズ結果 (.arf2): {len(deserialized_data):,} スポット x {len(df.columns)} 項目 =====")
        print(f"抽出項目: {list(df.columns)}")
        with pd.option_context("display.max_columns", None, "display.width", 240):
            print(df.head(20).to_string(index=False))
        print(f"... 全 {len(deserialized_data):,} スポット")

        # .arf2 はサンプル別強度を持たないため PCA は実行不可（メタデータ要約のみ）
        print("\n" + generate_text_summary(deserialized_data))
    else:
        print("デシリアライズ結果が空です。")