import argparse
import glob
import io
import os
from enum import Enum
import sys
import re

import lz4.block
import msgpack
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')  # 非インタラクティブバックエンドを使用
import matplotlib.pyplot as plt
import pprint

from lipidmix.analysis.pca import run_pca  # 後方互換の再エクスポート。消すと ARF テストのモックが効かなくなる
from lipidmix.core.data_config import get_data_dir
from lipidmix.msdial.classes import get_sample_class_id
from lipidmix.msdial.tags import get_sample_peak_tag_info


class IonMode(Enum):
    Positive = 0
    Negative = 1
    Both = 2
    Unknown = -1

    @classmethod
    def from_int(cls, val):
        try:
            return cls(val).name
        except Exception:
            return cls.Unknown.name


def _convert_to_times(data: list) -> dict:
    result = {}
    for t in data[:4]:
        if t[0] == 1 and isinstance(t[1], (list, tuple)) and len(t[1]) > 0:
            result["rt"] = t[1][0]
        elif t[0] == 2 and isinstance(t[1], (list, tuple)) and len(t[1]) > 0:
            result["ri"] = t[1][0]
        elif t[0] == 3 and isinstance(t[1], (list, tuple)) and len(t[1]) > 0:
            result["m/z"] = t[1][0]
        elif t[0] == 4 and isinstance(t[1], (list, tuple)) and len(t[1]) > 0:
            result["dt"] = t[1][0]
    return result


def _decode(value):
    return value.decode('utf-8', errors='ignore') if isinstance(value, bytes) else value


def _is_arf_feature_container(datas) -> bool:
    if not isinstance(datas, list) or len(datas) < 1:
        return False
    first = datas[0]
    if not isinstance(first, list) or len(first) < 10:
        return False
    if not isinstance(first[0], int) or not isinstance(first[1], int):
        return False

    nested_groups = [item for item in first[2:] if isinstance(item, list) and len(item) >= 3 and all(isinstance(sub, list) for sub in item)]
    if nested_groups:
        sample_group = nested_groups[0]
        print(f"[INFO] Detected ARF format: peak groups have {len(sample_group)} elements", file=sys.stderr)

    return len(nested_groups) >= max(10, len(first[2:]) // 2)


def _convert_arf_peak_group(
    group: list,
    group_index: int | None = None,
    source_block_index: int | None = None,
    source_local_index: int | None = None,
) -> dict | None:
    if not isinstance(group, list) or len(group) < 3 or not all(isinstance(row, list) for row in group):
        return None

    representative = group[0]
    times = _convert_to_times(representative[15]) if len(representative) > 15 and isinstance(representative[15], list) else {}

    return {
        "MasterAlignmentID": group_index,
        "AlignmentID": group_index,
        "RT": times.get("rt", None),
        "MassCenter": representative[22] if len(representative) > 22 else None,
        "IonMode": IonMode.from_int(representative[23]) if len(representative) > 23 else IonMode.Unknown.name,
        "Name": _decode(representative[24]) if len(representative) > 24 else None,
        "HeightAverage": representative[18] if len(representative) > 18 else None,
        "AlignedPeakProperties": group,
        "SourceBlockIndex": source_block_index,
        "SourceLocalIndex": source_local_index,
    }

def make_loading_details(pca_result: dict, deserialized_list: list[dict], feature_names: list[str]) -> str:
    """
    PCAのloadings情報を詳細化して、書き出し用のJSON文字列を返す。

    戻り値: JSON文字列（基本部分はindent=2で整形、loadingsは各アイテムを1行で保持）
    """
    import json

    detailed_loadings = []
    for pc_idx, pc_loadings in enumerate(pca_result.get("loadings", [])):
        pc_loadings_list = []
        for feat_idx, loading_value in enumerate(pc_loadings):
            feat_name = feature_names[feat_idx]
            parts = feat_name.split("_")
            if len(parts) >= 3:
                try:
                    master_id = int(parts[1])
                except Exception:
                    master_id = None

                if master_id is not None and master_id < len(deserialized_list):
                    spot = deserialized_list[master_id]
                    annotation = spot.get("Name", "")
                    m_z = spot.get("MassCenter")
                    rt = spot.get("RT")
                else:
                    annotation = ""
                    m_z = None
                    rt = None

                pc_loadings_list.append({
                    "id": master_id,
                    "value": loading_value,
                    "annotation": annotation,
                    "m_z": m_z,
                    "rt": rt,
                })
            else:
                pc_loadings_list.append({
                    "id": None,
                    "value": loading_value,
                    "annotation": "",
                    "m_z": None,
                    "rt": None,
                })

        pc_loadings_list.sort(key=lambda x: x["value"], reverse=True)
        detailed_loadings.append(pc_loadings_list)

    # base JSON (loadings を除いたもの) を生成
    pca_copy = pca_result.copy()
    if "loadings" in pca_copy:
        pca_copy.pop("loadings")

    base_json = json.dumps(pca_copy, ensure_ascii=False, indent=2)

    # loadings 部分を独自フォーマットで構築
    loadings_json_lines = []
    loadings_json_lines.append('  "loadings": [')
    for i, pc_loadings_list in enumerate(detailed_loadings):
        loadings_json_lines.append('    [')
        for j, item in enumerate(pc_loadings_list):
            id_val = json.dumps(item["id"], ensure_ascii=False)
            val_val = json.dumps(item["value"], ensure_ascii=False)
            ann_val = json.dumps(item["annotation"], ensure_ascii=False)
            mz_val = json.dumps(item["m_z"], ensure_ascii=False)
            rt_val = json.dumps(item["rt"], ensure_ascii=False)

            item_str = f'{{ "id": {id_val}, "value": {val_val}, "m_z": {mz_val}, "rt": {rt_val}, "annotation": {ann_val} }}'
            comma = ',' if j < len(pc_loadings_list) - 1 else ''
            loadings_json_lines.append(f'      {item_str}{comma}')

        comma = ',' if i < len(detailed_loadings) - 1 else ''
        loadings_json_lines.append(f'    ]{comma}')
    loadings_json_lines.append('  ]')

    final_json = base_json[:-2] + ',\n' + '\n'.join(loadings_json_lines) + '\n}'
    return final_json


def _msgpack_stream(data: bytes):
    return msgpack.Unpacker(
        io.BytesIO(data),
        raw=False,
        # POSモードの .arf はヘッダ相当フィールドに不正なUTF-8バイトを含む
        # ことがある（NEGは含まない）。strict だと1バイトで全体が
        # UnicodeDecodeError で停止するため、既存の _decode(errors='ignore')
        # 規約に合わせて不正バイトを無視して読み進める。
        unicode_errors="ignore",
        strict_map_key=False,
        max_buffer_size=1024 * 1024 * 1024,
    )


def _decode_lz4_msgpack_payload(compressed_data: bytes) -> list:
    stream = io.BytesIO(compressed_data)
    unpacker = msgpack.Unpacker(
        stream,
        raw=False,
        unicode_errors="ignore",
        strict_map_key=False,
        max_buffer_size=1024 * 1024 * 1024,
    )
    size = next(unpacker)
    read_size = unpacker.tell()
    decompressed_data = lz4.block.decompress(
        compressed_data[read_size:],
        uncompressed_size=size,
    )

    return list(_msgpack_stream(decompressed_data))


def _payload_from_top_level_object(obj):
    if isinstance(obj, msgpack.ExtType):
        return obj.data

    if isinstance(obj, (list, tuple)) and len(obj) == 2:
        _header, payload = obj
        if isinstance(payload, msgpack.ExtType):
            return payload.data
        if isinstance(payload, (bytes, bytearray)):
            return bytes(payload)

    return None


def _iter_lz4_msgpack_blocks(data: bytes):
    for block_index, obj in enumerate(_msgpack_stream(data)):
        payload = _payload_from_top_level_object(obj)
        if payload is None:
            print(
                f"[WARNING] ARF top-level object {block_index} is not a supported LZ4 MessagePack block; skipped",
                file=sys.stderr,
            )
            continue
        yield block_index, _decode_lz4_msgpack_payload(payload)


def _is_arf_peak_group(item) -> bool:
    """item が 1スポット分の「グループ（= サンプル行の並び）」か判定する。

    グループは行(AlignmentChromPeakFeature)のリストなので、非空かつ全要素が
    リストになる。行そのものは大半がスカラー要素を持つため、この条件では
    グループと誤認されない。
    """
    return (
        isinstance(item, list)
        and len(item) > 0
        and all(isinstance(row, list) for row in item)
    )


def _is_arf_feature_container_object(item) -> bool:
    """item が特徴量コンテナ（グループを内包する）か判定する。

    ヘッダ長は ion mode で異なる（NEG は [id, id] の2要素、POS は先頭に
    プロジェクト文字列を含む99要素などになる）ため、位置を仮定せず
    『グループ要素を1つ以上含むか』で判定する。単独グループ自身は
    グループ要素を持たないためコンテナ扱いされない。
    """
    return isinstance(item, list) and any(_is_arf_peak_group(e) for e in item)


def _iter_arf_peak_groups(data: bytes):
    for block_index, items in _iter_lz4_msgpack_blocks(data):
        if not items:
            continue

        local_index = 0
        for item in items:
            if _is_arf_feature_container_object(item):
                # コンテナ: グループ要素のみを取り出す（先頭のヘッダ
                # スカラー群は自動的にスキップされる）。
                for element in item:
                    if _is_arf_peak_group(element):
                        yield block_index, local_index, element
                        local_index += 1
            elif _is_arf_peak_group(item):
                # 単独グループ（コンテナ外に置かれた末尾スポット等）。
                yield block_index, local_index, item
                local_index += 1
            # それ以外（グループを持たないメタデータ）はスキップ。


def deserialize_lz4_packed_msgpack(data: bytes) -> list:
    """LZ4圧縮されたMsgPackデータを解凍し、内部のリストを返します。"""
    all_spots = []
    for _block_index, items in _iter_lz4_msgpack_blocks(data):
        all_spots.extend(items)
    return all_spots


def extract_arf_data(
    data,
    group_index: int | None = None,
    source_block_index: int | None = None,
    source_local_index: int | None = None,
) -> dict | None:
    if isinstance(data, list) and len(data) >= 3 and all(isinstance(row, list) for row in data):
        return _convert_arf_peak_group(
            data,
            group_index=group_index,
            source_block_index=source_block_index,
            source_local_index=source_local_index,
        )

    # .arf ファイルのみを扱うため、他の形式はサポートしない
    return None


class ArfMemoryError(MemoryError):
    """`.arf` を読み切る前に RAM が尽きたことと、そこまでの規模を運ぶ。

    `str(MemoryError())` は空文字なので、素の `MemoryError` のままだと呼び出し側は
    理由なしのエラーしか返せない。400 検体規模で最初に当たる失敗がこれなので、
    どこまで読めたか（スポット数・行数）を付けて投げ直す。
    """

    def __init__(self, *, spots: int, rows: int):
        super().__init__(f"{spots} スポット / {rows} 行を読んだところで RAM が尽きた")
        self.spots = spots
        self.rows = rows


def deserialize(file_like_object) -> list[dict]:
    data = file_like_object.read()
    results = []
    rows = 0

    try:
        for group_index, (block_index, local_index, group) in enumerate(
            _iter_arf_peak_groups(data),
            start=0,
        ):
            formatted = extract_arf_data(
                group,
                group_index=group_index,
                source_block_index=block_index,
                source_local_index=local_index,
            )
            if formatted is not None:
                formatted["MasterAlignmentID"] = group_index
                formatted["AlignmentID"] = group_index
                results.append(formatted)
                rows += len(group)
    except MemoryError as exc:
        # 規模を運ぶ例外へ詰め替える。ここで確保するのは小さなオブジェクト 1 つで、
        # 失敗した確保（数 GB 単位）とは桁が違う。
        raise ArfMemoryError(spots=len(results), rows=rows) from exc
    return results


def _is_alignment_peak_row(data) -> bool:
    """`_convert_to_alignment_feature` が実経路に入る行かどうか。

    file_name だけ欲しい呼び出し（`exclusions._entry_file_name`）と、dict を組む
    本体とで判定が食い違わないよう、ここ 1 か所に置く。
    """
    return (isinstance(data, list) and len(data) >= 3
            and len(data) > 25 and isinstance(data[18], (int, float)))


def _file_name_from_row(data: list) -> str | None:
    """行が実経路のものだと分かっている前提で、先頭付近の文字列を拾う。

    呼び出し側が `_is_alignment_peak_row` を通したあとに使う（ここで再判定すると、
    spot × 注入ぶんの isinstance/len が二重に走る）。
    """
    for item in data[:10]:
        if isinstance(item, (bytes, str)):
            decoded = _decode(item)
            if len(decoded) > 0:
                return decoded
    return None


def file_name_of(data: list) -> str | None:
    """`AlignedPeakProperties` の 1 行から file_name だけを取り出す公開境界。

    `_convert_to_alignment_feature(data).get("file_name")` と必ず同じ答えを返す
    （両者が同じ導出を共有している）。除外判定・roster は file_name しか要らない
    のに、以前は 1 行ごとに feature dict を丸ごと組んでいた——spot × 注入ぶん
    呼ばれるので、393MB の `.arf` では 294 万回の無駄な dict 生成になっていた。
    """
    if not _is_alignment_peak_row(data):
        return None
    return _file_name_from_row(data)


#: `_convert_to_alignment_feature` がそのまま添字を引くだけのフィールド。
#: 実経路の行は len(data) > 25 が保証されているので、範囲判定は要らない。
_DIRECT_FIELD_INDEX = {
    "file_id": 0,            # 真の FileID は Key0
    "master_peak_id": 2,
    "peak_id": 3,
    "ms2_raw_id": 9,         # Key9
    "height": 18,
    "area": 20,
    "area_above_baseline": 21,  # Key21
    "m_z": 22,
}

#: 行に依らず定数のフィールド（同定側があとから書き換える前提の初期値）。
_CONSTANT_FIELDS = {"is_msms_matched": False, "is_matched": False}


def peak_field(data: list, name: str):
    """`_convert_to_alignment_feature(data).get(name)` と同じ値を dict 無しで返す。

    `build_pca_matrix` は file_name / is_gap_filled / 要求された props しか読まない
    のに、1 行ごとに 16 キーの dict を組んでいた（393MB の `.arf` で 294 万回。
    既定の `props=["height"]` では `rt` のための `_convert_to_times` まで毎行走る）。

    **同値性は `tests/test_arf_hot_paths.py` が全キー × 行の種類で固定している。**
    フィールドを増やすときは両方に足すこと。
    """
    if not _is_alignment_peak_row(data):
        return None
    index = _DIRECT_FIELD_INDEX.get(name)
    if index is not None:
        return data[index]
    if name == "file_name":
        return _file_name_from_row(data)
    if name == "is_gap_filled":
        # MasterPeakID(Key2)が負(-2)のサンプルはギャップフィル（未検出→補間値）
        return isinstance(data[2], int) and data[2] < 0
    if name == "is_msms":
        # MS2: Key10 MS2RawSpectrumID2CE が非空ならMS/MS取得済み
        return bool(data[10]) if isinstance(data[10], dict) else False
    if name == "rt":
        return _rt_of(data)
    if name in ("estimated_noise", "signal_to_noise"):
        # Key37 PeakShape = [EstimatedNoise, SignalToNoise, ...]（msgpack配列）
        peak_shape = data[37] if len(data) > 37 and isinstance(data[37], list) else []
        position = 0 if name == "estimated_noise" else 1
        return peak_shape[position] if len(peak_shape) > position else None
    return _CONSTANT_FIELDS.get(name)


def _rt_of(data: list):
    rt_value = None
    if len(data) > 15 and isinstance(data[15], list):
        rt_value = _convert_to_times(data[15]).get("rt")
    elif len(data) > 16 and isinstance(data[16], list):
        rt_value = _convert_to_times(data[16]).get("rt")
    return rt_value


def _convert_to_alignment_feature(data: list) -> dict:
    if not isinstance(data, list) or len(data) < 3:
        return {}

    # .arf の AlignmentChromPeakFeature row
    if _is_alignment_peak_row(data):
        rt_value = _rt_of(data)

        file_name = _file_name_from_row(data)

        # Key37 PeakShape = [EstimatedNoise, SignalToNoise, ...]（msgpack配列）
        peak_shape = data[37] if len(data) > 37 and isinstance(data[37], list) else []
        estimated_noise = peak_shape[0] if len(peak_shape) > 0 else None
        signal_to_noise = peak_shape[1] if len(peak_shape) > 1 else None

        # MasterPeakID(Key2)が負(-2)のサンプルはギャップフィル（未検出→補間値）
        master_peak_id = data[2] if len(data) > 2 else None
        is_gap_filled = isinstance(master_peak_id, int) and master_peak_id < 0

        # MS2: Key10 MS2RawSpectrumID2CE が非空ならMS/MS取得済み
        ms2_ce = data[10] if len(data) > 10 else None
        is_msms = bool(ms2_ce) if isinstance(ms2_ce, dict) else False

        return {
            "peak_id": data[3] if len(data) > 3 else None,
            "file_id": data[0] if len(data) > 0 else None,  # 修正: 真のFileIDはKey0
            "file_name": file_name, # 追加
            "master_peak_id": master_peak_id,
            "height": data[18] if len(data) > 18 else None,
            "area": data[20] if len(data) > 20 else None,
            "area_above_baseline": data[21] if len(data) > 21 else None,  # 追加(Key21)
            "m_z": data[22] if len(data) > 22 else None,
            "rt": rt_value,
            "signal_to_noise": signal_to_noise,    # 追加(Key37[1])
            "estimated_noise": estimated_noise,    # 追加(Key37[0])
            "ms2_raw_id": data[9] if len(data) > 9 else None,  # 追加(Key9)
            "is_msms": is_msms,                    # 追加: MS/MS取得有無
            "is_gap_filled": is_gap_filled,        # 追加: 補間値フラグ
            "is_msms_matched": False,
            "is_matched": False,
        }

    # .arf ファイルのみを扱うため、他の形式はサポートしない
    return {}


def readable_block_count(file_like_object) -> int:
    """`.arf` として解釈できたトップレベルブロックの数を返す（公開境界）。

    `deserialize()` は壊れたファイルでも例外を出さず空リストを返す（別形式の
    トップレベルオブジェクトは warning で読み飛ばす仕様）。呼び出し側が
    「読めなかった」と「読めたがスポットが無い」を区別できるように、
    ブロック数だけを数える経路をここに置く。0 なら `.arf` として読めていない。
    """
    return sum(1 for _ in _iter_lz4_msgpack_blocks(file_like_object.read()))


def alignment_feature_row(data: list) -> dict:
    """`AlignedPeakProperties` の1行（注入1件分）を辞書へ落とす公開境界。

    `.arf` の msgpack Key 番号を知っているのはこのモジュールだけ、という約束を
    層をまたいで保つための入口（`docs/schema/AlignmentChromPeakFeature.md` が
    Key 番号の正準）。`lipidmix/mztab/evidence.py` と
    `lipidmix/analysis/assay_evidence.py` は private 名を掴まず、ここを呼ぶ。

    読めない行（別形式・短すぎる配列）では空 dict を返す——例外にしないのは、
    `.arf` に非ピーク行が混ざるのが常態だから。呼び出し側は空を「この行からは
    証拠を取れない」として扱う。
    """
    return _convert_to_alignment_feature(data)


def extract_peak_properties(deserialized_list: list[dict]) -> pd.DataFrame:
    rows = []
    for spot in deserialized_list:
        aligned = spot.get("AlignedPeakProperties")
        if not isinstance(aligned, list):
            continue
        for sample_index, sample in enumerate(aligned):
            if not isinstance(sample, list) or len(sample) < 3:
                continue
            
            feature = _convert_to_alignment_feature(sample)
            
            # 【追加】抽出したファイル名を優先し、無い場合はSample_0などにフォールバック
            file_name = feature.get("file_name")
            sample_label = file_name if file_name else f"Sample_{sample_index}"
            sample_tag_info = get_sample_peak_tag_info(
                spot, feature.get("file_id"), file_name,
            )
            sample_class_id = get_sample_class_id(
                spot, feature.get("file_id"), file_name,
            )

            row = {
                "MasterAlignmentID": spot.get("MasterAlignmentID"),
                "AlignmentID": spot.get("AlignmentID"),
                "SpotRT": spot.get("RT"),
                "SpotMassCenter": spot.get("MassCenter"),
                "IonMode": spot.get("IonMode"),
                "CompoundName": spot.get("Name"),
                "AlignmentTags": list(spot.get("Tags") or []),
                "SampleIndex": sample_index,
                "FileName": sample_label, # CSVにもファイル名を追加
                "ClassID": sample_class_id,
                "PeakID": feature.get("peak_id"),
                "FileID": feature.get("file_id"),
                "MasterPeakID": feature.get("master_peak_id"),
                "PeakHeight": feature.get("height"),
                "PeakArea": feature.get("area"),
                "PeakAreaAboveBaseline": feature.get("area_above_baseline"),
                "PeakMZ": feature.get("m_z"),
                "PeakRT": feature.get("rt"),
                "SignalToNoise": feature.get("signal_to_noise"),
                "IsMsms": feature.get("is_msms"),
                "IsGapFilled": feature.get("is_gap_filled"),
                "PeakTags": list(sample_tag_info.get("Tags") or []),
            }
            rows.append(row)
    
    return pd.DataFrame(rows) if rows else pd.DataFrame()


def count_peak_property_rows(deserialized_list: list[dict]) -> int:
    """`extract_peak_properties` が返す行数だけを、DataFrame を組まずに数える。

    `arf_parser` は「抽出された総ピークレコード数」と「平均サンプル数/スポット」
    しか使っていないのに、そのために spot × 注入ぶんの dict を作って pandas の
    DataFrame へ積んでいた（393MB の `.arf` で 881 万行）。行の採否条件は
    `extract_peak_properties` と同一に保つ——ずれたら報告値が嘘になる。
    """
    total = 0
    for spot in deserialized_list:
        aligned = spot.get("AlignedPeakProperties")
        if not isinstance(aligned, list):
            continue
        for sample in aligned:
            if not isinstance(sample, list) or len(sample) < 3:
                continue
            total += 1
    return total


def build_pca_matrix(deserialized_list: list[dict], use_properties: list[str] = None,
                     min_detection_rate: float = 0.0) -> tuple[np.ndarray, list[str], list[str]]:
    """
    指定された複数プロパティから多変量PCA用行列を構築

    引数:
      min_detection_rate: 特徴量(スポット)の検出率による足切り。0.0〜1.0。
        各スポットを「実検出(非ギャップフィル)できたサンプルの割合」で評価し、
        この閾値未満のスポットを行列から除外する。既定 0.0 = フィルタ無効（従来動作）。
        例: 0.5 を指定すると「半数以上のサンプルで実検出されたスポット」のみ使用。
    """
    # デフォルトではHeightとAreaを使用（m_zやrtも追加可能）
    if use_properties is None:
        use_properties = ["height"]

    sample_data_dict = {}
    detection_count: dict[str, int] = {}  # col_name -> 実検出(非ギャップフィル)サンプル数

    for spot in deserialized_list:
        aligned = spot.get("AlignedPeakProperties")
        if not isinstance(aligned, list) or not aligned:
            continue

        master_id = spot.get("MasterAlignmentID")

        for sample_index, sample in enumerate(aligned):
            if not isinstance(sample, list):
                continue

            # feature dict は組まない——読むのは file_name / is_gap_filled /
            # 要求された props だけで、16 キーを作るのは spot × 注入ぶんの無駄。
            # peak_field は dict 版と同じ値を返す（同値性はテストが固定している）。

            # 【追加】行列のキーとしてファイル名を使用
            file_name = file_name_of(sample)
            sample_key = file_name if file_name else f"Sample_{sample_index}"

            if sample_key not in sample_data_dict:
                sample_data_dict[sample_key] = {}

            is_detected = not peak_field(sample, "is_gap_filled")
            # 指定された複数のプロパティを列として追加
            for prop in use_properties:
                val = peak_field(sample, prop)
                col_name = f"Spot_{master_id}_{prop}"
                # 値が None（フィールド欠落・読めない値）のセルは 0.0 にする。
                # gap-fill のセルは通常 MS-DIAL の補間値を持つので、値として残る。
                sample_data_dict[sample_key][col_name] = float(val) if val is not None else 0.0
                if is_detected:
                    detection_count[col_name] = detection_count.get(col_name, 0) + 1

    if not sample_data_dict:
        return np.empty((0, 0)), [], []

    # 行がサンプル、列が「スポット×プロパティ」のデータフレームを作成
    df = pd.DataFrame.from_dict(sample_data_dict, orient='index')

    # ある試料の行がそのスポットに無かったセル（上の None→0.0 とは別）は NaN になる。
    # それだけを「その列（特徴量）の平均値」で埋める（Mean Imputation）。
    df = df.fillna(df.mean())
    # ※全サンプルで欠損だった列は NaN のまま残るので、その場合のみ 0.0 で埋める
    df = df.fillna(0.0)

    # 任意: 検出率(非ギャップフィル)による特徴量の足切り（既定0.0=無効）
    if min_detection_rate > 0.0:
        n_samples = len(df.index)
        keep_cols = [
            c for c in df.columns
            if detection_count.get(c, 0) / n_samples >= min_detection_rate
        ]
        df = df[keep_cols]

    # 改善点3: 分散が0（全サンプルで同一の値）の列を事前に行列から除外する
    df = df.loc[:, df.var(numeric_only=True) > 0]

    return df.to_numpy(dtype=float), list(df.index), list(df.columns)


def plot_pca(pca_result: dict, sample_names: list[str], props: list[str], file_path: str, output_file: str = "pca_plot.png"):
    """
    PCAスコアプロットを画像ファイルとして保存
    """
    components = np.array(pca_result["components"])
    
    if components.shape[1] < 2:
        f"[WARNING] PCA成分が2つ未満のため、プロットできません。"
        return
    
    plt.figure(figsize=(10, 8))
    plt.scatter(components[:, 0], components[:, 1], alpha=0.7)
    
    # サンプル名をラベルとして表示
    for i, name in enumerate(sample_names):
        plt.annotate(name, (components[i, 0], components[i, 1]), fontsize=8, alpha=0.8)
    
    plt.xlabel(f"PC1 ({pca_result['explained_variance_ratio'][0]*100:.1f}%)")
    plt.ylabel(f"PC2 ({pca_result['explained_variance_ratio'][1]*100:.1f}%)")
    file_name = os.path.basename(file_path)
    plt.title(f"PCA Score Plot - {file_name} (Properties: {', '.join(props)})")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()
    f"[INFO] PCAプロットを保存しました: {output_file}"


def plot_pca_scores_by_sample(pca_result: dict, sample_names: list[str], output_file: str = "pca_sample_scores.png"):
    """
    PC1とPC2について、横軸にサンプル名、縦軸にスコア(強度)をとったプロットを作成して保存します。
    """
    components = np.array(pca_result.get("components", []))
    
    if components.size == 0 or components.shape[1] < 2:
        print("[WARNING] PCA成分が2つ未満のため、サンプル別スコアをプロットできません。", file=sys.stderr)
        return
        
    pc1_scores = components[:, 0]
    pc2_scores = components[:, 1]
    
    # 2段のサブプロットを作成（x軸を共有）
    fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    
    # -------------------------
    # PC1のプロット (上段)
    # -------------------------
    axes[0].bar(sample_names, pc1_scores, color='#4C72B0', edgecolor='black', alpha=0.8)
    axes[0].set_title(f"PC1 Scores by Sample ({pca_result['explained_variance_ratio'][0]*100:.1f}%)")
    axes[0].set_ylabel("PC1 Score")
    axes[0].axhline(0, color='black', linewidth=1.0, linestyle='-') # 基準線（0）
    axes[0].grid(axis='y', linestyle='--', alpha=0.7)
    
    # -------------------------
    # PC2のプロット (下段)
    # -------------------------
    axes[1].bar(sample_names, pc2_scores, color='#55A868', edgecolor='black', alpha=0.8)
    axes[1].set_title(f"PC2 Scores by Sample ({pca_result['explained_variance_ratio'][1]*100:.1f}%)")
    axes[1].set_ylabel("PC2 Score")
    axes[1].axhline(0, color='black', linewidth=1.0, linestyle='-') # 基準線（0）
    axes[1].grid(axis='y', linestyle='--', alpha=0.7)
    
    # X軸のサンプル名が長い場合に備えて斜めに傾ける
    plt.xticks(rotation=45, ha="right")
    
    # レイアウトを整えて保存
    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"[INFO] サンプル別PCAスコアプロットを保存しました: {output_file}", file=sys.stderr)


def plot_peak_height_distribution(peak_df: pd.DataFrame, filter_keyword: str = None, group_regex: str = None, output_file: str = "peak_height_dist.png"):
    """
    指定したキーワード（脂質クラス等）でフィルタリングし、
    横軸にサンプル（またはグループ）、縦軸にPeakHeightをとったジッタープロットを作成します。
    """
    if peak_df.empty:
        print("[WARNING] データフレームが空のためプロットできません。", file=sys.stderr)
        return

    # 1. フィルタリング処理
    if filter_keyword:
        # Warning回避のため .copy() を付与
        df_filtered = peak_df[peak_df['CompoundName'].fillna("").astype(str).str.contains(filter_keyword, case=False)].copy()
        title_text = f"Class/Name: {filter_keyword}"
    else:
        df_filtered = peak_df.copy()
        title_text = "All Annotated Lipids"

    if df_filtered.empty:
        print(f"[WARNING] キーワード '{filter_keyword}' に一致するデータがありません。", file=sys.stderr)
        return

    # 2. レプリケートのグループ化
    if group_regex:
        # 正規表現にマッチした部分（_1 など）を消去してグループ名を作る
        df_filtered['GroupLabel'] = df_filtered['FileName'].apply(lambda x: re.sub(group_regex, '', str(x)))
        title_text += " (Grouped by Replicates)"
    else:
        df_filtered['GroupLabel'] = df_filtered['FileName']

    # 元の順序を保ったままユニークなグループ名を取得
    groups = df_filtered['GroupLabel'].unique()
    
    # 3. 描画設定
    plt.figure(figsize=(max(6, len(groups) * 1.2), 6))
    colors = plt.cm.tab10(np.linspace(0, 1, 10))
    
    for i, group in enumerate(groups):
        # グループに属する全レプリケートのデータを抽出
        group_data = df_filtered[df_filtered['GroupLabel'] == group]
        
        y_vals = pd.to_numeric(group_data['PeakHeight'], errors='coerce').dropna().values
        if len(y_vals) == 0:
            continue
            
        # データ点を少し散らす（ジッター）
        x_vals = np.random.normal(i, 0.08, size=len(y_vals))
        
        # プロット（レプリケートがまとまって表示される）
        plt.scatter(x_vals, y_vals, color=colors[i % 10], alpha=0.7, edgecolors='black', linewidths=0.5, s=25)
        
        # 平均値のバーを描画
        mean_val = np.mean(y_vals)
        plt.hlines(mean_val, i - 0.3, i + 0.3, colors='black', linewidth=2.5)

    # 軸とタイトルの設定
    plt.xticks(range(len(groups)), groups, rotation=45, ha="right")
    plt.ylabel("Peak Height")
    plt.title(title_text)
    
    # スタイル調整（上・右の枠線を消す）
    ax = plt.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"[INFO] 脂質分布プロットを保存しました: {output_file} (データ件数: {len(df_filtered)}, グループ数: {len(groups)})", file=sys.stderr)


def get_pca_loading_features(pca_result: dict, deserialized_list: list[dict],
                             feature_names: list[str], top_n: int = 10,
                             n_pcs: int = 2) -> list[dict]:
    """各主成分について Loading 値の正負トップを構造化して返す（表示・整形は呼び出し側）。

    返り値: [{"pc": "PC1", "var_ratio": <%>, "positive": [item...], "negative": [item...]}, ...]
      item = {"id", "value", "annotation", "m_z", "rt"}（feature名 "Spot_<id>_<prop>" から id を取得し
      deserialized_list[id] のメタデータを付与）。
    lipidmix.arf.tools の arf_parser / arf_pca_preprocessed が共通で利用する。
    """
    loadings = pca_result.get("loadings", [])
    evr = pca_result.get("explained_variance_ratio", [])
    results = []
    for pc_idx in range(min(n_pcs, len(loadings))):
        pc_loadings = loadings[pc_idx]
        feats = []
        for feat_idx, loading_value in enumerate(pc_loadings):
            parts = feature_names[feat_idx].split("_")
            if len(parts) < 3:
                continue
            try:
                master_id = int(parts[1])
            except ValueError:
                continue
            if not (0 <= master_id < len(deserialized_list)):
                continue
            spot = deserialized_list[master_id]
            feats.append({
                "id": master_id,
                "value": loading_value,
                "annotation": spot.get("Name", ""),
                "m_z": spot.get("MassCenter"),
                "rt": spot.get("RT"),
            })
        feats.sort(key=lambda x: x["value"], reverse=True)
        results.append({
            "pc": f"PC{pc_idx + 1}",
            "var_ratio": evr[pc_idx] * 100 if pc_idx < len(evr) else 0.0,
            "positive": feats[:top_n],
            "negative": feats[-top_n:][::-1],
        })
    return results


def extract_top_loading_features(pca_result: dict, deserialized_list: list[dict], feature_names: list[str], n: int = 1):
    """
    PCAのLoadingsからPC1とPC2それぞれについて、
    正に寄与している上位n件と、負に寄与している上位n件を抽出し、対応するピーク情報を出力
    """
    loadings = np.array(pca_result["loadings"])
    
    if loadings.size == 0:
        f"[WARNING] Loadingsが空です。"
        return
    
    # 少なくともPC1とPC2（または存在する分だけ）を処理
    for pc_idx in range(min(2, len(loadings))):
        pc_loadings = loadings[pc_idx]
        pc_name = f"PC{pc_idx + 1}"
        var_ratio = pca_result["explained_variance_ratio"][pc_idx] * 100

        # スポット名とロード値のペアを実値で降順ソート
        sorted_by_value = sorted(
            zip(feature_names, pc_loadings), key=lambda x: x[1], reverse=True
        )
        positive_top_n = sorted_by_value[:n]                 # 正の寄与 上位
        negative_top_n = sorted_by_value[-n:][::-1]          # 負の寄与 上位

        print(f"\n--- {pc_name} (説明分散比 {var_ratio:.2f}%) ---")
        print(f"  [正の寄与 上位{n}件]")
        _print_feature_details(positive_top_n, deserialized_list)
        print(f"  [負の寄与 上位{n}件]")
        _print_feature_details(negative_top_n, deserialized_list)


def _print_feature_details(feature_pairs: list[tuple], deserialized_list: list[dict]):
    """各ピークの詳細情報を1行で出力する。"""
    for spot_name, loading_value in feature_pairs:
        try:
            spot_id = int(spot_name.split("_")[1])
        except (IndexError, ValueError):
            continue
        if spot_id >= len(deserialized_list):
            continue

        spot = deserialized_list[spot_id]
        name = spot.get("Name") or "Unknown"
        mz = spot.get("MassCenter")
        rt = spot.get("RT")
        mz_s = f"{mz:.4f}" if isinstance(mz, (int, float)) else "N/A"
        rt_s = f"{rt:.2f}" if isinstance(rt, (int, float)) else "N/A"
        print(f"    ID {spot_id:>4}  loading={loading_value:+.6f}  m/z={mz_s}  RT={rt_s} min  {name}")


def find_input_file(file_path: str | None = None, index: int = 0) -> str | None:
    if file_path and os.path.exists(file_path):
        return file_path

    # .arf ファイルのみを検索（探索先は data_config 経由で環境変数上書き可）
    candidates = glob.glob(os.path.join(str(get_data_dir()), "*.arf"))
    return candidates[index] if candidates else None


def _list_arf_files() -> list[str]:
    """データディレクトリ内の .arf ファイルを名前順で返す。"""
    return sorted(glob.glob(os.path.join(str(get_data_dir()), "*.arf")))


def _select_arf_file(file_path: str | None, index: int | None) -> str | None:
    """解析対象 .arf を決定する（非対話対応）。

    優先順位: --file > --index > 対話プロンプト(端末時のみ) > 先頭ファイル[0]。
    """
    if file_path:
        if os.path.exists(file_path):
            return file_path
        print(f"[WARNING] 指定ファイルが見つかりません: {file_path}", file=sys.stderr)
        return None

    candidates = _list_arf_files()
    if not candidates:
        return None

    if index is not None:
        if 0 <= index < len(candidates):
            return candidates[index]
        print(f"[WARNING] index {index} は範囲外です (0..{len(candidates) - 1})。", file=sys.stderr)
        return None

    # 番号未指定: 一覧を提示
    print("解析可能な .arf ファイル:")
    for i, c in enumerate(candidates):
        print(f"  [{i}] {os.path.basename(c)}")

    # 対話端末のときだけプロンプト。非対話(パイプ/リダイレクト)時は先頭を既定採用
    if sys.stdin and sys.stdin.isatty():
        try:
            choice = int(input("何番目のファイルを解析しますか？ "))
        except (ValueError, EOFError):
            print("[INFO] 無効な入力のため [0] を使用します。", file=sys.stderr)
            choice = 0
        if not (0 <= choice < len(candidates)):
            print("[INFO] 範囲外のため [0] を使用します。", file=sys.stderr)
            choice = 0
        return candidates[choice]

    print("[INFO] 非対話モードのため先頭ファイル [0] を使用します。(--file/--index で明示指定可)",
          file=sys.stderr)
    return candidates[0]

def summarize_arf_data(deserialized_list):
    """
    ARFデータの統計的要約を生成
    """
    if not deserialized_list:
        
        return "データがありません。"

    df = pd.DataFrame(deserialized_list)
    
    summary = {
        "total_peaks": len(df),
        "rt_range": (df['RT'].min(), df['RT'].max()) if 'RT' in df.columns else None,
        "mass_range": (df['MassCenter'].min(), df['MassCenter'].max()) if 'MassCenter' in df.columns else None,
        "height_average_mean": df['HeightAverage'].mean() if 'HeightAverage' in df.columns else None,
        "height_average_max": df['HeightAverage'].max() if 'HeightAverage' in df.columns else None,
        "ion_modes": df['IonMode'].value_counts().to_dict() if 'IonMode' in df.columns else {},
        "named_compounds": df['Name'].notna().sum() if 'Name' in df.columns else 0
    }
    
    return summary

def main():
    parser = argparse.ArgumentParser(description="MS-DIAL .arf ファイルを解析します。")
    parser.add_argument("--file", "-f", help="解析する .arf ファイルのパス")
    parser.add_argument("--export", "-e", help="CSV に書き出すパス")
    parser.add_argument("--pca", action="store_true", help="PCA を実行する")
    parser.add_argument("--output-pca", help="PCA 結果を JSON で保存するパス")
    parser.add_argument("--output-plot", default=None, help="PCA プロットを PNG 保存するパス。未指定なら PNG を生成しない")
    parser.add_argument("--output-sample-scores", help="サンプル別PCAスコアプロットを画像ファイルとして保存するパス")
    parser.add_argument("--group-replicates", action="store_true", help="ファイル名からレプリケート番号（_1, _2など）を自動で除外してグループ化する")
    parser.add_argument("--group-regex", type=str, help="特殊なファイル名の場合に、削除したい文字列の正規表現を直接指定")
    parser.add_argument("--plot-distribution", action="store_true", help="脂質クラスごとのPeakHeight分布プロットを作成する")
    parser.add_argument("--filter-name", type=str, help="分布プロットで抽出する脂質名やクラスのキーワード (例: 'EtherPE_P', 'TG')")
    parser.add_argument("--output-dist-plot", type=str, default=None, help="分布プロットの PNG 保存先。未指定なら PNG を生成しない")
    parser.add_argument("--top-features", "-t", type=int, default=0, help="PCA Loadingから抽出する上位・下位ピーク件数")
    parser.add_argument("--props", nargs="+", default=["height"], 
                        help="PCAに使用するプロパティ (例: height area). m_zやrtは量ではないため非推奨です。")
    parser.add_argument("--components", type=int, default=None,
                        help="計算する主成分の数 (デフォルト: 計算可能な最大数)")
    parser.add_argument("--index", "-i", type=int, default=None,
                        help="data/ 内の .arf 一覧から解析するファイルを番号で指定（非対話）。--file指定時は無視")
    parser.add_argument("--log-transform", action="store_true",
                        help="[任意] PCA前に log10 変換を適用する（強度の歪みを抑え条件分離が向上しやすい）")
    parser.add_argument("--min-detection-rate", type=float, default=0.0,
                        help="[任意] 特徴量(スポット)の実検出率による足切り 0.0-1.0（既定0=無効。例 0.5）")
    args = parser.parse_args()

    file_path = _select_arf_file(args.file, args.index)
    if not file_path:
        raise FileNotFoundError(".arf ファイルが見つかりません。")

    print(f"[INFO] 解析対象ファイル: {file_path}")
    with open(file_path, "rb") as f:
        deserialized = deserialize(io.BytesIO(f.read()))

    print(f"[INFO] 読み込んだピーク数（Spots）: {len(deserialized)}")

    peak_df = extract_peak_properties(deserialized)
    print(f"[INFO] 抽出した Peak Property レコード数: {len(peak_df)}")
    
    if len(peak_df) > 0:
        avg_samples = len(peak_df) / len(deserialized)
        print(f"[INFO] 平均サンプル数/スポット: {avg_samples:.2f} (抽出が成功していればサンプル数と一致します)")

    # ===== デシリアライズ結果（抽出した per-sample ピークプロパティ）をCLI表示 =====
    if len(peak_df) > 0:
        print("\n===== デシリアライズ結果: 抽出ピークプロパティ (先頭20行) =====")
        with pd.option_context("display.max_columns", None, "display.width", 220):
            print(peak_df.head(20).to_string(index=False))
        print(f"... 全 {len(peak_df):,} レコード ({len(deserialized):,} スポット x 平均{avg_samples:.0f}サンプル)")

    if args.export:
        peak_df.to_csv(args.export, index=False)
        print(f"[INFO] CSV に出力しました: {args.export}")


    if args.plot_distribution:
        if len(peak_df) > 0:
            # 正規表現パターンの決定
            regex_pattern = None
            if args.group_regex:
                regex_pattern = args.group_regex
            elif args.group_replicates:
                # デフォルトの優秀な正規表現: 「_数字」の後に「_」か「行末」が続く場合のみ削除
                # 例: "..._0h_1_NEG" -> "..._0h_NEG" に自動変換される。(_0hは削除されない)
                regex_pattern = r'_\d+(?=_|$)'
                
            # PNG は明示要求時のみ。CLI が黙って画像ファイルを作らないようにする。
            if args.output_dist_plot:
                plot_peak_height_distribution(
                    peak_df=peak_df,
                    filter_keyword=args.filter_name,
                    group_regex=regex_pattern,  # 追加
                    output_file=args.output_dist_plot
                )
            else:
                print("[INFO] --output-dist-plot 未指定のため分布プロット PNG は生成しません")
        else:
            print("[WARNING] 抽出されたPeak Propertyレコードが存在しないため、分布プロットをスキップします。", file=sys.stderr)
    

    if args.pca:
        matrix, sample_names, feature_names = build_pca_matrix(
            deserialized, use_properties=args.props,
            min_detection_rate=args.min_detection_rate,
        )

        if matrix.size == 0:
            print("PCA 用データを構築できませんでした。")
            return

        print(f"[INFO] PCA入力行列の形状: {matrix.shape} (サンプル数 x 特徴量数)")
        print(f"[INFO] 使用プロパティ: {args.props}")
        if args.min_detection_rate > 0.0:
            print(f"[INFO] 検出率フィルタ: >= {args.min_detection_rate} (特徴量 {len(feature_names)} 件に絞り込み)")
        if args.log_transform:
            print("[INFO] log10 変換: 有効")
        try:
            pca_result = run_pca(matrix, n_components=args.components,
                                 log_transform=args.log_transform)
            print(f"[INFO] PCA 実行完了 (計算された主成分数: {len(pca_result['explained_variance_ratio'])})")

            # ===== PCA結果をCLI表示 =====
            evr = pca_result["explained_variance_ratio"]
            print("\n===== PCA 説明分散比 =====")
            for i, r in enumerate(evr[: min(5, len(evr))], 1):
                print(f"  PC{i}: {r * 100:.2f}%")
            if len(evr) > 0:
                print(f"  (PC1-PC2 累積: {sum(evr[:2]) * 100:.2f}%)")

            comps = pca_result.get("components", [])
            if comps and len(comps[0]) >= 2:
                print("\n===== サンプル別 PCA スコア (PC1, PC2) =====")
                for name, c in zip(sample_names, comps):
                    print(f"  {name:<42} PC1={c[0]:+9.3f}  PC2={c[1]:+9.3f}")

            # 寄与上位ピーク（--top-features 未指定でも上位5件を表示）
            n_show = args.top_features if args.top_features else 5
            print(f"\n===== PCA Loadings 寄与上位/下位 {n_show}件 =====")
            extract_top_loading_features(pca_result, deserialized, feature_names, n=n_show)

            # PNG は明示要求時のみ。CLI が黙って画像ファイルを作らないようにする。
            if args.output_plot:
                plot_pca(pca_result, sample_names, args.props, file_path, args.output_plot)
            else:
                print("[INFO] --output-plot 未指定のため PCA プロット PNG は生成しません")
            
            if args.output_sample_scores:
                plot_pca_scores_by_sample(pca_result, sample_names, output_file=args.output_sample_scores)

            if args.output_pca:
                # pca_result にサンプル名を保存しているため、そのまま渡す
                pca_result["samples"] = sample_names
                final_json = make_loading_details(pca_result, deserialized, feature_names)
                with open(args.output_pca, "w", encoding="utf-8") as f:
                    f.write(final_json)
                print(f"[INFO] PCA 結果を保存しました: {args.output_pca}")
                
        except ValueError as e:
            print(f"[ERROR] {e}")

    summary = summarize_arf_data(deserialized)
    pprint.pprint(summary)


if __name__ == "__main__":
    main()
