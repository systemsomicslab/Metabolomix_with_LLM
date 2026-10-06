from typing import BinaryIO
from enum import Enum
import io
import msgpack
import lz4.block
import numpy as np


class IonMode(Enum):
    Positive = 0
    Negative = 1
    Both = 2


def deserialize_lz4_packed_msgpack(data: bytes) -> list | dict:
    """
    Deserialize LZ4 compressed MsgPack data.

    Args:
        data (bytes): The MessagePack formatted data which content is LZ4 compressed.

    Returns:
        list | dict: The deserialized data.
    """
    _header, compressed_data = msgpack.unpackb(data, raw=False)

    stream = io.BytesIO(compressed_data)
    unpacker = msgpack.Unpacker(stream, raw=False)
    size = next(unpacker)
    read_size = unpacker.tell()

    # Decompress the LZ4 data
    decompressed_data = lz4.block.decompress(compressed_data[read_size:], uncompressed_size=size)

    # Unpack the MsgPack data
    unpacked_data = msgpack.unpackb(decompressed_data, raw=False, strict_map_key=False)

    return unpacked_data


def _convert_to_times(data: list) -> dict:
    result = {}
    if not isinstance(data, list): return result # 安全策

    # for t in data[:4]:
    #     if t[0] == 1:
    #         if t[1][0] >= 0:
    #             result["rt"] = t[1][0]
    #     elif t[0] == 2:
    #         if t[1][0] >= 0:
    #             result["ri"] = t[1][0]
    #     elif t[0] == 3:
    #         if t[1][0] >= 0:
    #             result["m/z"] = t[1][0]
    #     elif t[0] == 4:
    #         if t[1][0] >= 0:
    #             result["dt"] = t[1][0]

    for t in data[:4]:
        # t がリストであり、かつ2番目の要素もリストであることを確認
        if isinstance(t, (list, tuple)) and len(t) >= 2:
            val_list = t[1]
            if isinstance(val_list, (list, tuple)) and len(val_list) > 0:
                val = val_list[0]
                if val >= 0:
                    if t[0] == 1: result["rt"] = val
                    elif t[0] == 2: result["ri"] = val
                    elif t[0] == 3: result["m/z"] = val
                    elif t[0] == 4: result["dt"] = val
    return result


def _convert_to_peakfeature(data: list) -> dict:
    # formula (26番目) の安全な取得
    formula_data = data[26]
    formula = formula_data[0] if isinstance(formula_data, (list, tuple)) and len(formula_data) > 0 else "Unknown"
    
    # adduct (30番目) の安全な取得
    adduct_data = data[30]
    adduct = adduct_data[2] if isinstance(adduct_data, (list, tuple)) and len(adduct_data) > 2 else "Unknown"

    # MS/MS情報 (Key18 MS2RawSpectrumID / Key19 MS2RawSpectrumID2CE / Key24 Spectrum)
    # 注: MS/MS実スペクトル本体は別の .dcl ファイルに格納されるため、.pai2 単体では
    #     Key24 は通常 空配列。取得有無(Key18/19)と衝突エネルギーは .pai2 から利用可能。
    ms2_id = data[18] if len(data) > 18 and isinstance(data[18], int) else -1
    ce_map = data[19] if len(data) > 19 and isinstance(data[19], dict) else {}
    spectrum = data[24] if len(data) > 24 and isinstance(data[24], list) else []
    collision_energies = sorted({float(v) for v in ce_map.values()}) if ce_map else []
    # SpectrumPeak は [mass, intensity, ...] 形式（本データでは空）
    msms_spectrum = [
        [p[0], p[1]] for p in spectrum
        if isinstance(p, (list, tuple)) and len(p) >= 2
    ]

    return {
        "time": _convert_to_times(data[4]),
        "time_left": _convert_to_times(data[3]),
        "time_right": _convert_to_times(data[5]),
        "peak_height": data[7],
        "peak_height_left": data[6],
        "peak_height_right": data[8],
        "peak_area": data[9],
        "peak_area_above_baseline": data[10],
        "m/z": data[43],
        "S/N": data[40][1],
        "id": data[11],
        "ion_mode": IonMode(data[22]),
        "name": data[25],
        "formula": formula,
        "ontology": data[27],
        "smiles": data[28],
        "inchikey": data[29],
        "adduct": adduct,
        "collision_cross_section": data[31],
        "comment": data[38],
        # --- MS/MS (T5で追加) ---
        "has_msms": bool(ce_map) or ms2_id >= 0,
        "ms2_raw_id": ms2_id,
        "collision_energies": collision_energies,
        "msms_peak_count": len(msms_spectrum),
        "msms_spectrum": msms_spectrum,
    }




def deserialize(file: BinaryIO):
    """
    Deserialize a PAI2 file.

    Args:
        file (file like object): PAI2 file object to deserialize.

    Returns:
        list: The deserialized data.
    """
    packed_data = file.read()

    datas = deserialize_lz4_packed_msgpack(packed_data)
    return [_convert_to_peakfeature(data) for data in datas]



def get_signal_to_noise(feature: dict):
    """PAI2データからS/N比を抽出できる場合に返す。"""
    for key in ["signal_to_noise", "sn", "S/N", "SignalToNoise", "signalToNoise"]:
        if key in feature and feature[key] is not None:
            try:
                return float(feature[key])
            except (TypeError, ValueError):
                pass
    return None


def filter_features_by_params(features: list, filter_params: dict | None = None):
    """指定された条件でピークを絞り込む。"""
    if filter_params is None:
        filter_params = {}

    min_intensity = filter_params.get("min_intensity", filter_params.get("min_height", 0))
    min_sn = filter_params.get("min_sn", 0)

    filtered = []
    for f in features:
        if min_intensity and f.get("peak_height", 0) < min_intensity:
            continue
        if min_sn:
            sn = get_signal_to_noise(f)
            if sn is None or sn < min_sn:
                continue
        filtered.append(f)

    return filtered


# --- 2. 重要な代謝物の抽出 ---

def _is_confident_annotation(name) -> bool:
    """PAI2 の Name が確定同定として使えるかの粗い判定。

    空文字・Unknown・`no MS2:`・`low score:` は「注釈あり」に数えない
    （output-format の Annotation 注意に準拠）。存在するだけで確定同定ではない。
    """
    if not isinstance(name, str):
        return False
    s = name.strip()
    if not s:
        return False
    low = s.lower()
    if low == "unknown":
        return False
    if low.startswith("no ms2") or low.startswith("low score"):
        return False
    return True


def summarize_pai2_inventory(features: list) -> dict:
    """1測定ファイル（PAI2）のピーク在庫を、PCA を介さず素直に要約する。

    PAI2 は単一サンプルのピーク一覧なので、サンプル間比較（オミクス PCA）は原理的に
    できない。ここでは注釈状況・m/z・RT・強度・S/N の分布と、強度上位ピークという
    生化学的に意味のある指標だけを返す。MS/MS は同名 .dcl（dcl_index がリスト順に対応）
    を参照する。
    """
    n = len(features)
    heights = [f.get("peak_height") for f in features if isinstance(f.get("peak_height"), (int, float))]
    mzs = [f.get("m/z") for f in features if isinstance(f.get("m/z"), (int, float))]
    rts = [f.get("time", {}).get("rt") for f in features
           if isinstance(f.get("time", {}).get("rt"), (int, float))]
    sns = [sn for sn in (get_signal_to_noise(f) for f in features) if sn is not None]
    annotated = sum(1 for f in features if _is_confident_annotation(f.get("name")))

    ion_modes: dict[str, int] = {}
    for f in features:
        mode = f.get("ion_mode")
        key = mode.name if isinstance(mode, IonMode) else str(mode)
        ion_modes[key] = ion_modes.get(key, 0) + 1

    def _stats(xs):
        if not xs:
            return None
        return {"min": round(float(min(xs)), 4),
                "median": round(float(np.median(xs)), 4),
                "max": round(float(max(xs)), 4)}

    top_by_height = []
    for f in sorted(features, key=lambda f: f.get("peak_height") or 0.0, reverse=True)[:10]:
        top_by_height.append({
            "id": f.get("id"),
            "name": f.get("name", "Unknown"),
            "m/z": round(f.get("m/z"), 4) if isinstance(f.get("m/z"), (int, float)) else None,
            "rt": round(f.get("time", {}).get("rt"), 3)
                  if isinstance(f.get("time", {}).get("rt"), (int, float)) else None,
            "height": f.get("peak_height"),
            "signal_to_noise": get_signal_to_noise(f),
        })

    return {
        "total_peaks": n,
        "annotated_peaks": annotated,
        "annotated_fraction": round(annotated / n, 3) if n else 0.0,
        "ion_mode_counts": ion_modes,
        "mz_range": (_stats(mzs) or {}).get("min") is not None and {
            "min": _stats(mzs)["min"], "max": _stats(mzs)["max"]} or None,
        "rt_range": (_stats(rts) or {}).get("min") is not None and {
            "min": _stats(rts)["min"], "max": _stats(rts)["max"]} or None,
        "height_summary": _stats(heights),
        "sn_summary": {
            "available_fraction": round(len(sns) / n, 3) if n else 0.0,
            **( _stats(sns) or {}),
        },
        "top_by_height": top_by_height,
        "note": "PAI2 は単一測定ファイルのピーク一覧です。サンプル間比較（オミクス PCA）は"
                "この単位では行えません（複数サンプルの比較は ARF/ARF2 を使用）。"
                "MS/MS は同名 .dcl（dcl_index がリスト順に対応）を参照します。",
    }


def inspect_peak_details(filtered_features: list, peak_id: str | None = None, peak_name: str | None = None):
    """
    指定したピークについて、強度・S/N・MS/MS情報などの詳細を返す。
    """
    if peak_id is None and peak_name is None:
        return {
            "status": "error",
            "message": "peak_id か peak_name のいずれかを指定してください。"
        }

    matches = []
    for feat in filtered_features:
        if peak_id is not None and str(feat.get("id")) == str(peak_id):
            matches.append(feat)
        elif peak_name is not None and isinstance(feat.get("name"), str) and peak_name.lower() in feat.get("name", "").lower():
            matches.append(feat)

    if not matches:
        return {
            "status": "not_found",
            "message": "指定されたピークがフィルタ済みデータ内に見つかりませんでした。"
        }

    results = []
    for feat in matches:
        sn = get_signal_to_noise(feat)
        has_msms = any(
            isinstance(k, str) and any(token in k.lower() for token in ["msms", "fragment", "spectrum"]) 
            for k in feat.keys()
        )
        results.append({
            "id": feat.get("id"),
            "name": feat.get("name", "Unknown"),
            "m/z": feat.get("m/z"),
            "rt": feat.get("time", {}).get("rt"),
            "height": feat.get("peak_height"),
            "area": feat.get("peak_area"),
            "signal_to_noise": sn,
            "has_msms_like_fields": has_msms,
            "formula": feat.get("formula"),
            "adduct": feat.get("adduct"),
            "comment": feat.get("comment"),
        })

    return {
        "status": "success",
        "matches": results,
        "note": "PAI2データに明示的な S/N フィールドがない場合は null になります。"
    }