"""
.dcl (MSDecResult) パーサ
========================================================================
MS-DIAL が出力する .dcl ファイルは、デコンボリューション済み MS/MS スペクトルを
格納するバイナリ。他の MS-DIAL 出力(.arf/.arf2/.pai2)と異なり msgpack/lz4 ではなく
独自バイナリ形式（BinaryWriter/BitConverter, リトルエンディアン）。

バイト構造（公式 MsdecResultsReader.ReadMSDecResultVer1 準拠）:
  ヘッダ(11B): "DC"(2) + version int32(4) + annotationフラグ bool(1) + 結果数 int32(4)
  seekpointer配列: 結果数 × int64
  各 MSDecResult (seekpointから):
    Scan(60B)   : SeekPoint i64, ScanID i32, RawSpectrumID i32, PrecursorMz f64,
                  IonMode i32, RT f64, RI f64, Drift f64, Mz f64
    Quant(40B)  : ModelPeakMz, ModelPeakHeight, ModelPeakArea, IntegratedHeight, IntegratedArea (f64×5)
    Scoring(20B): AmplitudeScore, ModelPeakPurity, ModelPeakQuality, SignalNoiseRatio, EstimatedNoise (f32×5)
    Counts(12B) : spectraNumber i32, datapointNumber i32, modelMassNumber i32
    Spectrum    : spectraNumber × (Mass f64, Intensity f64, PeakQuality i32) = 20B/peak
    （以降の Chromatogram / ModelMasses / Annotation ブロックは本パーサでは未使用）

各結果は seekpointer で独立参照できるため、Spectrum までを読めば MS/MS は得られる。
"""

import argparse
import glob
import os
import struct
import sys
from pathlib import Path

from metabolomix.core.data_config import get_data_dir

DCL_MAGIC = b"DC"


def deserialize_dcl(file_path: str, include_spectrum: bool = True,
                    top_n_peaks: int | None = None) -> list[dict]:
    """.dcl を解析して MSDecResult のリストを返す。

    引数:
      include_spectrum: False なら MS/MS ピーク列を読まず軽量に（メタのみ）。
      top_n_peaks: 各スペクトルを強度上位 N 本に間引く（LLMトークン節約用）。None で全件。
    """
    with open(file_path, "rb") as f:
        raw = f.read()

    if raw[:2] != DCL_MAGIC:
        raise ValueError(f"{file_path} は .dcl 形式ではありません（先頭が 'DC' でない）。")

    version = struct.unpack_from("<i", raw, 2)[0]
    has_annotation = bool(raw[6])
    count = struct.unpack_from("<i", raw, 7)[0]
    seekpoints = struct.unpack_from(f"<{count}q", raw, 11)

    results = []
    for idx, sp in enumerate(seekpoints):
        o = sp
        # --- Scan block (60B) ---
        scan_id = struct.unpack_from("<i", raw, o + 8)[0]
        raw_spec_id = struct.unpack_from("<i", raw, o + 12)[0]
        precursor_mz = struct.unpack_from("<d", raw, o + 16)[0]
        ion_mode = struct.unpack_from("<i", raw, o + 24)[0]
        rt = struct.unpack_from("<d", raw, o + 28)[0]
        o += 60
        # --- Quant block (40B) ---
        model_mz, model_h, model_a, int_h, int_a = struct.unpack_from("<5d", raw, o)
        o += 40
        # --- Scoring block (20B) ---
        amp, purity, quality, sn, est_noise = struct.unpack_from("<5f", raw, o)
        o += 20
        # --- Counts block (12B) ---
        spectra_n, datapoint_n, model_n = struct.unpack_from("<3i", raw, o)
        o += 12
        # --- Spectrum peaks (spectra_n × 20B) ---
        spectrum = []
        if include_spectrum and spectra_n > 0:
            for _ in range(spectra_n):
                mass, inten = struct.unpack_from("<2d", raw, o)
                o += 20  # double mass + double intensity + int32 peakquality
                spectrum.append([mass, inten])
            if top_n_peaks is not None and len(spectrum) > top_n_peaks:
                spectrum = sorted(spectrum, key=lambda p: p[1], reverse=True)[:top_n_peaks]
                spectrum.sort(key=lambda p: p[0])  # m/z昇順に戻す

        results.append({
            "dcl_index": idx,          # .dcl 内の順序（= MSDecResultID, .pai2ピークと対応）
            "scan_id": scan_id,
            "raw_spec_id": raw_spec_id,
            "precursor_mz": precursor_mz,
            "ion_mode": ion_mode,
            "rt": rt,
            "model_peak_height": model_h,
            "signal_to_noise": sn,
            "estimated_noise": est_noise,
            "n_msms_peaks": spectra_n,
            "msms_spectrum": spectrum,
        })

    return results


def summarize_dcl(results: list[dict]) -> dict:
    """MSDecResult リストの要約統計を返す。"""
    if not results:
        return {"error": "データがありません。"}
    with_msms = [r for r in results if r["n_msms_peaks"] > 0]
    counts = sorted(r["n_msms_peaks"] for r in with_msms)
    mzs = [r["precursor_mz"] for r in results if r["precursor_mz"] > 0]
    rts = [r["rt"] for r in results]
    median = counts[len(counts) // 2] if counts else 0
    return {
        "total_results": len(results),
        "with_msms": len(with_msms),
        "msms_rate_pct": round(len(with_msms) / len(results) * 100, 1) if results else 0.0,
        "msms_peak_count_median": median,
        "msms_peak_count_max": max(counts) if counts else 0,
        "precursor_mz_range": (round(min(mzs), 4), round(max(mzs), 4)) if mzs else None,
        "rt_range": (round(min(rts), 2), round(max(rts), 2)) if rts else None,
    }


def get_msms_by_precursor(results: list[dict], precursor_mz: float,
                          tol: float = 0.01, rt: float | None = None,
                          rt_tol: float = 0.2) -> list[dict]:
    """指定 precursor m/z（任意で RT）に一致する MS/MS を返す。"""
    hits = []
    for r in results:
        if abs(r["precursor_mz"] - precursor_mz) > tol:
            continue
        if rt is not None and abs(r["rt"] - rt) > rt_tol:
            continue
        if r["n_msms_peaks"] > 0:
            hits.append(r)
    return hits


def find_dcl_for_pai2(pai2_path: str) -> str | None:
    """.pai2 と同じ basename の .dcl を同フォルダから探す。"""
    base = pai2_path[:-len(".pai2")] if pai2_path.endswith(".pai2") else pai2_path
    cand = base + ".dcl"
    return cand if os.path.exists(cand) else None


def attach_msms_to_features(features: list[dict], dcl_results: list[dict],
                            mz_tol: float = 0.01) -> int:
    """.pai2 の各ピーク(feature)に .dcl の MS/MS スペクトルを索引対応で付与する。

    .dcl の dcl_index は .pai2 のピーク順(= MasterPeakID)と一致するため索引で対応づける。
    安全のため precursor m/z が許容差を超える場合は付与をスキップする。
    feature に `msms_spectrum` / `n_msms_peaks` を書き込み、付与できた件数を返す。
    """
    attached = 0
    n = min(len(features), len(dcl_results))
    for i in range(n):
        feat, res = features[i], dcl_results[i]
        feat_mz = feat.get("m/z")
        if isinstance(feat_mz, (int, float)) and abs(feat_mz - res["precursor_mz"]) > mz_tol:
            continue  # 索引対応が崩れている可能性 → スキップ
        feat["msms_spectrum"] = res["msms_spectrum"]
        feat["n_msms_peaks"] = res["n_msms_peaks"]
        feat["msms_peak_count"] = res["n_msms_peaks"]
        if res["n_msms_peaks"] > 0:
            attached += 1
    return attached


def _select_dcl_file(file_path: str | None, index: int | None) -> str | None:
    """解析対象 .dcl を決定する（非対話対応）。--file > --index > 先頭[0]。"""
    if file_path:
        if os.path.exists(file_path):
            return file_path
        print(f"[WARNING] 指定ファイルが見つかりません: {file_path}", file=sys.stderr)
        return None
    candidates = sorted(glob.glob(os.path.join(str(get_data_dir()), "*.dcl")))
    if not candidates:
        return None
    if index is not None:
        if 0 <= index < len(candidates):
            return candidates[index]
        print(f"[WARNING] index {index} は範囲外です (0..{len(candidates) - 1})。", file=sys.stderr)
        return None
    print("解析可能な .dcl ファイル:")
    for i, c in enumerate(candidates):
        print(f"  [{i}] {os.path.basename(c)}")
    print("[INFO] 番号未指定のため先頭[0]を使用します。(--file/--index で指定可)", file=sys.stderr)
    return candidates[0]


def main():
    parser = argparse.ArgumentParser(description="MS-DIAL .dcl (MSDecResult/MS-MS) ファイルを解析します。")
    parser.add_argument("--file", "-f", help="解析する .dcl ファイルのパス")
    parser.add_argument("--index", "-i", type=int, default=None,
                        help="データディレクトリ内の .dcl 一覧から番号で指定（非対話）")
    parser.add_argument("--top-peaks", type=int, default=None,
                        help="各 MS/MS を強度上位 N 本に間引く（既定: 全件）")
    parser.add_argument("--show", type=int, default=10,
                        help="先頭から表示する MSDecResult 件数（既定: 10）")
    args = parser.parse_args()

    file_path = _select_dcl_file(args.file, args.index)
    if not file_path:
        raise FileNotFoundError(".dcl ファイルが見つかりません。")

    print(f"[INFO] 解析対象ファイル: {file_path}")
    results = deserialize_dcl(file_path, include_spectrum=True, top_n_peaks=args.top_peaks)
    print(f"[INFO] MSDecResult 数: {len(results)}")

    # ===== デシリアライズ結果をCLI表示 =====
    print(f"\n===== デシリアライズ結果 (.dcl): 先頭{min(args.show, len(results))}件 =====")
    for r in results[:args.show]:
        top = sorted(r["msms_spectrum"], key=lambda p: p[1], reverse=True)[:5]
        top_s = ", ".join(f"{m:.4f}({it:.0f})" for m, it in top) if top else "(MS/MSなし)"
        print(f"  [#{r['dcl_index']:>4}] precursor m/z={r['precursor_mz']:9.4f}  RT={r['rt']:6.2f}  "
              f"IonMode={r['ion_mode']}  MS/MS={r['n_msms_peaks']:>4}本  上位: {top_s}")

    # ===== 要約 =====
    s = summarize_dcl(results)
    print("\n===== .dcl 要約 =====")
    print(f"  総MSDecResult: {s['total_results']}")
    print(f"  MS/MS保有: {s['with_msms']} ({s['msms_rate_pct']}%)")
    print(f"  MS/MSピーク数: 中央値={s['msms_peak_count_median']}, 最大={s['msms_peak_count_max']}")
    print(f"  precursor m/z範囲: {s['precursor_mz_range']}")
    print(f"  RT範囲: {s['rt_range']} min")


if __name__ == "__main__":
    main()
