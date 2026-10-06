import io
import msgpack
import os
import glob
import struct


def parse_eic_aef_css1(file_path, include_chromatogram=False):
    results = []

    with open(file_path, 'rb') as f:
        # 1. ヘッダー (10 bytes)
        magic = f.read(4)
        if magic != b'CSS1':
            raise ValueError(f"Not a CSS1 file. Magic: {magic}")
        f.seek(10)  # 10バイト目までスキップ

        # 2. スポット数 (4 bytes, int)
        num_spots = struct.unpack('<i', f.read(4))[0]

        # 3. シークポインタの読み込み (8 bytes * num_spots, long long)
        pointers = [struct.unpack('<q', f.read(8))[0] for _ in range(num_spots)]

        # 4. 各スポットデータの解析
        for i, ptr in enumerate(pointers):
            f.seek(ptr)

            # Spot Chunk 読み込み
            # RT(f), RI(f), Mass(f), Drift(f) -> 4 floats (16 bytes)
            # MainType(b) -> 1 byte
            # NumSamples(i) -> 4 bytes
            rt, ri, mass, drift = struct.unpack('<ffff', f.read(16))
            main_type = struct.unpack('<b', f.read(1))[0]
            num_samples = struct.unpack('<i', f.read(4))[0]

            samples_data = []
            for _ in range(num_samples):
                # Peak Info Chunk
                # FileID(i), NumPeaks(i), Top(f), Left(f), Right(f) -> 20 bytes
                file_id, num_peaks, top, left, right = struct.unpack('<iifff', f.read(20))

                total_intensity = 0.0
                max_intensity = 0.0
                chrom_points = []
                for _ in range(num_peaks):
                    # Horiz(f), Intensity(f) -> 8 bytes
                    h_val, intensity = struct.unpack('<ff', f.read(8))
                    total_intensity += intensity
                    max_intensity = max(max_intensity, intensity)
                    if include_chromatogram:
                        chrom_points.append((h_val, intensity))

                mean_intensity = total_intensity / num_peaks if num_peaks else 0.0
                sample_info = {
                    "file_id": file_id,
                    "peak_top": top,
                    "num_peaks": num_peaks,
                    "mean_intensity": mean_intensity,
                    "max_intensity": max_intensity,
                }
                if include_chromatogram:
                    sample_info["chromatogram"] = chrom_points

                samples_data.append(sample_info)

            results.append({
                "spot_id": i,
                "rt": rt,
                "ri": ri,
                "mz": mass,
                "drift": drift,
                "main_type": main_type,
                "num_samples": num_samples,
                "samples": samples_data,
            })

    return results


def _read_exact(stream, size, context):
    data = stream.read(size)
    if len(data) != size:
        raise ValueError(
            f"Unexpected end of EIC/AEF while reading {context}: "
            f"expected {size} bytes, got {len(data)}"
        )
    return data


def read_eic_spots_css1(
    file_path,
    spot_ids,
    file_ids=None,
    *,
    max_traces=12,
    max_total_points=200_000,
    strict=False,
):
    """Read selected chromatogram traces for several CSS1 alignment spots.

    The CSS1 pointer table permits direct access, so the file is opened once and
    every requested spot is reached with ``seek``. Unselected sample point arrays
    are skipped so plotting a few traces per spot does not expand every
    chromatogram in the file into memory.

    ``strict=True`` reproduces the single-spot contract: an out-of-range
    ``spot_id`` and a missing requested FileID both raise. ``strict=False`` omits
    out-of-range spots from the result and returns spots whose requested FileID
    is absent with an empty ``samples`` list, so callers can tell the two cases
    apart.
    """
    requested_spots = list(spot_ids)
    if any(
        isinstance(value, bool) or not isinstance(value, int) or value < 0
        for value in requested_spots
    ):
        raise ValueError("spot_ids must contain non-negative integers only")
    if max_traces < 1:
        raise ValueError("max_traces must be at least 1")
    if max_total_points < 1:
        raise ValueError("max_total_points must be at least 1")

    requested = None
    requested_set = None
    if file_ids is not None:
        requested = list(file_ids)
        if not requested:
            raise ValueError("file_ids must contain at least one FileID")
        if any(isinstance(value, bool) or not isinstance(value, int) for value in requested):
            raise ValueError("file_ids must contain integers only")
        if len(set(requested)) != len(requested):
            raise ValueError("file_ids must not contain duplicates")
        if len(requested) > max_traces:
            raise ValueError(
                f"At most {max_traces} EIC traces can be plotted at once; "
                f"requested {len(requested)}"
            )
        requested_set = set(requested)

    spots = []
    total_points = 0
    with open(file_path, "rb") as stream:
        file_size = os.fstat(stream.fileno()).st_size
        magic = _read_exact(stream, 4, "CSS1 magic")
        if magic != b"CSS1":
            raise ValueError(f"Not a CSS1 file. Magic: {magic}")
        stream.seek(10)
        num_spots = struct.unpack("<i", _read_exact(stream, 4, "spot count"))[0]
        if num_spots < 0:
            raise ValueError(f"Invalid negative EIC spot count: {num_spots}")
        pointer_table_end = 14 + 8 * num_spots

        for spot_id in sorted(set(requested_spots)):
            if spot_id >= num_spots:
                if strict:
                    raise ValueError(
                        f"spot_id {spot_id} is out of range for {num_spots} EIC spots"
                    )
                continue

            stream.seek(14 + 8 * spot_id)
            pointer = struct.unpack("<q", _read_exact(stream, 8, "spot pointer"))[0]
            if pointer < pointer_table_end or pointer >= file_size:
                raise ValueError(
                    f"Invalid EIC spot pointer for spot_id {spot_id}: {pointer}"
                )
            stream.seek(pointer)

            spot_header = _read_exact(stream, 21, f"spot {spot_id} header")
            rt, ri, mass, drift, main_type, num_samples = struct.unpack(
                "<ffffbi", spot_header
            )
            if num_samples < 0:
                raise ValueError(
                    f"Invalid negative sample count for spot_id {spot_id}: {num_samples}"
                )
            if requested_set is None and num_samples > max_traces:
                raise ValueError(
                    f"spot_id {spot_id} contains {num_samples} samples. "
                    f"Specify file_ids (maximum {max_traces}) instead of silently truncating."
                )

            samples = []
            selected_points = 0
            found_ids = set()
            for sample_index in range(num_samples):
                sample_header = _read_exact(
                    stream, 20, f"spot {spot_id} sample {sample_index} header"
                )
                file_id, num_points, peak_top, peak_left, peak_right = struct.unpack(
                    "<iifff", sample_header
                )
                if num_points < 0:
                    raise ValueError(
                        f"Invalid negative chromatogram point count for FileID {file_id}: "
                        f"{num_points}"
                    )
                point_bytes = 8 * num_points
                if stream.tell() + point_bytes > file_size:
                    raise ValueError(
                        f"Chromatogram payload exceeds file size for spot_id {spot_id}, "
                        f"FileID {file_id}"
                    )

                if requested_set is not None and file_id not in requested_set:
                    stream.seek(point_bytes, os.SEEK_CUR)
                    continue

                if total_points + num_points > max_total_points:
                    raise ValueError(
                        f"Selected EIC traces exceed the {max_total_points} point safety "
                        f"limit; narrow the compound query (names / ontologies) or "
                        f"select fewer samples (file_ids)"
                    )
                raw_points = _read_exact(
                    stream, point_bytes, f"spot {spot_id} FileID {file_id} chromatogram"
                )
                chromatogram = [
                    [float(x), float(intensity)]
                    for x, intensity in struct.iter_unpack("<ff", raw_points)
                ]
                intensities = [point[1] for point in chromatogram]
                total_points += num_points
                selected_points += num_points
                found_ids.add(file_id)
                samples.append({
                    "file_id": file_id,
                    "peak_left": float(peak_left),
                    "peak_top": float(peak_top),
                    "peak_right": float(peak_right),
                    "num_points": num_points,
                    "mean_intensity": (
                        float(sum(intensities) / len(intensities)) if intensities else 0.0
                    ),
                    "max_intensity": float(max(intensities)) if intensities else 0.0,
                    "chromatogram": chromatogram,
                })

            if strict and requested_set is not None:
                missing = [file_id for file_id in requested if file_id not in found_ids]
                if missing:
                    raise ValueError(
                        f"Requested FileID values were not found in spot_id {spot_id}: "
                        f"{missing}"
                    )

            spots.append({
                "spot_id": spot_id,
                "rt": float(rt),
                "ri": float(ri),
                "mz": float(mass),
                "drift": float(drift),
                "main_type": int(main_type),
                "num_samples": num_samples,
                "selected_samples": len(samples),
                "selected_points": selected_points,
                "samples": samples,
            })

    return spots


def read_eic_spot_css1(
    file_path,
    spot_id,
    file_ids=None,
    *,
    max_traces=12,
    max_total_points=200_000,
):
    """Read selected chromatogram traces for one CSS1 alignment spot.

    Thin wrapper over :func:`read_eic_spots_css1` with the strict single-spot
    contract: an out-of-range ``spot_id`` or a missing requested FileID raises.
    """
    if isinstance(spot_id, bool) or not isinstance(spot_id, int) or spot_id < 0:
        raise ValueError("spot_id must be a non-negative integer")
    spots = read_eic_spots_css1(
        file_path,
        [spot_id],
        file_ids,
        max_traces=max_traces,
        max_total_points=max_total_points,
        strict=True,
    )
    return spots[0]


def summarize_eic_data(results):
    if not results:
        return {
            "total_spots": 0,
            "rt_range": None,
            "mz_range": None,
            "total_samples": 0,
            "total_peaks": 0,
            "peak_top_mean": 0.0,
            "peak_top_max": 0.0,
            "unique_file_ids": [],
        }

    rt_values = [spot["rt"] for spot in results]
    mz_values = [spot["mz"] for spot in results]
    total_samples = sum(spot["num_samples"] for spot in results)
    sample_entries = [sample for spot in results for sample in spot["samples"]]
    total_peaks = sum(sample["num_peaks"] for sample in sample_entries)
    peak_tops = [sample["peak_top"] for sample in sample_entries if sample.get("peak_top") is not None]
    unique_file_ids = sorted({sample["file_id"] for sample in sample_entries})

    return {
        "total_spots": len(results),
        "rt_range": (min(rt_values), max(rt_values)) if rt_values else None,
        "mz_range": (min(mz_values), max(mz_values)) if mz_values else None,
        "total_samples": total_samples,
        "total_peaks": total_peaks,
        "peak_top_mean": sum(peak_tops) / len(peak_tops) if peak_tops else 0.0,
        "peak_top_max": max(peak_tops) if peak_tops else 0.0,
        "unique_file_ids": unique_file_ids,
    }


def search_eic_by_mz_range(results, min_mz, max_mz):
    return [spot for spot in results if min_mz <= spot.get("mz", 0.0) <= max_mz]


def search_eic_by_rt_range(results, min_rt, max_rt):
    return [spot for spot in results if min_rt <= spot.get("rt", 0.0) <= max_rt]


def top_eic_spots_by_max_intensity(results, top_n=20):
    """EIC スポットを「サンプル間で最大のクロマトグラム最大強度」で降順に並べる。

    旧 ``top_eic_spots_by_peak_top`` は ``peak_top``（ピーク頂点の横軸=RT座標。強度ではない）
    で並べており、実質「遅い RT のスポット一覧」で強度上位ではなかった。強度で上位抽出する
    には各 sample の ``max_intensity`` を使う必要がある。
    """
    spots_with_max = []
    for spot in results:
        sample_max = max((sample.get("max_intensity", 0.0) for sample in spot["samples"]), default=0.0)
        spots_with_max.append({
            "spot_id": spot["spot_id"],
            "rt": spot["rt"],
            "mz": spot["mz"],
            "num_samples": spot["num_samples"],
            "max_intensity": sample_max,
        })
    return sorted(spots_with_max, key=lambda x: x["max_intensity"], reverse=True)[:top_n]


# --- 実行 ---
if __name__ == "__main__":
    from metabolomix.core.data_config import get_data_dir
    aef_files = glob.glob(os.path.join(str(get_data_dir()), "*.aef"))  # 探索先は環境変数で上書き可
    path = aef_files[0] if aef_files else None

    if os.path.exists(path):
        results = parse_eic_aef_css1(path)
        print(f"解析完了: {len(results)} 件のスポットを取得")
        if results:
            import pprint
            pprint.pprint(results[0])
