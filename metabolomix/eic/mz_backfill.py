"""m/z を欠いた `.EIC.aef` スポットを、同じバッチの ARF2 の `MassCenter` で補う。

上流 `DataObjConverter` は 2023-11-20（MsdialWorkbench 20dfa60b4 "Fix a bug of
alignment (#262)"）まで、EIC に書くスポット中心 m/z を `ChromXsTop.Mz` の平均で
求めていた。LC-MS のピークではこの値が未設定の -1 なので、その頃の `.EIC.aef` は
全スポットで m/z=-1 になる（RT は正しい）。m/z 照合・m/z 検索・図の見出しが
すべてこの -1 を実値として扱って壊れるため、読み込み直後にここで補う。

借りるのは同じ AlignmentID の ARF2 レコードで、RT が `RT_TOLERANCE` 内で一致する
ときだけ（バッチ取り違えの安全弁）。補えなかったスポットは -1 のまま残す。

deps: identity_map（RT 許容差）のみ。ARF2 の読み込みは呼び出し側が行う。
"""

from __future__ import annotations

from pathlib import Path

from metabolomix.eic.identity_map import RT_TOLERANCE

EIC_SUFFIX = ".EIC.aef"


def is_mz_missing(spot: dict) -> bool:
    """MS-DIAL は未設定の座標を -1 で書く。m/z は正なので負なら欠損。"""
    return float(spot["mz"]) < 0


def has_missing_mz(spots: list[dict]) -> bool:
    return any(is_mz_missing(spot) for spot in spots)


def sibling_arf2_path(eic_path: str | Path) -> Path | None:
    """`<stem>.EIC.aef` と同じバッチの `<stem>.arf2`。無ければ None。"""
    path = Path(eic_path)
    name = path.name
    if not name.lower().endswith(EIC_SUFFIX.lower()):
        return None
    candidate = path.with_name(name[: -len(EIC_SUFFIX)] + ".arf2")
    return candidate if candidate.is_file() else None


def backfill_missing_mz(
    spots: list[dict],
    records: list[dict],
    *,
    rt_tolerance: float = RT_TOLERANCE,
) -> dict[str, int]:
    """欠損 m/z を ARF2 `MassCenter` で上書きし、`{"missing", "filled"}` を返す。

    `spots` はその場で書き換える。補ったスポットには `mz_source="arf2"` を付ける。
    """
    by_id = {}
    for record in records:
        spot_id = record.get("AlignmentID")
        if isinstance(spot_id, int) and not isinstance(spot_id, bool):
            by_id.setdefault(spot_id, record)

    missing = filled = 0
    for spot in spots:
        if not is_mz_missing(spot):
            continue
        missing += 1
        record = by_id.get(int(spot["spot_id"]))
        if record is None:
            continue
        mass = record.get("MassCenter")
        rt = record.get("RT")
        if mass is None or rt is None or float(mass) <= 0:
            continue
        if abs(float(rt) - float(spot["rt"])) > rt_tolerance:
            continue
        spot["mz"] = float(mass)
        spot["mz_source"] = "arf2"
        filled += 1
    return {"missing": missing, "filled": filled}


def backfill_caveat(counts: dict[str, int], arf2_name: str | None) -> str | None:
    """補完結果を 1 行の注記にする。欠損が無ければ None。"""
    missing, filled = counts["missing"], counts["filled"]
    if not missing:
        return None
    note = (
        f"この .EIC.aef は {missing} 件のスポットで m/z を持たない"
        "（MS-DIAL 2023-11-20 より前の版が -1 を書いた既知の不具合）。"
    )
    if filled:
        note += (
            f"うち {filled} 件は同じバッチの ARF2（{arf2_name}）の MassCenter で補った"
            "（同じ AlignmentID で RT 一致を確認済み）。"
        )
    unresolved = missing - filled
    if unresolved:
        reason = "ARF2 が見つからない" if arf2_name is None else "ARF2 で RT が一致しない"
        note += f"残り {unresolved} 件は{reason}ため m/z=-1 のままで、m/z 検索に掛からない。"
    return note
