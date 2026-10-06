"""MS-DIAL の解析 param ファイル（`<Dataset>_param_<ts>.txt`）から、候補付けに要る値だけを読む。

読む行: `Ion mode:`、`Searched adduct ions:`（その解析の極性の一覧）、`MS1 tolerance for centroid:`、
`Retention time tolerance for alignment:`、`Minimum peak height:`。ファイルが無ければ極性ごとの既定値に落とし、どちらを
使ったかを `source` で返す。deps: stdlib だけ。
"""
from __future__ import annotations

from pathlib import Path

DEFAULT_ADDUCTS = {
    "Negative": ["[M-H]-", "[M+HCOO]-", "[M+CH3COO]-", "[M+Cl]-", "[M-H2O-H]-", "[2M-H]-", "[M-2H]2-"],
    "Positive": ["[M+H]+", "[M+NH4]+", "[M+Na]+", "[M+K]+", "[M+H-H2O]+", "[2M+H]+", "[M+2H]2+"],
}
DEFAULT_RT_WINDOW = 0.1
_KEYS = {"Ion mode": "ion_mode", "Searched adduct ions": "searched_adducts",
         "MS1 tolerance for centroid": "ms1_tolerance",
         "Retention time tolerance for alignment": "rt_tolerance_alignment",
         "Minimum peak height": "min_peak_height"}


def find_param_file(arf2_path) -> Path | None:
    """arf2_path と同じフォルダ内の param ファイルを探す。複数存在する場合は最新を返す。"""
    folder = Path(arf2_path).parent
    files = sorted(folder.glob("*_param_*.txt"), key=lambda p: p.name)
    return files[-1] if files else None


def _float(text: str) -> float | None:
    """テキストを float に変換する。失敗時は None。"""
    try:
        return float(text)
    except ValueError:
        return None


def read_analysis_params(path) -> dict:
    """param ファイルを読み、必要な値を抽出する。

    Returns:
        dict with keys: ion_mode, searched_adducts, rt_tolerance_alignment, ms1_tolerance, min_peak_height
    """
    found = {"ion_mode": None, "searched_adducts": [], "rt_tolerance_alignment": None,
             "ms1_tolerance": None, "min_peak_height": None}
    for line in Path(path).read_text(encoding="utf-8", errors="replace").splitlines():
        key, sep, value = line.partition(":")
        name = _KEYS.get(key.strip())
        if not sep or name is None:
            continue
        value = value.strip()
        if name == "searched_adducts":
            found[name] = [v.strip() for v in value.split(",") if v.strip()]
        elif name == "ion_mode":
            found[name] = value or None
        else:
            found[name] = _float(value)
    return found


def resolve_analysis_params(arf2_path, ion_mode: str | None) -> dict:
    """arf2_path の param ファイルを読み、見つからなければ極性ごとの既定値を返す。

    Returns:
        dict with keys: searched_adducts, rt_window, source ("param_file" or "default"), path
    """
    path = find_param_file(arf2_path)
    params = read_analysis_params(path) if path else None
    mode = (params or {}).get("ion_mode") or ion_mode or "Negative"
    adducts = (params or {}).get("searched_adducts") or DEFAULT_ADDUCTS.get(mode, DEFAULT_ADDUCTS["Negative"])
    rt_window = (params or {}).get("rt_tolerance_alignment") or DEFAULT_RT_WINDOW
    return {"searched_adducts": list(adducts), "rt_window": rt_window,
            "source": "param_file" if params else "default", "path": str(path) if path else None}
