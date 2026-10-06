"""LLM へ返す文字列の直列化。依存グラフの leaf（stdlib のみ）。

このサーバの戻り値は最終的に LLM の文脈を占める。`json.dumps(..., indent=2)` は
1 ネストにつき改行＋空白 2 文字を全キーに掛けるため、実測で戻り値の 15〜57% が
空白だった（`arf_list_classes` 3,513→1,502 字、`dcl_find_msms` 953→441 字）。
人間向けの整形はクライアント側の責務なので、サーバは最小形で返す。

tools_* / server / mcp_core を import しない（mcp_errors からも使うため）。
"""
import json

# m/z は 4 桁、RT は 3 桁（クロマトのピーク幅は秒オーダー）、それ以外の指標も
# 6 桁あれば分析的な情報は落ちない。float の既定 repr は 17 桁まで出すため、
# 座標配列では桁の大半が無意味な尾数になる。
DEFAULT_FLOAT_DIGITS = 6


def json_payload(obj) -> str:
    """ツールの戻り値を最小トークンの JSON 文字列にする。"""
    return json.dumps(obj, ensure_ascii=False, separators=(",", ":"))


def round_floats(obj, digits: int = DEFAULT_FLOAT_DIGITS):
    """入れ子の dict/list を辿って float を丸める（非破壊・型は保つ）。

    座標点列のように同じ形の数値が数百〜数千個並ぶ payload で効く。bool は int の
    派生だが数値として丸めてはいけないので素通しする。
    """
    if isinstance(obj, float):
        return round(obj, digits)
    if isinstance(obj, list):
        return [round_floats(item, digits) for item in obj]
    if isinstance(obj, tuple):
        return tuple(round_floats(item, digits) for item in obj)
    if isinstance(obj, dict):
        return {key: round_floats(value, digits) for key, value in obj.items()}
    return obj
