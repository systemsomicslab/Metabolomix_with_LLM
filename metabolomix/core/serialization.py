"""LLM へ返す文字列の直列化。依存グラフの leaf（stdlib のみ）。

このサーバの戻り値は最終的に LLM の文脈を占める。`json.dumps(..., indent=2)` は
1 ネストにつき改行＋空白 2 文字を全キーに掛けるため、実測で戻り値の 15〜57% が
空白だった（`dataset_status` 3,513→1,502 字、`dcl_find_msms` 953→441 字）。
人間向けの整形はクライアント側の責務なので、サーバは最小形で返す。

tools_* / server / mcp_core を import しない（mcp_errors からも使うため）。
numpy も import しない（numpy の値は `dtype` と `tolist` を持つかで見分ける）。
"""
import json
import math

# m/z は 4 桁、RT は 3 桁（クロマトのピーク幅は秒オーダー）、それ以外の指標も
# 6 桁あれば分析的な情報は落ちない。float の既定 repr は 17 桁まで出すため、
# 座標配列では桁の大半が無意味な尾数になる。
DEFAULT_FLOAT_DIGITS = 6

#: `truncate_lists` が既定で切る一覧のキー。統計エンジンの警告・前処理の段の記録に出る、
#: 個体・試料・特徴の数だけ伸びる一覧（`examples` は既に切られていることもある）。
TRUNCATE_KEYS = ("unit_ids", "samples", "sample_ids", "examples", "excluded_by_include",
                 "dropped_ids", "subjects", "factors", "factor_samples", "unscaled_samples",
                 "qc_ids", "unknown_samples", "unknown_features", "columns", "kept_features")


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


def _safe_float(value: float):
    if math.isnan(value):
        return None
    if math.isinf(value):
        return "inf" if value > 0 else "-inf"
    return float(value)


def _safe_key(key):
    if isinstance(key, str):
        return key
    key = json_safe(key)
    return key if key is None or isinstance(key, (str, int, float, bool)) else str(key)


def json_safe(obj):
    """厳密な JSON（NaN・Infinity を含まない）に直せる形へ写す（非破壊）。

    - float: NaN → None、+inf → `"inf"`、-inf → `"-inf"`（`json.dumps` の既定は `NaN` /
      `Infinity` を書き、厳密な JSON の読み手はそこで止まる。統計の `df_prior` は無限大を取りうる）。
    - numpy のスカラー・配列 → Python の値・入れ子の list（中の float も上の規則で）。
    - set / frozenset → 並べた list（並べられない混在は型名と repr の順。決定的にする）。
    - tuple → list、dict のキーは JSON が受ける型へ。それ以外の型はそのまま返す
      （直列化できないものを黙って文字列にしない）。
    bool は int の派生なので先に素通しする。
    """
    if obj is None or isinstance(obj, (bool, str)):
        return obj
    if hasattr(obj, "dtype") and hasattr(obj, "tolist"):   # numpy のスカラー・配列
        return json_safe(obj.tolist())
    if isinstance(obj, float):
        return _safe_float(obj)
    if isinstance(obj, int):
        return int(obj)
    if isinstance(obj, dict):
        return {_safe_key(key): json_safe(value) for key, value in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [json_safe(item) for item in obj]
    if isinstance(obj, (set, frozenset)):
        items = [json_safe(item) for item in obj]
        try:
            return sorted(items)
        except TypeError:
            return sorted(items, key=lambda item: (type(item).__name__, repr(item)))
    return obj


def truncate_lists(obj, limit: int = 20, keys=TRUNCATE_KEYS):
    """指定のキーの list を先頭 `limit` 件に切り、隣（直後）に `<key>_n_more` を置く（非破壊）。

    戻り値を肥大させないため。入れ子の dict / list も辿る。`<key>_n_more` が既にあれば
    （エンジンが先に切った一覧）、切った件数をそこへ足す（数を失わない）。
    """
    if isinstance(limit, bool) or not isinstance(limit, int) or limit < 0:
        raise ValueError(f"limit は 0 以上の整数: {limit!r}")
    wanted = frozenset(keys)

    def walk(value):
        if isinstance(value, dict):
            cut = {k for k, v in value.items()
                   if k in wanted and isinstance(v, (list, tuple)) and len(v) > limit}
            out = {}
            for k, v in value.items():
                if isinstance(k, str) and k.endswith("_n_more") and k[:-len("_n_more")] in cut:
                    continue                              # 切った一覧の直後に置き直す
                if k in cut:
                    out[k] = [walk(item) for item in list(v)[:limit]]
                    prev = value.get(f"{k}_n_more")
                    prev = prev if isinstance(prev, int) and not isinstance(prev, bool) else 0
                    out[f"{k}_n_more"] = prev + len(v) - limit
                else:
                    out[k] = walk(v)
            return out
        if isinstance(value, (list, tuple)):
            return [walk(item) for item in value]
        return value

    return walk(obj)
