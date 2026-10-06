"""参照ライブラリ照合の既定許容幅（stdlib のみの leaf）。

store の `search_params` が無い（`.msp` 由来など）ときのフォールバックとして、
`metabolomix.library.tools`（MCP 面）と `metabolomix.msdial.peak_verification`
（`verify_peak_annotation` の判断材料）の**両方**がここから読む。

以前は 2 箇所に同じ 3 つのリテラルを複写していた。上流の既定が変わったとき、
MCP に面していて目に付きやすい `library/tools.py` だけが更新され、
`peak_verification.py` 側がテストに気付かれず古くなる恐れがあったため、
値の置き場をここ 1 箇所に一本化した（Task 11 レビュー Important 1）。

`library.store`（sqlite3 / msgpack / reader 群に依存）には置かない——定数 3 つの
ためだけに `peak_verification.py`（下位レイヤ）へその依存を引き込むのは割に合わない。

**出所の訂正（最終レビュー Important 6）**: 以前の docstring は 3 定数すべてを
「`MsRefSearchParameterBase` の既定値」と説明していたが誤り。実際は 2 つの
別の窓が混じっている:

- `DEFAULT_MZ_TOL` (0.01) / `DEFAULT_MS2_TOL` (0.025) は **ライブラリ照合**用
  （`.dbs` の `MsRefSearchParameterBase`。`ms1_tolerance` / `ms2_tolerance`）。
  `.msp` 由来で `search_params` が無いときのフォールバックに使う。
- `DEFAULT_RT_TOL` (0.2) は上流 `MsRefSearchParameterBase.RtTolerance` の既定
  （`docs/schema/molecule_ms_reference.md` によれば **100.0**）**ではない**。
  これは `metabolomix/dcl/reader.py` / `dcl/tools.py` の `.dcl` 検索（precursor m/z
  から測定 MS/MS を引くときの RT 窓）の既定値で、値そのもの（0.2 分）は妥当だが
  出所が違う。`library.tools._measured_spectrum` が `.dcl` から測定スペクトルを
  引く窓としてここから読む（Important 7: ライブラリ照合用の許容幅と混ぜてはいけない
  ——`.dbs` の `RtTolerance`（実質 RT フィルタを無効化する 100.0）をこの窓に流用すると、
  同じ precursor m/z の別ピークを黙って拾う）。
"""
from __future__ import annotations

DEFAULT_MZ_TOL = 0.01
DEFAULT_MS2_TOL = 0.025
DEFAULT_RT_TOL = 0.2


def pick_tol(explicit: float | None, search_params: dict, key: str, default: float) -> float:
    """許容幅の決定を 1 箇所に一本化する: 明示指定 > store の `search_params` > 既定値。

    `search_params.get(key)` が `0.0` のような偽値でも正しく採用されるよう、
    `or` ではなく `is None` で判定する。

    最終レビュー Important 5: `metabolomix.library.tools._pick_tol` が正しくこの
    判定をしていた一方、`metabolomix.msdial.peak_verification._spectral_match_for_feature`
    は `search_params.get(key) or default` という `or` 判定を独自に持っていて、
    `0.0`（"足切りなし"のような正当な値）を黙って既定値に差し替えていた。
    2 箇所が同じロジックを別々に書くと再びずれるので、ここに一本化して両方から
    import する。
    """
    if explicit is not None:
        return explicit
    value = search_params.get(key)
    if value is not None:
        return value
    return default
