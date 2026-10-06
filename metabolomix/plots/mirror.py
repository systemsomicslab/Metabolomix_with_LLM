"""測定と参照の対向プロット（mirror plot）。

MS-DIAL の peak spot 画面が出す `Deconvolution vs. Reference` と同じ形——
上段が測定（上向き）、下段が参照（下向き。強度に `-1` を掛ける）、横軸 m/z 共通。

`metabolomix.analysis.spectral_match.match_spectrum` が返す `alignment`
（参照グリッドの窓ごとの `{"mz", "measured", "reference", "matched"}`）を
どの m/z が一致したかの注釈として使う。数値スコアだけでは分からない
「どのフラグメントが合っていて、どれが欠けているか」を人が見て判断できる
ようにするのがこの層の役割。

`metabolomix.plots.volcano` / `eic` と同じ役割分担: payload の組み立て（純ロジック、
matplotlib 非依存）と、明示的に図が要るときだけ使う描画（`render_mirror`）を分ける。
matplotlib は関数内 import に留め、payload 組み立てだけを使う経路には載せない。

依存グラフ: このモジュールは `metabolomix.plots.render` にのみ依存する。
`mcp_core` / `session_state` / `metabolomix.tools.*` / `metabolomix.library.*` /
`metabolomix.analysis.*` は import しない（「渡されたものを描く」だけの層）。

payload の `measured` / `reference` は**生のスペクトル**（正規化前）をそのまま
保持する。正規化は描画時に行う。座標をセッションへ保持できる形にしておくのが
狙いなので、ここでは payload を返すだけに留める（セッションへの保持は別レイヤの
責務）。
"""
from __future__ import annotations

import math
from typing import TypedDict

#: v2 で `unscored_mz` / `scored_peak_count` / `unscored_peak_count` を足した
#: （採点対象外のピークを図の上で区別するため）。読み手（Use-LLLM）との契約なので
#: フィールドを足したら版を上げる。
MIRROR_PLOT_SCHEMA = "lipidmix.mirror.v2"

DEFAULT_TOP_LABELS = 8


class MirrorLabel(TypedDict):
    mz: float
    intensity: float
    side: str  # "measured" | "reference"


class MirrorPayload(TypedDict):
    plot_schema: str
    plot_type: str
    title: str
    measured: list[list[float]]
    reference: list[list[float]]
    matched_mz: list[float]
    matched_measured_mz: list[float]
    unscored_mz: list[float]
    scored_peak_count: int
    unscored_peak_count: int
    labels: list[MirrorLabel]


def _ranked_labels(
    measured_points: list[list[float]], reference_points: list[list[float]],
) -> list[MirrorLabel]:
    """測定・参照を**それぞれ自分の最大値で正規化した上で**強度順に並べる。

    生の強度のまま混ぜてランク付けすると、絶対値の大きいほうの側（例えば装置の
    生カウントである測定側）がラベル枠を独占し、参照側のピークが一本も
    ラベル付けされないことがある（側ごとにスケールが異なるため——`spectral_match`
    が測定側だけ 0〜100 に再スケールするのと同じ理由）。ランク付けはそれぞれの
    側の最大値に対する相対強度で行い、`labels` に残す `intensity` 自体は生の値
    のまま返す（描画側の表示は生値ベースのままでよい）。
    """
    measured_max = max((intensity for _, intensity in measured_points), default=0.0)
    reference_max = max((intensity for _, intensity in reference_points), default=0.0)

    ranked = [
        (intensity / measured_max if measured_max > 0 else 0.0,
         {"mz": mz, "intensity": intensity, "side": "measured"})
        for mz, intensity in measured_points
    ] + [
        (intensity / reference_max if reference_max > 0 else 0.0,
         {"mz": mz, "intensity": intensity, "side": "reference"})
        for mz, intensity in reference_points
    ]
    ranked.sort(key=lambda item: item[0], reverse=True)
    return [label for _rank, label in ranked]


def build_mirror_payload(
    measured: list[list[float]],
    reference: list[list[float]],
    alignment: list[dict],
    *,
    title: str,
    top_labels: int = DEFAULT_TOP_LABELS,
    ms2_tol: float | None = None,
    unscored_mz: list[float] | None = None,
) -> MirrorPayload:
    """測定・参照スペクトル（生の値）と一致した m/z を対向プロット用にまとめる。

    `measured` / `reference` は `[[mz, intensity], ...]` の生のスペクトル
    （正規化前）をそのまま保持する。正規化は描画時（`render_mirror`）に行う。
    `alignment` は `match_spectrum` の戻り値の同名フィールド
    （参照グリッドの窓ごとの `{"mz", "measured", "reference", "matched"}`）で、
    `matched: True` の窓の `mz` だけを `matched_mz` として抜き出す。

    **`matched_mz` は常に参照側の m/z である**（`spectral_match._matched_peaks_walk`
    が参照グリッドだけを歩いて窓の中心を決めるため）。測定ピークは許容幅
    （`ms2_tol`）の中で一致するのであって、参照ピークとぴったり同じ m/z を
    持つことはまず無い。そのため測定側の「どのピークが一致したか」は
    `matched_mz` との厳密一致では判定できず、別に `matched_measured_mz` として
    計算する:

    - `ms2_tol` を渡したとき: `matched_mz` の各窓の中心から `ms2_tol` 以内に
      ある測定ピークを一致とみなす。
    - `ms2_tol` を渡さない（既定 `None`）とき: `matched_measured_mz` は空のまま
      にする——**許容幅を勝手に決め打ちしない**。実際の値はツール層が
      store の `search_params`（`ms2_tolerance`）から取って渡す。この場合、
      描画は参照側の色分けとガイド線だけに頼る（「一致しているように見えて
      実は許容幅が違う」を防ぐため、無根拠な厳密一致では色付けしない）。

    `unscored_mz` は**採点に入らなかった測定ピークの m/z**（`.dbs` の
    `relative_amp_cutoff` / `absolute_amp_cutoff` で `normalize_measured` が
    落とした分。判定は `spectral_match.cutoff_mask`）。`measured` は今までどおり
    全ピークを保ち、この一覧は「どれが採点対象外か」を示す印として持つ
    ——`matched_measured_mz` と同じ形で、座標を二重に持たない。渡さなければ空。
    足切りが 0 の run（既定）では常に空になり、図も payload も従来と変わらない。

    `labels` は測定・参照を**それぞれの最大値で正規化した相対強度**で
    ランク付けし、上位 `top_labels` 件（測定・参照の合計で打ち切り）だけを
    残す（`_ranked_labels` 参照）。ラベルが多すぎると重なって読めなくなるための
    上限で、`build_volcano_plot_payload` の `max_points` と同じ考え方。
    """
    if isinstance(top_labels, bool) or not isinstance(top_labels, int) or top_labels < 0:
        raise ValueError("top_labels must be a non-negative integer")

    measured_points = [[float(p[0]), float(p[1])] for p in (measured or [])]
    reference_points = [[float(p[0]), float(p[1])] for p in (reference or [])]
    matched_mz = [
        float(entry["mz"]) for entry in (alignment or []) if entry.get("matched")
    ]

    matched_measured_mz: list[float] = []
    if ms2_tol is not None and matched_mz:
        tol = float(ms2_tol)
        matched_measured_mz = [
            mz for mz, _intensity in measured_points
            if any(abs(mz - center) <= tol for center in matched_mz)
        ]

    # 側ごとに切る（上流は上下で別々の `Annotator`＝別枠。片側が枠を独占しない）。
    # 採点対象外のピークはラベル枠を取らない。m/z ラベルが付くと「一致候補として
    # 見た上で外れた」と読めてしまうが、実際には採点に入っていない。
    unscored = [float(mz) for mz in (unscored_mz or [])]
    unscored_set = set(unscored)
    scored_points = [point for point in measured_points if point[0] not in unscored_set]

    ranked = _ranked_labels(scored_points, reference_points)
    labels = [
        label for side in ("measured", "reference")
        for label in [item for item in ranked if item["side"] == side][:top_labels]
    ]

    return {
        "plot_schema": MIRROR_PLOT_SCHEMA,
        "plot_type": "mirror",
        "title": title,
        "measured": measured_points,
        "reference": reference_points,
        "matched_mz": matched_mz,
        "matched_measured_mz": matched_measured_mz,
        "unscored_mz": unscored,
        "scored_peak_count": len(scored_points),
        "unscored_peak_count": len(measured_points) - len(scored_points),
        "labels": labels,
    }


_MEASURED_COLOR = "#2471a3"  # 上段・測定（volcano の down と同系統の青）
_REFERENCE_COLOR = "#c0392b"  # 下段・参照（volcano の up と同系統の赤）
_MATCHED_COLOR = "#27ae60"  # 一致した m/z の印（緑）
_UNSCORED_COLOR = "#9aa0a6"  # 足切りで採点に入らなかったピーク（灰色・薄く背後に）


def _normalized_by_max(points: list[list[float]]) -> list[list[float]]:
    """自分自身の最大強度で正規化する（生の強度だと片方が潰れるため）。"""
    if not points:
        return []
    max_intensity = max(intensity for _, intensity in points)
    if max_intensity <= 0:
        return [[mz, 0.0] for mz, _ in points]
    return [[mz, intensity / max_intensity] for mz, intensity in points]


# 縦軸のスケール。MS-DIAL GUI は Relative / Absolute / Log10 / Sqrt を上下独立に
# 選ばせる（`ObservableMsSpectrum.CreateAxisPropertySelectors2`）。ここは対向プロット
# なので **Absolute は用意しない**——測定と参照は単位が違い、同じ図の上下に生の強度を
# 並べても比較にならない（上流は片側ずつ別の図として見られるので成立している）。
RELATIVE = "relative"
SQRT = "sqrt"
LOG10 = "log10"
SCALES = (RELATIVE, SQRT, LOG10)

# log10 軸の下限。これ未満は 0 と同じ高さに潰す（対数は 0 へ向けて発散するため、
# 下限を決めないと軸が引けない）。3 桁 = 0.1% は上流の Log10 軸が既定で見せる範囲。
_LOG10_FLOOR = 1e-3


def scale_intensity(value: float, scale: str = RELATIVE) -> float:
    """正規化済み強度（0〜1）を縦軸のスケールへ写す。戻り値も 0〜1。

    `relative` は恒等。`sqrt` と `log10` は小さいピークを持ち上げるためのもので、
    **precursor がベースピークのスペクトル**（脂質の [M-H]- など）で診断イオンが
    相対数 % に潰れて読めなくなる問題に効く。値の大小関係は保つので、
    「どちらが高いか」の読み取りは変わらない。
    """
    if scale == RELATIVE:
        return value
    if scale == SQRT:
        return math.sqrt(value) if value > 0 else 0.0
    if scale == LOG10:
        if value <= _LOG10_FLOOR:
            return 0.0
        decades = -math.log10(_LOG10_FLOOR)
        return (math.log10(value) + decades) / decades
    raise ValueError(f"scale は {SCALES} のいずれかを指定してください（受け取った値: {scale!r}）。")


# ラベルの置き方。上流 `Annotator.OnRender`（`Common/ChartDrawing/Chart/`）は
# 強度降順に走査し、**既に置いたラベルと重なるものを飛ばす**——本数上限ではなく
# 幾何で決める。MS2 ビューは `TopN` を指定しない（`MsSpectrumView.xaml`）。
AUTO = "auto"
MSDIAL = "msdial"
LABEL_POLICIES = (AUTO, MSDIAL)

# ラベルの文字サイズ（pt）。代表箱の幅もこのサイズで測る。
_LABEL_FONTSIZE = 6

# 忠実版が箱の大きさを測る代表文字列。上流は**ラベル自身の文字列ではなく**
# これ 1 つで全ラベルの箱を決める（`Annotator.OnRender` の `repText`）。
_MSDIAL_REPRESENTATIVE_LABEL = "1000.00000"

# 片側あたりの上限。上流の MS2 ビューには本数上限が無いので、これは暴走止めの
# 安全弁でしかない——`msdial` では水平棄却が先に効いて到達しない。
_MAX_LABELS_PER_SIDE = 25


# ラベル箱の高さは文字サイズ × この倍率（行送り。matplotlib の `text.linespacing`
# 既定と同じ 1.2）で決める。`get_window_extent` の高さは matplotlib の版・
# ヒンティング既定で変わる（3.10.9: 9.0px / 3.11.2: 8.33px、いずれも 6pt・dpi100）
# ので、縦方向の重なり判定に実測を使うと、密集ピークの縦間隔がその差に挟まったとき
# 版によって間引きが効かなくなる（GitHub #3）。幅だけは字形に依存するので実測する。
_LABEL_LINE_HEIGHT = 1.2


def _label_box_height_px(fig) -> float:
    """ラベル箱の高さ（ピクセル）。文字サイズと dpi だけで決まる決定的な値。"""
    return _LABEL_FONTSIZE * fig.dpi / 72.0 * _LABEL_LINE_HEIGHT


def _text_extent(ax, renderer, text: str) -> tuple[float, float]:
    """文字列の描画サイズ（ピクセル）。使い捨ての Text を置いて測って消す。

    呼び出し側は**幅だけ**を使う（高さは `_label_box_height_px`）。"""
    artist = ax.text(0, 0, text, fontsize=_LABEL_FONTSIZE)
    box = artist.get_window_extent(renderer=renderer)
    artist.remove()
    return box.width, box.height


def _is_overlap(policy: str, box_a, point_a, box_b, point_b) -> bool:
    """上流の `IOverlapMethod` に対応する判定（ピクセル座標）。

    `msdial` は **水平方向だけ**を見る。上流の MS2 ビューは
    `Overlap="Horizontal, Direct"` だが、この合成は OR で、しかも全ラベルが
    同じ代表箱を使うため `Direct`（水平 AND 垂直）は `Horizontal` の部分集合に
    なる——つまり実質 `Horizontal` 単独に縮退している。縦にどれだけ離れていても
    m/z が近ければ飛ばす。

    `auto` は水平と垂直の両方が近いときだけ飛ばす（上流の `Direct` 相当）。
    対向プロットは同じ側でもピークの高さが大きく違うので、2 次元で見るほうが
    読めるラベルを多く残せる。
    """
    horizontal = (box_a[0] + box_b[0]) / 2 > abs(point_a[0] - point_b[0])
    if policy == MSDIAL:
        return horizontal
    return horizontal and (box_a[1] + box_b[1]) / 2 > abs(point_a[1] - point_b[1])


def _draw_labels(ax, fig, points_by_side: dict, policy: str) -> None:
    """強度降順に走査し、重ならないものだけ描く（上流と同じ貪欲法）。

    候補は `payload["labels"]` ではなく**スペクトル全点**から採る。衝突判定は
    ピクセル座標が要るので描画時にしかできず、payload 側の `labels` は
    自前で描くクライアント向けの要約（件数で切った版）という役割分担になる。
    """
    if policy not in LABEL_POLICIES:
        raise ValueError(
            f"label_policy は {LABEL_POLICIES} のいずれかを指定してください"
            f"（受け取った値: {policy!r}）。"
        )

    fig.canvas.draw()          # transData を確定させてからピクセルへ写す
    renderer = fig.canvas.get_renderer()
    box_height = _label_box_height_px(fig)
    representative = (_text_extent(ax, renderer, _MSDIAL_REPRESENTATIVE_LABEL)[0], box_height)

    # 上流は上下で別々の Annotator を使う＝側をまたぐ衝突は起きない。
    for side, points in points_by_side.items():
        placed: list[tuple] = []
        sign = 1 if side == "measured" else -1
        for mz, height in sorted(points, key=lambda item: item[1], reverse=True):
            if len(placed) >= _MAX_LABELS_PER_SIDE:
                break
            text = f"{mz:.4f}"
            box = representative if policy == MSDIAL else (
                _text_extent(ax, renderer, text)[0], box_height)
            point = ax.transData.transform((mz, sign * height))
            if any(_is_overlap(policy, box, point, other_box, other_point)
                   for other_box, other_point in placed):
                continue
            placed.append((box, point))
            ax.annotate(
                text, xy=(mz, sign * height),
                xytext=(0, 4 if sign > 0 else -4), textcoords="offset points",
                fontsize=_LABEL_FONTSIZE, ha="center",
                va="bottom" if sign > 0 else "top",
            )


def render_mirror(payload: MirrorPayload, *, scale: str = RELATIVE,
                  label_policy: str = AUTO):
    """対向プロットを PNG バイト列にする（`metabolomix.plots.render.figure_to_png` 経由）。

    上段が測定（上向き）、下段が参照（下向き。正規化強度に `-1` を掛ける）。
    一致した m/z のガイド線は `payload["matched_mz"]`（参照側の窓の中心）に引く。
    ステムの色分けは側ごとに別の一致リストを使う——参照側は `matched_mz` との
    厳密一致でよい（`matched_mz` 自体が参照グリッドの m/z なので）が、測定側は
    許容幅の中で一致するため `matched_mz` とは値が揃わない。`build_mirror_payload`
    が計算済みの `payload["matched_measured_mz"]`（`ms2_tol` を渡さなかった場合は
    空）を使う。ラベルは `payload["labels"]` のみ（上位 `top_labels` 本）。

    `payload["unscored_mz"]`（足切りで採点に入らなかった測定ピーク）は灰色で薄く、
    測定パネルの背後に描く。消すと「取れていない」と読めてしまうので描いた上で
    区別する。ラベルの対象にはせず、凡例の項目も**該当ピークがあるときだけ**足す。

    **ピークは素のステムだけで、先端にマーカーを打たない。** 上流の
    `LineSpectrumControlSlim` も `DrawLine` だけで描いている。マーカーは
    m/z 軸上の見かけの太さを増やして近接ピークを潰すだけで、情報を足さない。

    `label_policy` はラベルの衝突回避（`_draw_labels` / `_is_overlap`）。
    `"auto"`（既定）は水平・垂直の両方を見る 2 次元判定、`"msdial"` は上流に
    忠実な**水平のみ**の判定で、箱の大きさも代表文字列 1 つで決める。
    画像のラベルは `payload["labels"]` ではなく**スペクトル全点**から選ぶ
    （衝突判定にピクセル座標が要るため描画時にしかできない）。

    `scale` は縦軸の写し方（`scale_intensity` 参照）。上下は**それぞれ自分の
    最大値で正規化**してから同じスケールを掛ける。軸ラベルに選んだスケールを
    書く——黙って `sqrt` で描くと、読む側が相対強度の比を誤読する。
    """
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from metabolomix.plots.render import figure_to_png

    # 未知の scale / label_policy はここで弾く（描き始めてから落ちないように）。
    scale_intensity(1.0, scale)
    if label_policy not in LABEL_POLICIES:
        raise ValueError(
            f"label_policy は {LABEL_POLICIES} のいずれかを指定してください"
            f"（受け取った値: {label_policy!r}）。"
        )
    measured = [[mz, scale_intensity(v, scale)] for mz, v in _normalized_by_max(payload["measured"])]
    reference = [[mz, scale_intensity(v, scale)] for mz, v in _normalized_by_max(payload["reference"])]
    matched_mz = set(payload.get("matched_mz") or [])
    matched_measured_mz = set(payload.get("matched_measured_mz") or [])
    unscored_mz = set(payload.get("unscored_mz") or [])

    # 採点対象外のピークは同じパネルに描くが別の層にする（消すと「取れていない」と
    # 読めてしまう）。正規化の分母は上と共通＝生の最大値なので、層を分けても
    # 採点されたピークの高さは変わらない。
    scored = [point for point in measured if point[0] not in unscored_mz]
    unscored = [point for point in measured if point[0] in unscored_mz]

    fig, ax = plt.subplots(figsize=(9, 5))

    for mz in sorted(matched_mz):
        ax.axvline(mz, color=_MATCHED_COLOR, alpha=0.25, linewidth=1.0, zorder=0)

    if unscored:
        ax.vlines(
            [mz for mz, _ in unscored], 0, [intensity for _, intensity in unscored],
            colors=[_UNSCORED_COLOR] * len(unscored), linewidth=1.4, alpha=0.45, zorder=1,
        )

    if scored:
        colors = [
            _MATCHED_COLOR if mz in matched_measured_mz else _MEASURED_COLOR
            for mz, _ in scored
        ]
        ax.vlines(
            [mz for mz, _ in scored], 0, [intensity for _, intensity in scored],
            colors=colors, linewidth=1.4,
        )

    if reference:
        colors = [_MATCHED_COLOR if mz in matched_mz else _REFERENCE_COLOR for mz, _ in reference]
        ax.vlines(
            [mz for mz, _ in reference], 0, [-intensity for _, intensity in reference],
            colors=colors, linewidth=1.4,
        )

    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_xlabel("m/z")
    ax.set_ylabel(f"Relative intensity [{scale}] (measured / reference)")
    ax.set_ylim(-1.15, 1.15)
    ax.set_title(payload.get("title") or "")

    from matplotlib.lines import Line2D

    legend_handles = [
        Line2D([0], [0], color=_MEASURED_COLOR, lw=1.4, label="Measured"),
        Line2D([0], [0], color=_REFERENCE_COLOR, lw=1.4, label="Reference"),
        Line2D([0], [0], color=_MATCHED_COLOR, lw=1.4, label="Matched"),
    ]
    # 採点対象外の層は**あるときだけ**凡例に出す。足切り 0 の run（既定）で
    # 常に出すと、無いものを探させることになる。
    if unscored:
        legend_handles.append(
            Line2D([0], [0], color=_UNSCORED_COLOR, lw=1.4, alpha=0.45,
                   label="Below cutoff (not scored)"))
    ax.legend(handles=legend_handles, fontsize=8, loc="best")

    # ラベルは軸範囲が確定してから置く（ピクセル座標で重なりを判定するため）。
    # 採点対象外のピークはラベルの対象にしない（payload 側の `labels` と揃える）。
    _draw_labels(ax, fig, {"measured": scored, "reference": reference}, label_policy)

    return figure_to_png(fig)
