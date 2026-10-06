"""対向プロット。座標は payload に持ち、描画は PNG で返す。"""
import pytest

from metabolomix.plots import mirror


def test_the_payload_keeps_both_spectra_and_the_matches():
    measured = [[100.0, 999.0], [200.0, 500.0]]
    reference = [[100.0, 999.0], [300.0, 250.0]]
    alignment = [{"mz": 100.0, "measured": 1.0, "reference": 1.0, "matched": True},
                 {"mz": 300.0, "measured": 0.0, "reference": 0.25, "matched": False}]

    payload = mirror.build_mirror_payload(measured, reference, alignment, title="GABA")

    assert payload["title"] == "GABA"
    assert payload["measured"] == [[100.0, 999.0], [200.0, 500.0]]
    assert payload["reference"] == [[100.0, 999.0], [300.0, 250.0]]
    assert payload["matched_mz"] == [100.0]


def test_the_labels_are_capped():
    """`top_labels` は**側ごと**の上限（上流は上下で別々の `Annotator`）。"""
    measured = [[float(i), float(1000 - i)] for i in range(50)]
    payload = mirror.build_mirror_payload(measured, measured, [], title="x", top_labels=3)
    assert len(payload["labels"]) <= 3 * 2


def test_render_returns_a_png():
    payload = mirror.build_mirror_payload([[100.0, 999.0]], [[100.0, 999.0]],
                                          [{"mz": 100.0, "measured": 1.0,
                                            "reference": 1.0, "matched": True}], title="t")
    png = mirror.render_mirror(payload)
    assert png[:8] == b"\x89PNG\r\n\x1a\n"


def test_measured_side_is_matched_within_tolerance():
    # alignment["mz"] is always a reference-grid m/z (see spectral_match._matched_peaks_walk),
    # so an exact-equality check against measured m/z would basically never fire on real data
    # (the measured peak sits within ms2_tol of the window center, not exactly on it).
    measured = [[100.021, 999.0]]
    reference = [[100.019, 950.0]]
    alignment = [{"mz": 100.019, "measured": 1.0, "reference": 1.0, "matched": True}]

    payload = mirror.build_mirror_payload(measured, reference, alignment, title="t", ms2_tol=0.01)

    assert payload["matched_mz"] == [100.019]
    assert payload["matched_measured_mz"] == [100.021]


def test_measured_side_is_not_matched_without_a_tolerance():
    # Without an explicit ms2_tol we must not guess one, so the measured side stays unmarked.
    measured = [[100.021, 999.0]]
    reference = [[100.019, 950.0]]
    alignment = [{"mz": 100.019, "measured": 1.0, "reference": 1.0, "matched": True}]

    payload = mirror.build_mirror_payload(measured, reference, alignment, title="t")

    assert payload["matched_measured_mz"] == []


def test_labels_are_ranked_per_side_not_by_raw_scale():
    # Measured intensities are raw instrument counts (tens of thousands); reference
    # intensities are library-relative (hundreds). Ranking by raw intensity would let
    # measured peaks crowd out every reference label.
    measured = [[100.0, 45000.0], [150.0, 30000.0], [200.0, 12000.0], [250.0, 8000.0]]
    reference = [[100.0, 999.0], [110.0, 800.0], [120.0, 650.0], [130.0, 500.0]]

    payload = mirror.build_mirror_payload(measured, reference, [], title="x", top_labels=4)

    sides = {label["side"] for label in payload["labels"]}
    assert sides == {"measured", "reference"}


# --------------------------------------------------------------------------
# 縦軸のスケール。MS-DIAL GUI は Relative / Absolute / Log10 / Sqrt を
# 上下独立に選ばせる（`ObservableMsSpectrum.CreateAxisPropertySelectors2`）。
# precursor がベースピークのスペクトルは Relative だと診断イオンが潰れる。
# --------------------------------------------------------------------------
def _payload():
    return mirror.build_mirror_payload(
        [[100.0, 20.0], [900.0, 999.0]], [[100.0, 25.0], [900.0, 999.0]],
        [{"mz": 100.0, "measured": 0.02, "reference": 0.025, "matched": True}], title="t")


def test_render_accepts_the_upstream_axis_scales():
    for scale in ("relative", "sqrt", "log10"):
        png = mirror.render_mirror(_payload(), scale=scale)
        assert png[:8] == b"\x89PNG\r\n\x1a\n"


def test_render_defaults_to_relative():
    assert mirror.render_mirror(_payload()) == mirror.render_mirror(_payload(), scale="relative")


def test_render_rejects_an_unknown_scale():
    """黙って relative に落とすと「Sqrt で見ている」と誤認したまま読むことになる。"""
    with pytest.raises(ValueError):
        mirror.render_mirror(_payload(), scale="ln")


def test_sqrt_lifts_a_peak_that_relative_flattens():
    """precursor 優勢のスペクトルで小さな診断イオンが見える高さになること。"""
    assert mirror.scale_intensity(0.02, "relative") == pytest.approx(0.02)
    assert mirror.scale_intensity(0.02, "sqrt") == pytest.approx(0.1414213562)
    assert mirror.scale_intensity(0.02, "log10") > 0.02


def test_log10_floors_at_the_documented_decade_and_never_goes_negative():
    """0 と「0.1% 未満」は同じ高さ 0 に落とす（対数は下へ発散するため）。"""
    assert mirror.scale_intensity(0.0, "log10") == 0.0
    assert mirror.scale_intensity(1e-9, "log10") == 0.0
    assert mirror.scale_intensity(1.0, "log10") == pytest.approx(1.0)


def test_peaks_are_bare_stems_without_tip_markers(monkeypatch):
    """ピーク先端にドットを打たない。上流 `LineSpectrumControlSlim` も
    `DrawLine` だけで描いており、マーカーは m/z 軸上の位置を太らせて
    近接ピークを潰すだけで情報を足さない。"""
    captured = {}

    from metabolomix.plots import render as plot_render
    real = plot_render.figure_to_png

    def spy(fig):
        captured["axes"] = fig.axes[0]
        return real(fig)

    monkeypatch.setattr(plot_render, "figure_to_png", spy)
    mirror.render_mirror(_payload())

    from matplotlib.collections import LineCollection, PathCollection
    collections = captured["axes"].collections
    assert any(isinstance(c, LineCollection) for c in collections)   # ステムは描く
    assert not any(isinstance(c, PathCollection) for c in collections)  # ドットは描かない


# --------------------------------------------------------------------------
# ラベルの衝突回避。上流 `Annotator.OnRender` は強度降順に走査し、既に置いた
# ラベルと重なるものを飛ばす（本数上限ではなく幾何で決める）。
# --------------------------------------------------------------------------
def _labelled(payload, **kwargs):
    """描画した図から、実際に置かれたラベルの文字列を拾う。"""
    captured = {}
    from metabolomix.plots import render as plot_render
    real = plot_render.figure_to_png

    def spy(fig):
        captured["texts"] = [t.get_text() for t in fig.axes[0].texts]
        return real(fig)

    import pytest as _pytest
    monkeypatch = _pytest.MonkeyPatch()
    try:
        monkeypatch.setattr(plot_render, "figure_to_png", spy)
        mirror.render_mirror(payload, **kwargs)
    finally:
        monkeypatch.undo()
    return captured["texts"]


def _crowded_payload():
    """m/z 軸が広い中に、狭い範囲へ密集したピークを置く（衝突を必ず起こす形）。

    衝突判定は**ピクセル距離**なので、m/z の差が小さいだけでは足りない——
    軸の幅に対して相対的に近くないと重ならない。そこで遠方に 1 本置いて軸を
    広げたうえで、0.1 m/z 幅に 12 本を詰める。
    """
    cluster = [[500.0 + i * 0.01, 1000.0 - i * 50] for i in range(12)]
    spread = [[100.0, 200.0], [1000.0, 200.0]]
    return mirror.build_mirror_payload(
        spread + cluster, spread + cluster, [], title="crowded")


def test_labels_no_longer_overlap_by_default():
    """既定（auto）では密集した 12 本を全部は描かない（以前は重ねて全部描いていた）。"""
    texts = _labelled(_crowded_payload())
    assert 0 < len(texts) < 2 * 14


def test_label_box_height_does_not_depend_on_the_measured_text_height(monkeypatch):
    """箱の高さを描画器の実測に頼らない（GitHub #3）。

    `get_window_extent` の高さはフォント・ヒンティングの既定が違う matplotlib の版で
    変わる（3.10.9 は 9.0px、3.11.2 は 8.33px）。密集ピークの縦間隔（8.37px）が
    その間に挟まると、ある版では重なり判定が成り立たず全ラベルが残る。高さは
    文字サイズから決める決定的な値にし、実測が小さく返ってもラベルは間引かれる。
    """
    monkeypatch.setattr(mirror, "_text_extent", lambda ax, renderer, text: (38.0, 8.0))
    texts = _labelled(_crowded_payload())
    assert 0 < len(texts) < 2 * 14


def test_the_msdial_policy_ignores_vertical_separation():
    """忠実版は水平距離だけを見る（`Overlap="Horizontal, Direct"` は全ラベルが
    同じ代表箱なので Direct ⊆ Horizontal に縮退する）。高さが違っても
    m/z が近ければ飛ばすので、2 次元で見る auto より必ずラベルが少ない。"""
    payload = _crowded_payload()
    assert len(_labelled(payload, label_policy="msdial")) \
        < len(_labelled(payload, label_policy="auto"))


def test_both_policies_label_well_separated_peaks():
    """m/z が十分離れていればどちらの方式でも全部に付く。"""
    payload = mirror.build_mirror_payload(
        [[100.0, 900.0], [400.0, 800.0], [700.0, 700.0]],
        [[100.0, 900.0], [400.0, 800.0], [700.0, 700.0]], [], title="sparse")
    for policy in ("auto", "msdial"):
        assert len(_labelled(payload, label_policy=policy)) == 6   # 上下 3 本ずつ


def test_the_two_sides_do_not_collide_with_each_other():
    """上流は上下で別々の `Annotator` を使うので側をまたぐ衝突は起きない。
    同じ m/z のピークが上下にあっても両方ラベルが付くこと。"""
    payload = mirror.build_mirror_payload(
        [[300.0, 900.0]], [[300.0, 900.0]], [], title="both sides")
    assert len(_labelled(payload, label_policy="msdial")) == 2


def test_render_rejects_an_unknown_label_policy():
    with pytest.raises(ValueError):
        mirror.render_mirror(_payload(), label_policy="clever")


def test_the_payload_label_cap_is_per_side():
    """上流は上下で別の `Annotator`＝別枠。片側が枠を独占しないよう側ごとに切る。"""
    measured = [[float(i), 1000.0] for i in range(20)]
    reference = [[float(i), 1.0] for i in range(20)]   # 相対では測定側と同順位
    payload = mirror.build_mirror_payload(measured, reference, [], title="x", top_labels=3)
    sides = [label["side"] for label in payload["labels"]]
    assert sides.count("measured") == 3
    assert sides.count("reference") == 3


# --------------------------------------------------------------------------
# 採点対象外のピーク（`.dbs` の amp cutoff で `normalize_measured` が落とした分）。
# 図には出るのに採点には入らない——今までこれが見分けられなかった。
# --------------------------------------------------------------------------
def _figure(payload, **kwargs):
    """描画した図の Axes を返す（`figure_to_png` を spy して横取りする）。"""
    captured = {}
    from metabolomix.plots import render as plot_render
    real = plot_render.figure_to_png

    def spy(fig):
        captured["axes"] = fig.axes[0]
        return real(fig)

    import pytest as _pytest
    monkeypatch = _pytest.MonkeyPatch()
    try:
        monkeypatch.setattr(plot_render, "figure_to_png", spy)
        mirror.render_mirror(payload, **kwargs)
    finally:
        monkeypatch.undo()
    return captured["axes"]


def _unscored_payload():
    return mirror.build_mirror_payload(
        [[100.0, 999.0], [200.0, 5.0]], [[100.0, 999.0]], [],
        title="t", unscored_mz=[200.0])


def test_the_payload_separates_scored_peaks_from_the_ones_the_cutoff_dropped():
    payload = _unscored_payload()

    assert payload["unscored_mz"] == [200.0]
    assert payload["scored_peak_count"] == 1
    assert payload["unscored_peak_count"] == 1
    # 生のスペクトルは今までどおり全ピークを保つ（描画も座標も壊さない）。
    assert payload["measured"] == [[100.0, 999.0], [200.0, 5.0]]


def test_a_dropped_peak_does_not_take_a_label_slot():
    """ラベル枠は採点に入ったピークのものを優先する。採点外のピークに
    m/z ラベルが付くと「一致候補として見た」と読めてしまう。"""
    payload = mirror.build_mirror_payload(
        [[100.0, 10.0], [200.0, 999.0]], [], [], title="t", unscored_mz=[200.0])

    assert [label["mz"] for label in payload["labels"]] == [100.0]


def test_the_schema_version_announces_the_added_layer():
    """追加フィールドは読み手（Use-LLLM）との契約。黙って足さない。"""
    assert mirror.MIRROR_PLOT_SCHEMA == "lipidmix.mirror.v2"
    assert _unscored_payload()["plot_schema"] == "lipidmix.mirror.v2"


def test_the_tool_description_names_the_current_schema():
    """ツール説明は LLM が payload の読み方を決める根拠。古い版名を載せない。"""
    import asyncio

    from metabolomix.core.mcp_core import mcp
    import metabolomix.library.tools  # noqa: F401  (登録の副作用)

    tools = {t.name: t for t in asyncio.run(mcp.list_tools())}
    description = tools["library_plot_mirror"].description
    assert mirror.MIRROR_PLOT_SCHEMA in description
    assert "lipidmix.mirror.v1" not in description


def test_a_payload_without_dropped_peaks_reports_an_empty_layer():
    payload = _payload()

    assert payload["unscored_mz"] == []
    assert payload["unscored_peak_count"] == 0


def test_the_legend_names_the_dropped_layer_only_when_there_is_one():
    plain = [t.get_text() for t in _figure(_payload()).get_legend().get_texts()]
    with_dropped = [t.get_text() for t in _figure(_unscored_payload()).get_legend().get_texts()]

    assert not any("cutoff" in text for text in plain)
    assert any("cutoff" in text for text in with_dropped)


def test_a_dropped_peak_is_still_drawn_but_in_its_own_colour():
    """消してしまうと「取れていない」と読めてしまう。描いた上で区別する。"""
    axes = _figure(_unscored_payload())

    colours = {}
    for collection in axes.collections:
        for segment, colour in zip(collection.get_segments(), collection.get_colors()):
            colours[round(float(segment[0][0]), 4)] = tuple(colour)

    assert 200.0 in colours                      # 描かれている
    assert colours[200.0] != colours[100.0]      # 採点されたピークと見分けられる
