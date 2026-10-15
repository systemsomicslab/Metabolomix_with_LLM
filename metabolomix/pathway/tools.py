# metabolomix/pathway/tools.py
"""MCP ツール: MSLipidMapper を通したパスウェイ投影（spec 2026-10-14）。

`metabolomix/pathway/` のうち session_state / mcp_core を import してよいのはこのモジュールだけ。
handoff / jsonld / compare / bridge は純ロジックのまま保つ。
"""
from __future__ import annotations

import json
import re
import uuid
from datetime import datetime, timezone
from pathlib import Path

from mcp.server.fastmcp import Image
from mcp.types import ToolAnnotations

from metabolomix.core import mcp_core, mcp_errors, session_state
from metabolomix.core.mcp_core import mcp
from metabolomix.core.serialization import json_payload
from metabolomix.pathway import bridge, compare, handoff, jsonld
from metabolomix.plots import pdf_raster

__all__ = ["pathway_project", "pathway_plot"]

_NODE_COLUMNS = ["node", "scope", "n_features", "mean_a", "mean_b", "log2fc", "q", "call"]
_PAIR_COLUMNS = ["node_1", "call_1", "log2fc_1", "node_2", "call_2", "log2fc_2", "via"]
_ADJACENCY_NOTE = (
    "隣接は MSLipidMapper のパスウェイ図のつながり（接合ノードだけを挟むもの）で、"
    "上流は反応だとは主張していない（not curated reactions）。逆向きの組は観察の提示であり、"
    "基質・生成物の解釈は文献で裏を取ること。"
    "掲載条件は両端とも q < 閾値かつ向きが逆。片方が ns/NA の組は件数だけ示す。")
#: pathway_project が採番する run_id の形。これ以外は reports/pathway の外を指しうるので受け付けない。
_RUN_ID_RE = re.compile(r"pw_\d{8}_\d{6}_[0-9a-f]{6}")
_FIGURE_NOTE = ("図のノード内の p 値ラベルは MSLipidMapper（Shiny と同じ計算）のもので、"
                "pathway_project の q 値とは別物。")


def _backend_text(result: dict) -> str:
    """どちらの MSLipidMapper で動いたか（サーバからローカルへ回ったならその理由も）の一文。"""
    backend = result.get("backend") or {}
    if backend.get("kind") == "remote":
        text = f"研究室サーバ {backend.get('api_url')}（job {backend.get('job_id')}）"
    else:
        text = "ローカルの R"
    fallback = result.get("remote_fallback")
    if fallback:
        text += (f"（サーバ {fallback.get('api_url')} が使えずローカルへ回した: "
                 f"{fallback.get('code')} {fallback.get('message')}）")
    return text


def _arf2_catalog():
    """兄弟 .arf2 のスポット一覧。無ければ (None, missing_state 封筒)。"""
    from metabolomix.arf.tools import _require_arf_with_arf2
    from metabolomix.arf2.reader import load_catalog
    found = _require_arf_with_arf2()
    if isinstance(found, str):
        return None, found
    return load_catalog(found), None


def _run_dir(run_id: str) -> Path:
    return mcp_core._resolve_report_dir() / "pathway" / run_id


def _load_manifest(run_id: str | None):
    """run の manifest（dict）か、run が無いときの missing_state 封筒（str）。"""
    run_id = run_id or session_state.session.pathway.last_run_id
    if not run_id:
        return mcp_errors.missing_state("pathway_run", ["pathway_project"],
                                        "先に pathway_project を実行してください。")
    if not _RUN_ID_RE.fullmatch(str(run_id)):
        return mcp_errors.pathway_error(
            "PATHWAY_BAD_REQUEST",
            f"run_id の形が不正です: {run_id!r}（pathway_project が返す pw_YYYYmmdd_HHMMSS_xxxxxx）。")
    path = _run_dir(run_id) / "manifest.json"
    if not path.is_file():
        return mcp_errors.missing_state("pathway_run", ["pathway_project"],
                                        f"run {run_id} の manifest.json がありません。")
    try:
        manifest = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        return mcp_errors.pathway_error("MSLM_OUTPUT_INVALID",
                                        f"run {run_id} の manifest.json を読めません: {exc}",
                                        {"manifest": str(path)})
    if not isinstance(manifest, dict):
        return mcp_errors.pathway_error("MSLM_OUTPUT_INVALID",
                                        f"run {run_id} の manifest.json が JSON オブジェクトではありません。",
                                        {"manifest": str(path)})
    return manifest


def _pick_source(source: str) -> str:
    s = session_state.session
    if source != "auto":
        return source
    if s.arf.feature_matrix is not None:
        return "arf"
    if s.dataset is not None:
        return "dataset"
    return "arf"


def _recipe_text(recipe: dict) -> str:
    parts = [f"{k}={v}" for k, v in sorted((recipe or {}).items()) if v not in (None, "", [])]
    return "; ".join(parts) or "none"


def _build_arf(group_a, group_b):
    from metabolomix.arf.tools import _pool_group_labels
    a = session_state.session.arf
    if a.feature_matrix is None:
        return None, mcp_errors.missing_state(
            "preprocessed_matrix", ["arf_preprocess"],
            "先に arf_preprocess を実行してください（前処理後行列が必要）。")
    catalog, err = _arf2_catalog()
    if err:
        return None, err
    labels = [None if (a.sample_meta.get(n) or {}).get("role", "sample") != "sample"
              else (a.sample_meta.get(n) or {}).get("group") for n in a.pp_sample_names]
    try:
        _, _, resolved = _pool_group_labels(a.pp_sample_names, labels, group_a, group_b)
    except ValueError as exc:
        raise handoff.HandoffError("PATHWAY_BAD_REQUEST", str(exc)) from exc
    table = handoff.build_from_arf(
        matrix=a.feature_matrix, sample_names=a.pp_sample_names,
        feature_names=a.pp_feature_names, sample_meta=a.sample_meta,
        excluded_samples=a.excluded_samples, excluded_spots=a.excluded_spots,
        catalog={int(s["MasterAlignmentID"]): s for s in catalog},
        props=(a.preprocessing_recipe or {}).get("props"),
        group_a=str(group_a), group_b=str(group_b), compared=resolved)
    return (table, _recipe_text(a.preprocessing_recipe), a.dataset_id,
            str(group_a), str(group_b)), None


def _build_dataset(group_a, group_b, label_a, label_b):
    ds = session_state.session.dataset
    if ds is None:
        return None, mcp_errors.missing_state("dataset", ["dataset_load"],
                                              "先に dataset_load を実行してください。")
    if ds.pp_matrix is None:
        return None, mcp_errors.missing_state("preprocessed_matrix", ["dataset_preprocess"],
                                              "先に dataset_preprocess を実行してください。")
    as_list = (lambda g: [g] if isinstance(g, str) else list(g or []))
    table = handoff.build_from_dataset(
        pp_matrix=ds.pp_matrix, pp_sample_names=ds.pp_sample_names,
        pp_feature_names=ds.pp_feature_names, roles=ds.roles,
        sample_names=ds.sample_names, sample_ids=ds.sample_ids,
        metadata_rows=ds.sample_metadata_rows, annotations=ds.feature_annotations,
        group_a=as_list(group_a), group_b=as_list(group_b),
        group_a_label=label_a, group_b_label=label_b)
    return (table, _recipe_text(ds.preprocessing_recipe), ds.dataset_id,
            label_a, label_b), None


# 呼ぶたびに新しい run（reports/pathway/<run_id>/）を作る → 冪等ではない。
# [mslipidmapper] api_url があれば研究室サーバへ中間 mzTab-M を送る（ローカルの R はネットワークに出ない）。
@mcp.tool(annotations=ToolAnnotations(
    readOnlyHint=False, destructiveHint=False, idempotentHint=False, openWorldHint=True),
    structured_output=False)
def pathway_project(source: str = "auto", group_a: str | list[str] | None = None,
                    group_b: str | list[str] | None = None,
                    group_a_label: str = "group_a", group_b_label: str = "group_b",
                    pathways: list[str] | None = None, acyl_chains: list[str] | None = None,
                    aggregation: str = "sum", q_threshold: float = 0.05) -> str:
    """前処理済みの脂質データを MSLipidMapper のパスウェイへ投影し、2 群を比べる。

    返すのは小さな表だけ（ノード表と、逆向きに動く隣接ペア）。JSON-LD・中間 mzTab-M は
    解析フォルダの reports/pathway/<run_id>/ に残る。図は pathway_plot。

    - source: auto | arf | dataset。arf は arf_preprocess 後、dataset は dataset_preprocess 後。
      auto は「前処理済み ARF があれば arf、無ければ dataset_load 済みなら dataset、
      どちらも無ければ arf（→ arf_preprocess を求める missing_state）」の順で選ぶ。
    - group_a / group_b: arf は Class ID か因子トークン（arf_differential と同じ）、
      dataset はサンプル名のリスト（group_a_label / group_b_label が群名になる）。
    - pathways: remodeling / ceramide / global（既定は 3 つ全部）。
    - acyl_chains: 例 ["18:0","20:4"]。指定すると鎖サブセットのスコープも比べる。
    - aggregation: sum | mean | median（サブクラス内の分子種の集計）。
    検定は log2(x+1) の Welch、BH はパスウェイ×スコープ種別ごと。log2fc は正なら group_b が高い。
    脂質アッセイ（record_objective の assay_kind=lipid）でだけ使える。
    MSLipidMapper は [mslipidmapper] api_url の研究室サーバを先に使い、使えなければローカルの R
    （[mslipidmapper] rscript）で動かす。どちらで動いたかは戻り値の backend。
    """
    s = session_state.session
    if s.assay_kind != "lipid":
        return mcp_errors.missing_state(
            "assay_kind", ["record_objective"],
            "MSLipidMapper は脂質専用です。record_objective(assay_kind='lipid') で種別を確定してください。")
    pathways = list(pathways or bridge.PATHWAY_NETWORKS)
    unknown = [p for p in pathways if p not in bridge.PATHWAY_NETWORKS]
    if unknown:
        return mcp_errors.pathway_error(
            "PATHWAY_BAD_REQUEST", f"未知のパスウェイ: {', '.join(unknown)}",
            {"valid": list(bridge.PATHWAY_NETWORKS)})
    if aggregation not in ("sum", "mean", "median"):
        return mcp_errors.pathway_error("PATHWAY_BAD_REQUEST",
                                        "aggregation は sum / mean / median のどれかです。")
    if source not in ("auto", "arf", "dataset"):
        return mcp_errors.pathway_error("PATHWAY_BAD_REQUEST",
                                        "source は auto / arf / dataset のどれかです。")
    if group_a is None or group_b is None:
        return mcp_errors.pathway_error("PATHWAY_BAD_REQUEST", "group_a と group_b を指定してください。")

    src = _pick_source(source)
    as_set = (lambda g: {str(g)} if isinstance(g, str) else {str(x) for x in g})
    if (as_set(group_a) == as_set(group_b) if src == "arf"
            else str(group_a_label) == str(group_b_label)):
        return mcp_errors.pathway_error(
            "PATHWAY_BAD_REQUEST",
            "group_a と group_b（dataset では group_a_label と group_b_label）が同じです。")
    try:
        built, err = (_build_arf(group_a, group_b) if src == "arf"
                      else _build_dataset(group_a, group_b, group_a_label, group_b_label))
    except handoff.HandoffError as exc:
        return mcp_errors.pathway_error(exc.code, exc.message, exc.details)
    if err:
        return err
    table, normalization, dataset_id, label_a, label_b = built
    # 群名は mzTab-M に書いた形（# → _）で JSON-LD に戻ってくる
    label_a, label_b = handoff.clean_text(label_a), handoff.clean_text(label_b)

    run_id = f"pw_{datetime.now():%Y%m%d_%H%M%S}_{uuid.uuid4().hex[:6]}"
    run_dir = _run_dir(run_id)
    run_dir.mkdir(parents=True, exist_ok=True)
    mztab_path = run_dir / "input.mztab"
    handoff.write_mztab(table, mztab_path, normalization=normalization)
    request = {"mode": "measure", "mztab": str(mztab_path), "out_dir": str(run_dir),
               "pathways": pathways,
               "networks": {p: bridge.PATHWAY_NETWORKS[p] for p in pathways},
               "acyl_chains": list(acyl_chains or []), "aggregation": aggregation,
               "normalization": normalization, "unit": "arbitrary", "dataset_id": run_id}
    try:
        result = bridge.run_bridge(request, run_dir)
        outputs = result.get("outputs") or {}
        absent = [p for p in pathways if not outputs.get(p)]
        if absent:
            raise jsonld.JsonldError(f"bridge-result.json に出力が無いパスウェイ: {', '.join(absent)}")
        graphs = {p: jsonld.load_graph(Path(outputs[p])) for p in pathways}
    except bridge.BridgeError as exc:
        return mcp_errors.pathway_error(exc.code, exc.message, exc.details)
    except jsonld.JsonldError as exc:
        return mcp_errors.pathway_error(exc.code, exc.message)

    # 上流は指定鎖ごとに全ノードへ空スコープを出すので、鎖がデータに無くても失敗しない。
    # 分子種の入った鎖スコープがどのパスウェイにも無い鎖は指定の誤りとして返す（spec §5）
    present = {c for g in graphs.values() for sc in g.scopes.values()
               if sc.scope_type != "subclass" and sc.n_features > 0 for c in sc.chains}
    missing_chains = [c for c in (acyl_chains or []) if c not in present]
    if missing_chains:
        return mcp_errors.pathway_error(
            "PATHWAY_BAD_REQUEST",
            f"指定した鎖がデータ（パスウェイ上のノード）にありません: {', '.join(missing_chains)}",
            {"missing_chains": missing_chains, "run_dir": str(run_dir)})

    per_pathway, all_rows = {}, []
    for p, graph in graphs.items():
        rows = compare.compare_scopes(graph, label_a, label_b, q_threshold)
        pairs, undeterminable = compare.opposite_pairs(graph, rows)
        # 値の付いた（分子種の入った whole スコープを持つ）脂質ノード / 脂質ノード。
        # PathwayEntity（CDP-DAG 等の化合物）と接合ノードは分母に入れない
        lipid_nodes = sum(1 for k in graph.node_kinds.values() if k == "lipid")
        valued = len({r["node_id"] for r in rows if r["scope"] == "whole"})
        per_pathway[p] = {
            "mapped_nodes": f"{valued}/{lipid_nodes}",
            "unmatched_subclasses": graph.unmatched_subclasses,
            "nodes_tsv": compare.rows_to_tsv(rows, _NODE_COLUMNS),
            "opposite_pairs_tsv": compare.rows_to_tsv(pairs, _PAIR_COLUMNS),
            "undeterminable_pairs": undeterminable,
        }
        all_rows += [dict(r, pathway=p) for r in rows]
    (run_dir / "nodes.tsv").write_text(
        compare.rows_to_tsv(all_rows, ["pathway", *_NODE_COLUMNS]), encoding="utf-8")

    groups: dict[str, int] = {}
    for g in map(handoff.clean_text, table.groups):
        groups[g] = groups.get(g, 0) + 1
    manifest = {"run_id": run_id, "schema_version": jsonld.SCHEMA_VERSION, "source": src,
                "created_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
                "group_a": label_a, "group_b": label_b, "groups": groups,
                "pathways": pathways, "acyl_chains": list(acyl_chains or []),
                "aggregation": aggregation, "q_threshold": q_threshold,
                "normalization": normalization, "dataset_id": dataset_id,
                "backend": result.get("backend"),
                "remote_fallback": result.get("remote_fallback"),
                "upstream_commit": result.get("package_sha"),
                "package_version": result.get("package_version"),
                "report": table.report,
                "mapped_nodes": {p: v["mapped_nodes"] for p, v in per_pathway.items()},
                "files": {"mztab": "input.mztab", "nodes": "nodes.tsv",
                          **{p: Path(outputs[p]).name for p in pathways}}}
    (run_dir / "manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=1),
                                           encoding="utf-8")
    s.pathway.last_run_id = run_id
    s.pathway.last_run_dir = str(run_dir)

    return json_payload({
        "run_id": run_id, "source": src, "run_dir": str(run_dir),
        "groups": groups, "group_a": label_a, "group_b": label_b,
        "normalization": normalization,
        "backend": _backend_text(result),
        "handoff": table.report,
        "pathways": per_pathway,
        "notes": {"adjacency": _ADJACENCY_NOTE,
                  "values": "MSLipidMapper は正規化しない。値はこちらの前処理後行列をサブクラスごとに集計したもの。",
                  "figure": "図は pathway_plot(run_id) で表示できる。"},
    })


# run ディレクトリの <pathway>.pdf を同じパスへ上書きする（同じ引数なら同じ結果 → 冪等）。
# 描画は研究室サーバか、ローカルなら unpkg.com から JS を取るのでネットワークに出る。
@mcp.tool(annotations=ToolAnnotations(
    readOnlyHint=False, destructiveHint=False, idempotentHint=True, openWorldHint=True),
    structured_output=False)
def pathway_plot(run_id: str | None = None, pathways: list[str] | None = None,
                 plot_type: str = "dot"):
    """pathway_project の run をパスウェイ図（PDF）に描き、チャットに PNG で表示する。

    - run_id: 省略時は直近の run。
    - pathways: その run に含まれるものから選ぶ（既定は全部）。1 枚約 2,400 画像トークン。
    - plot_type: dot | box | violin（ノード内のグラフ）。横軸の群は group_a, group_b, その他の群の順
      （上流は config の plot.x_order を軸順と配色の両方に使う）。
    [mslipidmapper] api_url の研究室サーバを先に使い（サーバは 3 パスウェイ全部を描く）、
    使えなければローカルの R で描く。ローカルではヘッドレス Chromium（[mslipidmapper] browser）と
    unpkg.com への接続が要る。
    """
    manifest = _load_manifest(run_id)
    if isinstance(manifest, str):
        return manifest
    chosen = list(pathways or manifest["pathways"])
    missing = [p for p in chosen if p not in manifest["pathways"]]
    if missing:
        return mcp_errors.pathway_error(
            "PATHWAY_BAD_REQUEST", f"この run に無いパスウェイ: {', '.join(missing)}",
            {"run_pathways": manifest["pathways"]})
    if plot_type not in ("dot", "box", "violin"):
        return mcp_errors.pathway_error("PATHWAY_BAD_REQUEST",
                                        "plot_type は dot / box / violin のどれかです。")
    run_dir = _run_dir(manifest["run_id"])
    others = [g for g in manifest["groups"] if g not in (manifest["group_a"], manifest["group_b"])]
    request = {"mode": "render", "mztab": str(run_dir / manifest["files"]["mztab"]),
               "out_dir": str(run_dir), "pathways": chosen,
               "networks": {p: bridge.PATHWAY_NETWORKS[p] for p in chosen},
               "aggregation": manifest["aggregation"], "browser": bridge.browser_path(),
               "plot_type": plot_type,
               "x_order": [manifest["group_a"], manifest["group_b"], *others],
               "timeout_seconds": 90}
    try:
        result = bridge.run_bridge(request, run_dir, timeout_s=120.0 * len(chosen))
    except bridge.BridgeError as exc:
        return mcp_errors.pathway_error(exc.code, exc.message, exc.details)

    outputs = result.get("outputs") or {}
    images, lines = [], [f"run {manifest['run_id']}（{manifest['group_a']} 対 {manifest['group_b']}）",
                         f"描画: {_backend_text(result)}"]
    for p in chosen:
        if not outputs.get(p):
            return mcp_errors.pathway_error(
                "PATHWAY_RENDER_FAILED", f"{p} の PDF が出力されませんでした。",
                {"pathway": p, "log": str(run_dir / "bridge.log")})
        pdf = Path(outputs[p])
        try:
            png, w, h = pdf_raster.pdf_first_page_png(pdf)
        except Exception as exc:  # 0 ページ・壊れた PDF など。PyMuPDF の例外型は版で揺れる
            return mcp_errors.pathway_error(
                "PATHWAY_RENDER_FAILED", f"{p} の PDF を画像にできません: {exc}",
                {"pathway": p, "pdf": str(pdf)})
        images.append(Image(data=png, format="png"))
        mapped = (manifest.get("mapped_nodes") or {}).get(p)
        mapped_text = f"、対応ノード {mapped}" if mapped else ""
        lines.append(f"- {p}: {pdf}（{w}×{h}px{mapped_text}）")
    lines.append(_FIGURE_NOTE)
    return ["\n".join(lines), *images]
