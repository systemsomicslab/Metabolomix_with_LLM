"""PAI2（1測定ファイルのピーク一覧）ツール群と単一ピーク検証。

pai2_parser（在庫要約）, pai2_inspect_peak, verify_peak_annotation。
PAI2 は単一サンプルなのでサンプル間比較（オミクス PCA）は行わない（複数サンプルの
多変量比較は ARF/ARF2 を使う）。deps: mcp_core / session_state / path_resolvers /
tool_helpers / pai2_reader / knowledge_store。tools_* / server は import しない。
"""
import io
from pathlib import Path

from metabolomix.corpus import knowledge_store
from metabolomix.core import mcp_core
from metabolomix.core import mcp_errors
from metabolomix.core import session_state
from mcp.types import ToolAnnotations
from metabolomix.core.mcp_core import mcp
from metabolomix.core.serialization import json_payload
from metabolomix.core.path_resolvers import resolve_pai2_file_path
from metabolomix.core.tool_helpers import _build_verification_dossier
from metabolomix.pai2.reader import inspect_peak_details, summarize_pai2_inventory

__all__ = [
    "pai2_parser",
    "pai2_inspect_peak",
    "verify_peak_annotation",
]


def _attach_sibling_msms(pai2_path: str, features: list[dict]) -> dict:
    """同名 `.dcl` の MS/MS を PAI2 ピークへ索引対応で付与する（ベストエフォート）。

    PAI2 の `has_msms` は取得参照の有無を示すだけで実スペクトルを持たない。実体は
    同名 `.dcl`（MSDecResult）にあり、`dcl_index` が PAI2 のピーク順に対応する。
    ここで付けておくと `verify_peak_annotation` が実フラグメントを根拠にできる。

    .dcl が無い/壊れている場合でも PAI2 の解析自体は続行させる（MS/MS は付加情報で
    あり、欠けても在庫要約は成立するため）。状況は caveat として返し、無言で
    「MS/MS 無し」と誤認させない。
    """
    from metabolomix.dcl.reader import attach_msms_to_features, deserialize_dcl, find_dcl_for_pai2

    dcl_path = find_dcl_for_pai2(pai2_path)
    if not dcl_path:
        return {"attached": 0, "dcl_file": None, "caveat": (
            "同名の .dcl が隣接していないため MS/MS を付与できませんでした。"
            "has_msms フラグだけでは実スペクトルの有無を判定できない点に注意。")}
    try:
        results = deserialize_dcl(dcl_path, include_spectrum=True, top_n_peaks=10)
        attached = attach_msms_to_features(features, results)
    except Exception as exc:  # noqa: BLE001 - MS/MS は付加情報。PAI2 解析は続行させる
        return {"attached": 0, "dcl_file": Path(dcl_path).name, "caveat": (
            f"隣接 .dcl の読み込みに失敗したため MS/MS は未付与です: {exc}")}
    report = {"attached": attached, "dcl_file": Path(dcl_path).name, "caveat": None}
    if attached == 0:
        report["caveat"] = (
            f"{Path(dcl_path).name} を読みましたが、MS/MS を付与できたピークは0件でした"
            "（precursor m/z が一致せず索引対応が崩れている可能性）。")
    return report


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=True), structured_output=False)
def pai2_parser(file_path: str, filter_threshold: float | None = None) -> str:
    """1つの .pai2（単一測定ファイル）を解析し、ピーク在庫の要約を返します。

    返すのは注釈状況・m/z・RT・強度・S/N の分布と、強度上位ピーク（生化学的に意味のある
    ランキング）です。PAI2 は単一サンプルなので、サンプル間比較（オミクス PCA）はこの単位
    では行えません（複数サンプルの多変量比較は ARF/ARF2 を使う）。個々のピークは
    pai2_inspect_peak / verify_peak_annotation で深掘りできます。MS/MS は
    同名 .dcl（dcl_index がリスト順に対応）を参照します。
    """
    file_path = resolve_pai2_file_path(file_path)
    if not file_path:
        return "データディレクトリに .pai2 ファイルが見つかりませんでした。"

    if filter_threshold is None:
        filter_threshold = 0.0

    from metabolomix.pai2.reader import deserialize

    try:
        with open(file_path, 'rb') as f:
            packed_data = f.read()

        file_like_object = io.BytesIO(packed_data)
        deserialized_and_formatted_data = deserialize(file_like_object)

        assert isinstance(deserialized_and_formatted_data, list)
        assert len(deserialized_and_formatted_data) > 0
        assert isinstance(deserialized_and_formatted_data[0], dict)

        # PAI2 専用スロットへ載せる。ARF（多サンプル）の解析基盤とは粒度が違い
        # 互いに代入不能なので、ARF の前処理行列・差次的結果・手動除外には触れない。
        # 以前はここで共有 features を上書きし reset_analysis_state() まで呼んでいたため、
        # ピークを1つ覗くだけで進行中の ARF 解析が無言で消えていた。
        session_state.session.pai2.load(
            file_path,
            deserialized_and_formatted_data,
            {"min_intensity": filter_threshold},
        )

        # 同名 .dcl の MS/MS を付けておく。verify_peak_annotation が「取得フラグ」ではなく
        # 実スペクトルを根拠に同定確度を語れるようにするため（欠けても解析は続行）。
        msms_report = _attach_sibling_msms(file_path, session_state.session.pai2.features)

        summary = summarize_pai2_inventory(session_state.session.pai2.filtered_features)
        summary["msms_attachment"] = msms_report
        text_report = (
            f"### PAI2 解析完了: {Path(file_path).name}\n"
            + json_payload(summary)
        )
        return session_state.session.maybe_prepend_caveat(text_report, topic="pai2")

    except Exception as e:
        return f"[ERROR] PAI2 解析に失敗しました: {str(e)}"


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=True), structured_output=False)
def pai2_inspect_peak(peak_id: str | None = None, peak_name: str | None = None) -> str:
    """特定のピークについて、強度・S/N・MS/MS相当の情報を返す。

    返り値には signal_to_noise フィールドが含まれます。peak_id か peak_name の
    いずれかを指定する。
    """
    if session_state.session.pai2.filtered_features is None:
        return mcp_errors.missing_state(
            "pai2_dataset", ["pai2_parser"],
            "先に pai2_parser を実行してデータを読み込んでください。")

    details = inspect_peak_details(
        session_state.session.pai2.filtered_features,
        peak_id=peak_id,
        peak_name=peak_name,
    )
    return json_payload(details)


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=True), structured_output=False)
def verify_peak_annotation(
    peak_id: str | None = None, peak_name: str | None = None
) -> str:
    """指定した1ピークのアノテーションが生化学的に妥当かを検証するドシエを返す。

    分析化学的な同定確度（精密質量誤差ppm・アダクト/イオンモード整合）を決定的に
    判定し、生物学的妥当性は関連 knowledge slug を添えて LLM の判断に委ねる。
    先に pai2_parser でデータを読み込むこと。peak_id か peak_name の
    いずれかを指定する。
    """
    if session_state.session.pai2.filtered_features is None:
        return mcp_errors.missing_state(
            "pai2_dataset", ["pai2_parser"],
            "先に pai2_parser を実行してデータを読み込んでください。")
    if peak_id is None and peak_name is None:
        return json_payload({"status": "error", "message": "peak_id か peak_name のいずれかを指定してください。"})

    matches = []
    for feat in session_state.session.pai2.filtered_features:
        if peak_id is not None and str(feat.get("id")) == str(peak_id):
            matches.append(feat)
        elif (
            peak_name is not None
            and isinstance(feat.get("name"), str)
            and peak_name.lower() in feat.get("name", "").lower()
        ):
            matches.append(feat)

    if not matches:
        return json_payload({"status": "not_found", "message": "指定された代謝物がフィルタ済みデータ内に見つかりませんでした。"})

    vocab = knowledge_store.load_vocab(mcp_core.KNOWLEDGE_DIR)
    dossiers = [_build_verification_dossier(feat, vocab) for feat in matches]
    payload = dossiers[0] if len(dossiers) == 1 else {"status": "success", "matches": dossiers}
    return json_payload(payload)
