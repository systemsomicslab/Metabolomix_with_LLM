"""解析セッションの状態: AnalysisSession とプロセス唯一の `session` シングルトン。

`session` は全ツールが共有する可変状態。ツールからの参照は必ず
`session_state.session`（module 修飾・動的）で行うこと。テストは
`session_state.session = AnalysisSession()` で丸ごと差し替えるため、
`from session_state import session` のようなスナップショット束縛を作ると
差し替えが見えなくなる。

**状態はパーサ単位に分離されている**（`session.arf` / `session.arf2` /
`session.pai2` / `session.eic`）。以前は `features` / `filtered_features` /
`current_file_path` を全パーサが共有していたため、単一サンプルの `pai2_parser` を
1回呼ぶだけで進行中の ARF 解析（前処理行列・差次的結果・手動除外）が無言で消え、
後続の `arf_pca_preprocessed` が「前処理後の行列がありません」で落ちていた。
LLM には破棄の事実が伝わらないため「手動除外が過度」という誤原因が報告される。
粒度の違う成果物（1スポット×1サンプル / 全サンプル統合 / 単一測定ファイル）は
互いに代入不能なので、スロットも分ける。新しいパーサを足すときも専用スロットを
作ること。無修飾の共有スロットへは絶対に書かない。

依存: 各種 reader / preprocessing（いずれも leaf）のみ。tools_* / server は
import しない（循環回避）。`arf_reader` は module オブジェクトのまま参照するので、
テストの `patch.object(server.arf_reader, ...)` が共有 module 経由でここにも効く。
"""
import io
import uuid
import os
import re

from metabolomix.arf import reader as arf_reader
from metabolomix.analysis import preprocessing
from metabolomix.arf2.reader import deserialize
from metabolomix.eic.reader import parse_eic_aef_css1
from metabolomix.pai2.reader import filter_features_by_params
from metabolomix.msdial.classes import (
    assign_sample_groups,
    attach_class_ids_to_spots,
    discover_arf_class_index,
)
from metabolomix.msdial.tags import (
    attach_tags_to_spots,
    discover_arf_tag_index,
    normalize_sample_name,
)

_BATCH_DATE_RE = re.compile(r"(\d{8})")

# output-format リソースを LLM が pull していないときに、解釈直結のパーサー出力の
# 先頭へ最大1回だけ前置する意味論ダイジェスト（自己完結・~7行）。round-trip 不要で
# ローカルLLM でも意味が届く。全文定義は docs/output_format/（共通核 core.md ＋
# パーサ別トピック）/ lipidmix://docs/output-format[/{topic}]。§2.1 は脂質名文法。
#
# **共通部（ここ）とアッセイ種別行（_ASSAY_CAVEAT_LINE）を分ける。** 脂質名文法は
# 脂質アッセイでしか成り立たない規則で、一般代謝物の出力に当てると「PC 34:1 形式で
# ないから未同定」のような誤読を生む。逆に一般代謝物では候補集合（アダクト・異性体・
# 順位）が同定の実体になる。種別は MS-DIAL のファイル形式からは決まらない
# （同じ .arf / .mzTab に両方が出る）ので、確定するまではどちらの規則も当てない。
SEMANTICS_CAVEAT = (
    "[意味論] 解釈前に lipidmix://docs/output-format を参照。要点:\n"
    "- 粒度: ARF行=1スポット×1サンプル / ARF2行=全サンプル統合スポット。\n"
    "- IsGapFilled=true は補間値（実測でない）。\n"
    "- Nameの存在≠確定同定。空/Unknown/no MS2:/low score: を区別。\n"
    "- EIC peak_top は横軸座標(RT)で強度でない。強度はmax_intensity。\n"
    "- PAI2は単一サンプル→PCA不能。多変量比較はARF/ARF2。"
)

#: アッセイ種別。`unknown` は「まだ確定していない」であって「脂質」ではない。
ASSAY_KINDS = ("lipid", "metabolite", "unknown")

#: 種別ごとの追加1行。行数は種別によらず1行（前置の総量を種別で変えない）。
_ASSAY_CAVEAT_LINE = {
    "lipid": "- 脂質名: 34:1(species)と 16:0/18:1(molecular)は別粒度。P-/O-は曖昧(§2.1)。",
    "metabolite": "- 一般代謝物: 名前1件は候補の代表。アダクト/異性体/候補順位を確認。脂質名文法(§2.1)は不適用。",
    "unknown": "- アッセイ種別が未確定。確定するまで脂質名文法(§2.1)も脂質クラス知識も当てない。",
}


def assay_digest(assay_kind: str) -> str:
    """共通部＋その種別の1行を組んだ意味論ダイジェストを返す。"""
    return f"{SEMANTICS_CAVEAT}\n{_ASSAY_CAVEAT_LINE[assay_kind]}"


def _build_sample_meta(sample_names, class_index):
    """サンプル名から role/group/run_order/batch を組み立てる。"""
    class_ids = {}
    orders = {}
    if class_index:
        by_name = {}
        for rec in class_index.get("records", []):
            key = normalize_sample_name(rec.get("file_name"), strip_processing_timestamp=False)
            if key:
                by_name[key] = rec
        for name in sample_names:
            rec = by_name.get(normalize_sample_name(name, strip_processing_timestamp=False))
            if rec:
                class_ids[name] = rec.get("class_id") or ""
                orders[name] = rec.get("analytical_order")
    roles = preprocessing.detect_sample_roles(sample_names, class_ids)
    groups = assign_sample_groups(sample_names, class_index, None)
    meta = {}
    for name in sample_names:
        m = _BATCH_DATE_RE.search(name)
        # バッチは MS-DIAL メタに専用項目が無いため、ファイル名中の8桁日付から推定する。
        # 出所を batch_source に明示し、交絡判定の脆さ（名前依存）を下流で開示できるようにする。
        meta[name] = {
            "role": roles.get(name, "sample"),
            "group": groups.get(name),
            "run_order": orders.get(name),
            "batch": m.group(1) if m else None,
            "batch_source": "filename_date" if m else None,
        }
    return meta


def _new_dataset_id() -> str:
    """ARF データセットの同一性 ID を発行する。"""
    return f"arf_{uuid.uuid4().hex}"


class ArfState:
    """ARF（1スポット×1サンプル）解析の状態。多変量解析の基盤はここだけが持つ。

    ここに載る派生状態（feature_matrix / last_differential / excluded_*）は
    ARF の行構造に強く紐づくので、他パーサのデータで置き換えることはできない。
    """

    def __init__(self):
        self.features = None            # デシリアライズ済みの全スポット
        self.filtered_features = None   # arf_parser の絞り込み後
        self.current_file_path = None
        self.current_tag_directory = None
        self.tag_index = None
        self.class_index = None
        self.pca_result = None          # 直近のPCA計算結果

        # --- P2a 前処理用の正準行列とサンプルメタ ---
        self.feature_matrix = None       # 前処理後のサンプル×特徴量行列
        self.pp_sample_names = None
        self.pp_feature_names = None
        self.sample_meta = {}            # {sample_name: {role, group, run_order, batch}}
        self.preprocessing_recipe = {}   # 直近適用した前処理レシピ（空=未適用）

        # --- 直近解析の描画用データ ---
        self.last_pca_plot = None        # save_pca_figure が参照
        self.last_differential = None    # arf_plot_volcano / save_volcano_figure が参照

        # --- 手動除外集合（PCA 外れサンプル / 特定ピークの可逆・非破壊除外） ---
        self.excluded_samples = set()   # 除外する file_name（サンプル）
        self.excluded_spots = set()     # 除外する MasterAlignmentID（スポット）

        # --- 同一性（図・出力が「どの結果か」を名指しできるようにする） ---
        # dataset_id: 今読んでいる ARF データセットの ID。ファイルを切り替えるたびに
        # 発行し直す（reset_analysis）。DatasetState 側の dataset_id と同じ役割で、
        # 「前のデータの結果を今のデータの図として保存する」取り違えを機械的に防ぐ。
        self.dataset_id = _new_dataset_id()

    def reset_analysis(self):
        """別 ARF データセットへ切り替える際に、前データ由来の解析成果を一括で破棄する。

        新ファイルの load 開始時に必ず呼ぶ。前回の feature_matrix / sample_meta /
        preprocessing_recipe / 差次的解析・PCA プロットなどが残ると、データ切替後に
        前回データの図や結果を「今のデータのもの」として保存・解釈してしまう。
        本メソッドは解析成果（派生状態）だけを消し、これから load される features や
        current_file_path は触らない。他パーサのスロットにも触れない（別データセットへの
        切替は ARF の話であって、PAI2 の在庫要約や EIC のクロマトグラムを消す理由はない）。
        """
        self.pca_result = None
        self.filtered_features = None
        self.last_pca_plot = None
        self.last_differential = None
        # 前処理由来の正準行列とサンプルメタ
        self.feature_matrix = None
        self.pp_sample_names = None
        self.pp_feature_names = None
        self.sample_meta = {}
        self.preprocessing_recipe = {}
        # 手動除外（別データに持ち越さない）
        self.excluded_samples = set()
        self.excluded_spots = set()
        # 別データセットになったので同一性も作り直す。ID を据え置くと、前データの
        # 結果 ID を指定した図の保存が「有効」に見えてしまう。
        self.dataset_id = _new_dataset_id()

    def load_data(self, file_path: str, tag_directory: str | None = None):
        """ファイルパスが前回と異なる場合のみデシリアライズを実行する"""
        if (
            self.current_file_path == file_path
            and self.features is not None
            and self.current_tag_directory == tag_directory
        ):
            if str(file_path).lower().endswith('.arf'):
                self.tag_index = discover_arf_tag_index(
                    file_path, self.features, tag_directory=tag_directory,
                )
                attach_tags_to_spots(self.features, self.tag_index)
                self.class_index = discover_arf_class_index(file_path)
                attach_class_ids_to_spots(self.features, self.class_index)
            return self.features

        # 別データセットへ切り替えるので、前データ由来の解析成果（差次・PCA・前処理行列・
        # 手動除外等）を open 前に一括破棄する。前回の結果を新データのものとして
        # 保存・解釈する取り違えを防ぐ。
        self.reset_analysis()

        with open(file_path, 'rb') as f:
            # 【修正点】ファイルの拡張子を見て正しいパーサーを呼び分ける
            file_ext = str(file_path).lower()
            if file_ext.endswith('.arf'):
                self.features = arf_reader.deserialize(io.BytesIO(f.read()))
                self.tag_index = discover_arf_tag_index(
                    file_path, self.features, tag_directory=tag_directory,
                )
                attach_tags_to_spots(self.features, self.tag_index)
                self.class_index = discover_arf_class_index(file_path)
                attach_class_ids_to_spots(self.features, self.class_index)
            else:
                self.features = deserialize(io.BytesIO(f.read())) # 元からインポートされている arf2_reader 用
                self.tag_index = None
                self.class_index = None

            self.current_file_path = file_path
            self.current_tag_directory = tag_directory

        return self.features


class Arf2State:
    """ARF2（全サンプル統合カタログ）の状態。

    サンプル別強度を持たないため多変量解析には使えない。ARF の解析基盤とは
    粒度が違うので、同じスロットに置かない。
    """

    def __init__(self):
        self.features = None
        self.current_file_path = None

    def load(self, file_path: str, features: list):
        self.features = features
        self.current_file_path = file_path
        return self.features


class Pai2State:
    """PAI2（単一測定ファイルのピーク一覧）の状態。

    単一サンプルなのでサンプル間比較（PCA・差次的解析）は原理的にできない。
    ARF の解析基盤とは独立で、ピーク検証のために往復しても ARF 側は壊れない。
    """

    def __init__(self):
        self.features = None
        self.filtered_features = None
        self.filter_params = {}
        self.current_file_path = None

    def load(self, file_path: str, features: list, filter_params: dict | None = None):
        self.features = features
        self.current_file_path = file_path
        self.apply_filter(filter_params)
        return self.features

    def apply_filter(self, filter_params: dict | None = None):
        """現データに対して動的にフィルタを適用する。"""
        if filter_params is None:
            filter_params = {}
        self.filter_params = filter_params
        if self.features is None:
            self.filtered_features = None
            return None

        self.filtered_features = filter_features_by_params(self.features, filter_params)
        return self.filtered_features


class EicState:
    """EIC/AEF（クロマトグラム）の状態。パスで自己検証するので再読込は安全。"""

    def __init__(self):
        self.features = None
        self.current_file_path = None
        self.last_plot = None  # 直近EICプロット情報（明示的なPNG保存時のみ参照）

    def load_data(self, file_path: str):
        """ファイルパスが前回と異なる場合のみEICデータを解析する"""
        if self.current_file_path == file_path and self.features is not None:
            return self.features

        self.features = parse_eic_aef_css1(file_path, include_chromatogram=False)
        self.current_file_path = file_path
        return self.features


class LibraryState:
    """参照ライブラリのスロット。他スロットとは共有しない。

    `store` は開いたままの `LibraryStore`。`last_match` は直近の
    `library_match_feature` の結果（候補ごとのスペクトルとアラインメントを含む）で、
    `library_plot_mirror` がここから座標を読む。**payload には載せない**。
    """

    def __init__(self):
        self.store = None
        self.source_path: str | None = None
        self.last_match: dict | None = None


class CurationState:
    """キュレーションのスロット。レビューの正準はディスク上の review-<id>.json で、
    ここは review_id からそのフォルダを引くためのポインタだけを持つ。"""

    def __init__(self):
        self.last_review_id: str | None = None
        self.review_dirs: dict[str, str] = {}


class AnalysisSession:
    """パーサ別スロットと、パーサ横断の意味論ガードを束ねる。

    スロット間に暗黙の共有は無い。あるパーサの結果が別パーサの状態を消す必要が
    あるときは、その旨をツールの返り値で明示的に開示すること（無言の破棄が
    「セッションが状態を保持できていない」という誤診につながる）。
    """

    def __init__(self):
        self.arf = ArfState()
        self.arf2 = Arf2State()
        self.pai2 = Pai2State()
        self.eic = EicState()

        # --- 意味論 caveat ガード（output-format 未 pull 時に1回だけ前置） ---
        # プロセス内で真に1回だけ発火させる。output-format リソースが読まれたら
        # output_format_seen=True になり以後は前置しない。データ切替では
        # リセットしない（データ切替のたびに再注入しないため）。
        self.output_format_seen = False
        self.caveat_emitted = False
        # どの種別で前置したか。種別が確定したとき（unknown→lipid/metabolite）に
        # 1回だけ再注入するために持つ。これが無いと、load_dataset 時点の unknown で
        # 1回出て終わり、後から種別が確定しても解釈規則が二度と届かない。
        self.caveat_emitted_kind: str | None = None
        # 既読の output-format トピック（core/arf/eic/…）。トピック別リソースが
        # 読まれるたびに増える。未読トピックのツール出力にだけ誘導1行を足す。
        self.sections_seen: set[str] = set()

        # --- アッセイ種別（解釈規則の選択軸。データ切替では消さない） ---
        # 確定は利用者との合意事項（record_objective / update_objective）で、
        # パーサ出力からは決まらない。既定の unknown は「脂質」ではない。
        self.assay_kind = "unknown"

        # --- mzTab-M / DatasetState スロット（session.arf とは独立） ---
        self.dataset = None  # DatasetState | None

        # --- 参照ライブラリスロット（MS/MS スペクトル照合。他スロットとは独立） ---
        self.library = LibraryState()

        # --- キュレーションスロット（アラインメント注釈の一覧確認・機械判別・フラグ） ---
        self.curation = CurationState()

        # --- Console ジョブスロット ---
        # ジョブ状態の正準はディスク上の analysis-job.json。ここはそのパスへのポインタ。
        # MCP ツールが job_path 省略で呼ばれたときのフォールバックにのみ使う。
        self.current_job_path: str | None = None

    def section_hint(self, topic: str) -> str | None:
        """未読トピックなら、該当セクションを引くよう促す1行を返す。既読なら None。

        core を読んでも各パーサの節を読んだことにはならないのでトピック単位で持つ。
        全文（約700行）を毎回前置する代わりに、1行の誘導＋LLM 側の pull で済ませる。
        """
        if topic in self.sections_seen:
            return None
        # 脂質名文法は脂質アッセイでしか成り立たないので、種別が脂質のときだけ挙げる。
        shared = "共通の粒度・脂質名文法" if self.assay_kind == "lipid" else "共通の粒度"
        return (
            f"[意味論] この出力の定義は `lipidmix://docs/output-format/{topic}` にある。"
            f"解釈前に参照すること（{shared}は `lipidmix://docs/output-format`）。"
        )

    def set_assay_kind(self, assay_kind: str) -> str:
        """アッセイ種別を確定する。正規化した値を返し、未知の値は拒否する。

        黙って `unknown` へ落とさない——落とすと、利用者が種別を伝えたつもりのまま
        サーバは未確定として振る舞い、どちらの規則も当たらない状態が続く。
        """
        normalized = str(assay_kind).strip().lower()
        if normalized not in ASSAY_KINDS:
            raise ValueError(
                f"assay_kind は {'/'.join(ASSAY_KINDS)} のいずれかです: {assay_kind!r}")
        self.assay_kind = normalized
        return normalized

    def maybe_prepend_caveat(self, text: str, topic: str | None = None) -> str:
        """解釈直結パーサー出力へ意味論ダイジェスト（最大1回）と誘導1行を付す。

        - 意味論ダイジェスト: output-format リソースが未 fetch（output_format_seen=False）
          のとき、その時点の `assay_kind` で1回だけ先頭へ前置する。種別が確定して
          変わったときだけもう1回出す（合計で最大2回。unknown のまま解析が進むと
          種別固有の規則が一度も届かないため）。
        - topic 誘導: そのトピックが未読のあいだ、毎回末尾に1行だけ付す（安いので
          既読になるまで出し続ける。これがオンデマンド参照の起点になる）。

        環境変数 LIPIDMIX_CAVEAT_MODE=off で両方とも無効化（ローカルの操作ナビゲータ
        専用デプロイ向け。既定 digest）。エラー文字列など解釈材料でない出力には呼ばない。
        """
        mode = os.getenv("LIPIDMIX_CAVEAT_MODE", "digest").strip().lower()
        if mode == "off":
            return text
        out = text
        if topic:
            hint = self.section_hint(topic)
            if hint:
                out = f"{out}\n\n{hint}"
        if self.output_format_seen:
            return out
        if self.caveat_emitted and self.caveat_emitted_kind == self.assay_kind:
            return out
        self.caveat_emitted = True
        self.caveat_emitted_kind = self.assay_kind
        return f"{assay_digest(self.assay_kind)}\n\n{out}"


# インスタンスを1つ作成（サーバー起動中に保持される）
session = AnalysisSession()
