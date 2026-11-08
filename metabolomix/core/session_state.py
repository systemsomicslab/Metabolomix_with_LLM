"""解析セッションの状態: AnalysisSession とプロセス唯一の `session` シングルトン。

`session` は全ツールが共有する可変状態。ツールからの参照は必ず
`session_state.session`（module 修飾・動的）で行うこと。テストは
`session_state.session = AnalysisSession()` で丸ごと差し替えるため、
`from session_state import session` のようなスナップショット束縛を作ると
差し替えが見えなくなる。

**パーサ別のスロット（旧 `session.arf` / `arf2` / `pai2` / `eic`）は無い**（B-2b で消した）。
読み込んだものはすべて `DatasetState` としてハンドルの登録簿に住む。粒度の違う成果物を
共有の無修飾スロットへ書くと、別の読み込みが進行中の解析を無言で消す——新しい入口を足すときも
共有スロットを作らない（`tests/test_session_namespaces.py` が固定している）。

**データセットとライブラリはハンドルの登録簿に住む**（spec 2026-11-06 §5）。
`session.dataset` / `session.library` は直近のハンドル（`current_handle_id`）を指す property で、
`session.dataset = ds` は `ds` を新しいハンドルとして登録簿に入れて直近にする。閉じた直近の
ハンドルは `None`（`missing_state`）になり、開いている別のハンドルへ黙って落ちない。

依存: `user_config`（stdlib だけの leaf）のみ。tools_* / server は import しない（循環回避）。
`core.handles` は property の中で遅延 import する（handles がこの module を import するため）。
"""
import contextlib
import contextvars
import threading

from metabolomix.core import user_config

# output-format リソースを LLM が pull していないときに、解釈の材料になる出力へ最大 1 回
# （種別が確定したらもう 1 回）届ける意味論ダイジェスト（自己完結・~6 行）。入口の
# `load_dataset` が JSON の `caveat` として運び、`record_objective(assay_kind=...)` が確定後の
# 1 回を前置する。round-trip 不要でローカル LLM でも意味が届く。全文定義は docs/output_format/
# （共通核 core.md ＋ トピック別）/ lipidmix://docs/output-format[/{topic}]。§2.1 は脂質名文法。
# 行は `DatasetState` の形（特徴×試料）で語る——旧パーサ別の粒度（ARF / ARF2 / PAI2）は無い。
#
# **共通部（ここ）とアッセイ種別行（_ASSAY_CAVEAT_LINE）を分ける。** 脂質名文法は
# 脂質アッセイでしか成り立たない規則で、一般代謝物の出力に当てると「PC 34:1 形式で
# ないから未同定」のような誤読を生む。逆に一般代謝物では候補集合（アダクト・異性体・
# 順位）が同定の実体になる。種別は MS-DIAL のファイル形式からは決まらない
# （同じ .arf / .mzTab に両方が出る）ので、確定するまではどちらの規則も当てない。
SEMANTICS_CAVEAT = (
    "[意味論] 解釈前に lipidmix://docs/output-format を参照。要点:\n"
    "- 粒度: データセットの 1 行 = 特徴×試料（特徴＝アラインメントスポット）。"
    "名前・m/z・RT は特徴ごとに 1 つで、試料別の値ではない。\n"
    "- detected_mask が False のセルは gap-fill（補間値）で、実測のピークではない。"
    "検出数・存在判定には使わない。\n"
    "- Nameの存在≠確定同定。空/Unknown/no MS2:/low score: を区別。\n"
    "- MS/MS の実スペクトルは .dcl にある。has_msms 系の印は取得の参照があるだけで、"
    "スペクトルの裏付けではない。\n"
    "- EIC peak_top は横軸座標(RT)で強度でない。強度は max_intensity。"
)

#: アッセイ種別。`unknown` は「まだ確定していない」であって「脂質」ではない。
ASSAY_KINDS = ("lipid", "metabolite", "unknown")


def normalize_assay_kind(assay_kind: str) -> str:
    """アッセイ種別を正規化する（状態は変えない）。未知の値は ValueError。"""
    normalized = str(assay_kind).strip().lower()
    if normalized not in ASSAY_KINDS:
        raise ValueError(
            f"assay_kind は {'/'.join(ASSAY_KINDS)} のいずれかです: {assay_kind!r}")
    return normalized

#: 種別ごとの追加1行。行数は種別によらず1行（前置の総量を種別で変えない）。
_ASSAY_CAVEAT_LINE = {
    "lipid": "- 脂質名: 34:1(species)と 16:0/18:1(molecular)は別粒度。P-/O-は曖昧(§2.1)。",
    "metabolite": "- 一般代謝物: 名前1件は候補の代表。アダクト/異性体/候補順位を確認。脂質名文法(§2.1)は不適用。",
    "unknown": "- アッセイ種別が未確定。確定するまで脂質名文法(§2.1)も脂質クラス知識も当てない。",
}


def assay_digest(assay_kind: str) -> str:
    """共通部＋その種別の1行を組んだ意味論ダイジェストを返す。"""
    return f"{SEMANTICS_CAVEAT}\n{_ASSAY_CAVEAT_LINE[assay_kind]}"


#: `dataset_handle` を省いた呼び出しが、ロックを取ったハンドルの ID（`core.handle_dispatch` が
#: ワーカースレッドの中で置く）。置かれている間、`session.dataset` / `session.library` は
#: 直近のハンドルでなくこのハンドルを読む——呼び出しの途中で `load_dataset` が直近を
#: 切り替えても、ロックを取っていない別のハンドルの状態を触らないため。
_PINNED_HANDLE_ID: contextvars.ContextVar[str | None] = contextvars.ContextVar(
    "pinned_handle_id", default=None)


@contextlib.contextmanager
def pinned_handle(handle_id: str):
    """この文脈の中の `session.dataset` / `session.library` を `handle_id` のハンドルに固定する。

    そのハンドルが閉じられていれば `dataset` は `None`（直近や別のハンドルへは落ちない）。
    """
    token = _PINNED_HANDLE_ID.set(handle_id)
    try:
        yield
    finally:
        _PINNED_HANDLE_ID.reset(token)


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


class PathwayState:
    """パスウェイ投影（pathway_project / pathway_plot）の直近 run へのポインタ。

    run の正準は `reports/pathway/<run_id>/manifest.json`。ここは run_id 省略時の
    既定を決めるためだけにある。他パーサのスロットに触らない。
    """

    def __init__(self):
        self.last_run_id: str | None = None
        self.last_run_dir: str | None = None


class AnalysisSession:
    """パーサ別スロットと、パーサ横断の意味論ガードを束ねる。

    スロット間に暗黙の共有は無い。あるパーサの結果が別パーサの状態を消す必要が
    あるときは、その旨をツールの返り値で明示的に開示すること（無言の破棄が
    「セッションが状態を保持できていない」という誤診につながる）。
    """

    def __init__(self):
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
        # `load_dataset` はワーカースレッドで走り、同時に呼ばれうる。出したかの判定と記録を
        # 1 つの区間にして、ダイジェストが 2 回出ないようにする。
        self._caveat_lock = threading.Lock()

        # --- アッセイ種別（解釈規則の選択軸。データ切替では消さない） ---
        # 確定は利用者との合意事項（record_objective）で、
        # パーサ出力からは決まらない。既定の unknown は「脂質」ではない。
        self.assay_kind = "unknown"

        # --- DatasetState（ハンドルの登録簿が持つ。session はその直近を指す） ---
        # データセットはハンドルの登録簿（core.handles.registry）に住む。ここが持つのは
        # 直近のハンドルの ID だけで、`dataset` / `library` はそれを読む property（下）。
        self.current_handle_id: str | None = None
        # 直近のハンドルの切り替えと、セッション自身のライブラリの移動・書き込みを守る。
        # 再入できる必要がある（ロックの下で `library` / `_handle` を読み、`_handle` が
        # 閉じられた直近を外すために同じロックを取り直す）。
        self._current_lock = threading.RLock()

        # --- 参照ライブラリ（MS/MS スペクトル照合。他スロットとは独立） ---
        # 直近のハンドルが無いときに `library` が指すセッション自身の LibraryState。
        self._library = LibraryState()

        # --- キュレーションスロット（アラインメント注釈の一覧確認・機械判別・フラグ） ---
        self.curation = CurationState()

        # --- パスウェイ投影スロット（MSLipidMapper。他スロットとは独立） ---
        self.pathway = PathwayState()

        # --- Console ジョブスロット ---
        # ジョブ状態の正準はディスク上の analysis-job.json。ここはそのパスへのポインタ。
        # MCP ツールが job_path 省略で呼ばれたときのフォールバックにのみ使う。
        self.current_job_path: str | None = None

    # ---------- データセットとライブラリ（直近のハンドルを指す） ----------

    def _handle(self):
        """`dataset` / `library` を読む先のハンドル（無ければ None）。

        `pinned_handle` の中ならそのハンドル、そうでなければ直近のハンドル。直近のハンドルが
        閉じられていたら `current_handle_id` を None に戻す——開いている別のハンドルへは
        落ちない。`core.handles.registry` は module 修飾で読む（テストが差し替える）。
        """
        from metabolomix.core import handles

        pinned = _PINNED_HANDLE_ID.get()
        if pinned is not None:
            return handles.registry.peek(pinned)
        return self._current_handle()

    def _current_handle(self):
        """固定を見ずに、直近のハンドル（無ければ None。閉じられていれば ID を None に戻す）。"""
        from metabolomix.core import handles

        handle_id = self.current_handle_id
        if handle_id is None:
            return None
        handle = handles.registry.peek(handle_id)
        if handle is None:
            with self._current_lock:
                if self.current_handle_id == handle_id:  # 別のスレッドが置いた直近は消さない
                    self.current_handle_id = None
        return handle

    @property
    def arf2_path(self):
        """直近のハンドルの `.arf2`（`DatasetHandle.arf2_path`。無い・ハンドルが無ければ None）。"""
        handle = self._handle()
        return handle.arf2_path if handle is not None else None

    @property
    def dcl_path(self):
        """直近のハンドルの `.dcl`（`DatasetHandle.dcl_path`。無い・ハンドルが無ければ None）。"""
        handle = self._handle()
        return handle.dcl_path if handle is not None else None

    @property
    def dataset(self):
        """直近のハンドルの DatasetState（無い・閉じられていれば None）。"""
        handle = self._handle()
        return handle.dataset if handle is not None else None

    @dataset.setter
    def dataset(self, value) -> None:
        """`None` なら直近のハンドルを外す（ハンドルは閉じない）。DatasetState なら新しい
        ハンドルとして登録簿に入れ、直近にする。

        既に直近のハンドルのデータセットと同じ実体なら何もしない（代入し直すたびにハンドルを
        増やして上限に当たらないため）。上限なら `DATASET_HANDLE_LIMIT`（DomainError）を投げ、
        直近もライブラリも動かさない——古いハンドルを黙って閉じると、閉じられた側の呼び出しが
        別のデータを読むため。

        直近のハンドルが無いあいだにセッションが持っていたライブラリ（store つき。
        `library_load` → `load_dataset` の順）は、新しいハンドルへ**移す**（写さない——2 か所が
        同じ store を閉じないため）。セッション自身には空の LibraryState が残るので、後で
        このハンドルを閉じても古いライブラリへ黙って落ちない。
        """
        if value is None:
            with self._current_lock:
                self.current_handle_id = None
            return
        current = self._handle()
        if current is not None and current.dataset is value:
            return
        from metabolomix.core import handles
        from metabolomix.mztab.polarity import dataset_polarity

        polarity, polarity_source, _notes = dataset_polarity(value)
        # 移すかの決定・ハンドルの登録・入れ替えは 1 つのロックの下で行う。`library_load` が
        # セッション自身のライブラリへ書く区間（`library_write_lock`）と排他にして、移す決定の
        # あとに書かれたストアが、見えない側のライブラリに取り残されないようにする。
        # ロックの順は `_current_lock` → 登録簿のロック（登録簿は session のロックを取らない）。
        with self._current_lock:
            current = self._handle()
            if current is not None and current.dataset is value:
                return
            own = self._library
            move = current is None and _PINNED_HANDLE_ID.get() is None and own.store is not None
            handle = handles.DatasetHandle(value, polarity=polarity,
                                           polarity_source=polarity_source,
                                           source={"kind": "session"},
                                           library=own if move else None)
            handles.registry.add(handle)  # 上限ならここで止まり、ライブラリはセッションに残る
            if move:
                self._library = LibraryState()
            self.current_handle_id = handle.handle_id

    def adopt_handle(self, handle, *, move_library: bool = False) -> None:
        """出来上がったハンドルを登録簿に入れ、**直近にする**（`load_dataset` が使う）。

        上限なら `DATASET_HANDLE_LIMIT`（DomainError）を投げ、直近もライブラリも動かさない。
        `move_library` が真で、ハンドルが自分のライブラリ（store）を持たず、直近のハンドルが
        無いあいだにセッションが store つきのライブラリを持っていたなら（`library_load` →
        `load_dataset` の順）、そのライブラリを新しいハンドルへ**移す**（`dataset` の setter と
        同じ規則）。登録・移動・直近の切り替えは 1 つのロックの下で行う。
        """
        from metabolomix.core import handles

        with self._current_lock:
            own = self._library
            move = (move_library and handle.library.store is None
                    and self._current_handle() is None
                    and _PINNED_HANDLE_ID.get() is None and own.store is not None)
            if move:
                handle.library = own
            try:
                handles.registry.add(handle)  # 上限ならここで止まり、ライブラリはセッションに残る
            except Exception:
                if move:
                    handle.library = LibraryState()
                raise
            if move:
                self._library = LibraryState()
            self.current_handle_id = handle.handle_id

    @property
    def library_write_lock(self) -> threading.RLock:
        """セッション自身のライブラリへ書く区間を守るロック（`dataset` の setter の移動と排他）。

        `library_load` が `dataset_handle` なしでストアを差し替える間、これを保持する。
        ロックの下で `library` を読み直せば、移動が済んだあとは新しいハンドルのライブラリへ
        書く。再入できる。
        """
        return self._current_lock

    @property
    def library(self) -> LibraryState:
        """直近のハンドルがあればその `library`、無ければセッション自身の LibraryState。

        `pinned_handle` の固定先が閉じられていたら、どこにも属さない空の LibraryState を返す
        （セッション自身の、古いかもしれないライブラリへは落ちない）。
        """
        handle = self._handle()
        if handle is not None:
            return handle.library
        if _PINNED_HANDLE_ID.get() is not None:
            return LibraryState()
        return self._library

    @library.setter
    def library(self, value: LibraryState) -> None:
        handle = self._handle()
        if handle is not None:
            handle.library = value
        elif _PINNED_HANDLE_ID.get() is None:
            self._library = value
        # 固定先が閉じられていれば置き先が無い（セッション自身へは書かない）。

    def set_assay_kind(self, assay_kind: str) -> str:
        """アッセイ種別を確定する。正規化した値を返し、未知の値は拒否する。

        黙って `unknown` へ落とさない——落とすと、利用者が種別を伝えたつもりのまま
        サーバは未確定として振る舞い、どちらの規則も当たらない状態が続く。
        """
        normalized = normalize_assay_kind(assay_kind)
        self.assay_kind = normalized
        return normalized

    def take_caveat(self) -> str | None:
        """いま届けるべき意味論ダイジェストを返す（無ければ None）。出したことを記録する。

        規則: `[server] caveat_mode = 'off'` で止める・output-format の核が読まれていれば止める・
        同じ種別では 1 回だけ・種別が変わったらもう 1 回（unknown のまま解析が進むと種別固有の
        規則が一度も届かないため）。JSON を返す入口（`load_dataset`）が最上位の `caveat` として載せ、
        `record_objective(assay_kind=...)` が確定後の 1 回を戻り値の先頭へ前置する。設定ファイルが
        読めない（`caveat_mode` の値が digest / off 以外、旧 `[msdial] lbm` が残っている、TOML の
        構文など）ときは**送出しない**: 既定の digest で続け、問題を 1 行だけ添える
        （状態を変えたあとで落ちないため。黙って既定に落とすのでもない）。
        """
        try:
            setting = user_config.get_setting("server.caveat_mode")
        except user_config.ConfigInvalidError as exc:
            digest = self._take_digest()
            note = _config_problem_note(exc)
            return f"{digest}\n\n{note}" if digest else note
        if setting is not None and setting.value == "off":
            return None
        return self._take_digest()

    def _take_digest(self) -> str | None:
        """ダイジェスト本体。出す条件を満たせば記録して返し、満たさなければ None。"""
        with self._caveat_lock:
            if self.output_format_seen:
                return None
            if self.caveat_emitted and self.caveat_emitted_kind == self.assay_kind:
                return None
            self.caveat_emitted = True
            self.caveat_emitted_kind = self.assay_kind
            return assay_digest(self.assay_kind)


def _config_problem_note(exc: "user_config.ConfigInvalidError") -> str:
    """設定ファイルが読めないことを出力の末尾に添える 1 行。"""
    name = exc.config_file.replace("\\", "/").rsplit("/", 1)[-1]
    return (f"⚠ 設定ファイル {name} が読めないので [server] caveat_mode を既定の digest として"
            f"扱った（CONFIG_INVALID: {exc.message}）。")


# インスタンスを1つ作成（サーバー起動中に保持される）
session = AnalysisSession()
