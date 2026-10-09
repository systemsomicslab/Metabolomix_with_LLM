# キュレーションのビューアから直接送信する（LLM を介さない送信経路） 設計

- 日付: 2026-10-09
- 発端: レビュー・候補付けのビューアの判断は「Copy submission text」→ チャットに貼る →
  LLM が `curation_submit(submission_text=...)` を呼ぶ、で記録している。ユーザーの懸念は
  「送信用 JSON を LLM に通すと文脈を圧迫し、余計な操作の余地を生む」。
- 関連: [アラインメントのキュレーション](2026-09-28-alignment-curation-design.md)、
  [キュレーションの候補付け](2026-09-29-curation-annotation-suggestion-design.md)

## 1. 目的と完成条件

ビューアのボタンから、`flags.jsonl` への記録と `_tags.xml` への反映を**直接**行う。
LLM はこの往復に入らない。

懸念の中身（調査で確認したもの）:

- 送信用テキストはプリセットの Wrong（kidney で数百件）と判定根拠のメモを含み数十 KB になる。
  貼った時点で入力文脈を占め、さらに LLM は**同じ全文を tool 引数として出力し直す**しかない。
- `validate_entries` は形と spot_id しか見ないので、LLM が書き写しで行を落とす・メモを
  言い換えるといった改変は検出されずに記録される。
- 判断はユーザーのもので、LLM は中継しかしていない。同意の取り方としても、本人がボタンを
  押すほうが「チャットに貼った」より明確。

完成条件:

1. `curation_review` / `curation_suggest` が作った HTML をブラウザで開くと、Submit と Finish が出る。
2. Submit で内訳つきの確認を経て、`curation_submit` と同じ記録・`_tags.xml` 反映が行われ、結果が
   ビューアに表示される。何度でも分けて送れる。
3. Finish か時間切れで受け口の登録が外れ、登録が 0 件になれば受け口は止まる。MCP サーバが
   終了すれば受け口も消える。
4. 受け口に届かないとき（時間切れ・サーバ再起動・古い HTML）は、従来の Copy 経路にそのまま戻る。
5. 既存の `curation_submit`（貼り付け・`flags` 直接・`source="llm"`）の挙動は変わらない。

## 2. 合意したユーザー判断（2026-10-09）

- 受け口は `curation_review` / `curation_suggest` のときに起動し、操作が終われば閉じる。
- 「操作の終わり」は **Finish ボタン＋時間切れ**。Submit では閉じない（分けて送れる）。
- 候補付けビューアにも同じ仕組みを入れる。
- 受け口はプロセスに 1 つを**共有**し、レビューごとにトークンを登録する。
- Submit の前に**内訳を出して確認**する。

## 3. 方式の選択

| 案 | 評価 |
|---|---|
| **A. 標準ライブラリ `http.server.ThreadingHTTPServer` を daemon スレッドで動かす** | **採用**。依存の追加なし、FastMCP のイベントループと分離できる |
| B. uvicorn/Starlette を FastMCP の asyncio ループに載せる | sync ツールから動いているループへ安全に立てて止める必要があり、stdio トランスポートと同じループなので詰まったときの影響が大きい |
| C. 別プロセスの補助サーバ | MCP サーバが落ちても残り「サーバと一緒に閉じる」と食い違う。後片付けも要る |

## 4. 構成

### 4.1 `metabolomix/curation/submit_server.py`（新規）

受け口と登録表を持つ。登録は `{token, kind, review_id, arf2_path, last_seen}`。

- `register(kind, review_id, arf2_path) -> {"port", "token", "idle_timeout_min"}`:
  受け口が無ければ `127.0.0.1` の空きポート（`port=0`）で起動し、`secrets.token_urlsafe(32)` を発行。
  `kind` は `"review"` か `"suggest"`。
- `unregister(token)`: 登録を外す。0 件になれば受け口を止める。
- 時間切れ: 最後の通信（ping / submit）から `IDLE_TIMEOUT_MIN = 30` 分で外す。
  掃除は daemon スレッドが 60 秒ごとに行う。時計は差し替え可能にする（テスト用）。
- **トークンに `review_id` と `arf2_path` を登録時に結び付ける**。ブラウザから送るのは flags
  だけで、書き込み先は登録側の値しか使わない（ページ側から差し替えられない）。
- 受け口の停止（`shutdown()`）は、要求を処理しているスレッドではなく別スレッドで、応答を返した後に行う。

### 4.2 送信本体の切り出し

[`curation_submit`](../../../metabolomix/tools/curation_tools.py) の「レビューの解決 → 検証
（`validate_entries` / `expand_entries`）→ 同一 spot_id の重複検出 → アラインメントの sha 照合 →
`FlagStore.append` → `sync_tags` → 件数集計」を 1 つの関数にし、MCP ツールと HTTP の両方が
これを呼ぶ。グローバルな session には依存せず、保存済みレビュー（または `arf2_path` と
`review_id`）を明示で受け取る。

`flags.jsonl` と `_tags.xml` に書く区間はモジュール単位の `threading.Lock` で直列化する
（MCP ツール呼び出しと HTTP スレッドが同時に書きうるため）。

### 4.3 HTML への埋め込み

`save_review` / `save_suggestion` が描画するとき、`{port, token}` を **HTML にだけ**入れる。
レビュー・候補付けの JSON には保存しない。MCP Apps 用の空テンプレート（`render_html(None)`）には
入れない（その経路は従来の `rpc` の Send を使う）。

## 5. HTTP の仕様

待受けは `127.0.0.1` のみ。ページは `http://127.0.0.1:<port>` を直接叩く（`localhost` を使わない）。

| パス | 本文 | 応答 |
|---|---|---|
| `POST /v1/ping` | `{token}` | `{status:"ok", kind, review_id, expires_at}`。期限を延ばす（heartbeat） |
| `POST /v1/submit` | `{token, review_id, flags}` | 送信本体の戻り値（`curation_submit` と同じ形。`tags_xml` を含む）。期限を延ばす |
| `POST /v1/finish` | `{token}` | `{status:"ok"}`。登録を外す |

`source` は `"user"` 固定（ボタンを押したのは本人）。

### 5.1 安全対策

- トークンは `hmac.compare_digest` で比べる。本文の `review_id` が登録と違えば拒否する
  （別の HTML の取り違え）。
- `Host` ヘッダが `127.0.0.1:<port>` でなければ拒否する（DNS rebinding）。
- 本文の上限は 8 MB（`MAX_SPOTS` = 3000 件にメモが付いても収まる）。JSON 以外・未知のパス・
  GET は拒否する（GET を返さないので、他のサイトから中身を覗けない）。
- CORS は `Access-Control-Allow-Origin: *`。クッキーは使わない。権限はトークンだけが持つ。
- **stdout に何も書かない**（stdio の MCP では stdout がプロトコルの通り道）。`http.server` の
  既定のアクセスログ（`log_message`）も止める。

### 5.2 エラー

| 状況 | 応答 | ページの扱い |
|---|---|---|
| トークンが不明・期限切れ | 401 | 送信用テキストを入れた Copy 欄を自動で開く。未送信の変更は失わない |
| `review_id` 不一致・`Host` 不正 | 403 | メッセージを表示 |
| 本文が上限超過 | 413 | メッセージを表示 |
| 検証エラー（不正な flag 等）・JSON でない | 400 | メッセージを表示し未送信の変更を残す |
| レビューの後でアラインメントが変わった | 409 | 同上（`curation_review` のやり直しを案内） |
| `_tags.xml` への反映だけ失敗 | 200（従来どおり記録は残し `tags_xml.error`） | 結果と一緒に失敗を表示 |

エラー本文は既存の `_error` と同じ `{status:"error", message, ...}`。

### 5.3 確かめておく危険

Chrome / Edge の Local Network Access が `file://` のページから `127.0.0.1` への要求を止めないか。
**実装計画の最初の工程で Edge と Chrome に当てて確かめる**。止められる場合は、許可の
ダイアログを受け入れるか、受け口が HTML も配る形（`http://127.0.0.1:<port>/view/<token>`）に
切り替えるかを、その時点でユーザーに諮る。`127.0.0.1` だけで待つので Windows ファイアウォールの
確認は出ない想定だが、同じ工程で確かめる。

## 6. ビューアの UI

送信クライアントは 2 つのビューアが共有する `viewer_common.js` に置く。

- 開いたときに ping。応答があればヘッダに **Submit**（主ボタン）と **Finish** を出し、
  **Copy submission text** は残して控えめにする。応答が無ければ従来どおり Copy だけ。
- **Submit** → 確認ダイアログ（`confirm()`）:
  - レビュー: `Wrong n / Suspect n / Confirmed n / clear n を記録し、_tags.xml を更新します。
    MS-DIAL でこのプロジェクトを開いているなら先に閉じてください。`
  - 候補付け: `assign n / redundant n / clear n`。`assign` / `redundant` は `_tags.xml` を変えないが、
    `clear`（元の注釈に戻す）は Misannotation と Confirmed を外すので、**`clear` が 1 件以上
    あるときだけ** MS-DIAL の注意を出す。
- 成功したら各カードを確定状態に更新し（MCP Apps の Send と同じ処理を共通関数にする）、未送信の
  変更を空にする。ステータスに記録件数、`_tags.xml` の変化（Misannotation の付与・除去、
  Confirmed の変化）、`tags_xml.note` を出す。何度送ってもよい。
- **Finish** → 未送信の変更があれば「n 件の未送信の変更があります。終了してよいですか
  （後から Copy で送れます）」と確認。終了後は Submit / Finish を消し Copy だけにする。
- 5 分おきと、タブが表に戻ったとき（`visibilitychange`）に heartbeat を送る。401 なら Copy の表示へ
  切り替え、理由（期限切れ・サーバ再起動）を出す。
- MCP Apps 経路（iframe 内の Send）はそのまま残す。
- 候補付けビューアは現状日本語なので、足す文言も日本語にそろえる（英語化は task.md の別項目）。

## 7. LLM への案内

- `curation_review` / `curation_suggest` の戻り値に `submit: {"via": "viewer", "idle_timeout_min": 30}` を足す。
- 両ツールの docstring を「ユーザーにはビューアの Submit で送ってもらう。送信用テキストを貼るよう
  頼まない。送った結果は `curation_flags` で確かめられる。Copy のテキストが貼られたときは
  従来どおり `curation_submit` に渡す」に改める。
- `curation_submit` は残す（貼り付け・LLM の提案にユーザーが同意した `source="llm"`）。
  「同意なしに呼ばない」も維持する。
- 同じ作業の中で直す文書: `MCP_INSTRUCTIONS`（curation の記述があれば）、`USAGE.md`、
  `docs/output_format/curation.md`、`docs/workflow/curation.md`（新モジュールと送信本体の関数）、
  vault の流れ図（利用者から見た送信の手順が変わるため）。

## 8. 入れないもの

- 受け口を無効にする設定、期限を変える引数・環境変数。
- 閉じた後に受け口を開き直すツール（閉じた後は Copy か `curation_review` のやり直し）。
- streamable-http モードでの同一ポート相乗り（stdio と同じく別ポートで動く）。

## 9. テスト

### 9.1 受け口と登録表（`tests/test_curation_submit_server.py`、新規）

- 登録で起動、2 件目の登録で同じポートを共有、最後の `unregister` で停止（接続できなくなる）。
- 時間切れ: 時計を差し替え、heartbeat で延び、来なければ外れて停止する。
- 実際の受け口を空きポートで立て `http.client` で叩く:
  - ping / submit / finish の正常系。submit 後に `flags.jsonl` と `_tags.xml` が変わる
    （fixture はテスト自身が tmp に作る）。
  - 拒否: トークン不明・期限切れ（401）、`review_id` 不一致・`Host` 不正（403）、上限超過（413）、
    JSON でない（400）、未知のパス（404）、GET（405）。
  - finish の後の ping が 401。
- 同時に 2 本送っても両方の行が記録され `_tags.xml` が壊れない（ロック）。
- 要求の処理中に stdout へ 1 文字も出ない。

### 9.2 既存経路の回帰（`tests/test_curation_tools.py`）

- 送信本体の切り出し後も `curation_submit` の既存テストが全て緑。
- `curation_review` / `curation_suggest` の戻り値に `submit` があり、HTML にポートとトークンが
  埋め込まれ、レビュー・候補付けの JSON にはトークンが入らない。

### 9.3 ビューアの純関数（既存の node 実行テストの流儀）

- 確認文の内訳（レビュー用・候補付け用。候補付けで `clear` があるときだけ MS-DIAL の注意）。
- 送信成功後のカード状態の更新（MCP Apps の Send と共通の関数）。
- ping の結果による表示の切り替え。

### 9.4 実地確認（自動化しない）

- 最初の工程: Edge / Chrome の `file://` から `127.0.0.1` へ届くか（§5.3）。
- 最後: kidney の実データで Submit → `_tags.xml` の変化 → MS-DIAL で開き直して見える、Finish、
  時間切れ後に Copy へ切り替わる、を通しで確かめる。`curation/` の既存の記録は消さない。
