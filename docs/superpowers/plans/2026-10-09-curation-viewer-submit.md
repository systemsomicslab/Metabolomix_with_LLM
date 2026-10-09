# キュレーションのビューアから直接送信する 実装計画

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** `curation_review` / `curation_suggest` のビューアの Submit ボタンから、LLM を介さずに `flags.jsonl` への記録と `_tags.xml` への反映を行う。

**Architecture:** 送信本体（検証 → 重複検出 → sha 照合 → 追記 → `_tags.xml` 反映）を `metabolomix/curation/submission.py` に切り出し、MCP ツール `curation_submit` と、新設の `metabolomix/curation/submit_server.py`（標準ライブラリ `ThreadingHTTPServer` を daemon スレッドで `127.0.0.1` の空きポートに立てる受け口）の両方がこれを呼ぶ。受け口はプロセスに 1 つを共有し、レビューごとにトークンを登録する。トークンとポートは HTML にだけ埋め込み、ビューア共通の `viewer_common.js` が ping / submit / finish を叩く。

**Tech Stack:** Python 3.14 標準ライブラリ（`http.server` `threading` `secrets` `hmac`）、素の JS（`fetch`）、pytest、node（ビューアの純関数テスト）。

**Spec:** [docs/superpowers/specs/2026-10-09-curation-viewer-submit-design.md](../specs/2026-10-09-curation-viewer-submit-design.md)（executor は spec と本計画の両方を読む）

## Global Constraints

- 作業はすべて worktree `C:\Users\yuu18\Metabolomix_with_LLM\.worktrees\feat-curation-viewer-submit`（ブランチ `feat/curation-viewer-submit`）で行う。`git stash` を使わない。
- Python は `C:/Python314/python.exe`。テストはリポジトリルート（＝worktree ルート）から `C:/Python314/python.exe -m pytest ...`。
- 依存を追加しない（標準ライブラリだけで受け口を作る。spec §3 案 A）。
- 待受けは `127.0.0.1` のみ。ページは `http://127.0.0.1:<port>` を叩く（`localhost` を使わない）。
- 受け口のコードは **stdout に何も書かない**（stdio の MCP では stdout がプロトコルの通り道）。`log_message` を止める。
- `IDLE_TIMEOUT_MIN = 30`、掃除は 60 秒ごと、本文上限 8 MB、トークンは `secrets.token_urlsafe(32)`、比較は `hmac.compare_digest`。
- HTTP 経路の `source` は `"user"` 固定。書き込み先（`review_id` / `arf2_path`）は登録側の値だけを使う。
- `{port, token}` は HTML にだけ入れる。レビュー・候補付けの JSON に保存しない。`render_html(None)`（MCP Apps 用）には入れない。
- 既存の `curation_submit`（貼り付け・`flags` 直接・`source="llm"`）の挙動とエラー文言を変えない。
- 受け口を無効にする設定・期限を変える引数・環境変数・開き直すツールは作らない（spec §8）。
- 文言: レビュービューアは英語、候補付けビューアは日本語にそろえる（spec §6 最終項）。
- ツール戻り値は `json_payload()`、`structured_output=False` を維持（新ツールは足さない）。
- `docs/workflow/curation.md` の呼び出し連鎖に書く関数は実在させる（`tests/test_workflow_docs.py` が AST で検証する）。
- コミットごとに pre-commit で全テストが走り数分かかる。`git commit` はバックグラウンド実行にし、出力をファイルへ落として `FAILED` を後から grep する。実行中に作業ツリーを編集しない。
- `docs/HISTRY.md` / `docs/task.md` は **main ツリー側**（`C:\Users\yuu18\Metabolomix_with_LLM\docs\`）に書く。task.md は末尾追記のみ。

## Review Focus

1. **MCP サーバ再起動後に古い HTML を開いたまま Submit** — 受け口に届かない（status 0）か 401 になり、未送信の変更を失わずに Copy 欄が自動で開くべき。→ Task 6 の `stateAfter` の純関数テストで「live から 0/401 で lost」「起動時の 0 は unavailable」を固定する。
2. **本文が JSON だが dict でない（配列・文字列）・`token` が文字列でない・非 ASCII のトークン** — 500 や接続断ではなく 400 / 401 を返すべき。→ Task 3 にテストを置く。
3. **登録後にレビューの JSON が消された（手で消した・別フォルダへ移した）** — 例外で接続が落ちるのではなく 410 とメッセージを返すべき。→ Task 4 にテストを置く。
4. **送信処理中の想定外の例外** — 500 の JSON を返し、受け口は生き続けて次の ping に 200 を返すべき。→ Task 4 にテストを置く。
5. **未送信の変更が 0 件で Submit** — 要求を送らず「No unsent changes.」を出すべき。→ Task 6 / 7 でクリックハンドラのソースに対する検査を置く。

---

## ファイル構成

| ファイル | 役割 |
|---|---|
| `metabolomix/curation/submission.py`（新規） | 送信本体 `submit_flags()`、保存済みレビューの読み分け `load_saved()`、例外 `SubmissionError`、書き込みロック |
| `metabolomix/curation/submit_server.py`（新規） | 受け口（HTTP ハンドラ）と登録表・時間切れ・起動停止 |
| `metabolomix/tools/curation_tools.py` | `curation_submit` を `submission` 経由に、`curation_review` / `curation_suggest` で受け口に登録し `submit` を返す |
| `metabolomix/curation/viewer.py` | `render_html` / `render_suggest_html` に `submit_endpoint` を足す |
| `metabolomix/curation/review.py` / `suggest.py` | `save_review` / `save_suggestion` に `submit_endpoint` を足す |
| `metabolomix/curation/viewer_common.js` | `SUBMIT_ENDPOINT`、送信クライアントの純関数と `fetch` 部 |
| `metabolomix/curation/viewer.html` / `suggest_viewer.html` | Submit / Finish ボタン・結果表示・配線 |
| `tests/curation_fixtures.py` | `build_saved_review()` を足す |
| `tests/conftest.py` | 受け口を毎テスト後に止める autouse fixture |
| `tests/test_curation_submission.py`（新規） | 送信本体の単体テスト |
| `tests/test_curation_submit_server.py`（新規） | 受け口と登録表のテスト |
| `tests/test_curation_tools.py` / `test_curation_review.py` / `test_curation_suggest_viewer.py` | 回帰と埋め込み・純関数のテスト |
| `USAGE.md` / `docs/output_format/curation.md` / `docs/workflow/curation.md` | 文書 |

---

### Task 1: ブラウザの `file://` から `127.0.0.1` へ届くかの実地確認（spec §5.3。コミットなし）

**Files:**
- Create（スクラッチ、コミットしない）: `<scratchpad>/lna-spike/spike_server.py`, `<scratchpad>/lna-spike/spike.html`

**Interfaces:**
- Consumes: なし
- Produces: 判定結果（「届く」/「許可ダイアログが出る」/「止められる」）を Task 6 以降の前提として記録する。止められる場合は**ここで作業を止めてユーザーに諮る**（許可ダイアログを受け入れるか、受け口が HTML を配る `http://127.0.0.1:<port>/view/<token>` 方式に切り替えるか）。

- [ ] **Step 1: 確認用サーバを書く**

```python
# spike_server.py — 実験用。CORS と Private Network Access の応答を返し、POST の本文をそのまま返す。
import json, sys
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

class H(BaseHTTPRequestHandler):
    def log_message(self, *a): pass
    def _cors(self):
        self.send_header("Access-Control-Allow-Origin", "*")
        self.send_header("Access-Control-Allow-Methods", "POST, OPTIONS")
        self.send_header("Access-Control-Allow-Headers", "Content-Type")
        self.send_header("Access-Control-Allow-Private-Network", "true")
    def do_OPTIONS(self):
        sys.stderr.write(f"OPTIONS {dict(self.headers)}\n")
        self.send_response(204); self._cors(); self.send_header("Content-Length", "0"); self.end_headers()
    def do_POST(self):
        n = int(self.headers.get("Content-Length") or 0)
        body = self.rfile.read(n)
        sys.stderr.write(f"POST Host={self.headers.get('Host')} Origin={self.headers.get('Origin')} body={body!r}\n")
        data = json.dumps({"status": "ok", "echo": body.decode()}).encode()
        self.send_response(200); self._cors()
        self.send_header("Content-Type", "application/json"); self.send_header("Content-Length", str(len(data)))
        self.end_headers(); self.wfile.write(data)

srv = ThreadingHTTPServer(("127.0.0.1", 0), H)
sys.stderr.write(f"PORT {srv.server_address[1]}\n")
srv.serve_forever()
```

- [ ] **Step 2: 確認用ページを書く**（`PORT` を Step 3 で出た番号に置き換える）

```html
<!doctype html><meta charset="utf-8"><title>LNA spike</title>
<pre id="out">…</pre>
<script>
const PORT = 0; // ← 置き換える
const out = document.getElementById("out");
(async () => {
  for (const ct of ["text/plain;charset=UTF-8", "application/json"]) {
    try {
      const r = await fetch(`http://127.0.0.1:${PORT}/v1/ping`, {method: "POST", headers: {"Content-Type": ct},
                                                                  body: JSON.stringify({token: "x"})});
      out.textContent += `\n${ct}: ${r.status} ${await r.text()}`;
    } catch (e) { out.textContent += `\n${ct}: FAILED ${e}`; }
  }
})();
</script>
```

- [ ] **Step 3: サーバを起動する**

Run: `C:/Python314/python.exe <scratchpad>/lna-spike/spike_server.py`（バックグラウンド）
Expected: stderr に `PORT <番号>`。Windows ファイアウォールのダイアログが**出ない**こと。

- [ ] **Step 4: Edge と Chrome で `spike.html` を `file:///...` として開く**

Chrome は claude-in-chrome で開いてよい（新しいタブ）。Edge はユーザーに開いてもらう。
Expected: 両方の Content-Type で `200 {"status":"ok",...}`。許可ダイアログが出た場合はその文言とスクリーンショットを残す。サーバの stderr で `Origin: null`、`Host: 127.0.0.1:<port>` であることを確かめる（Task 4 の Host 検査の前提）。

- [ ] **Step 5: 判定する**

- 両ブラウザで届く → 続行。ページは `Content-Type: text/plain;charset=UTF-8`（事前確認なしの単純要求）で送る前提を確定し、サーバは念のため `OPTIONS` にも応える（Task 3）。
- 許可ダイアログが出る／止められる → **作業を止め、結果を添えてユーザーに諮る**（spec §5.3）。
- 結果は main ツリーの `docs/HISTRY.md` に Task 8 でまとめて書くので、ここでは数行のメモを手元に残す。サーバのプロセスは止める。

---

### Task 2: 送信本体を `submission.py` に切り出す

**Files:**
- Create: `metabolomix/curation/submission.py`
- Modify: `metabolomix/tools/curation_tools.py`（`_load_any` と `curation_submit` の後半）
- Modify: `tests/curation_fixtures.py`（`build_saved_review` を足す）
- Modify: `docs/workflow/curation.md`（`## curation_submit` の連鎖と mermaid）
- Test: `tests/test_curation_submission.py`（新規）

**Interfaces:**
- Consumes: `flag_log.validate_entries` / `alignment_key` / `curation_dir` / `FlagStore` / `FlagFileError`、`suggest.expand_entries` / `is_valid_suggestion_id` / `load_suggestion`、`review.load_review`、`msdial_writeback.sync_tags`
- Produces:
  - `submission.SubmissionError(ValueError)`: 属性 `kind: str`（`"invalid"` / `"alignment_changed"` / `"flag_file"`）と `details: dict`
  - `submission.load_saved(directory, review_id: str) -> dict`（無ければ `FileNotFoundError`、形が不正なら `ValueError`）
  - `submission.submit_flags(saved: dict, entries, *, review_id: str, source: str) -> dict` — 戻り値は従来の `curation_submit` の JSON と同じ dict（`status, recorded, review_id, n_wrong, n_suspect, n_confirmed, n_assign, n_redundant, tags_xml`）
  - `tests.curation_fixtures.build_saved_review(folder: Path) -> dict` — `{"paths": <write_alignment_set の戻り値>, "review": <保存済みレビュー dict>}`

- [ ] **Step 1: fixture ヘルパを足す**（`tests/curation_fixtures.py` の末尾）

```python
def build_saved_review(folder: Path) -> dict:
    """`write_alignment_set` → `review.run_review` → `review.save_review` まで済ませる（送信系のテスト用）。
    `LIBRARY_CACHE_ENV` は呼び出し側が tmp に向けておく。"""
    from metabolomix.arf2.reader import load_catalog
    from metabolomix.curation import evidence, judge, review
    from metabolomix.library import store as library_store

    paths = write_alignment_set(folder)
    store = library_store.open_store(paths["msp"])
    try:
        spots = evidence.select_spots(load_catalog(paths["arf2"]), ontology=None, name_contains=None)
        result = review.run_review(paths["arf2"], spots, store=store, ms2_tol=0.025,
                                   th=judge.resolve_thresholds(None), file_ids=None, max_traces=12,
                                   selection={"kind": "annotated"})
    finally:
        store.close()
    review.save_review(result)
    return {"paths": paths, "review": result}
```

- [ ] **Step 2: 失敗するテストを書く**（`tests/test_curation_submission.py`）

```python
import pytest

from metabolomix.curation import flags, submission
from metabolomix.library import store as library_store
from metabolomix.msdial import tags as msdial_tags
from tests.curation_fixtures import build_saved_review


@pytest.fixture()
def saved(tmp_path, monkeypatch):
    monkeypatch.setenv(library_store.LIBRARY_CACHE_ENV, str(tmp_path / "cache"))
    return build_saved_review(tmp_path / "neg")


def _load(saved):
    r = saved["review"]
    return submission.load_saved(flags.curation_dir(r["arf2_path"]), r["review_id"])


def test_submit_flags_records_and_updates_the_tags_file(saved):
    r = saved["review"]
    body = submission.submit_flags(_load(saved), [{"spot_id": 1, "flag": "wrong", "note": "x"}],
                                   review_id=r["review_id"], source="user")
    assert body["status"] == "ok" and body["recorded"] == 1 and body["n_wrong"] == 1
    assert body["tags_xml"]["added"] == [1]
    rows = flags.FlagStore(flags.curation_dir(r["arf2_path"])).rows()
    assert rows[-1]["source"] == "user" and rows[-1]["review_id"] == r["review_id"]
    tag_path = msdial_tags.alignment_tag_path(saved["paths"]["arf2"])
    assert msdial_tags.parse_tag_file(tag_path)["peaks"] == {1: frozenset({3})}


def test_invalid_entries_raise_invalid_without_writing(saved):
    r = saved["review"]
    with pytest.raises(submission.SubmissionError) as info:
        submission.submit_flags(_load(saved), [{"spot_id": 1, "flag": "nope"}],
                                review_id=r["review_id"], source="user")
    assert info.value.kind == "invalid"
    assert not flags.FlagStore(flags.curation_dir(r["arf2_path"])).exists()


def test_duplicate_spot_ids_are_invalid(saved):
    r = saved["review"]
    with pytest.raises(submission.SubmissionError) as info:
        submission.submit_flags(_load(saved), [{"spot_id": 0, "flag": "wrong"}, {"spot_id": 0, "flag": "suspect"}],
                                review_id=r["review_id"], source="user")
    assert info.value.kind == "invalid" and "[0]" in str(info.value)


def test_alignment_change_is_its_own_kind(saved):
    r = saved["review"]
    with open(saved["paths"]["arf2"], "ab") as handle:
        handle.write(b"\x00")
    with pytest.raises(submission.SubmissionError) as info:
        submission.submit_flags(_load(saved), [{"spot_id": 0, "flag": "wrong"}],
                                review_id=r["review_id"], source="user")
    assert info.value.kind == "alignment_changed"


def test_a_corrupt_flags_file_is_flag_file_with_details(saved):
    r = saved["review"]
    path = flags.curation_dir(r["arf2_path"]) / flags.FLAGS_FILENAME
    path.write_text("not json\n", encoding="utf-8")
    with pytest.raises(submission.SubmissionError) as info:
        submission.submit_flags(_load(saved), [{"spot_id": 0, "flag": "wrong"}],
                                review_id=r["review_id"], source="user")
    assert info.value.kind == "flag_file"
    assert info.value.details == {"flags_file": str(path), "line": 1}
    assert path.read_text(encoding="utf-8") == "not json\n"


def test_load_saved_raises_file_not_found_for_a_missing_review(saved):
    r = saved["review"]
    with pytest.raises(FileNotFoundError):
        submission.load_saved(flags.curation_dir(r["arf2_path"]), "cr-20000101-000000-abcd")
```

- [ ] **Step 3: 失敗を確かめる**

Run: `C:/Python314/python.exe -m pytest tests/test_curation_submission.py -q`
Expected: FAIL（`ImportError: cannot import name 'submission'`）

- [ ] **Step 4: `submission.py` を書く**

```python
"""キュレーションの送信本体: 保存済みレビュー（cr-…）か候補付け（cs-…）に対するフラグを検証し、
`flags.jsonl` へ追記して `_tags.xml` に反映する。MCP ツール `curation_submit` とビューアの受け口
（`submit_server`）の両方がここを呼ぶ。spec 2026-10-09。

グローバルな session に依存しない（受け口の HTTP スレッドからも呼ぶため）。`flags.jsonl` と
`_tags.xml` に書く区間は `_WRITE_LOCK` で直列化する（MCP ツールと HTTP スレッドが同時に書きうる）。
deps: curation.flags / review / suggest / msdial_writeback。
"""
from __future__ import annotations

import threading

from metabolomix.curation import flags as flag_log
from metabolomix.curation import msdial_writeback, review, suggest

_WRITE_LOCK = threading.Lock()


class SubmissionError(ValueError):
    """送信を記録しなかった理由。`kind` は invalid（形・候補・重複）/ alignment_changed /
    flag_file（`flags.jsonl` が読めない。`details` に `flags_file` と `line`）。"""

    def __init__(self, message: str, *, kind: str, details: dict | None = None):
        super().__init__(message)
        self.kind = kind
        self.details = details or {}


def load_saved(directory, review_id: str) -> dict:
    """`review_id` の接頭辞で読み分ける。無ければ FileNotFoundError、形が不正なら ValueError。"""
    if suggest.is_valid_suggestion_id(review_id):
        return suggest.load_suggestion(directory, review_id)
    return review.load_review(directory, review_id)


def _count(effective: dict, flag: str) -> int:
    return sum(1 for r in effective.values() if r["flag"] == flag)


def submit_flags(saved: dict, entries, *, review_id: str, source: str) -> dict:
    """検証 → 同一 spot_id の重複検出 → アラインメントの sha 照合 → 追記 → `_tags.xml` 反映。
    不正があれば何も書かずに `SubmissionError`。`_tags.xml` の失敗は記録を残して `tags_xml.error`。"""
    try:
        if suggest.is_valid_suggestion_id(review_id):
            cleaned = suggest.expand_entries(entries, saved)
        else:
            cleaned = flag_log.validate_entries(
                entries, allowed_spot_ids={s["spot_id"] for s in saved["spots"]})
    except ValueError as exc:
        raise SubmissionError(str(exc), kind="invalid") from exc
    seen_spot_ids = set()
    duplicated = sorted({e["spot_id"] for e in cleaned if e["spot_id"] in seen_spot_ids
                         or seen_spot_ids.add(e["spot_id"])})
    if duplicated:
        raise SubmissionError(f"同じ spot_id を 1 回の送信で複数回指定しています: {duplicated}",
                              kind="invalid")
    with _WRITE_LOCK:
        current = flag_log.alignment_key(saved["arf2_path"])
        if current["alignment_sha256"] != saved["alignment"]["alignment_sha256"]:
            raise SubmissionError(
                "レビューの後でアラインメント（.arf2）が変わっています。curation_review（候補付けなら "
                "curation_suggest）をやり直してください。", kind="alignment_changed")
        store = flag_log.FlagStore(flag_log.curation_dir(saved["arf2_path"]))
        try:
            store.rows()                   # 壊れた記録に追記しない（先に読めるか確かめる）
        except flag_log.FlagFileError as exc:
            raise SubmissionError(str(exc), kind="flag_file", details=exc.details()) from exc
        n = store.append(cleaned, alignment=current, review_id=review_id, source=source)
        effective = store.effective(current["alignment_sha256"])
        tags_xml = msdial_writeback.sync_tags(saved["arf2_path"], cleaned)
    return {"status": "ok", "recorded": n, "review_id": review_id,
            "n_wrong": _count(effective, "wrong"), "n_suspect": _count(effective, "suspect"),
            "n_confirmed": _count(effective, "confirmed"), "n_assign": _count(effective, "assign"),
            "n_redundant": _count(effective, "redundant"), "tags_xml": tags_xml}
```

注意: alignment_changed のメッセージは現行 `curation_tools.py` の文字列と**完全に同じ**にする（上の 2 行連結で同一になる。改行位置で空白を入れない）。

- [ ] **Step 5: `curation_tools.py` を `submission` 経由にする**

import に `submission` を足す:

```python
from metabolomix.curation import evidence, judge, msdial_writeback, review, submission, suggest, viewer
```

（`msdial_writeback` が他で使われていなければ import から外す。`grep -n msdial_writeback metabolomix/tools/curation_tools.py` で確かめる。）

`_load_any` を置き換える:

```python
def _load_any(directory, review_id: str) -> dict:
    return submission.load_saved(directory, review_id)
```

`curation_submit` の `saved, searched = _find_review(...)` / `_not_found` の後ろ（`try: if suggest.is_valid_suggestion_id(review_id): ...` から `return json_payload({...})` まで）を次に置き換える:

```python
    try:
        result = submission.submit_flags(saved, entries, review_id=review_id, source=source)
    except submission.SubmissionError as exc:
        return _error(str(exc), **exc.details)
    return json_payload(result)
```

- [ ] **Step 6: テストを走らせる**

Run: `C:/Python314/python.exe -m pytest tests/test_curation_submission.py tests/test_curation_tools.py -q`
Expected: 全 PASS（既存の `curation_submit` のテストも含む）

- [ ] **Step 7: `docs/workflow/curation.md` の `## curation_submit` を直す**

mermaid の `CS --> VAL` `CS --> APP` `CS --> WB` を、`CS --> SUB[curation.submission.submit_flags]` と `SUB --> VAL` `SUB --> APP` `SUB --> WB` に、`CS --> LS` `CS --> EX` を `SUB --> EX` と `CS --> LD[curation.submission.load_saved]` `LD --> LS` に付け替える。ついでに古いノード名 `curation.msdial_writeback.sync_misannotation` を `curation.msdial_writeback.sync_tags` に直す。連鎖（番号付きの行）は次に置き換える:

```markdown
1. metabolomix/curation/flags.py  parse_submission_text()（`submission_text` のとき）
2. metabolomix/curation/review.py  is_valid_review_id()（形が違えばパスに使う前にエラー。`cs-…` は metabolomix/curation/suggest.py  is_valid_suggestion_id() も許す）
3. metabolomix/tools/curation_tools.py  _find_review()
4. └─ metabolomix/tools/curation_tools.py  _load_any()
5.    └─ metabolomix/curation/submission.py  load_saved()（ID の接頭辞で読み分ける）
6.       ├─ [cr-…] metabolomix/curation/review.py  load_review()（候補フォルダを順に）
7.       └─ [cs-…] metabolomix/curation/suggest.py  load_suggestion()（候補フォルダを順に）
8. metabolomix/curation/submission.py  submit_flags()（ビューアの受け口と共有。`SubmissionError` → エラー payload）
9. ├─ [cr-…] metabolomix/curation/flags.py  validate_entries()
10. ├─ [cs-…] metabolomix/curation/suggest.py  expand_entries()（候補 ID を記録行へ展開。偽の候補 ID は書く前に拒否）
11. └─ 以下は `_WRITE_LOCK` の中（MCP ツールと HTTP スレッドの書き込みを直列化）
12.    ├─ metabolomix/curation/flags.py  alignment_key()（レビュー時の sha256 と一致しなければ拒否）
13.    ├─ metabolomix/curation/flags.py  FlagStore.rows()（読めない行があれば追記せずにエラー）
14.    ├─ metabolomix/curation/flags.py  FlagStore.append()
15.    └─ metabolomix/curation/msdial_writeback.py  sync_tags()（失敗しても記録は残し `tags_xml.error`。assign / redundant は触らない）
16.       └─ metabolomix/msdial/tags.py  update_alignment_tags()（Misannotation と Confirmed を 1 回で書く）
17.          └─ metabolomix/core/atomic_io.py  atomic_write_bytes()
```

注意: 11 行目はファイル名を含まないので `CHAIN_RE` に掛からず検証対象外（説明行として許される）。

- [ ] **Step 8: 文書テストを走らせる**

Run: `C:/Python314/python.exe -m pytest tests/test_workflow_docs.py -q`
Expected: PASS

- [ ] **Step 9: コミット**（バックグラウンドで、出力をファイルへ）

```bash
git add metabolomix/curation/submission.py metabolomix/tools/curation_tools.py tests/curation_fixtures.py tests/test_curation_submission.py docs/workflow/curation.md
git commit -m "refactor(curation): 送信本体を submission.submit_flags に切り出し書き込みをロックで直列化する

Co-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>" > <scratchpad>/commit-task2.log 2>&1
```

完了後 `grep -n "FAILED\|passed\|failed" <scratchpad>/commit-task2.log` で確かめる。

---

### Task 3: 受け口の登録表・起動停止・時間切れ・ping

**Files:**
- Create: `metabolomix/curation/submit_server.py`
- Modify: `tests/conftest.py`（autouse で受け口を止める）
- Test: `tests/test_curation_submit_server.py`（新規）

**Interfaces:**
- Consumes: なし（Task 4 で `submission` を使う）
- Produces（Task 4〜5 が使う）:
  - 定数 `IDLE_TIMEOUT_MIN = 30`, `SWEEP_INTERVAL_S = 60`, `MAX_BODY_BYTES = 8 * 1024 * 1024`, `KINDS = ("review", "suggest")`
  - `register(kind: str, review_id: str, arf2_path) -> {"port": int, "token": str, "idle_timeout_min": int}`（`kind` が不正なら `ValueError`、bind 失敗は `OSError` をそのまま上げる）
  - `unregister(token: str, *, wait: bool = True) -> None`
  - `sweep() -> None`（時間切れの登録を外し、0 件なら止める）
  - `is_running() -> bool`, `current_port() -> int | None`, `reset() -> None`（テスト用）
  - モジュール属性 `_clock`（既定 `time.monotonic`。テストが差し替える）
  - ハンドラ内部: `_Handler._reply(code, body)`、`_Handler._read_json()`、ルート表 `_ROUTES: dict[str, callable]`（`(handler, token, entry, body) -> (code, payload)`）

- [ ] **Step 1: conftest に autouse fixture を足す**（`tests/conftest.py` の末尾。`import sys` を先頭の import に足す）

```python
@pytest.fixture(autouse=True)
def _stop_curation_submit_server():
    """curation_review / curation_suggest はビューアの受け口（daemon スレッドの HTTP サーバ）を
    立てる。テストをまたいで登録と待受けを残さない。未 import なら何もしない（重い import を避ける）。"""
    yield
    module = sys.modules.get("metabolomix.curation.submit_server")
    if module is not None:
        module.reset()
```

- [ ] **Step 2: 失敗するテストを書く**（`tests/test_curation_submit_server.py`）

```python
import http.client
import json

import pytest

from metabolomix.curation import submit_server


class FakeClock:
    def __init__(self):
        self.now = 1000.0

    def __call__(self):
        return self.now

    def advance(self, minutes: float):
        self.now += minutes * 60


@pytest.fixture()
def clock(monkeypatch):
    fake = FakeClock()
    monkeypatch.setattr(submit_server, "_clock", fake)
    return fake


def _post(port, path, body=None, *, raw=None, headers=None, method="POST"):
    conn = http.client.HTTPConnection("127.0.0.1", port, timeout=10)
    try:
        data = raw if raw is not None else (None if body is None else json.dumps(body).encode("utf-8"))
        conn.request(method, path, body=data,
                     headers={"Content-Type": "text/plain;charset=UTF-8", **(headers or {})})
        res = conn.getresponse()
        text = res.read().decode("utf-8")
        return res.status, (json.loads(text) if text else None)
    finally:
        conn.close()


def test_register_starts_the_endpoint_and_later_registrations_share_it():
    a = submit_server.register("review", "cr-20260101-000000-aaaa", "C:/x.arf2")
    b = submit_server.register("suggest", "cs-20260101-000000-bbbb", "C:/x.arf2")
    assert a["port"] == b["port"] == submit_server.current_port()
    assert a["token"] != b["token"] and len(a["token"]) >= 40
    assert a["idle_timeout_min"] == 30
    status, body = _post(a["port"], "/v1/ping", {"token": b["token"]})
    assert status == 200
    assert body["status"] == "ok" and body["kind"] == "suggest"
    assert body["review_id"] == "cs-20260101-000000-bbbb" and body["expires_at"]


def test_the_last_unregister_stops_the_endpoint():
    a = submit_server.register("review", "cr-20260101-000000-aaaa", "C:/x.arf2")
    b = submit_server.register("review", "cr-20260101-000000-bbbb", "C:/x.arf2")
    submit_server.unregister(a["token"])
    assert submit_server.is_running()
    submit_server.unregister(b["token"])
    assert not submit_server.is_running() and submit_server.current_port() is None
    with pytest.raises(OSError):
        _post(a["port"], "/v1/ping", {"token": b["token"]})


def test_register_rejects_an_unknown_kind():
    with pytest.raises(ValueError):
        submit_server.register("other", "cr-20260101-000000-aaaa", "C:/x.arf2")


def test_ping_extends_the_deadline_and_idle_registrations_are_swept(clock):
    a = submit_server.register("review", "cr-20260101-000000-aaaa", "C:/x.arf2")
    clock.advance(29)
    assert _post(a["port"], "/v1/ping", {"token": a["token"]})[0] == 200     # heartbeat で延びる
    clock.advance(29)
    submit_server.sweep()
    assert submit_server.is_running()
    clock.advance(2)
    submit_server.sweep()
    assert not submit_server.is_running()


def test_an_expired_token_is_refused_before_the_sweep(clock):
    a = submit_server.register("review", "cr-20260101-000000-aaaa", "C:/x.arf2")
    clock.advance(31)
    assert _post(a["port"], "/v1/ping", {"token": a["token"]})[0] == 401


@pytest.mark.parametrize("raw", [b"[1, 2]", b'"token"', b'{"token": 5}', '{"token": "ｔｏｋｅｎ"}'.encode("utf-8")])
def test_non_object_bodies_and_bad_tokens_are_400_or_401(raw):
    a = submit_server.register("review", "cr-20260101-000000-aaaa", "C:/x.arf2")
    status, body = _post(a["port"], "/v1/ping", raw=raw)
    assert status in (400, 401) and body["status"] == "error"
    assert _post(a["port"], "/v1/ping", {"token": a["token"]})[0] == 200    # 受け口は生きている


def test_options_preflight_allows_private_network_access():
    a = submit_server.register("review", "cr-20260101-000000-aaaa", "C:/x.arf2")
    conn = http.client.HTTPConnection("127.0.0.1", a["port"], timeout=10)
    try:
        conn.request("OPTIONS", "/v1/submit", headers={"Access-Control-Request-Private-Network": "true"})
        res = conn.getresponse(); res.read()
    finally:
        conn.close()
    assert res.status == 204
    assert res.getheader("Access-Control-Allow-Origin") == "*"
    assert res.getheader("Access-Control-Allow-Private-Network") == "true"
```

- [ ] **Step 3: 失敗を確かめる**

Run: `C:/Python314/python.exe -m pytest tests/test_curation_submit_server.py -q`
Expected: FAIL（`ImportError: cannot import name 'submit_server'`）

- [ ] **Step 4: `submit_server.py` を書く**（ping だけのルート表。submit / finish・Host 検査・上限は Task 4）

```python
"""キュレーションのビューアが直接送信する受け口（127.0.0.1 の HTTP）。spec 2026-10-09。

プロセスに 1 つを共有し、curation_review / curation_suggest のたびにトークンを登録する。
登録は `{kind, review_id, arf2_path, last_seen}`。書き込み先は登録側の値だけを使い、
ページからは flags しか受け取らない。最後の通信から `IDLE_TIMEOUT_MIN` 分で登録を外し、
0 件になれば待受けを止める。daemon スレッドなので MCP サーバが終われば一緒に消える。

**stdout に何も書かない**（stdio の MCP ではプロトコルの通り道）。`log_message` も止める。
deps: curation.submission、標準ライブラリ。session を import しない（HTTP スレッドから呼ぶ）。
"""
from __future__ import annotations

import hmac
import json
import secrets
import socketserver
import threading
import time
from datetime import datetime, timedelta, timezone
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

IDLE_TIMEOUT_MIN = 30
SWEEP_INTERVAL_S = 60
MAX_BODY_BYTES = 8 * 1024 * 1024        # MAX_SPOTS（3000 件）にメモが付いても収まる
KINDS = ("review", "suggest")

_clock = time.monotonic                 # テストが差し替える
_lock = threading.Lock()
_entries: dict[str, dict] = {}          # token -> {kind, review_id, arf2_path, last_seen}
_server: "_Server | None" = None
_sweeper_stop: threading.Event | None = None


def _error_body(message: str, **details) -> dict:
    return {"status": "error", "message": message, **details}


class _Server(ThreadingHTTPServer):
    allow_reuse_address = False         # Windows の SO_REUSEADDR は同じポートの横取りを許す

    def server_bind(self):
        # HTTPServer.server_bind は socket.getfqdn() を呼び、Windows では逆引きで数秒止まることがある。
        socketserver.TCPServer.server_bind(self)
        self.server_name, self.server_port = self.server_address[:2]


class _Handler(BaseHTTPRequestHandler):
    server_version = "ms-data-parser-curation"
    sys_version = ""

    def log_message(self, format, *args):       # 既定のアクセスログ（stderr）を出さない
        pass

    def _cors(self):
        self.send_header("Access-Control-Allow-Origin", "*")
        self.send_header("Access-Control-Allow-Methods", "POST, OPTIONS")
        self.send_header("Access-Control-Allow-Headers", "Content-Type")
        self.send_header("Access-Control-Allow-Private-Network", "true")

    def _reply(self, code: int, body: dict):
        data = json.dumps(body, ensure_ascii=False).encode("utf-8")
        self.send_response(code)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(data)))
        self._cors()
        self.end_headers()
        self.wfile.write(data)

    def do_OPTIONS(self):
        self.send_response(204)
        self._cors()
        self.send_header("Content-Length", "0")
        self.end_headers()

    def _read_json(self):
        """本文を dict として読む。読めなければ (None, (code, payload))。"""
        try:
            length = int(self.headers.get("Content-Length") or "")
        except ValueError:
            return None, (400, _error_body("Content-Length がありません。"))
        try:
            body = json.loads(self.rfile.read(length).decode("utf-8"))
        except (UnicodeDecodeError, ValueError):
            return None, (400, _error_body("本文が JSON ではありません。"))
        if not isinstance(body, dict):
            return None, (400, _error_body("本文は JSON オブジェクトにしてください。"))
        return body, None

    def do_POST(self):
        route = _ROUTES.get(self.path)
        if route is None:
            return self._reply(404, _error_body(f"不明なパスです: {self.path}"))
        body, problem = self._read_json()
        if problem:
            return self._reply(*problem)
        token = body.get("token")
        entry = _touch(token) if isinstance(token, str) else None
        if entry is None:
            return self._reply(401, _error_body("トークンが不明か期限切れです。送信用テキストをコピーしてチャットに貼ってください。"))
        self._reply(*route(self, token, entry, body))


def _ping(handler, token, entry, body):
    return 200, {"status": "ok", "kind": entry["kind"], "review_id": entry["review_id"],
                 "expires_at": _expires_at(entry)}


_ROUTES = {"/v1/ping": _ping}


def _expired(entry: dict, now: float) -> bool:
    return now - entry["last_seen"] >= IDLE_TIMEOUT_MIN * 60


def _expires_at(entry: dict) -> str:
    remaining = IDLE_TIMEOUT_MIN * 60 - (_clock() - entry["last_seen"])
    return (datetime.now(timezone.utc) + timedelta(seconds=remaining)).isoformat(timespec="seconds")


def _touch(token: str) -> dict | None:
    """トークンを定数時間で照合し、生きていれば期限を延ばして登録の写しを返す。"""
    now = _clock()
    given = token.encode("utf-8")
    with _lock:
        for key, entry in _entries.items():
            if hmac.compare_digest(key.encode("utf-8"), given):
                if _expired(entry, now):
                    return None
                entry["last_seen"] = now
                return dict(entry)
    return None


def _ensure_running() -> None:
    global _server, _sweeper_stop
    if _server is not None:
        return
    server = _Server(("127.0.0.1", 0), _Handler)
    threading.Thread(target=server.serve_forever, name="curation-submit", daemon=True).start()
    stop = threading.Event()
    threading.Thread(target=_sweep_loop, args=(stop,), name="curation-submit-sweep", daemon=True).start()
    _server, _sweeper_stop = server, stop


def _sweep_loop(stop: threading.Event) -> None:
    while not stop.wait(SWEEP_INTERVAL_S):
        sweep()


def _detach_if_idle():
    """（_lock の中で呼ぶ）登録が 0 件なら待受けを切り離して返す。止めるのは呼び出し側が鍵の外で。"""
    global _server, _sweeper_stop
    if _entries or _server is None:
        return None
    detached = (_server, _sweeper_stop)
    _server = _sweeper_stop = None
    return detached


def _shutdown(detached, *, wait: bool) -> None:
    server, stop = detached
    stop.set()

    def run():
        server.shutdown()
        server.server_close()
    if wait:
        run()
    else:       # 要求を処理しているスレッドからは、応答を返した後に別スレッドで止める
        threading.Thread(target=run, name="curation-submit-stop", daemon=True).start()


def register(kind: str, review_id: str, arf2_path) -> dict:
    if kind not in KINDS:
        raise ValueError(f"kind は {KINDS} のいずれか: {kind!r}")
    with _lock:
        _ensure_running()
        token = secrets.token_urlsafe(32)
        _entries[token] = {"kind": kind, "review_id": review_id, "arf2_path": str(arf2_path),
                           "last_seen": _clock()}
        port = _server.server_address[1]
    return {"port": port, "token": token, "idle_timeout_min": IDLE_TIMEOUT_MIN}


def unregister(token: str, *, wait: bool = True) -> None:
    with _lock:
        _entries.pop(token, None)
        detached = _detach_if_idle()
    if detached:
        _shutdown(detached, wait=wait)


def sweep() -> None:
    now = _clock()
    with _lock:
        for token in [t for t, e in _entries.items() if _expired(e, now)]:
            del _entries[token]
        detached = _detach_if_idle()
    if detached:
        _shutdown(detached, wait=True)


def is_running() -> bool:
    with _lock:
        return _server is not None


def current_port() -> int | None:
    with _lock:
        return None if _server is None else _server.server_address[1]


def reset() -> None:
    """テスト用: 登録を全部外して待受けを止める。"""
    with _lock:
        _entries.clear()
        detached = _detach_if_idle()
    if detached:
        _shutdown(detached, wait=True)
```

- [ ] **Step 5: テストを走らせる**

Run: `C:/Python314/python.exe -m pytest tests/test_curation_submit_server.py -q`
Expected: 全 PASS（停止後の接続テストは Windows で 2 秒ほどかかることがある）

- [ ] **Step 6: コミット**（バックグラウンド、ログを `<scratchpad>/commit-task3.log`）

```bash
git add metabolomix/curation/submit_server.py tests/test_curation_submit_server.py tests/conftest.py
git commit -m "feat(curation): ビューアの送信を受ける 127.0.0.1 の受け口と登録表・時間切れを足す

Co-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>"
```

---

### Task 4: submit / finish と安全対策・エラー応答

**Files:**
- Modify: `metabolomix/curation/submit_server.py`
- Test: `tests/test_curation_submit_server.py`（追記）

**Interfaces:**
- Consumes: `submission.load_saved` / `submission.submit_flags` / `submission.SubmissionError`（Task 2）、`flag_log.curation_dir`、Task 3 の `_ROUTES` / `_touch` / `unregister`
- Produces: HTTP 仕様（spec §5）。`POST /v1/submit` 本文 `{token, review_id, flags}` → 200 で `submit_flags` の戻り値。状態コード: 401 / 403（`review_id` 不一致・`Host` 不正）/ 404 / 405（GET）/ 409（alignment_changed）/ 410（保存済みレビューが無い）/ 413 / 400（invalid・JSON でない）/ 500（flag_file・想定外の例外）。`POST /v1/finish` → `{status:"ok"}` の後に登録を外す。

- [ ] **Step 1: 失敗するテストを追記する**

```python
import threading

from metabolomix.curation import flags, submission
from metabolomix.library import store as library_store
from metabolomix.msdial import tags as msdial_tags
from tests.curation_fixtures import build_saved_review


@pytest.fixture()
def saved(tmp_path, monkeypatch):
    monkeypatch.setenv(library_store.LIBRARY_CACHE_ENV, str(tmp_path / "cache"))
    built = build_saved_review(tmp_path / "neg")
    r = built["review"]
    built["endpoint"] = submit_server.register("review", r["review_id"], r["arf2_path"])
    return built


def _submit(saved, flag_entries, **overrides):
    ep, r = saved["endpoint"], saved["review"]
    body = {"token": ep["token"], "review_id": r["review_id"], "flags": flag_entries, **overrides}
    return _post(ep["port"], "/v1/submit", body)


def _n_rows(saved):
    return len(flags.FlagStore(flags.curation_dir(saved["review"]["arf2_path"])).rows())


def test_submit_records_as_the_user_and_updates_the_tags_file(saved):
    status, body = _submit(saved, [{"spot_id": 1, "flag": "wrong", "note": "x"}])
    assert status == 200 and body["recorded"] == 1 and body["tags_xml"]["added"] == [1]
    rows = flags.FlagStore(flags.curation_dir(saved["review"]["arf2_path"])).rows()
    assert rows[-1]["source"] == "user"
    tag_path = msdial_tags.alignment_tag_path(saved["paths"]["arf2"])
    assert msdial_tags.parse_tag_file(tag_path)["peaks"] == {1: frozenset({3})}


def test_submit_can_be_sent_in_several_parts(saved):
    assert _submit(saved, [{"spot_id": 1, "flag": "wrong", "note": ""}])[0] == 200
    assert _submit(saved, [{"spot_id": 0, "flag": "confirmed", "note": ""}])[0] == 200
    assert _n_rows(saved) == 2


def test_a_mismatched_review_id_is_403_and_writes_nothing(saved):
    status, body = _submit(saved, [{"spot_id": 1, "flag": "wrong"}], review_id="cr-20000101-000000-abcd")
    assert status == 403 and body["status"] == "error"
    assert _n_rows(saved) == 0


def test_a_foreign_host_header_is_403(saved):
    ep = saved["endpoint"]
    status, _ = _post(ep["port"], "/v1/ping", {"token": ep["token"]},
                      headers={"Host": f"evil.example:{ep['port']}"})
    assert status == 403


def test_an_oversized_body_is_413(saved, monkeypatch):
    monkeypatch.setattr(submit_server, "MAX_BODY_BYTES", 64)
    status, _ = _submit(saved, [{"spot_id": 1, "flag": "wrong", "note": "x" * 100}])
    assert status == 413
    assert _n_rows(saved) == 0


def test_non_json_is_400_unknown_path_404_and_get_405(saved):
    ep = saved["endpoint"]
    assert _post(ep["port"], "/v1/submit", raw=b"not json")[0] == 400
    assert _post(ep["port"], "/v1/nope", {"token": ep["token"]})[0] == 404
    assert _post(ep["port"], "/v1/ping", method="GET")[0] == 405


def test_an_invalid_flag_is_400_and_writes_nothing(saved):
    status, body = _submit(saved, [{"spot_id": 1, "flag": "nope"}])
    assert status == 400 and "flag" in body["message"]
    assert _n_rows(saved) == 0


def test_an_alignment_change_is_409(saved):
    with open(saved["paths"]["arf2"], "ab") as handle:
        handle.write(b"\x00")
    status, body = _submit(saved, [{"spot_id": 1, "flag": "wrong"}])
    assert status == 409 and "curation_review" in body["message"]


def test_a_deleted_review_is_410(saved):
    r = saved["review"]
    (flags.curation_dir(r["arf2_path"]) / f"review-{r['review_id']}.json").unlink()
    status, body = _submit(saved, [{"spot_id": 1, "flag": "wrong"}])
    assert status == 410 and body["status"] == "error"


def test_an_unexpected_error_is_500_and_the_endpoint_survives(saved, monkeypatch):
    def boom(*args, **kwargs):
        raise RuntimeError("boom")
    monkeypatch.setattr(submission, "submit_flags", boom)
    status, body = _submit(saved, [{"spot_id": 1, "flag": "wrong"}])
    assert status == 500 and body["status"] == "error"
    ep = saved["endpoint"]
    assert _post(ep["port"], "/v1/ping", {"token": ep["token"]})[0] == 200


def test_finish_unregisters_the_token(saved):
    other = submit_server.register("review", "cr-20260101-000000-zzzz", "C:/x.arf2")   # 待受けを残す
    ep = saved["endpoint"]
    status, body = _post(ep["port"], "/v1/finish", {"token": ep["token"]})
    assert status == 200 and body == {"status": "ok"}
    assert _post(ep["port"], "/v1/ping", {"token": ep["token"]})[0] == 401
    assert _post(ep["port"], "/v1/ping", {"token": other["token"]})[0] == 200


def test_concurrent_submits_are_both_recorded(saved):
    results = []
    threads = [threading.Thread(target=lambda sid=sid: results.append(
        _submit(saved, [{"spot_id": sid, "flag": "wrong", "note": ""}])[0])) for sid in (0, 1)]
    for t in threads:
        t.start()
    for t in threads:
        t.join(timeout=60)
    assert results == [200, 200]
    assert _n_rows(saved) == 2
    tag_path = msdial_tags.alignment_tag_path(saved["paths"]["arf2"])
    assert msdial_tags.parse_tag_file(tag_path)["peaks"] == {0: frozenset({3}), 1: frozenset({3})}


def test_requests_write_nothing_to_stdout(saved, capfd):
    capfd.readouterr()
    ep = saved["endpoint"]
    _post(ep["port"], "/v1/ping", {"token": ep["token"]})
    _post(ep["port"], "/v1/nope", {"token": ep["token"]})
    _submit(saved, [{"spot_id": 1, "flag": "nope"}])
    _submit(saved, [{"spot_id": 1, "flag": "wrong"}])
    assert capfd.readouterr().out == ""
```

- [ ] **Step 2: 失敗を確かめる**

Run: `C:/Python314/python.exe -m pytest tests/test_curation_submit_server.py -q`
Expected: 新しいテストが FAIL（404 / Host 未検査 など）。Task 3 のテストは PASS のまま。

- [ ] **Step 3: 実装する**

import に足す:

```python
from metabolomix.curation import flags as flag_log
from metabolomix.curation import submission
```

`_SUBMISSION_STATUS` とルートを足し、`_ROUTES` を置き換える:

```python
_SUBMISSION_STATUS = {"invalid": 400, "alignment_changed": 409, "flag_file": 500}


def _submit(handler, token, entry, body):
    if body.get("review_id") != entry["review_id"]:
        return 403, _error_body("review_id がこのビューアの登録と違います（別の HTML の取り違え）。")
    try:
        saved = submission.load_saved(flag_log.curation_dir(entry["arf2_path"]), entry["review_id"])
    except FileNotFoundError:
        return 410, _error_body(f"review_id={entry['review_id']} の保存済みレビューが見つかりません。"
                                "curation_review（候補付けなら curation_suggest）をやり直してください。")
    try:
        return 200, submission.submit_flags(saved, body.get("flags"), review_id=entry["review_id"],
                                            source="user")
    except submission.SubmissionError as exc:
        return _SUBMISSION_STATUS[exc.kind], _error_body(str(exc), **exc.details)


def _finish(handler, token, entry, body):
    return 200, {"status": "ok"}


_ROUTES = {"/v1/ping": _ping, "/v1/submit": _submit, "/v1/finish": _finish}
```

`_Handler` に Host 検査・GET・上限・例外の受け止め・finish 後の登録解除を入れる（`do_OPTIONS` と `do_POST` と `_read_json` を置き換え、`do_GET` と `_host_ok` を足す）:

```python
    def _host_ok(self) -> bool:      # DNS rebinding 対策: 127.0.0.1:<port> 以外の Host は拒む
        return self.headers.get("Host") == f"127.0.0.1:{self.server.server_address[1]}"

    def do_OPTIONS(self):
        if not self._host_ok():
            return self._reply(403, _error_body("Host が不正です。"))
        self.send_response(204)
        self._cors()
        self.send_header("Content-Length", "0")
        self.end_headers()

    def do_GET(self):               # GET を返さないので、他のサイトから中身を覗けない
        self._reply(405, _error_body("POST だけを受け付けます。"))

    def _read_json(self):
        """本文を dict として読む。読めなければ (None, (code, payload))。"""
        try:
            length = int(self.headers.get("Content-Length") or "")
        except ValueError:
            return None, (400, _error_body("Content-Length がありません。"))
        if length > MAX_BODY_BYTES:
            # 読まずに閉じると Windows では RST になり、ページに 413 が届かない。上限の 4 倍までは
            # 読み捨ててから答え、それを超える申告は読まずに閉じる。
            self.close_connection = True
            if length <= 4 * MAX_BODY_BYTES:
                remaining = length
                while remaining > 0:
                    chunk = self.rfile.read(min(remaining, 65536))
                    if not chunk:
                        break
                    remaining -= len(chunk)
            return None, (413, _error_body(f"本文が上限 {MAX_BODY_BYTES} バイトを超えています。"))
        try:
            body = json.loads(self.rfile.read(length).decode("utf-8"))
        except (UnicodeDecodeError, ValueError):
            return None, (400, _error_body("本文が JSON ではありません。"))
        if not isinstance(body, dict):
            return None, (400, _error_body("本文は JSON オブジェクトにしてください。"))
        return body, None

    def do_POST(self):
        if not self._host_ok():
            return self._reply(403, _error_body("Host が不正です。"))
        route = _ROUTES.get(self.path)
        if route is None:
            return self._reply(404, _error_body(f"不明なパスです: {self.path}"))
        body, problem = self._read_json()
        if problem:
            return self._reply(*problem)
        token = body.get("token")
        entry = _touch(token) if isinstance(token, str) else None
        if entry is None:
            return self._reply(401, _error_body("トークンが不明か期限切れです。送信用テキストをコピーしてチャットに貼ってください。"))
        try:
            code, payload = route(self, token, entry, body)
        except Exception as exc:     # noqa: BLE001 — 接続を落とさず理由をページへ返す
            code, payload = 500, _error_body(f"{type(exc).__name__}: {exc}")
        self._reply(code, payload)
        if route is _finish:         # 応答を返した後で外す（最後の 1 件なら別スレッドで止まる）
            unregister(token, wait=False)
```

注意: `_submit` は `submission.submit_flags` を**モジュール属性として**引く（`from ... import submit_flags` にしない）。500 のテストが `monkeypatch.setattr(submission, "submit_flags", ...)` で差し替えるため。

- [ ] **Step 4: テストを走らせる**

Run: `C:/Python314/python.exe -m pytest tests/test_curation_submit_server.py tests/test_curation_submission.py -q`
Expected: 全 PASS

- [ ] **Step 5: コミット**（バックグラウンド、ログを `<scratchpad>/commit-task4.log`）

```bash
git add metabolomix/curation/submit_server.py tests/test_curation_submit_server.py
git commit -m "feat(curation): 受け口に submit / finish と Host・上限・トークンの検査を足す

Co-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>"
```

---

### Task 5: レビュー・候補付けで受け口に登録し、HTML にだけ埋め込む

**Files:**
- Modify: `metabolomix/curation/viewer.py`
- Modify: `metabolomix/curation/viewer_common.js`（先頭に `SUBMIT_ENDPOINT` の宣言だけ）
- Modify: `metabolomix/curation/review.py`（`save_review`）、`metabolomix/curation/suggest.py`（`save_suggestion`）
- Modify: `metabolomix/tools/curation_tools.py`（`curation_review` / `curation_suggest` / docstring / `_open_submit`）
- Test: `tests/test_curation_tools.py`、`tests/test_curation_review.py`

**Interfaces:**
- Consumes: `submit_server.register`（Task 3）
- Produces:
  - `viewer.render_html(review, submit_endpoint=None)`、`viewer.render_suggest_html(suggestion, submit_endpoint=None)`（`submit_endpoint` は `{"port", "token", ...}`。埋め込むのは `port` と `token` だけ）
  - `review.save_review(review, submit_endpoint=None)`、`suggest.save_suggestion(s, submit_endpoint=None)`
  - JS グローバル `const SUBMIT_ENDPOINT = {port, token} | null;`（`viewer_common.js` の先頭。Task 6 / 7 が使う）
  - ツール戻り値の `submit`: `{"via": "viewer", "idle_timeout_min": 30}`、受け口を立てられなかったときは `{"via": "copy", "reason": "..."}`

- [ ] **Step 1: 失敗するテストを書く**

`tests/test_curation_tools.py` の末尾に:

```python
import re
from pathlib import Path

from metabolomix.curation import submit_server


def _embedded_endpoint(html_path):
    html = Path(html_path).read_text(encoding="utf-8")
    match = re.search(r"const SUBMIT_ENDPOINT = (\{[^;]*\});", html)
    assert match, "SUBMIT_ENDPOINT が埋め込まれていない"
    return json.loads(match.group(1))


def test_review_registers_a_submit_endpoint_and_embeds_it_only_in_the_html(ready):
    body = json.loads(curation_tools.curation_review())
    assert body["submit"] == {"via": "viewer", "idle_timeout_min": 30}
    endpoint = _embedded_endpoint(body["html_path"])
    assert set(endpoint) == {"port", "token"}
    assert endpoint["port"] == submit_server.current_port()
    json_text = Path(body["html_path"]).with_suffix(".json").read_text(encoding="utf-8")
    assert endpoint["token"] not in json_text
    assert endpoint["token"] not in json.dumps(body)          # LLM の文脈にも出さない


def test_suggest_registers_its_own_token_on_the_shared_endpoint(suggest_env):
    arf2 = str(suggest_env["arf2"])
    review_ep = _embedded_endpoint(json.loads(curation_tools.curation_review(file_path=arf2))["html_path"])
    out = json.loads(curation_tools.curation_suggest(file_path=arf2))
    assert out["submit"]["via"] == "viewer"
    suggest_ep = _embedded_endpoint(out["html_path"])
    assert suggest_ep["port"] == review_ep["port"] and suggest_ep["token"] != review_ep["token"]


def test_review_falls_back_to_copy_when_the_endpoint_cannot_start(ready, monkeypatch):
    def refuse(*args, **kwargs):
        raise OSError("no port")
    monkeypatch.setattr(submit_server, "register", refuse)
    body = json.loads(curation_tools.curation_review())
    assert body["submit"]["via"] == "copy" and "no port" in body["submit"]["reason"]
    assert "const SUBMIT_ENDPOINT = null;" in Path(body["html_path"]).read_text(encoding="utf-8")
```

`tests/test_curation_review.py` の `test_app_template_has_no_embedded_data` の下に:

```python
def test_app_template_has_no_submit_endpoint():
    assert "const SUBMIT_ENDPOINT = null;" in viewer.render_html(None)
    assert "const SUBMIT_ENDPOINT = null;" in viewer.render_suggest_html(None)


def test_submit_endpoint_is_embedded_before_the_data(built):
    _, result = built
    html = viewer.render_html({**result, "warnings": ["/*__SUBMIT_ENDPOINT__*/null"]},
                              {"port": 5, "token": "t", "idle_timeout_min": 30})
    assert 'const SUBMIT_ENDPOINT = {"port":5,"token":"t"};' in html
    assert "/*__SUBMIT_ENDPOINT__*/null" in html            # データの中の同じ文字列は書き換えない
```

- [ ] **Step 2: 失敗を確かめる**

Run: `C:/Python314/python.exe -m pytest tests/test_curation_tools.py tests/test_curation_review.py -q -k "submit_endpoint or registers or falls_back"`
Expected: FAIL（`KeyError: 'submit'` / `SUBMIT_ENDPOINT` が無い）

- [ ] **Step 3: `viewer_common.js` の先頭に宣言を足す**

```js
// ビューアの受け口（curation_review / curation_suggest が立てる 127.0.0.1 の HTTP）。MCP Apps と
// 受け口を立てられなかったときは null（Copy / Send だけ）。spec 2026-10-09。
const SUBMIT_ENDPOINT = /*__SUBMIT_ENDPOINT__*/null;
```

- [ ] **Step 4: `viewer.py` を直す**

モジュール docstring の置換順序の段落に「`__SUBMIT_ENDPOINT__` もデータより**前**に入れる」を足し、次のように変える:

```python
_SUBMIT_PLACEHOLDER = "/*__SUBMIT_ENDPOINT__*/null"


def _endpoint(submit_endpoint) -> dict | None:
    """HTML に入れるのは port と token だけ（期限などは入れない）。"""
    if submit_endpoint is None:
        return None
    return {"port": int(submit_endpoint["port"]), "token": str(submit_endpoint["token"])}


def _render(template_path: Path, data, submit_endpoint=None) -> str:
    template = template_path.read_text(encoding="utf-8")
    template = template.replace(_COMMON_PLACEHOLDER, _COMMON_JS.read_text(encoding="utf-8"), 1)
    template = template.replace("__SUBMISSION_PREFIX__", SUBMISSION_PREFIX.strip())
    template = template.replace(_SUBMIT_PLACEHOLDER, _embed(_endpoint(submit_endpoint)), 1)
    return template.replace(_PLACEHOLDER, _embed(data), 1)


def render_html(review: dict | None, submit_endpoint: dict | None = None) -> str:
    return _render(_TEMPLATE, review, submit_endpoint)


def render_suggest_html(suggestion: dict | None, submit_endpoint: dict | None = None) -> str:
    return _render(_SUGGEST_TEMPLATE, suggestion, submit_endpoint)
```

- [ ] **Step 5: `save_review` / `save_suggestion` に引数を通す**

```python
def save_review(review: dict, submit_endpoint: dict | None = None) -> dict:
    ...
    html_path.write_text(viewer.render_html(review, submit_endpoint), encoding="utf-8")
```

```python
def save_suggestion(s: dict, submit_endpoint: dict | None = None) -> dict:
    ...
    html_path.write_text(viewer.render_suggest_html(s, submit_endpoint), encoding="utf-8")
```

- [ ] **Step 6: `curation_tools.py` で登録する**

import に `submit_server` を足す（`from metabolomix.curation import ..., submission, submit_server, suggest, viewer`）。ヘルパを足す:

```python
def _open_submit(kind: str, review_id: str, arf2_path) -> tuple[dict | None, dict]:
    """ビューアの受け口に登録する。立てられなければ Copy の経路だけ（記録は curation_submit）。"""
    try:
        endpoint = submit_server.register(kind, review_id, arf2_path)
    except OSError as exc:
        return None, {"via": "copy", "reason": f"受け口を立てられませんでした: {exc}"}
    return endpoint, {"via": "viewer", "idle_timeout_min": endpoint["idle_timeout_min"]}
```

`curation_review` の `saved = review.save_review(result)` を:

```python
    endpoint, submit_info = _open_submit("review", result["review_id"], result["arf2_path"])
    saved = review.save_review(result, submit_endpoint=endpoint)
```

戻り値の dict に `"submit": submit_info,` を `"html_path"` の次に足す。`curation_suggest` も同様に `_open_submit("suggest", result["suggestion_id"], result["arf2_path"])` → `suggest.save_suggestion(result, submit_endpoint=endpoint)`、戻り値に `"submit": submit_info`。

docstring（spec §7）:
- `curation_review` の「ユーザーには `html_path` をブラウザで開いてもらい、ビューアで付けたフラグを「送信用テキストをコピー」→ チャットに貼ってもらう。貼られたら curation_submit(submission_text=...) に渡す。」を次に置き換える:
  「ユーザーには `html_path` をブラウザで開いてもらい、ビューアの **Submit** で送ってもらう（記録と `_tags.xml` の反映はビューアが直接行う。`submit.via == "viewer"`）。送信用テキストを貼るよう頼まない。送った結果は curation_flags で確かめられる。Copy のテキストが貼られたとき（受け口の時間切れ・`submit.via == "copy"`）は従来どおり curation_submit(submission_text=...) に渡す。」
- `curation_suggest` の対応する 2 文も同じ趣旨に置き換える（「選んだ内容をビューアの **送信** で送ってもらう …」）。「**ユーザーの同意なしに curation_submit を呼ばない。**」は残す。

- [ ] **Step 7: テストを走らせる**

Run: `C:/Python314/python.exe -m pytest tests/test_curation_tools.py tests/test_curation_review.py tests/test_curation_suggest.py tests/test_curation_suggest_viewer.py -q`
Expected: 全 PASS

- [ ] **Step 8: コミット**（バックグラウンド、ログを `<scratchpad>/commit-task5.log`）

```bash
git add metabolomix/curation/viewer.py metabolomix/curation/viewer_common.js metabolomix/curation/review.py metabolomix/curation/suggest.py metabolomix/tools/curation_tools.py tests/test_curation_tools.py tests/test_curation_review.py
git commit -m "feat(curation): レビューと候補付けで受け口に登録し、ポートとトークンを HTML にだけ埋め込む

Co-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>"
```

---

### Task 6: 送信クライアント（共通）とレビュービューアの Submit / Finish

**Files:**
- Modify: `metabolomix/curation/viewer_common.js`
- Modify: `metabolomix/curation/viewer.html`
- Test: `tests/test_curation_review.py`

**Interfaces:**
- Consumes: `SUBMIT_ENDPOINT`（Task 5）、HTTP 仕様（Task 4）
- Produces（`viewer_common.js`。Task 7 が使う）:
  - 純関数ブロック `// --- submit client (pure) ---` 〜 `// --- end submit client ---`:
    `countFlags(flags) -> {flag: n}`、`reviewConfirmText(flags) -> string`、`suggestConfirmText(flags) -> string`、
    `stateAfter(current, httpStatus) -> "unavailable" | "live" | "lost" | "finished"`、
    `submitControls(state) -> {submit: bool, finish: bool, copyPrimary: bool}`、
    `tagsResultText(body, lang) -> string`、`lostText(lang)`、`finishConfirmText(n, lang)`、`finishedText(lang)`（`lang` は `"en"` / `"ja"`）
  - 非純関数: `postEndpoint(path, extra) -> Promise<{status, body}>`（届かなければ `status: 0`）、
    `startSubmitClient({onState}) -> {state(), submit(payload), finish()}`、`showSubmitControls(state)`、`openFallback(text, message)`
  - `viewer.html` の純関数ブロック `// --- recorded (pure) ---`: `applyRecorded(spots, edits)`
  - 両ビューアのヘッダに `#submit` `#finish` `#result` が要る（Task 7 も同じ ID を使う）

- [ ] **Step 1: 失敗するテストを書く**（`tests/test_curation_review.py` の末尾。既存の `_run_block` を使う）

```python
def test_review_confirm_text_counts_each_flag(tmp_path):
    text = _run_block(tmp_path, "submit client", "", "reviewConfirmText([{flag:'wrong'},{flag:'wrong'},{flag:'confirmed'},{flag:'clear'}])")
    assert "Wrong 2 / Suspect 0 / Confirmed 1 / clear 1" in text and "MS-DIAL" in text


def test_suggest_confirm_text_warns_about_ms_dial_only_with_clear(tmp_path):
    without = _run_block(tmp_path, "submit client", "", "suggestConfirmText([{flag:'assign'},{flag:'redundant'}])")
    assert "assign 1 / redundant 1 / clear 0" in without and "MS-DIAL" not in without
    with_clear = _run_block(tmp_path, "submit client", "", "suggestConfirmText([{flag:'clear'}])")
    assert "clear 1" in with_clear and "MS-DIAL" in with_clear


def test_state_after_switches_to_copy_when_the_endpoint_is_lost(tmp_path):
    cases = _run_block(tmp_path, "submit client", "", "[" + ",".join([
        "stateAfter('unavailable', 200)", "stateAfter('unavailable', 0)", "stateAfter('live', 401)",
        "stateAfter('live', 0)", "stateAfter('live', 400)", "stateAfter('live', 409)",
        "stateAfter('lost', 200)", "stateAfter('finished', 200)"]) + "]")
    assert cases == ["live", "unavailable", "lost", "lost", "live", "live", "lost", "finished"]


def test_submit_controls_show_submit_and_finish_only_when_live(tmp_path):
    out = _run_block(tmp_path, "submit client", "", "['live','unavailable','lost','finished'].map(submitControls)")
    assert out[0] == {"submit": True, "finish": True, "copyPrimary": False}
    assert all(c == {"submit": False, "finish": False, "copyPrimary": True} for c in out[1:])


def test_tags_result_text_reports_changes_and_failures(tmp_path):
    ok = _run_block(tmp_path, "submit client", "", "tagsResultText({recorded:2, tags_xml:{added:[1], removed:[], confirmed:{added:[0], removed:[]}, note:'N'}}, 'en')")
    assert ok.startswith("Recorded 2.") and "Misannotation +1 / -0" in ok and "Confirmed +1 / -0" in ok and ok.endswith("N")
    failed = _run_block(tmp_path, "submit client", "", "tagsResultText({recorded:1, tags_xml:{error:'OSError: x'}}, 'ja')")
    assert "1 件を記録しました" in failed and "OSError: x" in failed


def test_apply_recorded_moves_edits_onto_spots(tmp_path):
    spots = json.dumps([{"spot_id": 0, "flag": None}, {"spot_id": 1, "flag": "wrong", "flag_note": "a"}])
    out = _run_block(tmp_path, "recorded", f"const SPOTS = {spots};",
                     "(applyRecorded(SPOTS, new Map([[0, {flag:'confirmed', note:'ok'}], [1, {flag:'', note:''}]])), SPOTS)")
    assert out[0] == {"spot_id": 0, "flag": "confirmed", "flag_note": "ok", "flag_cleared": False, "confirmed": True}
    assert out[1] == {"spot_id": 1, "flag": None, "flag_note": None, "flag_cleared": True, "confirmed": False}


def test_review_submit_button_skips_empty_submissions_and_confirms_first():
    script = _script()
    handler = script[script.index('getElementById("submit").addEventListener'):]
    handler = handler[:handler.index("\n});")]
    assert handler.index("if (!edits.size)") < handler.index("confirm(reviewConfirmText(")
    assert handler.index("confirm(reviewConfirmText(") < handler.index("submitClient.submit(")
    assert "applyRecorded" in script[script.index("function afterRecorded"):]
```

既存の `test_app_loader_shows_error_payloads_and_send_updates_spots` を、Send が共通関数を使う形に直す:

```python
def test_app_loader_shows_error_payloads_and_send_updates_spots():
    script = _script()
    assert "showError" in script
    send = script[script.index('getElementById("send").addEventListener'):]
    send = send[:send.index("\n});")]
    assert "afterRecorded(" in send
    after = script[script.index("function afterRecorded"):]
    after = after[:after.index("\n}")]
    assert after.index("applyRecorded(REVIEW.spots, edits)") < after.index("edits.clear()")
```

- [ ] **Step 2: 失敗を確かめる**

Run: `C:/Python314/python.exe -m pytest tests/test_curation_review.py -q`
Expected: 新しいテストが FAIL（`// --- submit client (pure) ---` が無い等）

- [ ] **Step 3: `viewer_common.js` の末尾に送信クライアントを足す**

```js
// --- submit client (pure) ---
// ビューアの Submit / Finish（spec 2026-10-09）。受け口の状態: unavailable（受け口なし・未接続）/
// live（ping が通った）/ lost（live の後で 401 か届かない＝時間切れ・サーバ再起動）/ finished。
function countFlags(flags) { const n = {}; for (const f of flags) n[f.flag] = (n[f.flag] || 0) + 1; return n; }
function reviewConfirmText(flags) {
  const n = countFlags(flags);
  return `Record Wrong ${n.wrong || 0} / Suspect ${n.suspect || 0} / Confirmed ${n.confirmed || 0} / clear ${n.clear || 0} and update _tags.xml.\n` +
    "If this project is open in MS-DIAL, close it first.";
}
// assign / redundant は _tags.xml を変えない。clear（元の注釈に戻す）だけが Misannotation と Confirmed を外す。
function suggestConfirmText(flags) {
  const n = countFlags(flags);
  const head = `assign ${n.assign || 0} / redundant ${n.redundant || 0} / clear ${n.clear || 0} を記録します。`;
  return n.clear ? head + "\nclear は _tags.xml の Misannotation と Confirmed を外します。" +
    "MS-DIAL でこのプロジェクトを開いているなら先に閉じてください。" : head;
}
function stateAfter(current, httpStatus) {
  if (current === "finished" || current === "lost") return current;
  if (httpStatus === 401 || httpStatus === 0) return current === "live" ? "lost" : "unavailable";
  if (httpStatus === 200) return "live";
  return current;             // 400 / 409 などは接続の状態を変えない
}
function submitControls(state) {
  const live = state === "live";
  return {submit: live, finish: live, copyPrimary: !live};
}
function tagsResultText(body, lang) {
  const t = body.tags_xml || {}, c = t.confirmed || {}, n = body.recorded, len = a => (a || []).length;
  if (t.error) return lang === "ja" ? `${n} 件を記録しました。_tags.xml への反映に失敗: ${t.error}`
                                    : `Recorded ${n}. Updating _tags.xml failed: ${t.error}`;
  const tags = `Misannotation +${len(t.added)} / -${len(t.removed)}, Confirmed +${len(c.added)} / -${len(c.removed)}`;
  return (lang === "ja" ? `${n} 件を記録しました。${tags}。` : `Recorded ${n}. ${tags}. `) + (t.note || "");
}
function lostText(lang) {
  return lang === "ja" ? "送信の受け口がありません（時間切れかサーバの再起動）。下の欄をコピーしてチャットに貼ってください。未送信の選択は残っています。"
    : "The submit endpoint is gone (timed out or the server restarted). Copy the text below and paste it into the chat; your unsent changes are kept.";
}
function finishConfirmText(n, lang) {
  return lang === "ja" ? `${n} 件の未送信の変更があります。終了してよいですか（後から「送信用テキストをコピー」で送れます）。`
    : `${n} unsent change(s). Finish anyway? (You can still send them later with Copy submission text.)`;
}
function finishedText(lang) {
  return lang === "ja" ? "終了しました。以降は「送信用テキストをコピー」で送れます。"
    : "Finished. Use Copy submission text to send anything else.";
}
// --- end submit client ---

const HEARTBEAT_MS = 5 * 60 * 1000;

// 本文は text/plain で送る（事前確認の要らない単純要求。中身は JSON）。届かなければ status 0。
async function postEndpoint(path, extra) {
  if (!SUBMIT_ENDPOINT) return {status: 0, body: null};
  try {
    const res = await fetch(`http://127.0.0.1:${SUBMIT_ENDPOINT.port}${path}`, {method: "POST",
      headers: {"Content-Type": "text/plain;charset=UTF-8"},
      body: JSON.stringify({token: SUBMIT_ENDPOINT.token, ...extra})});
    let body = null; try { body = await res.json(); } catch { /* 本文なし */ }
    return {status: res.status, body};
  } catch { return {status: 0, body: null}; }
}

function startSubmitClient(hooks) {
  let state = "unavailable";
  const update = status => { const next = stateAfter(state, status);
    if (next !== state) { state = next; hooks.onState(state); } };
  const ping = async () => { if (state === "lost" || state === "finished") return;
    update((await postEndpoint("/v1/ping", {})).status); };
  if (SUBMIT_ENDPOINT) {
    ping(); setInterval(ping, HEARTBEAT_MS);
    document.addEventListener("visibilitychange", () => { if (document.visibilityState === "visible") ping(); });
  }
  return {
    state: () => state,
    async submit(payload) { const r = await postEndpoint("/v1/submit", payload); update(r.status); return r; },
    async finish() { const r = await postEndpoint("/v1/finish", {}); state = "finished"; hooks.onState(state); return r; },
  };
}

function showSubmitControls(state) {
  const c = submitControls(state);
  document.getElementById("submit").hidden = !c.submit;
  document.getElementById("finish").hidden = !c.finish;
  const copy = document.getElementById("copy");
  copy.classList.toggle("primary", c.copyPrimary); copy.classList.toggle("quiet", !c.copyPrimary);
}

function openFallback(text, message) {
  const area = document.getElementById("fallback");
  area.hidden = false; area.value = text; area.select();
  document.getElementById("status").textContent = message;
}
```

- [ ] **Step 4: `viewer.html` を直す**

CSS（`button.primary` の行の下）:

```css
button.quiet { opacity:.7; }
#result { color:var(--muted); overflow-wrap:anywhere; }
```

ヘッダのボタン列を次にする:

```html
  <button id="send" class="primary" hidden>Send</button>
  <button id="submit" class="primary" hidden>Submit</button>
  <button id="finish" hidden>Finish</button>
  <button id="copy">Copy submission text</button>
  <span id="status"></span>
  <span id="result"></span>
```

`function submission()` の直前に純関数ブロックを足す:

```js
// --- recorded (pure) ---
// 送信が記録されたら、未送信の変更をスポットへ写す（MCP Apps の Send とビューアの Submit で共通）。
function applyRecorded(spots, edits) {
  for (const [spotId, e] of edits) {
    const spot = spots.find(s => s.spot_id === spotId); if (!spot) continue;
    spot.flag = e.flag || null; spot.flag_note = e.flag ? (e.note || "") : null; spot.flag_cleared = !e.flag;
    spot.confirmed = e.flag === "confirmed";
  }
}
// --- end recorded ---
```

Copy ハンドラの直後に、結果表示と Submit / Finish の配線を足す（`render()` は状態欄を後から書き直すので、結果は `#result` に出す）:

```js
function afterRecorded(text) {
  applyRecorded(REVIEW.spots, edits);
  edits.clear(); presetCount = 0; presetCauses = {}; render();
  document.getElementById("result").textContent = text;
}

const submitClient = startSubmitClient({onState: state => {
  showSubmitControls(state);
  if (state === "lost") openFallback(PREFIX + " " + JSON.stringify(submission()), lostText("en"));
}});

document.getElementById("submit").addEventListener("click", async () => {
  if (!edits.size) { document.getElementById("status").textContent = "No unsent changes."; return; }
  const payload = submission();
  if (!confirm(reviewConfirmText(payload.flags))) return;
  document.getElementById("result").textContent = "Submitting…";
  const r = await submitClient.submit(payload);
  if (r.status === 200) { afterRecorded(tagsResultText(r.body, "en")); return; }
  document.getElementById("result").textContent = "";
  if (r.status !== 401 && r.status !== 0) showError(errorMessage(r.body) || `HTTP ${r.status}`);   // 未送信の変更は残す
});

document.getElementById("finish").addEventListener("click", async () => {
  if (edits.size && !confirm(finishConfirmText(edits.size, "en"))) return;
  await submitClient.finish();
  document.getElementById("result").textContent = finishedText("en");
});
```

MCP Apps の Send ハンドラの末尾（`for (const [spotId, e] of edits) { ... }` から `status` の書き込みまで）を次に置き換える:

```js
  afterRecorded(text.slice(0, 200));
```

注意: `showError` / `errorMessage` は下の MCP Apps 節で関数宣言されている（巻き上げで参照できる）。`submitClient` の宣言は `if (REVIEW) start();` より前に置く。

- [ ] **Step 5: テストを走らせる**

Run: `C:/Python314/python.exe -m pytest tests/test_curation_review.py tests/test_curation_suggest_viewer.py -q`
Expected: 全 PASS（`test_extracted_script_is_valid_javascript` も含む）

- [ ] **Step 6: コミット**（バックグラウンド、ログを `<scratchpad>/commit-task6.log`）

```bash
git add metabolomix/curation/viewer_common.js metabolomix/curation/viewer.html tests/test_curation_review.py
git commit -m "feat(curation): レビューのビューアに Submit / Finish と受け口の状態表示を足す

Co-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>"
```

---

### Task 7: 候補付けビューアの 送信 / 終了

**Files:**
- Modify: `metabolomix/curation/suggest_viewer.html`
- Test: `tests/test_curation_suggest_viewer.py`

**Interfaces:**
- Consumes: Task 6 の `startSubmitClient` / `showSubmitControls` / `openFallback` / `suggestConfirmText` / `tagsResultText` / `lostText` / `finishConfirmText` / `finishedText`、`#submit` `#finish` `#result` の ID
- Produces: `// --- suggest (pure) ---` ブロックに `applySuggestRecorded(spots, choiceMap)` と `recordedText(rec) -> string`

- [ ] **Step 1: 失敗するテストを書く**（`tests/test_curation_suggest_viewer.py` の末尾）

```python
def test_apply_suggest_recorded_marks_spots(tmp_path):
    out = _run(tmp_path, "(applySuggestRecorded(SPOTS, new Map([[1, {candidate: 'L1', flag: 'assign', level: 'species', note: ''}], [2, {flag: 'clear', note: ''}]])), SPOTS.map(s => s.recorded))",
               prelude="const SPOTS = [{spot_id: 1}, {spot_id: 2}, {spot_id: 3}];")
    assert out == [{"flag": "assign", "candidate": "L1", "level": "species"}, {"flag": "clear"}, None]


def test_recorded_text_names_the_decision(tmp_path):
    out = _run(tmp_path, "[recordedText({flag:'assign', candidate:'L1', level:'sum'}), recordedText({flag:'redundant', candidate:'R1'}), recordedText({flag:'clear'})]")
    assert out[0].startswith("記録済み") and "L1" in out[0]
    assert "R1" in out[1] and "元の注釈" in out[2]


def _script():
    html = viewer.render_suggest_html(None)
    return html[html.index("<script>"):html.index("</script>")]


def test_suggest_viewer_has_submit_and_finish_wired():
    html = viewer.render_suggest_html(None)
    for element in ('id="submit"', 'id="finish"', 'id="result"'):
        assert element in html
    script = _script()
    handler = script[script.index('getElementById("submit").addEventListener'):]
    handler = handler[:handler.index("\n});")]
    assert handler.index("if (!choices.size)") < handler.index("confirm(suggestConfirmText(")
    assert handler.index("confirm(suggestConfirmText(") < handler.index("submitClient.submit(")


def test_suggest_extracted_script_is_valid_javascript(tmp_path):
    node = shutil.which("node")
    if not node:
        pytest.skip("node が無いので構文チェックを省略")
    path = tmp_path / "suggest.js"
    path.write_text(_script()[len("<script>"):], encoding="utf-8")
    result = subprocess.run([node, "--check", str(path)], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
```

注意: `_run` は `undefined` を JSON にできないので、`SPOTS.map(s => s.recorded)` で `recorded` の無いスポットは `JSON.stringify` で `null` になる（配列要素の `undefined` は `null`）。

- [ ] **Step 2: 失敗を確かめる**

Run: `C:/Python314/python.exe -m pytest tests/test_curation_suggest_viewer.py -q`
Expected: 新しいテストが FAIL

- [ ] **Step 3: `suggest_viewer.html` を直す**

CSS（`button.primary` の行の下）:

```css
button.quiet { opacity:.7; }
#result { color:var(--muted); overflow-wrap:anywhere; }
```

ヘッダ:

```html
  <button id="submit" class="primary" hidden>送信</button>
  <button id="finish" hidden>終了</button>
  <button id="copy" class="primary">送信用テキストをコピー</button>
  <span id="status"></span>
  <span id="result"></span>
```

`// --- end suggest ---` の直前に足す:

```js
// 送信が記録されたら、選択をスポットへ写す（カードに「記録済み」を出す）。
function applySuggestRecorded(spots, choiceMap) {
  for (const [id, c] of choiceMap) {
    const spot = spots.find(s => s.spot_id === id); if (!spot) continue;
    spot.recorded = c.flag === "clear" ? {flag: "clear"}
      : {flag: c.flag, candidate: c.candidate, level: c.flag === "assign" ? (c.level || "sum") : undefined};
    if (spot.recorded.level === undefined) delete spot.recorded.level;
  }
}
function recordedText(rec) {
  if (rec.flag === "clear") return "記録済み: 元の注釈に戻した";
  if (rec.flag === "assign") return `記録済み: 候補 ${rec.candidate} を採用（${rec.level === "species" ? "分子種" : "和組成"}）`;
  return `記録済み: 別イオン ${rec.candidate}`;
}
```

`spotCard` の最初の `card.appendChild(el("div", {class: "meta"}, ...))`（`#${spot.spot_id} ...` の行）の直後に:

```js
  if (spot.recorded) card.appendChild(el("div", {class: "meta"}, recordedText(spot.recorded)));
```

Copy ハンドラの直後に:

```js
function afterRecorded(text) {
  applySuggestRecorded(DATA.spots, choices); choices.clear(); render();
  document.getElementById("result").textContent = text;
}

const submitClient = startSubmitClient({onState: state => {
  showSubmitControls(state);
  if (state === "lost") openFallback(PREFIX + " " + JSON.stringify(suggestSubmission(DATA, choices)), lostText("ja"));
}});

document.getElementById("submit").addEventListener("click", async () => {
  if (!choices.size) { document.getElementById("status").textContent = "未送信の選択はありません。"; return; }
  const payload = suggestSubmission(DATA, choices);
  if (!confirm(suggestConfirmText(payload.flags))) return;
  document.getElementById("result").textContent = "送信中…";
  const r = await submitClient.submit(payload);
  if (r.status === 200) { afterRecorded(tagsResultText(r.body, "ja")); return; }
  document.getElementById("result").textContent = "";
  if (r.status !== 401 && r.status !== 0)          // 未送信の選択は残す
    document.getElementById("status").textContent = "エラー: " + String((r.body && r.body.message) || `HTTP ${r.status}`).slice(0, 300);
});

document.getElementById("finish").addEventListener("click", async () => {
  if (choices.size && !confirm(finishConfirmText(choices.size, "ja"))) return;
  await submitClient.finish();
  document.getElementById("result").textContent = finishedText("ja");
});
```

注意: `submitClient` の宣言は末尾の `if (DATA) start(); ...` より前に置く。

- [ ] **Step 4: テストを走らせる**

Run: `C:/Python314/python.exe -m pytest tests/test_curation_suggest_viewer.py tests/test_curation_review.py tests/test_curation_tools.py -q`
Expected: 全 PASS

- [ ] **Step 5: コミット**（バックグラウンド、ログを `<scratchpad>/commit-task7.log`）

```bash
git add metabolomix/curation/suggest_viewer.html tests/test_curation_suggest_viewer.py
git commit -m "feat(curation): 候補付けのビューアにも 送信 / 終了 を足す

Co-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>"
```

---

### Task 8: 文書・記録・流れ図

**Files:**
- Modify: `USAGE.md`、`docs/output_format/curation.md`、`docs/workflow/curation.md`
- Modify（追跡外・main ツリー）: `C:\Users\yuu18\Metabolomix_with_LLM\docs\HISTRY.md`、`C:\Users\yuu18\Metabolomix_with_LLM\docs\task.md`
- Modify（リポジトリ外）: `C:\Users\yuu18\Documents\KnowledgeVault\30_Projects\ms-data-parser\ms-data-parser-flow.md`

**Interfaces:**
- Consumes: Task 2〜7 の関数名・戻り値（`submit`）・HTTP 仕様
- Produces: なし

- [ ] **Step 1: `USAGE.md`**

221 行目付近の「ユーザーがビューアでフラグを付けて「送信用テキストをコピー」→ チャットに貼る → `curation_submit`。」を次にする:
「ユーザーがビューアでフラグを付けて **Submit** を押すと、ビューアが直接記録し `_tags.xml` に反映する（LLM を介さない。受け口は `127.0.0.1` で、Finish か 30 分の無通信で閉じる）。受け口が無いとき（時間切れ・サーバ再起動）は「Copy submission text」→ チャットに貼る → `curation_submit`。」
`curation_review` / `curation_suggest` の行に「戻り値の `submit`（`via: viewer` ならビューアの Submit で送れる。`via: copy` なら貼り付け）」を足し、`curation_suggest` の「ビューアで選んで送信用テキストを貼る → `curation_submit`」を「ビューアで選んで **送信**（受け口が無ければ送信用テキストを貼る → `curation_submit`）」に直す。

- [ ] **Step 2: `docs/output_format/curation.md`**

- `curation_review` の戻り値表と `curation_suggest` の戻り値表に行を足す:
  `| submit | ビューアから直接送れるか。`{"via": "viewer", "idle_timeout_min": 30}` ならビューアの Submit（候補付けは「送信」）が記録と `_tags.xml` の反映を行う。`{"via": "copy", "reason"}` なら受け口を立てられず、送信用テキストを貼る経路だけ |`
- `html_path`（候補付け、232 行目付近）の「「送信用テキストをコピー」でチャットへ貼る」を「「送信」で直接送る（受け口が無いときは「送信用テキストをコピー」でチャットへ貼る）」に。
- `### curation_submit の戻り値の tags_xml` の節の後に、新しい節を足す:

```markdown
### ビューアからの直接送信（2026-10-09）

`curation_review` / `curation_suggest` は、`127.0.0.1` の空きポートに受け口を立て（プロセスで 1 つを共有）、
レビューごとのトークンを HTML にだけ埋め込む（レビュー・候補付けの JSON とツールの戻り値には入れない）。
ビューアの Submit は `curation_submit` と同じ検証・記録・`_tags.xml` 反映を行い、戻り値も同じ形。`source` は常に `user`。

- 受け口は Finish か、最後の通信（ping は 5 分おき）から 30 分で登録を外し、登録が 0 件になれば止まる。MCP サーバが終われば消える。
- 応答の状態コード: 401（トークン不明・期限切れ → ビューアは送信用テキストの欄を開く）/ 403（`review_id` 不一致・`Host` 不正）/
  400（検証エラー・JSON でない）/ 409（レビューの後でアラインメントが変わった）/ 410（保存済みレビューが無い）/
  413（本文が 8 MB 超）/ 500（`flags.jsonl` が読めない・想定外の例外）。本文は `{status: "error", message, ...}`。
- 送った結果は `curation_flags` で確かめられる。
```

- [ ] **Step 3: `docs/workflow/curation.md`**

`## curation_review` の連鎖の `save_review()` の行の直前に（番号は前後に合わせて振り直す）:

```markdown
N. metabolomix/tools/curation_tools.py  _open_submit()
N+1. └─ metabolomix/curation/submit_server.py  register()（受け口が無ければ 127.0.0.1 の空きポートで起動。失敗したら `submit.via = "copy"`）
```

`## curation_suggest` にも同じ 2 行を `save_suggestion()` の直前に足す。状態変更の段落に「受け口の登録表にトークンを足す（プロセス内）」を足す。`## curation_submit` の節の末尾に小見出しを足す:

```markdown
### ビューアからの直接送信（HTTP。MCP ツールではない）

`POST /v1/submit` は登録時の `review_id` / `arf2_path` だけを使い、`source="user"` で記録する。

1. metabolomix/curation/submit_server.py  _Handler.do_POST()（`Host` → パス → 本文の上限と JSON → トークン）
2. └─ metabolomix/curation/submit_server.py  _touch()（定数時間で照合し、期限を延ばす）
3. metabolomix/curation/submit_server.py  _submit()（本文の `review_id` が登録と違えば 403）
4. ├─ metabolomix/curation/submission.py  load_saved()（無ければ 410）
5. └─ metabolomix/curation/submission.py  submit_flags()（`curation_submit` と同じ。`SubmissionError.kind` → 400 / 409 / 500）
6. metabolomix/curation/submit_server.py  unregister()（`/v1/finish` の応答の後。0 件なら別スレッドで止める）
7. metabolomix/curation/submit_server.py  sweep()（daemon スレッドが 60 秒ごと。30 分無通信の登録を外す）
```

mermaid に `CR --> REG[curation.submit_server.register]`、`CG --> REG`、`HTTP[viewer Submit] --> SUB` を足す。

- [ ] **Step 4: 文書テストを走らせる**

Run: `C:/Python314/python.exe -m pytest tests/test_workflow_docs.py tests/test_readme_links.py tests/test_output_format_sections.py -q`
Expected: 全 PASS

- [ ] **Step 5: コミット**（バックグラウンド、ログを `<scratchpad>/commit-task8.log`）

```bash
git add USAGE.md docs/output_format/curation.md docs/workflow/curation.md
git commit -m "docs(curation): ビューアからの直接送信の経路を文書に足す

Co-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>"
```

- [ ] **Step 6: vault の流れ図を直す**

`ms-data-parser-flow.md` のキュレーションの経路で、「ビューア → 送信用テキストをコピー → チャット → curation_submit」を
「ビューア → Submit（直接記録・`_tags.xml` 反映）／受け口が無いとき Copy → チャット → curation_submit」に分岐させ、
Finish と 30 分の時間切れで閉じることを書く。末尾の出典行の日付を 2026-10-09 に直す。

- [ ] **Step 7: main ツリーの記録**

`C:\Users\yuu18\Metabolomix_with_LLM\docs\HISTRY.md` の末尾に `## 2026-10-09 キュレーションのビューアから直接送信` の節を足す（Task 1 の実地確認の結果＝Edge / Chrome で `file://` から届いたか・ダイアログの有無・ファイアウォールの確認、方式の要点、新モジュール）。
`docs/task.md` は並列のエージェントがいれば末尾追記だけ。該当する未完了項目（ビューアからの送信に関するもの）があれば、並列がいないときに消す。

---

### Task 9: 実データでの通し確認（自動化しない）と全テスト

**Files:** なし（確認のみ）

**Interfaces:**
- Consumes: すべて
- Produces: 完了報告

- [ ] **Step 1: 全テスト**

Run: `C:/Python314/python.exe -m pytest tests -q`（バックグラウンド、出力をファイルへ）
Expected: 全 PASS（失敗があれば `FAILED` 行を grep して直す）

- [ ] **Step 2: worktree の MCP サーバで kidney の実データを通す**（ユーザーと一緒に）

`library_load` → `curation_review`（kidney、ontology を絞って少数）→ 戻り値に `submit.via == "viewer"`、トークンが戻り値に無いことを確かめる → `html_path` を Edge で開く → Submit と Finish が出る → 数件を Wrong / Confirmed にして Submit → 確認ダイアログの内訳が合う → 結果欄に記録件数と `_tags.xml` の変化 → `curation_flags` で `source=user` の行を確かめる → MS-DIAL でプロジェクトを開き直して Misannotation / Confirmed が見える → もう一度 Submit（分けて送れる）→ Finish（未送信があれば確認が出る）→ Submit が消え Copy だけになる。

`curation/` の既存の記録は消さない（`flags.jsonl` は追記のみ。試しに付けたフラグは clear で取り消す）。

- [ ] **Step 3: 時間切れ・再起動の確認**

別のレビュー HTML を開いたまま MCP サーバを再起動（または 30 分放置）→ タブを表に戻す → Copy 欄が開き、理由が出て、未送信の変更が残っていることを確かめる。`curation_suggest` のビューアでも Submit → `assign` が記録され `_tags.xml` が変わらない（`clear` を含むときだけ MS-DIAL の注意が出る）ことを確かめる。

- [ ] **Step 4: 結果を main ツリーの `docs/HISTRY.md` の同じ日付の節に追記し、ブランチの統合は superpowers:finishing-a-development-branch に従う**（main へは `--no-ff`、`Merge feat/curation-viewer-submit: キュレーションのビューアから直接送信する`。push は都度確認）。
