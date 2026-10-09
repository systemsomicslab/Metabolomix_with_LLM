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
            store.rows(tolerate_partial_tail=False)   # 壊れた記録（書きかけの末尾を含む）に追記しない
        except flag_log.FlagFileError as exc:
            raise SubmissionError(str(exc), kind="flag_file", details=exc.details()) from exc
        n = store.append(cleaned, alignment=current, review_id=review_id, source=source)
        effective = store.effective(current["alignment_sha256"])
        tags_xml = msdial_writeback.sync_tags(saved["arf2_path"], cleaned)
    return {"status": "ok", "recorded": n, "review_id": review_id,
            "n_wrong": _count(effective, "wrong"), "n_suspect": _count(effective, "suspect"),
            "n_confirmed": _count(effective, "confirmed"), "n_assign": _count(effective, "assign"),
            "n_redundant": _count(effective, "redundant"), "tags_xml": tags_xml}
