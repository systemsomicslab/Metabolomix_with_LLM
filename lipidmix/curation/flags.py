"""キュレーションのフラグ記録（追記専用 JSON Lines）と送信用テキスト。spec §6。

記録するのはユーザーが付けたフラグだけ（wrong / suspect、取り消しは clear）。
無印のスポットは「間違っていない」で何も書かない。キーはアラインメントファイルの
sha256 と MasterAlignmentID の組で、MS-DIAL を再実行して `.arf2` が作り直されたら
古いフラグは当たらない（新しい ID に黙って当てない）。
"""
from __future__ import annotations

import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path

FLAG_VALUES = ("wrong", "suspect", "clear")
SUBMISSION_PREFIX = "CURATION_SUBMIT "
FLAGS_FILENAME = "flags.jsonl"


def curation_dir(arf2_path) -> Path:
    return Path(arf2_path).parent / "curation"


def alignment_key(arf2_path) -> dict:
    digest = hashlib.sha256()
    with open(arf2_path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return {"alignment_file": Path(arf2_path).name, "alignment_sha256": digest.hexdigest()}


def validate_entries(entries, *, allowed_spot_ids: set[int] | None) -> list[dict]:
    if not isinstance(entries, list) or not entries:
        raise ValueError("flags は 1 件以上のリストで渡してください。")
    cleaned = []
    for index, entry in enumerate(entries):
        if not isinstance(entry, dict):
            raise ValueError(f"flags[{index}] が dict ではありません。")
        spot_id = entry.get("spot_id")
        if isinstance(spot_id, bool) or not isinstance(spot_id, int):
            raise ValueError(f"flags[{index}].spot_id が整数ではありません: {spot_id!r}")
        if allowed_spot_ids is not None and spot_id not in allowed_spot_ids:
            raise ValueError(f"flags[{index}].spot_id={spot_id} はこのレビューの対象外です。")
        flag = entry.get("flag")
        if flag not in FLAG_VALUES:
            raise ValueError(f"flags[{index}].flag={flag!r} は {FLAG_VALUES} のいずれかにしてください。")
        note = entry.get("note") or ""
        cleaned.append({"spot_id": spot_id, "flag": flag, "note": str(note)[:500]})
    return cleaned


class FlagStore:
    def __init__(self, directory):
        self.directory = Path(directory)
        self.path = self.directory / FLAGS_FILENAME

    def append(self, entries, *, alignment: dict, review_id: str, source: str) -> int:
        self.directory.mkdir(parents=True, exist_ok=True)
        now = datetime.now(timezone.utc).isoformat()
        lines = [json.dumps({"ts": now, **alignment, "review_id": review_id,
                             "source": source, **entry}, ensure_ascii=False)
                 for entry in entries]
        with open(self.path, "a", encoding="utf-8", newline="\n") as handle:
            handle.write("\n".join(lines) + "\n")
            handle.flush()
            os.fsync(handle.fileno())
        return len(lines)

    def _rows(self):
        if not self.path.exists():
            return []
        rows = []
        for line in self.path.read_text(encoding="utf-8").splitlines():
            if line.strip():
                rows.append(json.loads(line))
        return rows

    def effective(self, alignment_sha256: str) -> dict[int, dict]:
        latest: dict[int, dict] = {}
        for row in self._rows():
            if row.get("alignment_sha256") == alignment_sha256:
                latest[int(row["spot_id"])] = row
        return {spot: row for spot, row in latest.items() if row.get("flag") != "clear"}

    def digest(self, alignment_sha256: str) -> str:
        effective = self.effective(alignment_sha256)
        canonical = json.dumps(sorted((spot, row["flag"]) for spot, row in effective.items()))
        return hashlib.sha256(canonical.encode("utf-8")).hexdigest()


def build_submission_text(review_id: str, entries: list[dict]) -> str:
    return SUBMISSION_PREFIX + json.dumps({"review_id": review_id, "flags": entries},
                                          ensure_ascii=False, separators=(",", ":"))


def parse_submission_text(text: str) -> dict:
    for line in (text or "").splitlines():
        line = line.strip()
        if line.startswith(SUBMISSION_PREFIX):
            try:
                body = json.loads(line[len(SUBMISSION_PREFIX):])
            except json.JSONDecodeError as exc:
                raise ValueError(f"送信用テキストの JSON が壊れています: {exc}") from exc
            if not isinstance(body, dict) or "review_id" not in body or "flags" not in body:
                raise ValueError("送信用テキストに review_id と flags がありません。")
            return {"review_id": str(body["review_id"]), "flags": body["flags"]}
    raise ValueError(f"送信用テキストが見つかりません（`{SUBMISSION_PREFIX.strip()}` で始まる行）。")
