"""キュレーションのフラグ記録（追記専用 JSON Lines）と送信用テキスト。spec §6。

記録するのはユーザーの判断: レビューの wrong / suspect / confirmed、候補付けの assign / redundant、取り消しの clear。
confirmed（注釈が正しいと確かめた。2026-10-09）は wrong / suspect と排他で、スポットごとの最新 1 行が効く。
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

FLAG_VALUES = ("wrong", "suspect", "confirmed", "clear", "assign", "redundant")
REVIEW_FLAG_VALUES = ("wrong", "suspect", "confirmed", "clear")
SUGGEST_FLAG_VALUES = ("assign", "redundant", "clear")
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
        if flag not in REVIEW_FLAG_VALUES:
            raise ValueError(f"flags[{index}].flag={flag!r} は {REVIEW_FLAG_VALUES} のいずれかにしてください。")
        note = str(entry.get("note") or "")
        for control in ("\r", "\n", "\t"):     # TSV の 1 セル・JSONL の 1 行に収める
            note = note.replace(control, " ")
        cleaned.append({"spot_id": spot_id, "flag": flag, "note": note[:500]})
    return cleaned


class FlagFileError(ValueError):
    """`flags.jsonl` に読めない行がある。黙って読み飛ばすと、その行の wrong が
    エクスポートから消えずに残る（＝判断が静かに失われる）ので、名指しで止める。"""

    def __init__(self, path, line_no: int, reason: str):
        self.path = Path(path)
        self.line_no = line_no
        self.reason = reason
        super().__init__(f"フラグ記録 {self.path} の {line_no} 行目が読めません（{reason}）。"
                         "その行を直すか消してから、もう一度実行してください。")

    def details(self) -> dict:
        return {"flags_file": str(self.path), "line": self.line_no}


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

    def exists(self) -> bool:
        return self.path.exists()

    def rows(self, *, tolerate_partial_tail: bool = True) -> list[dict]:
        """全行。読めない行があれば `FlagFileError`（行番号は 1 始まり）。"""
        if not self.path.exists():
            return []
        rows = []
        text = self.path.read_text(encoding="utf-8")
        lines = text.splitlines()
        # 受け口の HTTP スレッドが追記している最中に読むと、最終行が書きかけ（改行なし）のことがある。
        # ロックを持たない読み手（既定 tolerate_partial_tail=True）はその最終行だけ読み飛ばす。
        # 改行で終わっている壊れた行は常に FlagFileError。
        # 書込み側の事前確認（submission.submit_flags）は _WRITE_LOCK の中で呼ぶので、他に追記中の
        # 書き手はいない。そこでの書きかけ末尾は本物の破損で、その後ろへ追記すると新しい行が
        # 癒着して判断が失われるため、tolerate_partial_tail=False で FlagFileError にする。
        if lines and not text.endswith(("\n", "\r")):
            if tolerate_partial_tail:
                lines.pop()
            elif lines[-1].strip():
                raise FlagFileError(self.path, len(lines), "最終行が改行で終わっていない（書きかけ）")
        for line_no, line in enumerate(lines, 1):
            if not line.strip():
                continue
            try:
                row = json.loads(line)
            except json.JSONDecodeError as exc:
                raise FlagFileError(self.path, line_no, f"JSON として読めない: {exc.msg}") from exc
            if not isinstance(row, dict):
                raise FlagFileError(self.path, line_no, "JSON オブジェクトではない")
            spot_id = row.get("spot_id")
            if isinstance(spot_id, bool) or not isinstance(spot_id, int):
                raise FlagFileError(self.path, line_no, f"spot_id が整数ではない: {spot_id!r}")
            if row.get("flag") not in FLAG_VALUES:
                raise FlagFileError(self.path, line_no, f"flag が不正: {row.get('flag')!r}")
            rows.append(row)
        return rows

    def effective(self, alignment_sha256: str, rows: list[dict] | None = None) -> dict[int, dict]:
        return effective_flags(self.rows() if rows is None else rows, alignment_sha256)

    def digest(self, alignment_sha256: str, rows: list[dict] | None = None) -> str:
        return flags_digest(self.effective(alignment_sha256, rows))


def effective_flags(rows: list[dict], alignment_sha256: str) -> dict[int, dict]:
    """その sha256 のフラグのうち、スポットごとの最新 1 行（clear 済みは除く）。"""
    latest: dict[int, dict] = {}
    for row in rows:
        if row.get("alignment_sha256") == alignment_sha256:
            latest[int(row["spot_id"])] = row
    return {spot: row for spot, row in latest.items() if row.get("flag") != "clear"}


def cleared_spots(rows: list[dict], alignment_sha256: str) -> set[int]:
    """その sha256 で最新の行が clear のスポット（＝人が明示的に取り消した）。"""
    latest: dict[int, str] = {}
    for row in rows:
        if row.get("alignment_sha256") == alignment_sha256:
            latest[int(row["spot_id"])] = row.get("flag")
    return {spot for spot, flag in latest.items() if flag == "clear"}


def _digest_item(spot: int, row: dict) -> list:
    flag = row["flag"]
    if flag == "assign":
        # inchikey も含める: 同じ和組成・同じ粒度でも別レコードへ付け替えるとエクスポートの InChIKey が変わる
        return [spot, flag, row.get("name"), row.get("level"), row.get("inchikey")]
    if flag == "redundant":
        return [spot, flag, row.get("of"), row.get("relation")]
    return [spot, flag]


def flags_digest(effective: dict[int, dict]) -> str:
    # wrong / suspect は以前と同じ [spot, flag]（既存の curation_flags_sha256 を変えない）
    canonical = json.dumps(sorted((_digest_item(spot, row) for spot, row in effective.items()),
                                  key=lambda item: item[0]))
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()


def split_decisions(effective: dict[int, dict]) -> dict:
    return {"wrong": {s for s, r in effective.items() if r["flag"] == "wrong"},
            "suspect": {s for s, r in effective.items() if r["flag"] == "suspect"},
            "confirmed": {s for s, r in effective.items() if r["flag"] == "confirmed"},
            "assign": {s: r for s, r in effective.items() if r["flag"] == "assign"},
            "redundant": {s: r for s, r in effective.items() if r["flag"] == "redundant"}}


def orphaned_count(rows: list[dict], alignment: dict) -> int:
    """同じファイル名で sha256 が違う（＝以前の版の .arf2 に付いた）有効フラグの件数。
    それらはスポット番号が今の版と対応する保証が無いので当てない（spec §6）。"""
    other_shas = {row.get("alignment_sha256") for row in rows
                  if row.get("alignment_file") == alignment["alignment_file"]
                  and row.get("alignment_sha256") != alignment["alignment_sha256"]}
    return sum(len(effective_flags(rows, sha)) for sha in other_shas)


def orphaned_warning(n: int) -> str:
    return (f"{n} flag(s) were recorded against an earlier version of this alignment (same file name, "
            "different .arf2 sha256) and were not applied.")


def build_submission_text(review_id: str, entries: list[dict], arf2_path: str | None = None) -> str:
    """ビューアの「送信用テキストをコピー」と同じ形。`arf2_path`（絶対パス）を運ぶので、
    サーバ再起動やデータフォルダの切り替えの後でもレビューを探し当てられる。"""
    body = {"review_id": review_id, "flags": entries}
    if arf2_path is not None:
        body["arf2_path"] = str(arf2_path)
    return SUBMISSION_PREFIX + json.dumps(body, ensure_ascii=False, separators=(",", ":"))


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
            arf2_path = body.get("arf2_path")
            return {"review_id": str(body["review_id"]), "flags": body["flags"],
                    "arf2_path": str(arf2_path) if arf2_path else None}
    raise ValueError(f"送信用テキストが見つかりません（`{SUBMISSION_PREFIX.strip()}` で始まる行）。")
