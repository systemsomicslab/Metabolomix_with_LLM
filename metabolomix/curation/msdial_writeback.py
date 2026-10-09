"""送信されたフラグを MS-DIAL のアラインメントの `_tags.xml` に反映する（Misannotation / Confirmed）。

規則（ユーザー決定 2026-09-29 / 2026-10-09）:
- wrong → Misannotation を付け、Confirmed を外す
- confirmed → Confirmed を付け、Misannotation を外す
- suspect → Confirmed を外す（Misannotation は触らない）
- clear → 両方外す
他のタグは残す。
正本は `flags.jsonl`。こちらは MS-DIAL で見るための写しで、失敗しても記録は巻き戻さない。
書く前に `<curation>/tags-backup/` へ既存ファイルの控えを取る。

MS-DIAL はアラインメントを保存するたびにメモリ上のタグで `_tags.xml` を丸ごと書き直し、
読むのはプロジェクトを開くときだけ（上流 `AlignmentResultContainer.Save` / `Load`）。
deps: msdial.tags、curation.flags（curation_dir）。
"""
from __future__ import annotations

import shutil
import xml.etree.ElementTree as ET
from datetime import datetime, timezone

from metabolomix.curation.flags import curation_dir
from metabolomix.msdial import tags as msdial_tags

BACKUP_DIR = "tags-backup"
MSDIAL_NOTE = ("MS-DIAL でこのプロジェクトを開いたままだと、GUI の保存で _tags.xml が上書きされ"
               "この反映は消えます。プロジェクトを閉じてから開き直すと Misannotation / Confirmed が見えます。")


def _ids(entries, *flag_values) -> list[int]:
    return [e["spot_id"] for e in entries if e["flag"] in flag_values]


def sync_tags(arf2_path, entries: list[dict]) -> dict:
    """`entries`（検証済みの `{spot_id, flag}`）を `_tags.xml` の Misannotation と Confirmed に反映する。
    戻り値の `added` / `removed` は Misannotation（従来どおり）、`confirmed` は Confirmed の変化。"""
    changes = {
        msdial_tags.MISANNOTATION_TAG_ID: {"add": _ids(entries, "wrong"),
                                           "remove": _ids(entries, "confirmed", "clear")},
        msdial_tags.CONFIRMED_TAG_ID: {"add": _ids(entries, "confirmed"),
                                       "remove": _ids(entries, "wrong", "suspect", "clear")},
    }
    tag_path = msdial_tags.alignment_tag_path(arf2_path)
    if not any(c["add"] or c["remove"] for c in changes.values()):
        return {"path": str(tag_path), "added": [], "removed": [], "confirmed": {"added": [], "removed": []},
                "backup": None, "note": MSDIAL_NOTE}
    backup = None
    try:
        if tag_path.is_file():
            folder = curation_dir(arf2_path) / BACKUP_DIR
            folder.mkdir(parents=True, exist_ok=True)
            stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
            backup = folder / f"{tag_path.name}.{stamp}"
            shutil.copy2(tag_path, backup)
        result = msdial_tags.update_alignment_tags(tag_path, changes)
    except (OSError, ValueError, ET.ParseError) as exc:   # ParseError: 壊れた XML / ValueError: 別形式
        return {"path": str(tag_path), "error": f"{type(exc).__name__}: {exc}",
                "backup": None if backup is None else str(backup),
                "message": "フラグは記録済みです。_tags.xml への反映だけが失敗しました。"}
    return {"path": result["path"], **result["tags"][msdial_tags.MISANNOTATION_TAG_ID],
            "confirmed": result["tags"][msdial_tags.CONFIRMED_TAG_ID], "created": result["created"],
            "backup": None if backup is None else str(backup), "note": MSDIAL_NOTE}
