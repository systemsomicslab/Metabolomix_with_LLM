"""送信されたフラグを MS-DIAL のアラインメントの `_tags.xml` に反映する（Misannotation）。

規則（ユーザー決定 2026-09-29）: 「間違い」(wrong) → 付ける、取消 (clear) → 外す、
「疑わしい」(suspect) → 触らない。他のタグ（Confirmed など）は残す。
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
               "この反映は消えます。プロジェクトを閉じてから開き直すと Misannotation が見えます。")


def sync_misannotation(arf2_path, entries: list[dict]) -> dict:
    """`entries`（検証済みの `{spot_id, flag}`）を `_tags.xml` の Misannotation に反映する。"""
    add = [e["spot_id"] for e in entries if e["flag"] == "wrong"]
    remove = [e["spot_id"] for e in entries if e["flag"] == "clear"]
    tag_path = msdial_tags.alignment_tag_path(arf2_path)
    if not add and not remove:
        return {"path": str(tag_path), "added": [], "removed": [], "backup": None, "note": MSDIAL_NOTE}
    backup = None
    try:
        if tag_path.is_file():
            folder = curation_dir(arf2_path) / BACKUP_DIR
            folder.mkdir(parents=True, exist_ok=True)
            stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
            backup = folder / f"{tag_path.name}.{stamp}"
            shutil.copy2(tag_path, backup)
        result = msdial_tags.update_alignment_tag(
            tag_path, tag_id=msdial_tags.MISANNOTATION_TAG_ID, add=add, remove=remove)
    except (OSError, ValueError, ET.ParseError) as exc:   # ParseError: 壊れた XML / ValueError: 別形式
        return {"path": str(tag_path), "error": f"{type(exc).__name__}: {exc}",
                "backup": None if backup is None else str(backup),
                "message": "フラグは記録済みです。_tags.xml への反映だけが失敗しました。"}
    return {**result, "backup": None if backup is None else str(backup), "note": MSDIAL_NOTE}
