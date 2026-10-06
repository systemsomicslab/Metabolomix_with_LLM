"""アラインメントの `_tags.xml` への書き戻し（MS-DIAL 5 の AlignmentResultContainer.Save / Load と
同じ形: `<Peaks><Peak Id="<MasterAlignmentID>"><Tag>3</Tag></Peak></Peaks>`）。"""
from metabolomix.msdial import tags

EXISTING = """\ufeff<?xml version="1.0" encoding="utf-8"?>
<PeakSpotTags>
  <Definitions>
    <Tag>
      <Id>1</Id>
      <Label>Confirmed</Label>
    </Tag>
    <Tag>
      <Id>3</Id>
      <Label>Misannotation</Label>
    </Tag>
  </Definitions>
  <Peaks>
    <Peak Id="0">
      <Tag>1</Tag>
    </Peak>
    <Peak Id="4">
      <Tag>1</Tag>
      <Tag>3</Tag>
    </Peak>
    <Peak Id="5">
      <Tag>3</Tag>
    </Peak>
  </Peaks>
</PeakSpotTags>
"""


def write_existing(tmp_path):
    path = tmp_path / "AlignmentResult_x_tags.xml"
    path.write_bytes(EXISTING.encode("utf-8"))
    return path


def test_adding_misannotation_keeps_other_tags_and_definitions(tmp_path):
    path = write_existing(tmp_path)
    result = tags.update_alignment_tag(path, tag_id=tags.MISANNOTATION_TAG_ID, add=[0, 7, 4], remove=[])
    parsed = tags.parse_tag_file(path)
    assert parsed["definitions"] == {1: "Confirmed", 3: "Misannotation"}
    assert parsed["peaks"] == {0: frozenset({1, 3}), 4: frozenset({1, 3}), 5: frozenset({3}),
                               7: frozenset({3})}
    assert result["added"] == [0, 7]          # 4 は既に付いていた
    assert result["removed"] == [] and result["created"] is False


def test_removing_misannotation_drops_peaks_left_without_tags(tmp_path):
    path = write_existing(tmp_path)
    result = tags.update_alignment_tag(path, tag_id=tags.MISANNOTATION_TAG_ID, add=[], remove=[4, 5, 9])
    assert tags.parse_tag_file(path)["peaks"] == {0: frozenset({1}), 4: frozenset({1})}
    assert result["removed"] == [4, 5]        # 9 は元から付いていない
    assert '<Peak Id="5">' not in path.read_text(encoding="utf-8-sig")


def test_missing_file_is_created_with_the_msdial_definitions(tmp_path):
    path = tmp_path / "AlignmentResult_x_tags.xml"
    result = tags.update_alignment_tag(path, tag_id=tags.MISANNOTATION_TAG_ID, add=[2], remove=[])
    parsed = tags.parse_tag_file(path)
    assert parsed["definitions"] == dict(tags.MSDIAL_TAG_DEFINITIONS)
    assert parsed["peaks"] == {2: frozenset({3})}
    assert result["created"] is True


def test_written_file_has_the_msdial_layout(tmp_path):
    path = write_existing(tmp_path)
    tags.update_alignment_tag(path, tag_id=tags.MISANNOTATION_TAG_ID, add=[0], remove=[])
    raw = path.read_bytes()
    assert raw.startswith(b"\xef\xbb\xbf<?xml version=\"1.0\" encoding=\"utf-8\"?>")   # XElement.Save と同じ BOM
    text = raw.decode("utf-8-sig")
    # XElement.Save（Windows）と同じ: 2 空白字下げ・CRLF・末尾改行なし
    assert '  <Peaks>\r\n    <Peak Id="0">\r\n      <Tag>1</Tag>\r\n      <Tag>3</Tag>\r\n    </Peak>' in text
    assert text.endswith("</PeakSpotTags>")


def test_alignment_tag_path_sits_next_to_the_arf2(tmp_path):
    arf2 = tmp_path / "AlignmentResult_2026_09_09_17_31_52.arf2"
    assert tags.alignment_tag_path(arf2) == tmp_path / "AlignmentResult_2026_09_09_17_31_52_tags.xml"
