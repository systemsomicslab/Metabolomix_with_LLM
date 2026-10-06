import tempfile
import unittest
from pathlib import Path

from metabolomix.msdial.tags import (
    attach_tags_to_spots,
    discover_arf_tag_index,
    filter_arf_by_tags,
    get_sample_peak_tag_info,
    normalize_sample_name,
    parse_tag_file,
)


TAG_XML = """<?xml version="1.0" encoding="utf-8"?>
<PeakSpotTags>
  <Definitions>
    <Tag><Id>1</Id><Label>Confirmed</Label></Tag>
    <Tag><Id>3</Id><Label>Misannotation</Label></Tag>
  </Definitions>
  <Peaks>{peaks}</Peaks>
</PeakSpotTags>
"""


def write_tags(path: Path, peaks: str) -> None:
    path.write_text(TAG_XML.format(peaks=peaks), encoding="utf-8")


class MsdialTagsTests(unittest.TestCase):
    def test_normalize_sample_name_removes_timestamp_and_extensions(self):
        expected = "sample_a"
        self.assertEqual(normalize_sample_name("sample_A"), expected)
        self.assertEqual(normalize_sample_name("sample_A.wiff"), expected)
        self.assertEqual(normalize_sample_name("sample_A_202605151012_tags.xml"), expected)
        self.assertEqual(
            normalize_sample_name(
                "sample_A_202605151012_tags.xml",
                strip_processing_timestamp=False,
            ),
            "sample_a_202605151012",
        )

    def test_normalize_sample_name_removes_every_msdial_raw_extension(self):
        """MS-DIAL は `GetFileNameWithoutExtension` でサンプル名を作る。

        剥がす拡張子の集合がそれより狭いと、`QC01.d` 由来の名前が `qc01.d` の
        まま残り、`qc01` を名乗るタグサイドカーや ARF の FileName と結合できない。
        """
        for ext in ("abf", "ibf", "cdf", "mzml", "wiff", "raw",
                    "d", "wiff2", "qgd", "lcd", "lrp", "imzml"):
            with self.subTest(extension=ext):
                self.assertEqual(normalize_sample_name(f"sample_A.{ext}"), "sample_a")

    def test_parse_tag_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "sample_tags.xml"
            write_tags(path, '<Peak Id="10"><Tag>1</Tag><Tag>3</Tag></Peak>')
            parsed = parse_tag_file(path)
            self.assertEqual(parsed["definitions"], {1: "Confirmed", 3: "Misannotation"})
            self.assertEqual(parsed["peaks"][10], frozenset({1, 3}))

    def test_discover_attach_and_filter_sample_peaks(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            arf_path = directory / "AlignmentResult_1_PeakProperties.arf"
            arf_path.touch()
            write_tags(
                directory / "sample_A_202605151012_tags.xml",
                '<Peak Id="10"><Tag>1</Tag></Peak>',
            )
            write_tags(
                directory / "sample_B_202605151012_tags",
                '<Peak Id="20"><Tag>3</Tag></Peak>',
            )
            write_tags(
                directory / "AlignmentResult_1_tags.xml",
                '<Peak Id="7"><Tag>1</Tag></Peak>',
            )
            features = [{
                "MasterAlignmentID": 7,
                "AlignedPeakProperties": [
                    [0, "sample_A", 10],
                    [1, "sample_B", 20],
                ],
            }]

            index = discover_arf_tag_index(arf_path, features)
            self.assertEqual(index["summary"]["sample_tag_files"], 2)
            self.assertEqual(index["summary"]["matched_arf_samples"], 2)
            attach_tags_to_spots(features, index)
            self.assertEqual(features[0]["Tags"], ["Confirmed"])
            self.assertEqual(features[0]["SamplePeakTags"][0]["Tags"], ["Confirmed"])

            filtered, stats = filter_arf_by_tags(
                features, index, ["Confirmed"], mode="any", scope="sample_peak",
            )
            self.assertEqual(filtered[0]["AlignedPeakProperties"], [[0, "sample_A", 10]])
            self.assertEqual(stats["after_sample_peaks"], 1)

            filtered, _ = filter_arf_by_tags(
                features, index, ["Misannotation"], mode="none", scope="sample_peak",
            )
            self.assertEqual(filtered[0]["AlignedPeakProperties"], [[0, "sample_A", 10]])

    def test_filter_alignment_spots_and_all_mode(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            arf_path = directory / "AlignmentResult_1_PeakProperties.arf"
            arf_path.touch()
            write_tags(
                directory / "AlignmentResult_1_tags.xml",
                '<Peak Id="1"><Tag>1</Tag><Tag>3</Tag></Peak>',
            )
            features = [
                {"MasterAlignmentID": 1, "AlignedPeakProperties": [[0, "sample_A", 10]]},
                {"MasterAlignmentID": 2, "AlignedPeakProperties": [[0, "sample_A", 11]]},
            ]
            index = discover_arf_tag_index(arf_path, features)
            filtered, _ = filter_arf_by_tags(
                features,
                index,
                ["Confirmed", "Misannotation"],
                mode="all",
                scope="alignment_spot",
            )
            self.assertEqual([spot["MasterAlignmentID"] for spot in filtered], [1])

    def test_variable_sample_count_is_not_fixed(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            arf_path = directory / "AlignmentResult_1_PeakProperties.arf"
            arf_path.touch()
            sample_count = 137
            rows = []
            for sample_id in range(sample_count):
                name = f"sample_{sample_id}"
                peak_id = 1000 + sample_id
                rows.append([sample_id, name, peak_id])
                write_tags(
                    directory / f"{name}_202605151012_tags.xml",
                    f'<Peak Id="{peak_id}"><Tag>1</Tag></Peak>',
                )
            features = [{"MasterAlignmentID": 1, "AlignedPeakProperties": rows}]

            index = discover_arf_tag_index(arf_path, features)
            filtered, stats = filter_arf_by_tags(
                features, index, ["Confirmed"], scope="sample_peak",
            )
            self.assertEqual(index["summary"]["matched_arf_samples"], sample_count)
            self.assertEqual(len(filtered[0]["AlignedPeakProperties"]), sample_count)
            self.assertEqual(stats["after_sample_peaks"], sample_count)

    def test_missing_sample_tag_file_requires_explicit_policy(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            arf_path = directory / "AlignmentResult_1_PeakProperties.arf"
            arf_path.touch()
            write_tags(
                directory / "sample_A_202605151012_tags.xml",
                '<Peak Id="10"><Tag>1</Tag></Peak>',
            )
            features = [{
                "MasterAlignmentID": 1,
                "AlignedPeakProperties": [[0, "sample_A", 10], [1, "sample_B", 20]],
            }]
            index = discover_arf_tag_index(arf_path, features)

            with self.assertRaisesRegex(ValueError, "sample_B"):
                filter_arf_by_tags(features, index, ["Confirmed"])

            excluded, _ = filter_arf_by_tags(
                features, index, ["Confirmed"], missing_sample_policy="exclude",
            )
            self.assertEqual(excluded[0]["AlignedPeakProperties"], [[0, "sample_A", 10]])

            untagged, _ = filter_arf_by_tags(
                features,
                index,
                ["Confirmed"],
                mode="none",
                missing_sample_policy="untagged",
            )
            self.assertEqual(untagged[0]["AlignedPeakProperties"], [[1, "sample_B", 20]])

    def test_exact_name_wins_and_ambiguous_relaxed_match_fails(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            arf_path = directory / "AlignmentResult_1_PeakProperties.arf"
            arf_path.touch()
            features = [{
                "AlignedPeakProperties": [
                    [0, "patient", 10],
                    [1, "patient_202605151012", 20],
                ],
            }]
            write_tags(
                directory / "patient_202605151012_tags.xml",
                '<Peak Id="20"><Tag>1</Tag></Peak>',
            )
            index = discover_arf_tag_index(arf_path, features)
            self.assertIn("patient_202605151012", index["sample_files"])

            write_tags(
                directory / "patient_202701010101_tags.xml",
                '<Peak Id="10"><Tag>1</Tag></Peak>',
            )
            with self.assertRaisesRegex(ValueError, "Ambiguous timestamp-stripped"):
                discover_arf_tag_index(arf_path, features)

    def test_duplicate_files_for_one_sample_uses_latest_timestamp(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            arf_path = directory / "AlignmentResult_1_PeakProperties.arf"
            arf_path.touch()
            features = [{"AlignedPeakProperties": [[0, "sample_A", 10]]}]
            write_tags(directory / "sample_A_202501010101_tags.xml", "")
            write_tags(
                directory / "sample_A_202605151012_tags.xml",
                '<Peak Id="10"><Tag>1</Tag></Peak>',
            )

            index = discover_arf_tag_index(arf_path, features)

            self.assertEqual(index["summary"]["duplicate_sample_tag_files"], 1)
            self.assertTrue(
                index["sample_files"]["sample_a"]["path"].endswith(
                    "sample_A_202605151012_tags.xml"
                )
            )
            attach_tags_to_spots(features, index)
            self.assertEqual(features[0]["SamplePeakTags"][0]["Tags"], ["Confirmed"])

    def test_custom_tag_directory_applies_to_alignment_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            arf_directory = root / "arf"
            tag_directory = root / "tags"
            arf_directory.mkdir()
            tag_directory.mkdir()
            arf_path = arf_directory / "AlignmentResult_1_PeakProperties.arf"
            arf_path.touch()
            alignment_path = tag_directory / "AlignmentResult_1_tags.xml"
            write_tags(alignment_path, '<Peak Id="1"><Tag>1</Tag></Peak>')

            index = discover_arf_tag_index(
                arf_path,
                [{"MasterAlignmentID": 1, "AlignedPeakProperties": []}],
                tag_directory=tag_directory,
            )
            self.assertEqual(index["alignment_file"], str(alignment_path))
            self.assertEqual(index["alignment_peaks"][1], frozenset({1}))

    def test_alignment_filter_requires_alignment_tag_file(self):
        index = {
            "definitions": {1: "Confirmed"},
            "alignment_file": None,
            "alignment_peaks": {},
        }
        features = [{"MasterAlignmentID": 1, "AlignedPeakProperties": []}]
        with self.assertRaisesRegex(ValueError, "No alignment-level"):
            filter_arf_by_tags(
                features,
                index,
                ["Confirmed"],
                scope="alignment_spot",
            )

    def test_sample_tag_lookup_falls_back_to_file_name(self):
        spot = {
            "SamplePeakTags": {
                "sample_a": {"Tags": ["Confirmed"]},
            },
        }
        info = get_sample_peak_tag_info(spot, None, "sample_A")
        self.assertEqual(info["Tags"], ["Confirmed"])


if __name__ == "__main__":
    unittest.main()
