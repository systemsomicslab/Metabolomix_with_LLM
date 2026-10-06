import io
import tempfile
import unittest
from pathlib import Path
import zipfile

import lz4.block
import msgpack

from metabolomix.msdial.classes import (
    assign_sample_groups,
    attach_class_ids_to_spots,
    decode_msdial_messagepack,
    discover_arf_class_index,
    expand_class_specs,
    filter_arf_by_class_ids,
    get_sample_class_id,
    parse_analysis_file_classes,
    resolve_mddata_path,
)
from metabolomix.msdial.tags import normalize_sample_name


def pack_msdial(value) -> bytes:
    unpacked = msgpack.packb(value, use_bin_type=True)
    size_header = msgpack.packb(len(unpacked))
    compressed = lz4.block.compress(unpacked, store_size=False)
    extension = msgpack.ExtType(99, size_header + compressed)
    return msgpack.packb(extension, use_bin_type=True)


def analysis_file_row(file_id: int, name: str, class_id: str) -> list:
    return [
        f"C:/data/{name}.wiff",
        name,
        0,
        class_id,
        file_id + 1,
        file_id,
        True,
    ]


class MsdialClassTests(unittest.TestCase):
    def test_decode_and_parse_analysis_file_classes(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "dataset.mddata"
            storage = [[
                analysis_file_row(0, "sample_A", "control"),
                analysis_file_row(1, "sample_B", "treated"),
            ], [], [], [], [], None, None, None, None]
            path.write_bytes(pack_msdial(storage))

            self.assertEqual(decode_msdial_messagepack(path.read_bytes()), storage)
            records = parse_analysis_file_classes(path)
            self.assertEqual([record["class_id"] for record in records], ["control", "treated"])
            self.assertEqual([record["file_id"] for record in records], [0, 1])

    def test_resolve_mddata_from_mdproject_prefers_current_directory(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            mddata = directory / "dataset.mddata"
            mddata.write_bytes(pack_msdial([[], [], [], [], [], None, None, None, None]))
            project_parameter = [None, None, "version", "C:/old/location", "dataset.mddata"]
            project = {"ProjectParameters": [project_parameter]}
            mdproject = directory / "project.mdproject"
            with zipfile.ZipFile(mdproject, "w") as archive:
                archive.writestr("Project", pack_msdial(project))

            self.assertEqual(resolve_mddata_path(mdproject), mddata.resolve())

    def test_resolve_mddata_uses_latest_timestamp_when_multiple_exist(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            older = directory / "Dataset_2024_06_13_03_04_16.mddata"
            newer = directory / "Dataset_2026_06_17_00_34_26.mddata"
            older.write_bytes(pack_msdial([[analysis_file_row(0, "sample_A", "old")], [], [], [], [], None, None, None, None]))
            newer.write_bytes(pack_msdial([[analysis_file_row(0, "sample_A", "new")], [], [], [], [], None, None, None, None]))

            self.assertEqual(resolve_mddata_path(directory), newer.resolve())
            index = discover_arf_class_index(directory / "result_PeakProperties.arf")
            self.assertEqual(index["mddata_path"], str(newer.resolve()))
            self.assertEqual(index["class_counts"], {"new": 1})

    def test_resolve_mddata_for_timestamped_arf_uses_latest_not_newer_dataset(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            older = directory / "Dataset_2024_06_13_03_04_16.mddata"
            newer = directory / "Dataset_2026_06_17_00_34_26.mddata"
            older.write_bytes(pack_msdial([[analysis_file_row(0, "sample_A", "matching")], [], [], [], [], None, None, None, None]))
            newer.write_bytes(pack_msdial([[analysis_file_row(0, "sample_A", "too_new")], [], [], [], [], None, None, None, None]))
            arf_path = directory / "AlignmentResult_2024_06_13_09_21_45_PeakProperties.arf"
            arf_path.touch()

            self.assertEqual(resolve_mddata_path(arf_path), older.resolve())
            index = discover_arf_class_index(arf_path)
            self.assertEqual(index["class_counts"], {"matching": 1})

    def test_attach_and_filter_variable_sample_count(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            arf_path = directory / "result_PeakProperties.arf"
            arf_path.touch()
            sample_count = 73
            rows = []
            analysis_files = []
            for file_id in range(sample_count):
                name = f"sample_{file_id}"
                class_id = "A" if file_id % 2 == 0 else "B"
                rows.append([file_id, name, 1000 + file_id])
                analysis_files.append(analysis_file_row(file_id, name, class_id))
            storage = [analysis_files, [], [], [], [], None, None, None, None]
            (directory / "dataset.mddata").write_bytes(pack_msdial(storage))
            features = [{"MasterAlignmentID": 1, "AlignedPeakProperties": rows}]

            index = discover_arf_class_index(arf_path)
            self.assertEqual(index["class_counts"], {"A": 37, "B": 36})
            attach_class_ids_to_spots(features, index)
            self.assertEqual(get_sample_class_id(features[0], 0, "sample_0"), "A")

            filtered, stats = filter_arf_by_class_ids(features, index, ["B"])
            self.assertEqual(len(filtered[0]["AlignedPeakProperties"]), 36)
            self.assertEqual(stats["after_sample_peaks"], 36)

    def test_file_id_and_file_name_mismatch_fails(self):
        class_index = {
            "by_file_id": {0: {"file_id": 0, "file_name": "sample_A", "class_id": "A"}},
            "by_file_name": {
                "sample_b": {"file_id": 1, "file_name": "sample_B", "class_id": "B"},
            },
            "records": [],
        }
        features = [{"AlignedPeakProperties": [[0, "sample_B", 10]]}]
        with self.assertRaisesRegex(ValueError, "different mddata samples"):
            attach_class_ids_to_spots(features, class_index)

    def test_unknown_class_id_fails(self):
        index = {
            "records": [{"class_id": "A"}],
            "by_file_id": {},
            "by_file_name": {},
        }
        with self.assertRaisesRegex(ValueError, "matched"):
            filter_arf_by_class_ids([], index, ["missing"])


CLASS_VOCAB = [
    "Cerebellum_gf_AIN",
    "Cerebellum_gf_HFD",
    "Cerebellum_spf_AIN",
    "Hippocampus_gf_AIN",
]


class ExpandClassSpecsTests(unittest.TestCase):
    def test_single_token_matches_all_with_token(self):
        result = expand_class_specs(["gf"], CLASS_VOCAB)
        self.assertEqual(
            set(result["gf"]),
            {"Cerebellum_gf_AIN", "Cerebellum_gf_HFD", "Hippocampus_gf_AIN"},
        )

    def test_multi_token_spec_is_and(self):
        result = expand_class_specs(["Cerebellum_gf"], CLASS_VOCAB)
        self.assertEqual(
            set(result["Cerebellum_gf"]),
            {"Cerebellum_gf_AIN", "Cerebellum_gf_HFD"},
        )

    def test_exact_full_id_matches_only_itself(self):
        result = expand_class_specs(["Cerebellum_gf_AIN"], CLASS_VOCAB)
        self.assertEqual(result["Cerebellum_gf_AIN"], ["Cerebellum_gf_AIN"])

    def test_case_insensitive(self):
        result = expand_class_specs(["CEREBELLUM_GF"], CLASS_VOCAB)
        self.assertEqual(
            set(result["CEREBELLUM_GF"]),
            {"Cerebellum_gf_AIN", "Cerebellum_gf_HFD"},
        )

    def test_zero_match_spec_raises(self):
        with self.assertRaisesRegex(ValueError, "matched"):
            expand_class_specs(["liver"], CLASS_VOCAB)


def _name_class_index(name_to_class: dict[str, str]) -> dict:
    records = []
    by_file_name = {}
    for file_id, (name, class_id) in enumerate(name_to_class.items()):
        record = {"file_id": file_id, "file_name": name, "class_id": class_id}
        records.append(record)
        by_file_name[normalize_sample_name(name, strip_processing_timestamp=False)] = record
    return {"records": records, "by_file_id": {}, "by_file_name": by_file_name}


class AssignSampleGroupsTests(unittest.TestCase):
    def setUp(self):
        self.index = _name_class_index({
            "s1": "Cerebellum_gf_AIN",
            "s2": "Cerebellum_spf_AIN",
            "s3": "Cerebellum_gf_HFD",
        })

    def test_default_group_is_full_class_id(self):
        groups = assign_sample_groups(["s1", "s2"], self.index)
        self.assertEqual(groups["s1"], "Cerebellum_gf_AIN")
        self.assertEqual(groups["s2"], "Cerebellum_spf_AIN")

    def test_group_levels_collapse_to_one_factor(self):
        groups = assign_sample_groups(["s1", "s2", "s3"], self.index, group_levels=["gf", "spf"])
        self.assertEqual(groups["s1"], "gf")
        self.assertEqual(groups["s3"], "gf")
        self.assertEqual(groups["s2"], "spf")

    def test_sample_matching_no_level_is_other(self):
        groups = assign_sample_groups(["s1"], self.index, group_levels=["liver"])
        self.assertEqual(groups["s1"], "other")

    def test_sample_matching_multiple_levels_raises(self):
        with self.assertRaisesRegex(ValueError, "multiple"):
            assign_sample_groups(["s1"], self.index, group_levels=["gf", "AIN"])

    def test_no_class_index_yields_none(self):
        groups = assign_sample_groups(["s1", "s2"], None)
        self.assertEqual(groups, {"s1": None, "s2": None})


class FilterTokenSubsetTests(unittest.TestCase):
    def test_partial_spec_filters_subset_and_reports_matched(self):
        index = _name_class_index({
            "s1": "Cerebellum_gf_AIN",
            "s2": "Cerebellum_spf_AIN",
            "s3": "Hippocampus_gf_AIN",
        })
        features = [{
            "MasterAlignmentID": 1,
            "AlignedPeakProperties": [
                [0, "s1", 100],
                [1, "s2", 200],
                [2, "s3", 300],
            ],
        }]
        filtered, stats = filter_arf_by_class_ids(features, index, ["gf"])
        kept = {row[1] for row in filtered[0]["AlignedPeakProperties"]}
        self.assertEqual(kept, {"s1", "s3"})
        self.assertEqual(
            set(stats["matched_class_ids"]),
            {"Cerebellum_gf_AIN", "Hippocampus_gf_AIN"},
        )


if __name__ == "__main__":
    unittest.main()
