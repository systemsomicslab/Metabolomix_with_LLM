import asyncio
import json
import builtins
import os
import struct
import tempfile
import unittest
import unittest.mock
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

from metabolomix.core import mcp_core
import server
from metabolomix.core import session_state
from metabolomix.eic.reader import read_eic_spot_css1, read_eic_spots_css1


def _write_css1(path: Path, spots: list[dict]) -> None:
    chunks = []
    for spot in spots:
        chunk = bytearray(struct.pack(
            "<ffffbi",
            spot.get("rt", 0.0), spot.get("ri", 0.0), spot.get("mz", 0.0),
            spot.get("drift", 0.0), spot.get("main_type", 0), len(spot["samples"]),
        ))
        for sample in spot["samples"]:
            points = sample["points"]
            chunk.extend(struct.pack(
                "<iifff", sample["file_id"], len(points), sample["peak_top"],
                sample["peak_left"], sample["peak_right"],
            ))
            for x, y in points:
                chunk.extend(struct.pack("<ff", x, y))
        chunks.append(bytes(chunk))

    header_size = 14 + 8 * len(chunks)
    offsets = []
    cursor = header_size
    for chunk in chunks:
        offsets.append(cursor)
        cursor += len(chunk)

    body = bytearray(b"CSS1" + b"\x00" * 6)
    body.extend(struct.pack("<i", len(chunks)))
    for offset in offsets:
        body.extend(struct.pack("<q", offset))
    for chunk in chunks:
        body.extend(chunk)
    path.write_bytes(body)


class EicRandomAccessTests(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.path = Path(self._tmp.name) / "alignment.EIC.aef"
        _write_css1(self.path, [
            {"rt": 1.0, "mz": 100.0, "samples": [{
                "file_id": 1, "peak_top": 1.0, "peak_left": 0.9,
                "peak_right": 1.1, "points": [(0.9, 10.0), (1.0, 20.0)],
            }]},
            {"rt": 2.0, "mz": 200.0, "samples": [
                {"file_id": 10, "peak_top": 2.0, "peak_left": 1.8,
                 "peak_right": 2.2,
                 "points": [(1.8, 2.0), (2.0, 8.0), (2.2, 3.0)]},
                {"file_id": 11, "peak_top": 2.1, "peak_left": 1.9,
                 "peak_right": 2.3, "points": [(1.9, 4.0), (2.1, 12.0)]},
            ]},
        ])

    def tearDown(self):
        self._tmp.cleanup()

    def test_reads_one_spot_and_selected_file_id(self):
        spot = read_eic_spot_css1(self.path, 1, file_ids=[11])
        self.assertEqual(spot["spot_id"], 1)
        self.assertAlmostEqual(spot["mz"], 200.0)
        self.assertEqual(spot["selected_samples"], 1)
        self.assertEqual(spot["selected_points"], 2)
        self.assertEqual(spot["samples"][0]["file_id"], 11)
        self.assertAlmostEqual(spot["samples"][0]["peak_top"], 2.1, places=5)
        self.assertAlmostEqual(spot["samples"][0]["chromatogram"][0][0], 1.9, places=5)
        self.assertAlmostEqual(spot["samples"][0]["chromatogram"][1][0], 2.1, places=5)
        self.assertEqual(
            [point[1] for point in spot["samples"][0]["chromatogram"]],
            [4.0, 12.0],
        )

    def test_rejects_missing_file_id_and_out_of_range_spot(self):
        with self.assertRaisesRegex(ValueError, "not found"):
            read_eic_spot_css1(self.path, 1, file_ids=[999])
        with self.assertRaisesRegex(ValueError, "out of range"):
            read_eic_spot_css1(self.path, 2)

    def test_requires_selection_when_trace_limit_is_exceeded(self):
        crowded = Path(self._tmp.name) / "crowded.EIC.aef"
        samples = [{
            "file_id": index, "peak_top": 1.0, "peak_left": 0.9,
            "peak_right": 1.1, "points": [(1.0, float(index))],
        } for index in range(13)]
        _write_css1(crowded, [{"rt": 1.0, "mz": 100.0, "samples": samples}])
        with self.assertRaisesRegex(ValueError, "Specify file_ids"):
            read_eic_spot_css1(crowded, 0)


class EicBatchReadTests(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.path = Path(self._tmp.name) / "batch.EIC.aef"
        _write_css1(self.path, [
            {"rt": 1.0, "mz": 100.0, "samples": [
                {"file_id": 1, "peak_top": 1.0, "peak_left": 0.9,
                 "peak_right": 1.1, "points": [(0.9, 10.0), (1.0, 20.0)]},
                {"file_id": 2, "peak_top": 1.0, "peak_left": 0.9,
                 "peak_right": 1.1, "points": [(0.9, 1.0), (1.0, 2.0)]},
            ]},
            {"rt": 2.0, "mz": 200.0, "samples": [
                {"file_id": 2, "peak_top": 2.0, "peak_left": 1.8,
                 "peak_right": 2.2, "points": [(1.8, 3.0), (2.0, 9.0)]},
            ]},
            {"rt": 3.0, "mz": 300.0, "samples": [
                {"file_id": 1, "peak_top": 3.0, "peak_left": 2.8,
                 "peak_right": 3.2, "points": [(2.8, 4.0), (3.0, 8.0)]},
            ]},
        ])

    def tearDown(self):
        self._tmp.cleanup()

    def test_reads_requested_spots_in_one_file_open(self):
        opened = []
        real_open = builtins.open

        def counting_open(*args, **kwargs):
            opened.append(args[0])
            return real_open(*args, **kwargs)

        with unittest.mock.patch("builtins.open", counting_open):
            spots = read_eic_spots_css1(self.path, [2, 0], file_ids=[1])

        self.assertEqual(len(opened), 1)
        self.assertEqual([spot["spot_id"] for spot in spots], [0, 2])
        self.assertEqual(spots[0]["samples"][0]["file_id"], 1)
        self.assertEqual(
            [point[1] for point in spots[1]["samples"][0]["chromatogram"]],
            [4.0, 8.0],
        )

    def test_out_of_range_spot_is_omitted_when_not_strict(self):
        spots = read_eic_spots_css1(self.path, [0, 99], file_ids=[1])
        self.assertEqual([spot["spot_id"] for spot in spots], [0])

    def test_spot_without_requested_file_id_is_returned_empty(self):
        spots = read_eic_spots_css1(self.path, [0, 1], file_ids=[1])
        self.assertEqual([spot["spot_id"] for spot in spots], [0, 1])
        self.assertEqual(spots[1]["samples"], [])
        self.assertEqual(spots[1]["selected_samples"], 0)

    def test_strict_mode_raises_for_out_of_range_and_missing_file_id(self):
        with self.assertRaisesRegex(ValueError, "out of range"):
            read_eic_spots_css1(self.path, [99], file_ids=[1], strict=True)
        with self.assertRaisesRegex(ValueError, "not found"):
            read_eic_spots_css1(self.path, [1], file_ids=[1], strict=True)

    def test_point_budget_is_shared_across_spots(self):
        with self.assertRaisesRegex(ValueError, "point safety limit"):
            read_eic_spots_css1(
                self.path, [0, 2], file_ids=[1], max_total_points=3,
            )


class EicPlotToolTests(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = Path(self._tmp.name)
        self.path = self.tmp / "alignment.EIC.aef"
        _write_css1(self.path, [{
            "rt": 3.25, "mz": 512.3456, "samples": [{
                "file_id": 7, "peak_top": 3.25, "peak_left": 3.1,
                "peak_right": 3.4,
                "points": [(3.1, 5.0), (3.25, 10.0), (3.4, 2.0)],
            }],
        }])
        self._saved_data_dir = mcp_core.DATA_DIR
        self._saved_reports = os.environ.get("LIPIDMIX_REPORTS_DIR")
        self._saved_plot = session_state.session.eic.last_plot
        mcp_core.DATA_DIR = self.tmp
        os.environ["LIPIDMIX_REPORTS_DIR"] = str(self.tmp / "fallback")
        session_state.session.eic.last_plot = None

    def tearDown(self):
        mcp_core.DATA_DIR = self._saved_data_dir
        session_state.session.eic.last_plot = self._saved_plot
        if self._saved_reports is None:
            os.environ.pop("LIPIDMIX_REPORTS_DIR", None)
        else:
            os.environ["LIPIDMIX_REPORTS_DIR"] = self._saved_reports
        self._tmp.cleanup()

    def test_returns_structured_plot_data_without_writing_png(self):
        payload = json.loads(server.eic_plot_chromatograms(
            0, file_path=str(self.path), file_ids=[7], normalize="per_trace_max",
        ))
        self.assertEqual(payload["plot_schema"], "lipidmix.eic.v1")
        self.assertEqual(payload["plot_type"], "line")
        self.assertEqual(payload["axes"]["x"]["unit"], "min")
        for actual, expected in zip(payload["series"][0]["x"], [3.1, 3.25, 3.4]):
            self.assertAlmostEqual(actual, expected, places=5)
        self.assertEqual(payload["series"][0]["y"], [0.5, 1.0, 0.2])
        self.assertEqual(session_state.session.eic.last_plot, payload)
        self.assertEqual(list(self.tmp.rglob("*.png")), [])

    def test_fastmcp_does_not_publish_an_output_schema(self):
        """outputSchema を publish しない（structured_output=False）。

        publish すると payload が content と structuredContent へ二重に載る。
        クライアントは content のテキストから `plot_schema` を読む契約。
        """
        tools = asyncio.run(server.mcp.list_tools())
        tool = next(item for item in tools if item.name == "eic_plot_chromatograms")
        self.assertIsNone(tool.outputSchema)

    def test_png_is_written_only_by_explicit_save_tool(self):
        server.eic_plot_chromatograms(0, file_path=str(self.path), file_ids=[7])
        message = server.save_eic_figure("sample-eic")
        png = self.tmp / "fallback" / "figures" / "sample-eic_eic.png"
        self.assertTrue(png.is_file())
        self.assertIn("figures/sample-eic_eic.png", message)

    def test_save_tool_guides_when_plot_is_missing(self):
        message = server.save_eic_figure("sample-eic")
        self.assertIn("eic_plot_chromatograms", message)
        self.assertEqual(list(self.tmp.rglob("*.png")), [])


if __name__ == "__main__":
    unittest.main()
