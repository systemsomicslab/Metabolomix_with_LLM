# Lipidmix with LLM

An MCP server (`ms-data-parser`) that parses MS-DIAL's binary lipidomics outputs, runs the analyses
server-side, and returns compact summaries, so an LLM client (Claude Desktop, Claude Code, a local WebUI)
can drive an analysis without seeing the raw matrices.

## What it does

- **Parses MS-DIAL binaries** — alignment results, per-measurement peak lists, deconvoluted MS/MS
  spectra, extracted ion chromatograms.
- **Runs the analysis server-side** — preprocessing and QC (normalization, blank subtraction, QC-RSD
  filtering, drift correction), PCA, differential analysis with BH-FDR.
- **Keeps results out of the context window** — plots are rendered to PNG on the server; full result
  tables stay in session state until exported.
- **Records the analysis** — objectives and reports go to `reports/` in the analysis folder (next to the
  raw data and MS-DIAL output, never under the clone; `data/` is untracked).
- **Can drive MS-DIAL itself** — two entries meet at an MS-DIAL analysis folder: read one MS-DIAL has
  written (e.g. in the GUI), or build one from raw data with MS-DIAL Console. Given a folder, the server
  uses existing MS-DIAL output if present; with raw files only, it runs Console per polarity
  (`pipeline_run`), checks each analysis folder against its mzTab-M, and opens the result with
  `load_dataset(pipeline_path=..., polarity=...)`.

## Supported inputs

| Format | Contents |
|---|---|
| `.arf` | Alignment result, per-sample peaks. Main source for cross-sample comparison and PCA. |
| `.arf2` | Alignment result, spot representatives. Dataset-wide catalog and annotations. |
| `.dcl` | Deconvoluted MS/MS spectra (MS-DIAL's own binary format, not MessagePack). |
| `.EIC.aef` | Extracted ion chromatograms. |
| mzTab-M 2.0 | Standard exchange format, produced by the MS-DIAL Console path. |

MS-DIAL sidecars (`*_tags.xml`, `.mddata`, `.mdproject`) are read for tags and sample class assignments
where present.

## Quick start

Run the setup script once from the repository root. It uses [uv](https://docs.astral.sh/uv/) when
available (uv also fetches Python) and otherwise a standard `venv` with Python 3.13+; dependency
versions are pinned by `uv.lock`.

```powershell
powershell -ExecutionPolicy Bypass -File setup.ps1          # Windows
```

```bash
sh setup.sh          # macOS / Linux
```

It creates `.venv/`, checks that every dependency imports, copies `lipidmix.local.toml` from the example
if missing, writes `.mcp.json` for Claude Code, and prints the snippet for Claude Desktop's
`claude_desktop_config.json`. Existing `.mcp.json` and `lipidmix.local.toml` are never overwritten;
differences are reported.

There is no default data directory or environment variable; pass the MS-DIAL output folder to the client:

> Load the dataset in `/path/to/your/msdial/output` and show me its samples and groups.

That calls `load_dataset`, which picks the latest alignment batch, opens it as the current dataset, and
reports its samples, groups and library.

### MS-DIAL Console and reference libraries (optional)

MS-DIAL Console from raw data (`pipeline_*`, with `console_method_*` to pick or make a method file) and
your own `.msp` library (`library_load`) need file locations. Copy `lipidmix.example.toml` to
`lipidmix.local.toml` next to `server.py` and fill in the paths (Windows paths in single quotes). The file
is untracked and read on every call (no restart needed), and these locations are read from it only; the
lipid library (`.lbm2`) goes under `[library] lbm`. It also holds where reports and the library cache go
(`[paths]`) and server behaviour such as the open dataset-handle limit (`[server]`). A tool's error names
the file and key of a missing setting.

## Repository map

`server.py` registers the tools by import side effect and stays at the repository root (MCP client configs
reference it by absolute path). The implementation is under `metabolomix/`; the package map is in
[CLAUDE.md](CLAUDE.md).

## Documentation

| To find out | Read |
|---|---|
| The whole system with diagrams (Japanese) | [ms-data-parser specification](https://yuukamegai.github.io/ms-data-parser-spec/) |
| What each MCP tool takes and does | [USAGE.md](USAGE.md) |
| What an output field *means* (also the MCP resource `lipidmix://docs/output-format`) | [docs/output_format/](docs/output_format/core.md) |
| Which functions a tool calls, in what order | [docs/workflow/](docs/workflow/index.md) |
| MessagePack key indices, transcribed from MS-DIAL's `[Key(N)]` attributes | [docs/schema/](docs/schema/AlignmentSpotProperty.md) |
| Running the parsers from the command line | [docs/cli.md](docs/cli.md) |
| Setting up the local MCP server on a member's machine | [DEPLOY.md](DEPLOY.md) |
| Working on this repository | [CLAUDE.md](CLAUDE.md) |

Tests run from the repository root and need `pytest` (the development group):

```powershell
powershell -ExecutionPolicy Bypass -File setup.ps1 -Dev
.venv/Scripts/python.exe -m pytest tests -q
```

## Status and caveats

A research prototype, developed against one lab's MS-DIAL outputs.

- **Version-coupled.** The binary readers follow MS-DIAL's internal class layout; a release that changes
  those classes needs the key indices in `docs/schema/` re-checked.
- **Identification confidence is reported, not assumed.** An MS/MS *flag* and an acquired spectrum are
  kept distinct; identification is reported as a confidence band.

## License

[MIT](LICENSE).

This project reads MS-DIAL's output files; it contains no MS-DIAL source code and is not a derivative of it.
The MessagePack key indices in [docs/schema/](docs/schema/AlignmentSpotProperty.md) are transcribed from the
`[Key(N)]` attributes of MS-DIAL's own classes, with the upstream repository and commit recorded in each
file. MS-DIAL ([systemsomicslab/MsdialWorkbench](https://github.com/systemsomicslab/MsdialWorkbench)) is
licensed under LGPL-3.0 by its authors.
