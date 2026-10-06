# Lipidmix with LLM

An MCP server (`ms-data-parser`) that reads MS-DIAL's binary lipidomics outputs and
exposes them to an LLM in a form it can actually reason over.

MS-DIAL writes its alignment results as compressed MessagePack and custom binary blobs.
They are unreadable outside the GUI, and far too large to paste into a chat. This server
parses them, runs the standard lipidomics analyses server-side, and returns compact
summaries — so an LLM client such as Claude Desktop, Claude Code, or a local WebUI can
drive a full analysis without ever seeing the raw matrices.

## What it does

- **Parses MS-DIAL binaries** — alignment results, per-measurement peak lists,
  deconvoluted MS/MS spectra, and extracted ion chromatograms.
- **Runs the analysis server-side** — preprocessing and QC (normalization, blank
  subtraction, QC-RSD filtering, drift correction), PCA, and two-group differential
  analysis with BH-FDR.
- **Keeps results out of the context window** — plots are rendered to PNG on the server
  rather than shipped as coordinate lists, and full result tables stay in session state
  until exported.
- **Accumulates knowledge** — literature notes and reusable analysis playbooks are served
  back to the LLM as MCP resources, so findings compound across sessions.
- **Can drive MS-DIAL itself (auxiliary)** — the main workflow reads output you have already
  processed in the MS-DIAL GUI. As a secondary path, the server can plan and run MS-DIAL
  Console over raw acquisition files and load the resulting mzTab-M. It only does this when
  you choose to; given raw data alone, it first suggests processing it in the GUI.

## Supported inputs

| Format | Contents |
|---|---|
| `.arf` | Alignment result, per-sample peaks. The main source for cross-sample comparison and PCA. |
| `.arf2` | Alignment result, spot representatives. Dataset-wide catalog and annotations. |
| `.pai2` | A single measurement's detected peaks. |
| `.dcl` | Deconvoluted MS/MS spectra. MS-DIAL's own binary format — not MessagePack. |
| `.EIC.aef` | Extracted ion chromatograms. |
| mzTab-M 2.0 | Standard exchange format, produced by the MS-DIAL Console path. |

MS-DIAL sidecars (`*_tags.xml`, `.mddata`, `.mdproject`) are read for tags and sample
class assignments where present.

## Quick start

Requires Python 3.13+ (developed on 3.14).

```bash
python -m pip install -r requirements.txt
```

Register the server with your MCP client — for example, in `.mcp.json`:

```json
{
  "mcpServers": {
    "ms-data-parser": {
      "type": "stdio",
      "command": "python",
      "args": ["/absolute/path/to/Lipidmix_with_LLM/server.py"]
    }
  }
}
```

Point it at your data and ask the client to start:

```bash
export LIPIDMIX_DATA_DIR=/path/to/your/msdial/output
```

> Load the dataset in `LIPIDMIX_DATA_DIR` and show me the PCA.

That calls `load_dataset`, the entry point: it picks the latest alignment batch, runs the
standard overview, and primes the session for everything that follows.

### MS-DIAL Console and reference libraries (optional)

Running MS-DIAL Console from raw data (`console_*`, `pipeline_*`) or matching against your own
`.msp` library (`library_load`) needs to know where those files live. Copy
`lipidmix.example.toml` to `lipidmix.local.toml` next to `server.py` and fill in the paths
(Windows paths in single quotes). The file is not tracked by git and is read on every call, so
edits take effect without restarting the server. The environment variables `MSDIAL_EXE`,
`MSDIAL_LBM`, `MSDIAL_MSP_POS` and `MSDIAL_MSP_NEG` still work and take precedence over the file.
If something is missing, the tool's error names the file and the key to fill in.

To run the parsers without an MCP client, see [docs/cli.md](docs/cli.md).

## Repository map

`server.py` is a thin facade that registers tools by import side effect; it must stay at
the repository root because MCP client configs reference it by absolute path. The
implementation lives under `metabolomix/`.

| Package | Responsibility |
|---|---|
| `metabolomix/core/` | FastMCP instance, configuration, session state, path resolution. Depends on nothing else in the tree. |
| `metabolomix/{arf,arf2,pai2,dcl,eic}/` | One `reader.py` (parser) and `tools.py` (MCP surface) per input format. |
| `metabolomix/mztab/` | mzTab-M reader and canonical `DatasetState`. |
| `metabolomix/console/` | MS-DIAL Console execution: job planning, running, output collection. |
| `metabolomix/analysis/` | Format-independent numerics: preprocessing/QC, PCA, differential analysis, the export contract. |
| `metabolomix/plots/` | Renderer-neutral plot payloads and matplotlib rendering. |
| `metabolomix/msdial/` | MS-DIAL-specific sidecars, identification, peak verification. |
| `metabolomix/corpus/` | Pure logic for the knowledge and playbook notes. |
| `metabolomix/tools/` | MCP tools not tied to a single input format. |

## Documentation

| To find out | Read |
|---|---|
| The whole system at a glance — every tool's inputs and outputs, the analysis routes, the statistics, and export formats, with diagrams (Japanese) | [ms-data-parser specification](https://yuukamegai.github.io/ms-data-parser-spec/) |
| What each MCP tool takes and does | [USAGE.md](USAGE.md) |
| What an output field *means* — row granularity, lipid-name grammar, caveats | [docs/output_format/](docs/output_format/core.md), also served as the MCP resource `lipidmix://docs/output-format` |
| Which functions a tool calls, in what order | [docs/workflow/](docs/workflow/index.md) |
| MessagePack key indices, transcribed from MS-DIAL's `[Key(N)]` attributes | [docs/schema/](docs/schema/AlignmentSpotProperty.md) |
| Running the parsers from the command line | [docs/cli.md](docs/cli.md) |
| Setting up the local MCP server on a member's machine | [DEPLOY.md](DEPLOY.md) |
| Working on this repository | [CLAUDE.md](CLAUDE.md) |

Tests run from the repository root. They need `pytest`, which is not part of the
runtime dependencies; install it through the development requirements first:

```bash
python -m pip install -r requirements-dev.txt
python -m pytest tests -q
```

## Status and caveats

This is a research prototype, developed against one lab's MS-DIAL outputs. Treat the
following as known limits rather than surprises:

- **Version-coupled.** The binary readers follow MS-DIAL's internal class layout. A
  MS-DIAL release that changes those classes will need the key indices in `docs/schema/`
  re-checked.
- **Multi-group ANOVA is deliberately unavailable.** MS-DIAL metadata carries no
  factor-to-level mapping, so the server asks you to carve out two groups explicitly
  instead of guessing a design.
- **Identification confidence is reported, not assumed.** An MS/MS *flag* and an actual
  acquired spectrum are distinguished throughout; check the reported band before claiming
  an MSI level.
- **Generated artifacts are not checked in.** `data/`, `analyses/`, and `reports/` are
  untracked; a clean checkout will not have them.

## License

[MIT](LICENSE).

This project reads MS-DIAL's output files; it contains no MS-DIAL source code and is not
a derivative of it. The MessagePack key indices in [docs/schema/](docs/schema/AlignmentSpotProperty.md)
are transcribed from the `[Key(N)]` attributes of MS-DIAL's own classes, with the upstream
repository and the commit they were checked against recorded in each file. MS-DIAL
([systemsomicslab/MsdialWorkbench](https://github.com/systemsomicslab/MsdialWorkbench)) is
licensed under LGPL-3.0 by its authors.
