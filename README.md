# EPIC Data Viewer

A coordinate-based viewer and reproducible preparation pipeline for the [EPIC promoter and transcription initiation challenge](https://epic.autosome.org/).

The pipeline reads the original EPIC Zenodo files and produces:

1. `data_fixed/`: coordinate-indexed reference DNA, eligibility intervals, sparse pooled counts, and multiresolution peak indexes for fast browsing.
2. `parquet/`: **seven Parquet files with one row per reference base**, without artificial 4,096-bp segments. All reference bases on each contig represented in a dataset's whitelist are retained, including excluded context with null targets. Parquet row groups are storage units, not biological or model boundaries.

No experimental data are included in this repository. EPIC's [data usage agreement](https://epic.autosome.org/docs/data-usage-agreement) restricts redistribution outside an authorized team without permission. The code license does not grant rights to the experimental data. Configure your own hosting access for your authorized audience.

## Install

Python 3.12+, sufficient disk space, and the extracted source files are required. Allow at least 30 GB beyond the compressed source archives for reference indexes, Parquet output, and intermediate files. Exporting several billion base rows takes time; it streams bounded batches rather than loading the genomes into RAM.

```bash
git clone https://github.com/jsomeara/EPIC-data-viewer.git
cd EPIC-data-viewer
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
```

## Obtain and extract the raw EPIC data

Download the datasets linked by the organizers at [Zenodo DOI 10.5281/zenodo.22285753](https://doi.org/10.5281/zenodo.22285753), including the Nematostella worked example. Follow the organizers' access and usage requirements. Extract the archives into `raw/`, preserving the assembly directories. Do not flatten or rename their contents.

Expected layout:

```text
raw/
  ASM119413v2/
  Ofas_2.0/
  cgigas_uk_roslin_v1/
  ilPloInte3.2/
  sSquSuc5.hap2/
  jaNemVect1.1/
```

Each assembly directory must contain:

```text
genome/genome.fa.gz
genome/genome.fa.fai
genome/assembly_report.txt
for_participants/whitelist.train.bed.gz
for_participants/template.train.bed.gz
for_participants/csRNA-train.txt.gz
```

The toy `jaNemVect1.1` directory additionally needs `for_participants/whitelist.test.bed.gz`, `for_participants/template.test.bed.gz`, and `ground_truth/csRNA-test.txt.gz`. The other species' test ground truth is withheld; this pipeline does not invent or export it. Control sRNA-seq and individual replicate files are not needed: the target is the **supplied pooled** csRNA vector. See the [official file descriptions](https://epic.autosome.org/docs/epic-data).

## Convert raw data into the seven segment-free files

Run these commands from the repository root:

```bash
# Build the coordinate-native data directly from FASTA, whitelist, template, and pooled vector.
.venv/bin/python prepare/from_raw.py --raw raw --output data_fixed

# Export and fully verify the seven Parquet files. Use fewer workers on smaller machines.
.venv/bin/python prepare/export_parquet.py --fixed data_fixed --output parquet --workers 3
```

For a smaller first run, add `--assembly jaNemVect1.1` to the first command. That builds both toy splits; the exporter processes the datasets present in `data_fixed/catalog.json`.

The output files are:

| File | Organism | Split |
|---|---|---|
| `california_two_spot_octopus_train.parquet` | Octopus bimaculoides | Train |
| `indianmeal_moth_train.parquet` | Plodia interpunctella | Train |
| `milkweed_bug_train.parquet` | Oncopeltus fasciatus | Train |
| `pacific_oyster_train.parquet` | Magallana gigas | Train |
| `puget_sound_dogfish_train.parquet` | Squalus suckleyi | Train |
| `starlet_sea_anemone_toy_train.parquet` | Nematostella vectensis | Train |
| `starlet_sea_anemone_toy_test.parquet` | Nematostella vectensis | Test |

`parquet/manifest.json` contains row counts, eligibility totals, non-zero target counts, count sums, file sizes, SHA-256 hashes, and read-back verification status. Per-file JSON records are also written. Incomplete files use `.partial`; they are renamed only after every row group passes verification. Completed outputs are reused on rerun. To rebuild after changing the input, use new output directories; do not mix an old index with changed raw files.

### What the conversion does

- Reads and uppercases the reference FASTA, preserving N and other ambiguous bases. Checks every sequence length against the supplied `.fai`.
- Retains the exact supplied whitelist and train/test separation.
- Validates the template's actual order against the whitelist, with the plus half followed by the minus half. A changed input format fails explicitly rather than silently reassigning counts.
- Reads the dense pooled vector in chunks, validating its length and nonnegative integer counts. Maps each value back to its reference position and strand using the template. Missing positive entries inside the whitelist are measured zeros; outside the whitelist targets are null.
- Verifies that indexing preserves the source's non-zero counts and total counts. Records source file hashes.
- Writes contiguous per-base rows across each whole included contig, with no cropping or circular shifts. The two strands use the same increasing reference coordinates.
- Reads back every Parquet row group with page-checksum verification and compares its logical values against the original arrays before publishing it.

No normalization, log transform, background subtraction, peak calling, read-count threshold, or zero downsampling is applied. The viewer's log scale and peak summaries affect display only.

## Parquet schema

| Column | Type | Meaning |
|---|---|---|
| `assembly`, `organism`, `split`, `sequence_id` | string / dictionary-encoded string | Reference assembly, species, split, chromosome/contig. |
| `position` | uint64 | Zero-based reference coordinate. |
| `base` | string / dictionary-encoded string | One uppercase reference nucleotide, including ambiguous letters. |
| `eligible` | bool | Whether the base belongs to the original whitelist. |
| `csrna_plus`, `csrna_minus` | nullable uint64 | Supplied pooled counts; zero is measured zero, null is excluded. |
| `label_plus`, `label_minus` | nullable uint8 | `count > 0` at eligible bases; otherwise null. |

Read regions lazily rather than loading an entire genome into pandas:

```python
import pyarrow.dataset as ds

data = ds.dataset('parquet/milkweed_bug_train.parquet', format='parquet')
region = ((ds.field('sequence_id') == 'YOUR_CONTIG_ID')
          & (ds.field('position') >= 100_000)
          & (ds.field('position') < 110_000))
for batch in data.scanner(filter=region, batch_size=4096).to_batches():
    # Choose your own model windows and retain the eligibility mask.
    print(batch)
```

## Run the viewer

```bash
.venv/bin/python viewer/app.py --port 8767
```

Open http://127.0.0.1:8767. The viewer reads `data_fixed/`; verified downloads come from `parquet/manifest.json`. You can browse while the export is running; only completed, manifested downloads appear. The interface has a GitHub link and a Downloads menu next to GitHub for the seven files and their checksum manifest.

The chart uses canvas, viewport-only queries and multiresolution indexes. Scroll zooms, drag selects a zoom range, Shift-drag temporarily swaps zoom/pan, and double-click resets. The overview supports pan/resize. Log10(count + 1) is enabled by default, with raw counts in hover/table views. Nucleotides appear at readable zoom. The exact base table paginates 100 rows and supports a non-zero filter. PNG and SVG export use the visible resolution.

### Docker

After preparing `data_fixed/` and optionally `parquet/`:

```bash
docker compose up -d --build
```

The default Compose configuration publishes only `127.0.0.1:8767`. Data mounts are read-only and containers restart unless explicitly stopped. Put your access-controlled reverse proxy in front of it for team hosting. For a URL prefix such as `/epic-data-viewer/`, redirect the non-slash path to the slash path and strip that prefix when proxying; frontend URLs are relative. See `deploy/nginx.conf` for an example, including efficient direct downloads with byte-range support.

### Deploy as a Rune

For an existing Runes installation, clone this repository into `~/runes/runes/epic-data-viewer/src/`, copy `deploy/rune.compose.yaml` to the parent `compose.yaml`, and place `data_fixed/` and `parquet/` under `~/runes/data/epic-data-viewer/`. Use `deploy/runes-route.conf` as the gateway route. The gateway needs a read-only mount of the Parquet directory at `/downloads/epic-data-viewer`. Then run:

```bash
~/runes/bin/runes deploy epic-data-viewer
```

The route strips the URL prefix for app requests and serves verified Parquet files directly through Nginx. The viewer also supports streaming/ranged downloads itself when run without Nginx.

## Tests

Synthetic tests need no EPIC data:

```bash
.venv/bin/python -m unittest discover -s tests -v
node tests/test_plot.cjs
```

The raw converter was additionally checked against the existing independently verified toy train/test data: every sparse signal/whitelist byte and the full reference DNA matched. The hosted export verifies all generated row groups. Tests cover raw-to-Parquet alignment, null versus zero, counts on both strands, download ranges, scroll zoom, Shift-drag, log scaling, and nucleotide rendering.
