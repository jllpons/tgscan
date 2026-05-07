# tgscan.nf

**tgscan.nf** is a pipeline that scans a 6-frame translated genome for sequences of interest using [HMMER](http://hmmer.org/).

## Contents

- [Overview](#overview)
- [Quick Start](#quick-start)
- [Input](#input)
- [Output](#output)
- [Parameters](#parameters)

## Overview

`tgscan.nf` takes a genome FASTA file and a set of HMM profiles, performs a 6-frame translation of the genome, and runs `hmmsearch` to identify sequences matching the provided profiles. Finally, the genomic coordinates of all hits meeting the specified score thresholds are back-mapped to the original genome and reported in GFF3 format.

## Quick Start

```bash
nextflow run main.nf \
   --input samplesheet.csv \
   --outdir results
```

## Input

The pipeline requires a CSV samplesheet passed via `--input`. The file must have the following columns:

| Column   | Description                                                                 |
|----------|-----------------------------------------------------------------------------|
| `sample` | Unique sample name (no spaces)                                              |
| `genome` | Path to the genome FASTA file (`.fa`, `.fna`, `.fasta`, or `.gz` variants) |
| `hmm`    | Path to the HMM profile file (`.hmm`, uncompressed)                        |

**Example `samplesheet.csv`:**

```csv
sample,genome,hmm
sample1,/path/to/genome1.fasta,/path/to/profiles.hmm
sample2,/path/to/genome2.fna.gz,/path/to/profiles.hmm
```

> **Note:** Multiple HMM profiles can be concatenated into a single `.hmm` file and passed as the `hmm` field.

## Output

All results are written to the directory specified by `--outdir` (default: `tgscan_results/`), organised per sample:

```
tgscan_results/
└── <sample>/
    ├── hmmsearch/
    │   ├── hmmsearch.<sample>.domtbl.gff3.gz
    │   ├── hmmsearch.<sample>.domtbl.gz
    │   ├── hmmsearch.<sample>.sto.gz
    │   ├── hmmsearch.<sample>.tbl.gz
    │   └── hmmsearch.<sample>.txt.gz
    └── translate/
        └── <genome>.6frameORFs.fasta.gz
```

### `hmmsearch/`

| File | Description |
|------|-------------|
| `*.domtbl.gff3.gz` | Per-domain hits from `hmmsearch` in GFF3 format, compressed |
| `*.domtbl.gz` | Per-domain hits table (`--domtblout`), compressed |
| `*.tbl.gz` | Per-sequence hits table (`--tblout`), compressed |
| `*.sto.gz` | Multiple sequence alignment of hits in Stockholm format (`-A`), compressed |
| `*.txt.gz` | Full `hmmsearch` text output, compressed |

### `translate/`

| File | Description |
|------|-------------|
| `*.6frameORFs.fasta.gz` | 6-frame translated ORFs of the input genome, compressed |

## Parameters

### Required

| Parameter | Description                          |
|-----------|--------------------------------------|
| `--input` | Path to the input samplesheet (CSV)  |

### Optional

| Parameter          | Default           | Description                                    |
|--------------------|-------------------|------------------------------------------------|
| `--outdir`         | `tgscan_results`  | Directory where results will be saved          |
| `--hmmsearch_T`    | `null`            | `hmmsearch` bitscore target cutoff (`-T`)      |
| `--hmmsearch_domT` | `null`            | `hmmsearch` per-domain bitscore cutoff (`--domT`) |

### CLI

| Parameter   | Description                    |
|-------------|--------------------------------|
| `--help`    | Print help message and exit    |
| `--version` | Print pipeline version and exit |

## Testing & Playground

If you want to explore or adjust how hits are back-mapped to genomic coordinates and the GFF3 file is generated, the logic lives in a AWK script:

- Script: `bin/hmmsearch_tblout2gff.awk`
- Test file: `bin/test/test.domtbl`

You can run it directly with gawk to experiment without needing to execute the full pipeline:
```bash
gawk -f bin/hmmsearch_tblout2gff.awk bin/test/test.domtbl
```

## Dependencies

The pipeline relies on the following tools, which should be available in the execution environment (e.g. via a container or module system):

| Tool | Description |
|------|-------------|
| [HMMER](http://hmmer.org/) | Provides `hmmsearch` for profile-based sequence searching |
| [Easel miniapplications](http://hmmer.org/) | Utility programs bundled with HMMER. `esl-translate` is used for 6-frame translation of the genome. |
| [gawk](https://www.gnu.org/software/gawk/) | GNU AWK, used for processing `hmmsearch` output and back-mapping coordinates to the genome. |

---

## License

Until publication, no license is assigned to this code. Please contact the authors for permission to use or modify the code.
