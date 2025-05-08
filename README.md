# SNAP-ATAC Nextflow Pipeline

This repository contains a Nextflow DSL2 pipeline for processing single-cell ATAC-seq data using SnapATAC2. The workflow comprises two main steps:

1. Per-fragment preprocessing to generate individual fragment-level AnnData files (.h5ad).
2. Peak calling and concatenation of all fragment-level datasets into a combined AnnData object and peak-count matrix.

---

## Table of Contents

* [Prerequisites](#prerequisites)
* [Directory Structure](#directory-structure)
* [Parameters](#parameters)
* [Inputs](#inputs)
* [Outputs](#outputs)
* [Usage](#usage)
* [Example Command](#example-command)
* [Scripts](#scripts)
* [License](#license)

---

## Prerequisites

* [Nextflow](https://www.nextflow.io/) (v22.10.6 or later)
* Python packages:

  * scanpy 1.10.4
  * snapatac2 2.8.0

---

## Directory Structure

```
├── main.nf                   # Nextflow DSL2 pipeline script
├── params.config            # (Optional) Default parameter file
├── scripts/                 # Helper Python scripts
│   ├── snapatac2_per_fragment.py
│   └── snapatac2_peak_calling.py
├── data/
│   ├── fragment_files/      # Published per-sample .h5ad files
│   └── snapatac2_files/     # Output of preprocess step (snapatac2 format)
└── atac_names.txt           # List of sample names (one per line)
```

---

## Parameters

All parameters can be overridden via the command line or a custom `params.config` file.

| Parameter           | Default                  | Description                                       |
| ------------------- | ------------------------ | ------------------------------------------------- |
| `--concated_path`   | `./concated.h5ads`       | Filename of the output combined AnnData dataset   |
| `--peakmat_path`    | `./peakmat.h5ad`         | Filename of the output peak-count matrix          |
| `--atac_names_file` | `./atac_names.txt`       | Text file listing sample names (one per line)     |
| `--fragments_dir`   | `./data/fragment_files`  | Directory containing input fragment files         |
| `--snapatac_dir`    | `./data/snapatac2_files` | Directory to publish snapATAC2 intermediate files |

---

## Inputs

1. **Sample Names**: A plaintext file (`atac_names.txt`) listing each sample name on its own line.
2. **Fragment Files**: A directory (`data/fragment_files`) containing per-sample fragment files (.tsv.gz) from Cellranger output

---

## Outputs

* **Per-sample AnnData files**: `data/snapatac2_files/{sample_name}.h5ad`
* **Combined AnnData**: `data/concated.h5ads`
* **Peak-count Matrix**: `data/peakmat.h5ad`

---

## Usage

Run the pipeline with default parameters:

```bash
nextflow run main.nf -profile lsf
```

---

## Pipeline Scripts

* **scripts/snapatac2\_per\_fragment.py**: Processes individual fragment files to generate per-sample `.h5ad` AnnData objects.
* **scripts/snapatac2\_peak\_calling.py**: Performs peak calling on all per-sample files, concatenates them into a single AnnData object, and outputs a peak-count matrix.

