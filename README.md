<p align="center">
  <img src="docs/logo.png" alt="AnResist logo" width="200"/>
</p>

# AnResist 🧬

> **A reproducible Nextflow DSL2 pipeline for cross-tool antimicrobial resistance (AMR) gene detection and systematic comparison across multiple bacterial genomes.**

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A522.10.0-brightgreen)](https://www.nextflow.io/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Conda](https://img.shields.io/badge/conda-mamba-green)](https://mamba.readthedocs.io/)

## Table of Contents

1. [Overview](#overview)
2. [Biological Motivation](#biological-motivation)
3. [Pipeline Architecture](#pipeline-architecture)
4. [Tools Included](#tools-included)
5. [Requirements](#requirements)
6. [Installation](#installation)
7. [Quick Start](#quick-start)
8. [Input Format](#input-format)
9. [Configuration](#configuration)
10. [Output Structure](#output-structure)
11. [Test Dataset](#test-dataset)
12. [Reproducibility](#reproducibility)
13. [Contact](#contact)

## Overview

AnResist runs four established AMR detection tools — **AMRFinderPlus**, **RGI (CARD)**, **ABRicate**, and **ResFinder** — simultaneously on a panel of bacterial genome assemblies. It then harmonizes their outputs into a unified schema and performs systematic cross-tool comparison, highlighting genes where tools agree (high confidence) and where they disagree (requiring attention) with clear visual interpretations.

The key novelty of AnResist is not running the tools — it is **quantifying and explaining the disagreement** between them, producing per-sample and global figures and tables that make tool discordance actionable for clinical and research settings.

## Biological Motivation

AMR is one of the greatest threats to global public health. While multiple bioinformatics tools exist for detecting resistance genes from genome assemblies, they use different databases (NCBI NDARO, CARD, ResFinder DB) and different detection algorithms. This means the same genome can yield different resistance profiles depending on which tool is used.

AnResist addresses the question: **which AMR genes are robustly detected across tools, and which are tool-specific?How can we interpret it?**

A gene detected by all four tools is high-confidence. A gene detected by only one tool may represent a novel variant, a database gap, or a false positive — all of which have clinical implications.


## Pipeline Architecture

```
Raw FASTA assemblies
        │
        ▼
   01. QUAST (QC)
        │
        ▼
   02. FILTER_SAMPLES (size/quality thresholds)
        │
        ├─────────────────────────────────────┐
        ▼                                     ▼
   DB SETUP (once)                    Tool execution (parallel)
   ├── AMRFINDER_SETUP                ├── RUN_AMRFINDER
   ├── RGI_SETUP                      ├── RUN_RGI
   └── RESFINDER_SETUP                ├── RUN_ABRICATE
                                      └── RUN_RESFINDER
                                             │
                                             ▼
                                      03. HARMONIZE
                                      (unified TSV per sample)
                                             │
                                             ▼
                                      04. COMPARE
                                      (cross-tool analysis +
                                       publication figures)

```

## Tools Included

| Tool | Database | Version | Purpose |
|------|----------|---------|---------|
|Quast|-|5.2.0|Filter sequences|
| AMRFinderPlus | NCBI NDARO | 2026-03-24.1 | Nucleotide + protein AMR detection |
| RGI | CARD | latest | Homology + SNP-based AMR prediction |
| ABRicate | NCBI | 1.0.1 | Fast sequence-based screening |
| ResFinder | ResFinder DB | latest | Acquired resistance gene detection |

---

## Requirements

- Linux or macOS (tested on Ubuntu 22.04, WSL2)
- [Nextflow](https://www.nextflow.io/) ≥ 22.10.0
- [Conda](https://docs.conda.io/) + [Mamba](https://mamba.readthedocs.io/) (recommended) 
- Git
- Internet access for first run (database downloads)

**Disk space:** ~5 GB for all databases and conda environments on first run.

**Memory:** Minimum 8 GB RAM recommended.

> [!WARNING]
> Install AnResist in your home directory (e.g. `~/AnResist`). Nextflow writes temporary files and conda environments to paths relative to the working directory — installing in system directories or network mounts will cause permission errors.
> 
> Install Nextflow inside a conda environment, not system-wide:
> ```bash
> conda create -n nextflow -c bioconda nextflow
> conda activate nextflow
> ```
> Always activate this environment before running the pipeline.


## Installation

Clone the repository

```bash
git clone https://github.com/aruna-venkat/AnResist.git
cd AnResist
```

## Quick Start

### Run with test data

```bash
nextflow run main.nf \
    -params-file config/params.yml \
    -profile conda
```

This will:
1. Download all required databases (first run only, ~10 min)
2. Run QC and filtering
3. Run all four AMR tools in parallel
4. Harmonize outputs
5. Generate comparison tables and figures

### Resume after interruption

```bash
nextflow run main.nf \
    -params-file config/params.yml \
    -profile conda \
    -resume
```

## Input Format

### To run with your own assembled fasta contigs-

AnResist requires a **samplesheet CSV** with two columns:

```csv
sample,fasta
S_aureus_MRSA,data/S_aureus_MRSA.fasta
K_pneumoniae,data/K_pneumoniae.fasta
E_faecium,data/E_faecium.fasta
```

- `sample` — unique identifier for the sample (no spaces)
- `fasta` — path to the genome assembly FASTA file (relative )

FASTA files should be assembled genome sequences (contigs or complete chromosomes). Raw reads are not supported.

Set the samplesheet path in `config/params.yml`:

```yaml
samplesheet: "config/samplesheet.csv"
```

## Configuration

All parameters are set in `config/params.yml`:

## Output Structure

```
results/
├── qc/                          # QUAST assembly statistics
├── amrfinder/{sample}/          # Raw AMRFinderPlus TSVs
├── rgi/{sample}/                # Raw RGI output
├── abricate/{sample}/           # Raw ABRicate TSVs
├── resfinder/{sample}/          # Raw ResFinder TSVs
├── harmonized/                  # Unified TSVs (one per sample)
│   └── {sample}_unified.tsv
└── comparison/                  # Cross-tool analysis
    ├── global/
    │   ├── all_hits.tsv                        # Every gene call, all tools, all samples
    │   ├── global_gene_summary.tsv             # Per-gene confidence across samples
    │   ├── systematically_discordant_genes.tsv # Genes always missed by some tools
    │   └── figures/
    │       ├── global_jaccard_heatmap.png      # Tool similarity matrix
    │       ├── gene_count_heatmap.png          # Genes detected per tool per sample
    │       └── global_consensus_discordance.png
    └── {sample}/
        ├── {sample}_all_hits.tsv               # All tool outputs for this sample
        ├── {sample}_gene_detection_summary.tsv # Per-gene confidence + explanation
        ├── {sample}_high_confidence_genes.tsv  # Genes found by ALL tools
        ├── {sample}_discordant_single_tool.tsv # Genes found by only ONE tool
        ├── {sample}_pairwise_overlap.tsv       # Tool pairs: shared and unique gene lists
        ├── {sample}_discordance_report.txt     # Human-readable plain English report
        └── figures/
            ├── {sample}_gene_detection_map.png     # Dot plot: gene × tool × identity
            ├── {sample}_upset.png                  # Set intersection plot
            ├── {sample}_drug_class_confidence.png  # Drug class breakdown
            ├── {sample}_identity_scatter_shared_genes.png
            └── {sample}_intersection_gene_lists.tsv
```

### Key output files explained

**`{sample}_gene_detection_summary.tsv`** — the most important per-sample file. One row per gene with:
- Which tools detected it
- Identity and coverage per tool
- Confidence classification (high / moderate / low)
- Plain-English note explaining discordance

**`{sample}_discordance_report.txt`** — human-readable report listing every gene by confidence tier with biological interpretation notes.

**`global_jaccard_heatmap.png`** — shows pairwise Jaccard similarity between all four tools across all samples. Lower values = more disagreement = more interesting biology.

## Test Dataset

A minimal test dataset is included in `test/` covering three clinically relevant organisms:

| Organism | Accession | Known resistance |
|----------|-----------|-----------------|
| S. aureus MRSA | NC_002952.2 | mecA, bleomycin, aminoglycoside |
| K. pneumoniae | NC_009648.1 | Beta-lactam, aminoglycoside |
| N. gonorrhoeae | NC_011035.1 | Fluoroquinolone, beta-lactam |

Expected outputs: 3 sample folders in `results/comparison/`, minimum 5 genes detected per sample, `global_jaccard_heatmap.png` showing similarity scores between 0.3 and 0.8.

Expected runtime: ~15 minutes on first run (database downloads), ~5 minutes on resume.

## Reproducibility

AnResist is designed for full reproducibility:

| Feature | Implementation |
|---------|---------------|
| Workflow management | Nextflow DSL2 with full dependency graph |
| Environment isolation | Per-process conda environments via mamba |
| Database versioning | `storeDir` caches databases across runs; pin versions in `nextflow.config` |
| Resume capability | `-resume` flag reuses cached task outputs |
| Portability | Runs identically on laptop, HPC, or cloud |
| Version control | Full git history with incremental commits |

## Contact

**Aruna** — M.S. Bioinformatics, Georgia Institute of Technology  
GitHub Issues: [github.com/aruna-venkat/AnResist/issues](https://github.com/venkat-aruna/AnResist/issues)

For bug reports please include:
- Your operating system and Nextflow version (`nextflow -version`)
- The `.nextflow.log` file
- The full error message from the failed process work directory (`.command.err`)

## Citation

If you use AnResist in your research, please cite:

> AnResist: A reproducible Nextflow pipeline for cross-tool AMR gene detection and systematic comparison. [2026]. GitHub: [https://github.com/aruna-venkat/AnResist](https://github.com/venkat-aruna/AnResist)


## License

MIT License — see [LICENSE](LICENSE) for details.
