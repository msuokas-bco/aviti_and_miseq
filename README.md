# DADA2 Pipeline for MiSeq and AVITI Reads

A Quarto-based R workflow for amplicon sequence variant (ASV) analysis from paired-end 16S rRNA sequencing data. The pipeline processes reads from two sequencing platforms — **Element Biosciences AVITI** and **Illumina MiSeq** — in parallel, then merges the results for a cross-platform comparison of community composition.

The workflow is adapted from the [DADA2 pipeline tutorial](https://benjjneb.github.io/dada2/tutorial.html) and the Bioconductor Workflow for Microbiome Data Analysis by Benjamin Callahan. The sequenced microbiome samples were collected in the Wood for Health project of the ERA-Net ForestValue program. The project’s lead partner was the Unit of Measurement Technology at the University of Oulu. Experimental design and laboratory work was performed by Ilse Ekman and Pekka Kilpeläinen and original data analysis by Oksana Suokas.

---

## Overview

```
Raw reads (AVITI + MiSeq)
        │
        ▼
  Primer trimming (cutadapt, pre-pipeline)
        │
        ▼
  Quality inspection → Filter & Trim
        │
        ▼
  Error learning → Denoising (DADA2)
        │
        ▼
  Read merging → Chimera removal
        │
        ▼
  Taxonomy assignment (SILVA)
        │
        ▼
  TreeSummarizedExperiment objects
        │
        ▼
  Merge platforms → Community analysis & comparison
```

Each platform is processed independently through the full DADA2 workflow. Results are stored in separate output directories before being merged for downstream comparative analyses.

---

## Prerequisites

### Software

- R (≥ 4.1 recommended)
- [Quarto](https://quarto.org/)
- [cutadapt](https://cutadapt.readthedocs.io/) — must be run on raw reads **before** this pipeline

### R Packages

| Package | Purpose |
|---|---|
| `dada2` | Core amplicon denoising workflow |
| `mia` | Microbiome data structures and analysis |
| `vegan` | Ordination and beta-diversity metrics |
| `Biostrings` | DNA sequence handling |
| `tidyverse` | Data wrangling and I/O |
| `ape` | Phylogenetic tree support |
| `kableExtra` | Formatted summary tables |
| `patchwork` | Multi-panel plot layouts |
| `ggthemes`, `ggrepel`, `ggpubr`, `ggsci`, `hrbrthemes` | Plotting utilities |
| `eulerr` | Euler/Venn diagrams for ASV overlap |
| `scater` | Single-cell/microbiome utilities |

```

### Reference Databases

Download the SILVA reference files and place them in `~/reference/`:

- `silva_nr99_v138.2_toGenus_trainset.fa.gz` — for genus-level taxonomy assignment
- `silva_v138.2_assignSpecies.fa.gz` — for species-level assignment

SILVA reference files are available from the [DADA2 taxonomy reference page](https://benjjneb.github.io/dada2/training.html).

---

## Project Structure

```
project/
├── dada2_pipeline.qmd          # Main analysis document
├── aviti_reads/                # Raw AVITI .fastq.gz files
│   └── filtered/               # Created automatically during pipeline
├── miseq_reads/                # Raw MiSeq .fastq.gz files
│   └── filtered/               # Created automatically during pipeline
├── aviti_metadata.tsv          # Sample metadata for AVITI run
├── miseq_metadata.tsv          # Sample metadata for MiSeq run
├── aviti_results/              # Output directory (created automatically)
└── miseq_results/              # Output directory (created automatically)
```
## Usage

### 1. Prepare reads

Trim primers from raw reads using cutadapt before running this pipeline. The pipeline assumes primers have already been removed.

### 2. Configure paths and parameters

At the top of each platform section in `dada2_pipeline.qmd`, update the path and parameter variables as needed:

```r
truncation <- c(260, 200)   # Forward, reverse truncation lengths
phi        <- TRUE           # Remove PhiX reads
```

Adjust `truncation` lengths based on the quality profiles visualised in the quality inspection step.

### 3. Render the document

```bash
quarto render dada2_pipeline.qmd
```

This produces an HTML report with all quality plots, summary tables, and comparative analyses.

> **Note:** Computationally intensive steps (quality profiling, filtering, error learning, denoising, merging, chimera removal, taxonomy assignment) are set to `eval = FALSE` to allow cached `.rds` intermediates to be reused on re-renders. Run these chunks interactively in RStudio on first use, or set `eval = TRUE` before rendering.

---

## Pipeline Steps

### Per-platform processing (AVITI and MiSeq)

1. **Quality inspection** — Aggregate quality profiles plotted across 200,000 reads; saved as `.rds` for re-use.
2. **Filter and trim** — `filterAndTrim()` with `maxN=0`, `maxEE=2`, `truncQ=2`; multithreaded.
3. **Error learning** — Error rate profiles learned separately for forward and reverse reads.
4. **Denoising** — `dada()` applied to filtered reads using the learned error model.
5. **Merge paired ends** — `mergePairs()` combines denoised forward and reverse reads.
6. **ASV table construction** — `makeSequenceTable()` builds the sample-by-sequence matrix.
7. **Chimera removal** — `removeBimeraDenovo()` with consensus method.
8. **Taxonomy assignment** — `assignTaxonomy()` to genus level + `addSpecies()` for exact species matches against SILVA v138.2.
9. **Object construction** — A `TreeSummarizedExperiment` (TSE) object is built incorporating counts, taxonomy, metadata, and representative sequences. Non-bacterial, mitochondrial, and chloroplastic ASVs are removed.
10. **Export** — TSE object, representative sequences (`.fasta`), taxonomy table (`.tsv`), ASV count table (`.tsv`), and metadata (`.tsv`) written to the results directory.

### Cross-platform comparison

After both platforms are processed, the pipeline:

- Merges AVITI and MiSeq TSE objects with a `Platform` metadata column.
- Performs **rarefaction** to normalise read depth across samples.
- Calculates **alpha diversity** (Shannon index, Observed ASVs) and correlates values between platforms.
- Performs **beta-diversity** analysis (Bray-Curtis PCoA) to assess community structure concordance.
- Identifies **shared and platform-specific ASVs** using Euler diagrams and abundance tables.
- Re-evaluates platform-specific ASVs after **filtering rare variants** (total abundance < 100), which removes the vast majority of platform-specific differences.

---

## Citation

If you use this workflow, please cite:

- Callahan BJ et al. (2016) DADA2: High-resolution sample inference from Illumina amplicon data. *Nature Methods*, 13:581–583.
- Callahan BJ et al. (2016) Bioconductor Workflow for Microbiome Data Analysis. *F1000Research*, 5:1492.
- Quandt LC et al. SILVA ribosomal RNA gene database project.
