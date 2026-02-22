DADA2 Pipeline for MiSeq and AVITI Reads
A Quarto-based R workflow for amplicon sequence variant (ASV) analysis from paired-end 16S rRNA sequencing data. The pipeline processes reads from two sequencing platforms — Element Biosciences AVITI and Illumina MiSeq — in parallel, then merges the results for a cross-platform comparison of community composition.
The workflow is adapted from the DADA2 pipeline tutorial and the Bioconductor Workflow for Microbiome Data Analysis by Benjamin Callahan.

Overview
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
Each platform is processed independently through the full DADA2 workflow. Results are stored in separate output directories before being merged for downstream comparative analyses.
