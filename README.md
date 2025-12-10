# epiTom-NatPop

_Pipeline for processing and analysing tomato natural population sequencing data._

This repository contains the full set of scripts and Makefiles used to process
raw sequencing data (`*.fq.gz`) and perform downstream analyses for the
**An epigenetic blueprint for tomato breeding** paper. The workflow covers:

- Raw read QC and organisation
- WGBS / bisulfite sequencing processing
- Pseudo–pan-genome and methylome comparison
- DMR calling and methylation VCF generation
- RNA-seq and phenotypic data processing
- Downstream analyses, figures and visualisation

The code is designed to use Makefiles to
parallelise jobs and keep sample-level processing under control.

---

## Repository structure

Top-level folders are numbered according to the main steps of the workflow:

- `01_raw-data/`  
  Organisation of raw sequencing files (`*.fq.gz`) and sample sheets.  
  No raw data are stored in this public repository; this folder contains
  only helper scripts and metadata.

- `02_get-pseudo-pan-genomes/`  
  Construction of pseudo–pan-genomes and reference-related preprocessing.

- `03_biseq-processing/`  
  Bisulfite / WGBS processing: trimming, mapping, deduplication, and
  methylation extraction.  
  Mostly driven by Makefiles to parallelise per sample / per chromosome.

- `04_methylome-comparison/`  
  Scripts to compare methylomes across samples / groups, e.g. differential
  methylation summaries, global patterns, QC.

- `05_DMR-processing/`  
  DMR calling and post-processing (filtering, merging, annotation helpers).

- `06_get-meth-vcf/`  
  Conversion of methylation calls into VCF-like formats and related utilities.

- `07_rnaseq-processing/`  
  Processing of RNA-seq data (alignment, quantification, normalisation, etc.).

- `08_fruit-processing/`  
  Scripts related to fruit phenotype data and related pre-processing.

- `09_KO-processing/`  
  Processing and comparison of knockout (KO) lines versus controls.

- `10_data-analysis/`  
  Downstream statistical analyses, GWAS/DMR integration, model fitting, etc.  
  This folder often contains R / Python scripts that assume inputs produced by
  steps 01–09.

- `11_figures/`  
  Code to generate the main and supplementary figures (R scripts, ggplot,
  etc.).

- `12_data-visualization/`  
  Additional visualisation utilities and exploratory plots.

Other useful directories:

- `bin/`  
  Helper scripts (bash, R, Python, etc.) used across multiple steps.

---

## Requirements
The pipeline was originally developed for a Linux cluster.
Typical dependencies include (non-exhaustive):

- **Core tools**: `bash`, `make`, `awk`, `sed`, `grep`
- **Sequencing tools**: e.g. `fastp`, `bismark`, `samtools`, `bedtools`,
  `bcftools`, `seqkit`, `plink`
- **R** (≥ 4.x) with packages such as `data.table`, `dplyr`, `tidyr`,
  `ggplot2`, `GenomicRanges`, `rtracklayer`, etc.
- **Python** (≥ 3.x) with standard scientific stack as needed

---

## Usage overview

1. **Clone the repository**

   ```bash
   git clone https://github.com/VeronicaNoe/epiTom-NatPop.git
   cd epiTom-NatPop
   ```

2. **Prepare your environment**
Load the necessary modules on your cluster, or
Activate the appropriate Conda environments.

3. **Organise raw data**
Place your raw *.fq.gz files and sample sheets according to the expected
layout in 01_raw-data/.

4.**Run the processing steps**
Each numbered directory contains one or more Makefiles. The typical
usage is:
  ```bash 
  cd 03_biseq-processing
  make -n          # dry-run: inspect the planned commands
  make -j 8        # run with 8 parallel jobs
  ```
Similar patterns apply to other steps (01–09), with targets defined per
sample / chromosome / context.

5. **Downstream analyses and figures**
Once all upstream processing is complete, use the scripts in:

- `10_data-analysis/` for statistical analyses and integration of methylation,
genetics, and phenotypes.

- `11_figures/` for reproducing the figures.

- `12_data-visualization/` for additional plots.

---

## Data availability

New WGBS data has been deposited in the European Nucleotide Archive (ENA) repository PRJEB80801. 
New WGS and RNAseq data has been deposited in the Short-Read Archive (SRA) BioProject ID PRJNA118967. 
The VCF file used for variant-based analyses has been deposited in the European Variation Archive (EVA) under accession number PRJEB82793.
Publicly available WGBS and RNAseq data reanalyzed in this study has been obtained from Genome Sequence Archive (GSA) at the Big Data Center project accession PRJCA009995 (CRA007189; http://bigd.big.ac.cn/gsa), 
the ENA project accession PRJNA397191 and PRJNA516166; and NCBI Short Read Archive (SRA) accessions SRA046092, SRA046132, SRA046131, SRA053345, and SRA046480.

---

## Citation

If you use this code or parts of this pipeline in your own work, please cite:
The associated epiTom-NatPop publication (once available), and
The relevant tools and methods used by the pipeline (bismark, samtools, etc.).
Citation details will be added here when the manuscript is published.

---

## Contact

For questions about the pipeline or requests related to the project, please
contact:

[Veronica N. Ibanez] – [veronicanoeibanez@gmail.com]

---
