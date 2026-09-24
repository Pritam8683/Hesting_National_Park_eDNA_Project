# Multi-Primer eDNA Metabarcoding Pipeline

This repository contains a reproducible workflow for processing paired-end amplicon sequencing data using **fastp → Cutadapt → DADA2 → BLAST → NCBI taxonomy/LCA**.

The pipeline is designed for four primer sets:

- **MCO** – mitochondrial COI marker
- **MiBird** – bird COI marker
- **MiMammal** – mammal COI marker
- **VERT** – vertebrate COI marker

## Workflow

```text
Raw paired-end FASTQ
        |
        v
      fastp
        |
        v
   Cutadapt
 barcode + primer matching / demultiplexing / primer removal
        |
        v
      DADA2
 filtering → error learning → denoising → paired-read merging
 → ASV table → chimera removal
        |
        v
 ASV FASTA / abundance tables
        |
        v
      BLASTn
 against NCBI nt database
        |
        v
 NCBI taxonomy + LCA
        |
        v
 Final ASV taxonomy tables
```

## Directory structure

```text
Spider_web_hs/
├── 00_data_raw/
│   └── sequences/              # Raw paired-end FASTQ files
├── 01_fastp/                   # fastp-filtered reads and QC reports
├── 02_demultiplexed/           # Cutadapt primer/barcode outputs
├── 04_dada2_output/
│   ├── MCO/
│   ├── MiBird/
│   ├── MiMammal/
│   └── VERT/
└── README.md
```

## Software requirements

The R pipeline requires:

```r
library(dada2)
library(Biostrings)
library(ShortRead)
library(tidyverse)
library(DECIPHER)
library(data.table)
library(taxizedb)
```

Command-line programs:

- `fastp`
- `cutadapt`
- `blastn`

The BLAST database used by the supplied script is:

```text
/data/databases/25_02_04_ncbi_database/nt
```

## Computational settings

Default number of threads:

```text
16
```

Working directory:

```text
/home/pritam/SERDP/Spider_web_hs
```

## Input naming convention

The raw FASTQ files are expected to contain `_R1_` and `_R2_` in their names. The pipeline identifies forward and reverse reads using these patterns.

Example:

```text
Sample01_R1_001.fastq
Sample01_R2_001.fastq
```

The sample name is extracted from the R1 filename by removing everything from `_R1_` onward.

## Barcode configuration

The current inline barcodes are:

| Barcode | Sequence |
|---|---|
| A | ATCACG |
| B | CGATGT |
| C | TTAGGC |

Cutadapt searches for each barcode immediately upstream of the primer sequence.

## Primer configuration

The primer sequences are defined in the main R script. IUPAC ambiguity codes are retained, allowing degenerate primer bases to be represented directly.

## Main outputs

For each primer, DADA2 generates:

```text
04_dada2_output/<PRIMER>/
├── filtered/
├── <PRIMER>_Quality_Forward.pdf
├── <PRIMER>_Quality_Reverse.pdf
├── <PRIMER>_Error_Forward.pdf
├── <PRIMER>_Error_Reverse.pdf
├── <PRIMER>_Read_Tracking.csv
├── <PRIMER>_ASV_abundance_table.csv
├── <PRIMER>_ASV_sequences.csv
└── <PRIMER>_ASV_table.rds
```

The downstream BLAST/LCA workflow additionally produces:

```text
<PRIMER>_blast.tsv
<PRIMER>_LCA_taxonomy.csv
```

## Important reproducibility notes

Record the versions of R, DADA2, fastp, Cutadapt, BLAST+, taxizedb, and the NCBI database used for every analysis. The NCBI taxonomy database and `nt` database are time-dependent resources, so the database date should be retained with the final results.

The current workflow uses primer-specific truncation lengths and a minimum read length of 150 bp. These parameters should be checked against the observed quality profiles and expected amplicon lengths before each sequencing run is analyzed.

## Recommended execution order

1. Place raw FASTQ files in `00_data_raw/sequences/`.
2. Run the fastp/Cutadapt/DADA2 R script.
3. Export or generate the ASV FASTA files required by BLAST.
4. Run the BLAST shell script.
5. Run the NCBI taxonomy/LCA R script.
6. Inspect read tracking, quality plots, error plots, BLAST results, and taxonomy assignments before downstream ecological analyses.
