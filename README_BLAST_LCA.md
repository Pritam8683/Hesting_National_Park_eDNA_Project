# BLAST and NCBI LCA Taxonomy README

## Purpose

This component assigns taxonomy to DADA2 ASVs using BLASTn against the NCBI `nt` database followed by an NCBI taxonomy lookup and a lowest-common-ancestor (LCA) procedure.

## BLAST database

The supplied shell script uses:

```text
/data/databases/25_02_04_ncbi_database/nt
```

The database date/version should be recorded because BLAST and taxonomy results can change as NCBI databases are updated.

## BLAST settings

The current workflow uses:

```text
-max_target_seqs 500
-evalue 1e-20
-perc_identity 80
-num_threads 16
```

The output contains:

```text
qseqid
sseqid
pident
length
mismatch
gapopen
qstart
qend
sstart
send
evalue
bitscore
staxids
```

## BLAST output

For each primer:

```text
04_dada2_output/<PRIMER>/<PRIMER>_blast.tsv
```

is generated from:

```text
04_dada2_output/<PRIMER>/<PRIMER>_zotu.fas
```

## Taxonomy extraction

The R script uses `taxizedb` and the NCBI taxonomy database to retrieve:

- Kingdom
- Phylum
- Class
- Order
- Family
- Genus
- Species

For every ASV, taxonomy is retrieved for the unique taxids present in the BLAST results.

## LCA logic

The current implementation compares taxonomic names across the retrieved BLAST taxids, beginning at Kingdom and proceeding downward. Once more than one taxonomic value is encountered at a rank, assignment stops at that rank and all lower ranks remain `NA`.

This provides a conservative LCA-style assignment across the taxids represented in the BLAST input.

## Best-hit statistics

The output also records:

- `Best_Identity` – maximum percent identity among the retained BLAST hits;
- `Best_Evalue` – minimum E-value;
- `Best_Bitscore` – maximum bit score.

These statistics describe the BLAST hit set but should not by themselves be treated as taxonomic confidence scores.

## Output

For each primer:

```text
<PRIMER>_LCA_taxonomy.csv
```

with columns:

```text
ASV_ID
Kingdom
Phylum
Class
Order
Family
Genus
Species
Best_Identity
Best_Evalue
Best_Bitscore
```

## Important interpretation note

The current LCA script uses all unique `staxids` present in the BLAST file for an ASV after the BLAST filtering parameters have been applied. For publication-quality taxonomic assignment, consider an explicit hit-selection strategy before LCA, such as filtering by alignment coverage, identity, E-value, bitscore, and/or a narrow top-hit window. Otherwise, biologically weak or divergent hits can influence the LCA.

Likewise, species-level assignment should be based on marker-specific validation rather than a single universal identity threshold. Short or highly conserved COI fragments may not provide enough resolution for reliable species identification.

## Recommended downstream QC

For each ASV, retain the original BLAST hits so that taxonomy can be audited. Review:

1. top BLAST matches;
2. alignment length and query coverage;
3. percent identity;
4. E-value and bitscore;
5. taxonomic agreement among strong hits;
6. unexpected host/non-target assignments;
7. negative controls and contaminant ASVs.

Taxonomy should be considered an inference from sequence similarity and database coverage, not an observation independent of the reference database.
