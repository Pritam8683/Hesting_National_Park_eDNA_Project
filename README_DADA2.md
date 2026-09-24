# DADA2 Processing README

## Purpose

This component processes paired-end amplicon reads after quality preprocessing and primer/barcode extraction. DADA2 is run independently for MCO, MiBird, MiMammal, and VERT.

## Processing steps

For each primer:

1. Locate primer-specific forward and reverse FASTQ files.
2. Filter reads with `filterAndTrim()`.
3. Generate quality-profile PDFs.
4. Learn forward and reverse error models.
5. Generate error-model PDFs.
6. Dereplicate reads.
7. Denoise forward and reverse reads.
8. Merge paired reads.
9. Construct the sequence table.
10. Remove chimeric sequences.
11. Generate read-tracking statistics.
12. Assign stable ASV IDs.
13. Export abundance and reference tables.

## Filtering parameters

Common settings:

```text
maxEE = c(2, 2)
truncQ = 2
maxN = 0
rm.phix = TRUE
multithread = 16
```

Primer-specific truncation:

| Primer | Forward | Reverse |
|---|---:|---:|
| MCO | 240 | 200 |
| VERT | 220 | 180 |
| MiMammal | 150 | 150 |
| MiBird | 150 | 150 |

These values should be validated using the generated quality-profile plots and expected amplicon length.

## Read tracking

`<PRIMER>_Read_Tracking.csv` contains:

- `Input`
- `Filtered`
- `DenoisedF`
- `DenoisedR`
- `Merged`
- `Nonchimera`

This table should be retained as the primary processing-QC record.

## ASV outputs

### ASV abundance table

`<PRIMER>_ASV_abundance_table.csv` contains one row per sample and one column per ASV.

### ASV reference table

`<PRIMER>_ASV_sequences.csv` contains:

```text
ASV_ID
Sequence
```

ASV identifiers are primer-specific, for example:

```text
MCO_ASV_1
MCO_ASV_2
...
```

### RDS object

`<PRIMER>_ASV_table.rds` stores the chimera-filtered DADA2 sequence table.

## Important QC checks

Before ecological analysis, inspect:

- fraction of reads retained after filtering;
- paired-read merging success;
- number of non-chimeric reads;
- number of ASVs generated;
- quality profiles;
- error-model plots;
- samples with unusually low read counts;
- blank/control samples.

Do not compare raw ASV counts between primers as though they represent directly comparable measures of biological abundance. Primer efficiency, amplification bias, copy number, and sequencing depth can differ among markers.
