# Cutadapt Demultiplexing and Primer Removal README

## Purpose

Cutadapt is used after fastp to identify inline barcode + primer combinations and produce primer-specific paired-end FASTQ files for DADA2.

## Barcode + primer strategy

For each primer and barcode combination, the forward matching sequence is constructed as:

```text
barcode + forward primer
```

The reverse primer is supplied independently for reverse reads.

The current command uses:

```text
--revcomp
-e 0.1
-O 3
--minimum-length 150
--discard-untrimmed
```

and paired-end input/output.

## Output naming

Files are written as:

```text
<sample>_<barcode>_<primer>_R1.fastq
<sample>_<barcode>_<primer>_R2.fastq
```

Examples:

```text
Sample01_A_MCO_R1.fastq
Sample01_A_MCO_R2.fastq
```

## Important considerations

### Barcode assignment

The current script loops over all barcode values for every sample. This means each sample/barcode/primer combination is tested independently. Samples with no matching sequence are discarded because `--discard-untrimmed` is used.

### Primer matching

Primer sequences contain IUPAC ambiguity codes. Cutadapt interprets these ambiguity codes during matching.

### Reverse-complement behavior

`--revcomp` allows Cutadapt to detect adapters/primers in reverse-complement orientation. Orientation should nevertheless be checked against the sequencing library design.

### Minimum length

Reads shorter than 150 bp after trimming are discarded. This value should be consistent with the expected amplicon size and DADA2 truncation settings.

## QC recommendations

Inspect Cutadapt summaries for each sample, primer, and barcode. In particular, record:

- number of input read pairs;
- number of reads with an accepted barcode/primer match;
- number discarded as untrimmed;
- read lengths after trimming;
- unexpected differences among barcode groups.

If the barcode is expected to occur only in a particular read orientation or at a specific position, the matching strategy should be checked carefully before interpreting barcode-specific results.
