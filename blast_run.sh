#!/bin/bash

THREADS=16

DB=/data/databases/25_02_04_ncbi_database/nt

BASE=04_dada2_output

for PRIMER in MCO MiBird MiMammal VERT
do

echo "=============================="
echo "Running BLAST for ${PRIMER}"
echo "=============================="

INPUT=${BASE}/${PRIMER}/${PRIMER}_zotu.fas
OUTDIR=${BASE}/${PRIMER}

mkdir -p "${OUTDIR}"

if [ ! -f "$INPUT" ]; then
    echo "ERROR: Input file not found: $INPUT"
    continue
fi

/usr/bin/blastn \
-query "$INPUT" \
-db "$DB" \
-out "${OUTDIR}/${PRIMER}_blast.tsv" \
-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids" \
-max_target_seqs 500 \
-evalue 1e-20 \
-perc_identity 80 \
-num_threads "$THREADS"

echo "Finished BLAST for ${PRIMER}"

done
