# =========================================================
# MULTI-PRIMER DADA2 METABARCODING PIPELINE
# fastp + Cutadapt + DADA2
# =========================================================

# =========================================================
# LOAD LIBRARIES
# =========================================================

library(dada2)
library(Biostrings)
library(ShortRead)
library(tidyverse)
library(DECIPHER)

# =========================================================
# SET WORKING DIRECTORY
# =========================================================

setwd("/home/pritam/SERDP/Spider_web_hs")

# =========================================================
# THREADS
# =========================================================

threads <-  16

# =========================================================
# DIRECTORY STRUCTURE
# =========================================================

raw_dir    <- "00_data_raw/sequences"
fastp_dir  <- "01_fastp"
demux_dir  <- "02_demultiplexed"
output_dir <- "04_dada2_output"

dir.create(fastp_dir,
           recursive = TRUE,
           showWarnings = FALSE)

dir.create(demux_dir,
           recursive = TRUE,
           showWarnings = FALSE)

dir.create(output_dir,
           recursive = TRUE,
           showWarnings = FALSE)

# =========================================================
# INLINE BARCODES
# =========================================================

barcodes <- c(
  A = "ATCACG",
  B = "CGATGT",
  C = "TTAGGC"
)

# =========================================================
# PRIMERS
# =========================================================

primer_fwd <- c(

  MCO = "GGWACWGGWTGAACWGTWTAYCCYCC",

  MiBird  = "GGGTTGGTAAATCTTGTGCCAGC",

MiMammal = "GGGTTGGTAAATTTCGTGCCAGC",

VERT = "TNTTYTCMACYAACCACAAAGA"
)

primer_rev <- c(

  MCO = "TANACYTCNGGRTGNCCRAARAAYCA",

  MiBird = "CATAGTGGGGTATCTAATCCCAGTTTG",

MiMammal = "CATAGTGGGGTATCTAATCCCAGTTTG",

VERT = "CARAAGCTYATGTTRTTYATDCG"

)

# =========================================================
# FIND RAW FASTQ FILES
# =========================================================

all_fastq <- list.files(
  raw_dir,
  pattern = "\\.fastq",
  full.names = TRUE
)

fnFs <- sort(
  all_fastq[
    grepl("_R1_", basename(all_fastq))
  ]
)

fnRs <- sort(
  all_fastq[
    grepl("_R2_", basename(all_fastq))
  ]
)

sample.names <- sub(
  "_R1_.*",
  "",
  basename(fnFs)
)

cat("Forward files:", length(fnFs), "\n")
cat("Reverse files:", length(fnRs), "\n")

# =========================================================
# STEP 1 — FASTP
# =========================================================

for(i in seq_along(fnFs)){

  sample_name <- sample.names[i]

  outF <- file.path(
    fastp_dir,
    paste0(sample_name,
           "_R1_fastp.fastq")
  )

  outR <- file.path(
    fastp_dir,
    paste0(sample_name,
           "_R2_fastp.fastq")
  )

  html_report <- file.path(
    fastp_dir,
    paste0(sample_name,
           "_fastp.html")
  )

  json_report <- file.path(
    fastp_dir,
    paste0(sample_name,
           "_fastp.json")
  )

  if(file.exists(outF) &&
     file.exists(outR)){

    cat("Skipping fastp:",
        sample_name,
        "\n")

    next
  }

  cmd <- paste(

    "fastp",

    "-i", fnFs[i],
    "-I", fnRs[i],

    "-o", outF,
    "-O", outR,

    "--detect_adapter_for_pe",

    "--cut_front",
    "--cut_tail",

    "--cut_window_size 4",

    "--cut_mean_quality 30",

    "--length_required 150",

    "--thread", threads,

    "--html", html_report,

    "--json", json_report
  )

  cat("Running fastp:",
      sample_name,
      "\n")

  system(cmd)

  cat("fastp completed:",
      sample_name,
      "\n")
}

# =========================================================
# STEP 2 — CUTADAPT
# DEMULTIPLEX + PRIMER REMOVAL
# =========================================================

fastpFs <- sort(list.files(
  fastp_dir,
  pattern = "_R1_fastp.fastq$",
  full.names = TRUE
))

fastpRs <- sort(list.files(
  fastp_dir,
  pattern = "_R2_fastp.fastq$",
  full.names = TRUE
))

for(primer_name in names(primer_fwd)){

  cat("\n")
  cat("====================================\n")
  cat("PROCESSING PRIMER:", primer_name, "\n")
  cat("====================================\n")

  for(i in seq_along(fastpFs)){

    sample_name <- sub(
      "_R1_fastp.fastq",
      "",
      basename(fastpFs[i])
    )

    for(bc in names(barcodes)){

      barcode_seq <- barcodes[bc]

      outF <- file.path(
        demux_dir,
        paste0(sample_name,
               "_",
               bc,
               "_",
               primer_name,
               "_R1.fastq")
      )

      outR <- file.path(
        demux_dir,
        paste0(sample_name,
               "_",
               bc,
               "_",
               primer_name,
               "_R2.fastq")
      )

      if(file.exists(outF) &&
         file.exists(outR)){

        cat("Skipping:",
            sample_name,
            bc,
            primer_name,
            "\n")

        next
      }

      FWD_FULL <- paste0(
        barcode_seq,
        primer_fwd[primer_name]
      )

      REV_FULL <- primer_rev[primer_name]

      cmd <- paste(

        "cutadapt",

        "-j", threads,

        "--revcomp",

        "-e 0.1",

        "-O 3",

        "--minimum-length 150",

        "--discard-untrimmed",

        "-g", FWD_FULL,

        "-G", REV_FULL,

        "-o", outF,

        "-p", outR,

        fastpFs[i],
        fastpRs[i]
      )

      cat("Running:",
          sample_name,
          bc,
          primer_name,
          "\n")

      system(cmd)

      cat("Completed:",
          sample_name,
          bc,
          primer_name,
          "\n")
    }
  }
}

# =========================================================
# STEP 3 — DADA2 PER PRIMER
# =========================================================

for(primer_name in names(primer_fwd)){

  cat("\n")
  cat("====================================\n")
  cat("RUNNING DADA2:", primer_name, "\n")
  cat("====================================\n")

  # =======================================================
  # CREATE PRIMER-SPECIFIC DIRECTORIES
  # =======================================================

  primer_filter_dir <- file.path(
    output_dir,
    primer_name,
    "filtered"
  )

  primer_output_dir <- file.path(
    output_dir,
    primer_name
  )

  dir.create(primer_filter_dir,
             recursive = TRUE,
             showWarnings = FALSE)

  dir.create(primer_output_dir,
             recursive = TRUE,
             showWarnings = FALSE)

  # =======================================================
  # GET PRIMER FILES
  # =======================================================

  demuxFs <- sort(list.files(

    demux_dir,

    pattern = paste0(
      "_",
      primer_name,
      "_R1.fastq$"
    ),

    full.names = TRUE
  ))

  demuxRs <- sort(list.files(

    demux_dir,

    pattern = paste0(
      "_",
      primer_name,
      "_R2.fastq$"
    ),

    full.names = TRUE
  ))

  if(length(demuxFs) == 0){

    cat("No files found for",
        primer_name,
        "\n")

    next
  }

  sample.names <- sub(
    "_R1.fastq$",
    "",
    basename(demuxFs)
  )

  cat("Samples found:",
      length(sample.names),
      "\n")

  # =======================================================
  # FILTER FILE PATHS
  # =======================================================

  filtFs <- file.path(
    primer_filter_dir,
    paste0(sample.names,
           "_F_filt.fastq")
  )

  filtRs <- file.path(
    primer_filter_dir,
    paste0(sample.names,
           "_R_filt.fastq")
  )
# =======================================================
# PRIMER-SPECIFIC FILTER SETTINGS
# =======================================================

if(primer_name == "MCO"){

  trunc_forward <- 240
  trunc_reverse <- 200

} else if(primer_name == "VERT"){

  trunc_forward <- 220
  trunc_reverse <- 180

} else if(primer_name == "MiMammal"){

  trunc_forward <- 150
  trunc_reverse <- 150

} else if(primer_name == "MiBird"){

  trunc_forward <- 150
  trunc_reverse <- 150
}

  # =======================================================
  # FILTERING
  # =======================================================

  filter_out <- filterAndTrim(

    demuxFs,
    filtFs,

    demuxRs,
    filtRs,

    truncLen = c(
  trunc_forward,
  trunc_reverse
),

    maxEE = c(2,2),

    truncQ = 2,

    maxN = 0,

    rm.phix = TRUE,

    compress = TRUE,

    multithread = threads
  )

  # =======================================================
  # QUALITY PLOTS
  # =======================================================

  pdf(file.path
    (primer_output_dir,
    paste0(primer_name,
           "_Quality_Forward.pdf")
  ))

  plotQualityProfile(filtFs[1:min(2,length(filtFs))])

  dev.off()

  pdf(file.path(
    primer_output_dir,
    paste0(primer_name,
           "_Quality_Reverse.pdf")
  ))

  plotQualityProfile(filtRs[1:min(2,length(filtRs))])

  dev.off()

  # =======================================================
  # LEARN ERRORS
  # =======================================================

  errF <- learnErrors(
    filtFs,
    multithread = threads
  )

  errR <- learnErrors(
    filtRs,
    multithread = threads
  )

  # =======================================================
  # SAVE ERROR PLOTS
  # =======================================================

  pdf(file.path(
    primer_output_dir,
    paste0(primer_name,
           "_Error_Forward.pdf")
  ))

  plotErrors(errF,
             nominalQ = TRUE)

  dev.off()

  pdf(file.path(
    primer_output_dir,
    paste0(primer_name,
           "_Error_Reverse.pdf")
  ))

  plotErrors(errR,
             nominalQ = TRUE)

  dev.off()

  # =======================================================
  # DEREPLICATION
  # =======================================================

  derepFs <- derepFastq(filtFs)
  derepRs <- derepFastq(filtRs)

  names(derepFs) <- sample.names
  names(derepRs) <- sample.names

  # =======================================================
  # DENOISING
  # =======================================================

  dadaFs <- dada(
    derepFs,
    err = errF,
    multithread = threads
  )

  dadaRs <- dada(
    derepRs,
    err = errR,
    multithread = threads
  )

  # =======================================================
  # MERGE PAIRED READS
  # =======================================================

  mergers <- mergePairs(
    dadaFs,
    derepFs,
    dadaRs,
    derepRs
  )

  # =======================================================
  # MAKE SEQUENCE TABLE
  # =======================================================

  seqtab <- makeSequenceTable(
    mergers
  )

  # =======================================================
  # REMOVE CHIMERAS
  # =======================================================

  seqtab.nochim <- removeBimeraDenovo(

    seqtab,

    method = "consensus",

    multithread = threads
  )

  # =======================================================
  # READ TRACKING
  # =======================================================

  getN <- function(x)
    sum(getUniques(x))

  track <- data.frame(

    SampleID = sample.names,

    Input = filter_out[,1],

    Filtered = filter_out[,2],

    DenoisedF = sapply(dadaFs, getN),

    DenoisedR = sapply(dadaRs, getN),

    Merged = sapply(mergers, getN),

    Nonchimera = rowSums(seqtab.nochim)
  )

  write.csv(
    track,
    file.path(
      primer_output_dir,
      paste0(primer_name,
             "_Read_Tracking.csv")
    ),
    row.names = FALSE
  )

  # =======================================================
  # CREATE ASV IDS
  # =======================================================

  asv_seqs <- colnames(
    seqtab.nochim
  )

  asv_headers <- paste0(
    primer_name,
    "_ASV_",
    seq_len(length(asv_seqs))
  )

  colnames(seqtab.nochim) <- asv_headers

  # =======================================================
  # ASV REFERENCE TABLE
  # =======================================================

  asv_reference <- data.frame(

    ASV_ID = asv_headers,

    Sequence = asv_seqs
  )

  # =======================================================
  # CLEAN ABUNDANCE TABLE
  # =======================================================

  asv_table <- as.data.frame(
    seqtab.nochim
  )

  asv_table <- cbind(
    SampleID = rownames(asv_table),
    asv_table
  )

  rownames(asv_table) <- NULL

  # =======================================================
  # EXPORT TABLES
  # =======================================================

  write.csv(

    asv_table,

    file.path(
      primer_output_dir,
      paste0(primer_name,
             "_ASV_abundance_table.csv")
    ),

    row.names = FALSE
  )

  write.csv(

    asv_reference,

    file.path(
      primer_output_dir,
      paste0(primer_name,
             "_ASV_sequences.csv")
    ),

    row.names = FALSE
  )

  saveRDS(

    seqtab.nochim,

    file.path(
      primer_output_dir,
      paste0(primer_name,
             "_ASV_table.rds")
    )
  )

  cat("DADA2 COMPLETE:",
      primer_name,
      "\n")
}

cat("\n")
cat("====================================\n")
cat("PIPELINE COMPLETED SUCCESSFULLY\n")
cat("====================================\n")
