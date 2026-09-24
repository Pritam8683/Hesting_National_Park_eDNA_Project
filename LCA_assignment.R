# =========================================================
# LCA TAXONOMY ASSIGNMENT PIPELINE
# BLAST + NCBI TAXONOMY + LCA
# =========================================================

# =========================================================
# LOAD LIBRARIES
# =========================================================

library(tidyverse)
library(data.table)
library(taxizedb)

# =========================================================
# SETTINGS
# =========================================================

base_dir <- "04_dada2_output"

primers <- c("MCO", "MiBird", "MiMammal", "VERT")

# =========================================================
# SAFE TAXONOMY FUNCTION
# =========================================================

get_taxonomy <- function(taxid){

  ranks_needed <- c(
    "kingdom",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species"
  )

  empty_return <- rep(NA, length(ranks_needed))

  # -------------------------------------------------------
  # CHECK EMPTY TAXID
  # -------------------------------------------------------

  if(is.na(taxid) || taxid == ""){
    return(empty_return)
  }

  taxid <- as.character(taxid)

  # -------------------------------------------------------
  # GET TAXONOMY
  # -------------------------------------------------------

  taxonomy <- tryCatch({

    classification(
      taxid,
      db = "ncbi"
    )[[1]]

  }, error = function(e) NULL)

  # -------------------------------------------------------
  # HANDLE FAILED TAXONOMY
  # -------------------------------------------------------

  if(is.null(taxonomy)){
    return(empty_return)
  }

  if(is.atomic(taxonomy)){
    return(empty_return)
  }

  if(!is.data.frame(taxonomy)){
    return(empty_return)
  }

  if(nrow(taxonomy) == 0){
    return(empty_return)
  }

  # -------------------------------------------------------
  # EXTRACT RANKS
  # -------------------------------------------------------

  out <- rep(NA, length(ranks_needed))

  for(i in seq_along(ranks_needed)){

    rank_match <- taxonomy[
      taxonomy$rank == ranks_needed[i],
    ]

    if(nrow(rank_match) > 0){

      out[i] <- tail(
        rank_match$name,
        1
      )
    }
  }

  return(out)
}

# =========================================================
# PROCESS EACH PRIMER
# =========================================================

for(primer in primers){

  cat("\n")
  cat("=====================================\n")
  cat("PROCESSING:", primer, "\n")
  cat("=====================================\n")

  # =======================================================
  # INPUT FILE
  # =======================================================

  blast_file <- file.path(
    base_dir,
    primer,
    paste0(primer, "_blast.tsv")
  )

  if(!file.exists(blast_file)){

    cat("BLAST FILE NOT FOUND\n")
    next
  }

  # =======================================================
  # READ BLAST
  # =======================================================

  blast <- fread(
    blast_file,
    header = FALSE
  )

  colnames(blast) <- c(
    "ASV_ID",
    "sseqid",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
    "staxids"
  )

  cat("Total BLAST hits:",
      nrow(blast),
      "\n")

  # =======================================================
  # UNIQUE ASVS
  # =======================================================

  asv_ids <- unique(
    blast$ASV_ID
  )
  taxonomy_results <- list()

  # =======================================================
  # PROCESS EACH ASV
  # =======================================================

  for(i in seq_along(asv_ids)){

    asv <- asv_ids[i]

    cat(
      "Processing ASV:",
      i,
      "/",
      length(asv_ids),
      "\r"
    )

    # -----------------------------------------------------
    # GET ASV HITS
    # -----------------------------------------------------

    sub <- blast %>%
      filter(ASV_ID == asv)

    # -----------------------------------------------------
    # UNIQUE TAXIDS
    # -----------------------------------------------------

    taxids <- unique(
      sub$staxids
    )

    # REMOVE EMPTY TAXIDS

    taxids <- taxids[
      !is.na(taxids)
    ]

    if(length(taxids) == 0){

      taxonomy_results[[i]] <- data.frame(

        ASV_ID = asv,

        Kingdom = NA,
        Phylum = NA,
        Class = NA,
        Order = NA,
        Family = NA,
        Genus = NA,
        Species = NA,

        Best_Identity = NA,
        Best_Evalue = NA,
        Best_Bitscore = NA
      )

      next
    }

    # -----------------------------------------------------
    # GET TAXONOMY FOR ALL HITS
    # -----------------------------------------------------

    tax_list <- lapply(
      taxids,
      get_taxonomy
    )

    # -----------------------------------------------------
    # COMBINE TAXONOMY
    # -----------------------------------------------------

    tax_matrix <- do.call(
      rbind,
      tax_list
    )

    # -----------------------------------------------------
    # LCA ASSIGNMENT
    # -----------------------------------------------------

    lca <- rep(NA, 7)

    for(j in 1:7){

      vals <- unique(
        na.omit(tax_matrix[,j])
      )

      if(length(vals) == 1){

        lca[j] <- vals

      } else {

        break
      }
    }

    # -----------------------------------------------------
    # BEST HIT STATS
    # -----------------------------------------------------

    best_identity <- max(
      sub$pident,
      na.rm = TRUE
    )

    best_evalue <- min(
      sub$evalue,
      na.rm = TRUE
    )

    best_bitscore <- max(
      sub$bitscore,
      na.rm = TRUE
    )

    # -----------------------------------------------------
    # STORE RESULT
    # -----------------------------------------------------

    taxonomy_results[[i]] <- data.frame(

      ASV_ID = asv,

      Kingdom = lca[1],
      Phylum = lca[2],
      Class = lca[3],
      Order = lca[4],
      Family = lca[5],
      Genus = lca[6],
      Species = lca[7],

      Best_Identity = best_identity,
      Best_Evalue = best_evalue,
      Best_Bitscore = best_bitscore
    )
  }

  # =======================================================
  # COMBINE RESULTS
  # =======================================================

  taxonomy_table <- bind_rows(
    taxonomy_results
  )

  # =======================================================
  # EXPORT
  # =======================================================

  output_file <- file.path(
    base_dir,
    primer,
    paste0(
      primer,
      "_LCA_taxonomy.csv"
    )
  )

  write.csv(

    taxonomy_table,

    output_file,

    row.names = FALSE
  )

  cat("\n")
  cat("LCA COMPLETE:", primer, "\n")
  cat("Output:", output_file, "\n")
}

# =========================================================
# FINISHED
# =========================================================

cat("\n")
cat("=====================================\n")
cat("ALL LCA ASSIGNMENTS COMPLETE\n")
cat("=====================================\n")
