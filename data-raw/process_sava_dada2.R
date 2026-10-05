# ============================================================
# data-raw/process_sava_dada2_single_end.R
# Single-end DADA2 workflow for Tele02 12S amplicons (PE256)
# Updated with correct Tele02 primer positions
# ============================================================

library(dada2)
library(Biostrings)

# --- Configuration ---
path <- "data-raw/sava_river/fastq"
filt_path <- file.path(path, "filtered")
if (!dir.exists(filt_path)) dir.create(filt_path, recursive = TRUE)

# Tele02 primer parameters (from Taberlet et al. 2018)
FWD_PRIMER <- "AAACTCGTGCCAGCCACC"
REV_PRIMER <- "GGGTATCTAATCCCAGTTTG"
TAG_LEN <- 8 # Sample identification tag
FWD_PRIMER_LEN <- as.numeric(nchar(FWD_PRIMER))
REV_PRIMER_LEN <- as.numeric(nchar(REV_PRIMER))

TRIM_LEFT_F <- TAG_LEN + FWD_PRIMER_LEN
TRIM_LEFT_R <- TAG_LEN + REV_PRIMER_LEN

# --- File organization ---
fnFs <- sort(list.files(path, pattern = "_1.fastq.gz$", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_2.fastq.gz$", full.names = TRUE))

if (length(fnFs) == 0) stop("No forward FASTQ files found.")
sample.names <- sub("_1.fastq.gz$", "", basename(fnFs))

# --- Filter and Trim ---
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(
  fnFs, filtFs, fnRs, filtRs,
  trimLeft = c(TRIM_LEFT_F, TRIM_LEFT_R),
  truncLen = c(176, 176),
  maxN = 0,
  maxEE = c(2,2),
  truncQ = 2,
  rm.phix = TRUE,
  compress = TRUE,
  multithread = TRUE
)
print(head(out))

# Check for samples that lost all reads
if (any(out[, 2] == 0)) {
  message("WARNING: Some samples had 0 reads after filtering:")
  print(rownames(out[out[, 2] == 0, ]))
}

# --- Learn Error Rates ---
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)

# --- Dereplication and Sample Inference ---
derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)
names(derepFs) <- sample.names
names(derepRs) <- sample.names

# Denoise 
dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
dadaRs <- dada(derepRs, err = errR, multithread = TRUE)

# Merge paired reads (this replaces FLASH)
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose = TRUE, minOverlap = 80)

# --- Construct Sequence Table ---
seqtab <- makeSequenceTable(mergers)

# --- Remove Chimeras ---
seqtab.nochim <- removeBimeraDenovo(
  seqtab, method = "consensus",
  multithread = TRUE, verbose = TRUE
)

# --- Apply Paper-Specific Filters ---

# Filter A: Length filter (150 bp to 190 bp)
# After trimming the 18bp primer, the remaining biological
# sequence should be around 167bp.
asv_seqs <- colnames(seqtab.nochim)
asv_lengths <- nchar(asv_seqs)
keep_len <- asv_lengths >= 150 & asv_lengths <= 190
seqtab.len <- seqtab.nochim[, keep_len, drop = FALSE]

# Filter B: Abundance and Prevalence filters
# Paper: "ASVs represented by fewer than 10 sequences or found
# in only a single sample were discarded"
multi_sample_mode <- (nrow(seqtab.len) > 1)
asv_totals <- colSums(seqtab.len)
keep_abund <- asv_totals >= 10

if (multi_sample_mode) {
  asv_prev <- colSums(seqtab.len > 0)
  keep_prev <- asv_prev > 1
  keep_final <- keep_abund & keep_prev
} else {
  message("Single sample detected. Skipping prevalence filter.")
  keep_final <- keep_abund
}

seqtab.final <- seqtab.len[, keep_final, drop = FALSE]
message(sprintf(
  "Final ASV table: %d samples x %d ASVs",
  nrow(seqtab.final), ncol(seqtab.final)
))

# --- Export Outputs ---
out_dir <- "data-raw/sava_river/processed"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Export ASV FASTA
final_asvs <- colnames(seqtab.final)
asv_ids <- paste0("ASV_", seq_along(final_asvs))
names(final_asvs) <- asv_ids

dna_set <- DNAStringSet(final_asvs)
writeXStringSet(
  dna_set,
  filepath = file.path(out_dir, "sava_asvs.fasta"),
  format = "fasta"
)

# Export Abundance Table (Rows = ASVs, Cols = Samples)
abund_df <- as.data.frame(t(seqtab.final))
abund_df$asv_id <- asv_ids
abund_df <- abund_df[, c("asv_id", setdiff(names(abund_df), "asv_id"))]

write.csv(
  abund_df,
  file.path(out_dir, "sava_abundance_table.csv"),
  row.names = FALSE
)

message("Exported sava_asvs.fasta and sava_abundance_table.csv to ", out_dir)
