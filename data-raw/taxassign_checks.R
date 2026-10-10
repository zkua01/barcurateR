# ============================================================
# evaluate_assignment_vs_article.R  (Part 1: audits)
# Input: RAW barcurateR two-pass assignment (pre-audit).
# Output: audited assignment + evidence tables for comparison.
# ============================================================
library(barcurateR)
library(dplyr)
library(tidyr)
library(readr)
library(purrr)
library(Biostrings)
library(DECIPHER)
library(ape)
library(dada2)

db_path    <- "data-raw/slovenian_fish/processed/slovenian_fish_refdb.sqlite"
asv_fasta  <- "data-raw/sava_river/processed/sava_asvs.fasta"
counts_csv <- "data-raw/sava_river/processed/sava_abundance_table.csv"
assign_csv <- "inst/extdata/sava_river_12s_assignments_two_pass.csv"  # RAW, pre-audit
site_metadata <- "data-raw/sava_site_metadata.csv"
spp_site_csv <- "inst/extdata/sava_river_community_matrix.csv"
audit_dir  <- "exports/evaluation"; dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)
FWD <- "AAACTCGTGCCAGCCACC"; RCREV <- "CAAACTGGGATTAGATACCC"
force_rebuild <- FALSE

sava_asvs <- readDNAStringSet(asv_fasta); asv_char <- as.character(sava_asvs)
sava_counts <- read_csv(counts_csv, show_col_types = FALSE)
assign <- read_csv(assign_csv, show_col_types = FALSE)   # RAW assignment
sava_site <- read_csv(site_metadata, show_col_types = FALSE)
sava_spp_site <- read_csv(spp_site_csv, show_col_types = FALSE)
site_names <- sava_site$sample_name[match(colnames(sava_spp_site), sava_site$sample_id)]
colnames(sava_spp_site) <- c(colnames(sava_spp_site)[1], site_names[-1])
con <- rb_connect(db_path)
all_refs <- rb_get_sequences(con, qc_flag = NULL)
all_12S  <- all_refs |> filter(seq_type == "12S")
sample_cols <- setdiff(names(sava_counts), "asv_id")
cnt_mat     <- as.matrix(sava_counts[, sample_cols, drop = FALSE])

abund <- tibble(
  asv_id    = sava_counts$asv_id,
  total     = rowSums(cnt_mat),
  n_samples = rowSums(cnt_mat > 0)
)

clean_for_align <- function(seqs) gsub("[-.]","",toupper(seqs))

pairdist <- function(seqs, ids) {
  aln <- DECIPHER::AlignSeqs(DNAStringSet(setNames(clean_for_align(seqs), ids)), verbose = FALSE)
  as.matrix(ape::dist.dna(as.DNAbin(aln), model = "raw", pairwise.deletion = TRUE))
}

# ---- 0. Provision overhang table (ASV vs 12S pass pool) ----
overhang_csv <- file.path(audit_dir, "asv_best_hits_overhangs.csv")
ref_fasta    <- file.path(audit_dir, "ref_12S_pass.fasta")
if (force_rebuild || !file.exists(overhang_csv) || !file.exists(ref_fasta)) {
  rp <- rb_get_sequences(con, marker = "12S", qc_flag = "pass")
  writeXStringSet(DNAStringSet(setNames(rp$sequence, rp$sequence_id)), ref_fasta)
  tsv <- file.path(audit_dir, "asv_vs_12S_coords.tsv")
  system2("blastn", c("-query ", shQuote(asv_fasta), "-subject ", shQuote(ref_fasta),
                      "-outfmt ", shQuote("6 qseqid sseqid pident length bitscore qstart qend sstart send qlen slen"),
                      "-max_target_seqs ", "5", "-evalue ", "1e-10", "-out ", shQuote(tsv)))
  h <- read.delim(tsv, header = FALSE, col.names = c("qseqid","sseqid","pident","length",
                                                     "bitscore","qstart","qend","sstart","send","qlen","slen"))
  write_csv(h |> group_by(qseqid) |> slice_max(bitscore, n = 1, with_ties = FALSE) |>
              ungroup() |> mutate(ov5 = qstart - 1, ov3 = qlen - qend), overhang_csv)
}
best <- read_csv(overhang_csv, show_col_types = FALSE)

# ---- 1. ASV-level audit ----
tot <- abund$total[match(names(asv_char), abund$asv_id)]; tot[is.na(tot)] <- 1L
bimera   <- isBimeraDenovo(setNames(as.integer(tot), asv_char))
chim_ids <- names(asv_char)[asv_char %in% names(bimera)[bimera]]

chimera_scan <- function(q, pool_seqs, pool_ids, min_seg = 60) {
  if (length(q) != 1 || is.na(q) || nchar(q) <= 2 * min_seg) return(tibble())
  keep <- !is.na(pool_seqs); ps <- pool_seqs[keep]; pid <- pool_ids[keep]
  hits <- list()
  for (i in min_seg:(nchar(q) - min_seg)) {
    L <- substr(q, 1, i); R <- substr(q, i + 1, nchar(q))
    lm <- substr(ps, 1, i) == L;            lm[is.na(lm)]  <- FALSE
    rm_ <- substr(ps, i + 1, nchar(q)) == R; rm_[is.na(rm_)] <- FALSE
    l_id <- pid[lm]; r_id <- pid[rm_]
    if (length(l_id) && length(r_id) && !is.na(l_id[1]) && !is.na(r_id[1]) && l_id[1] != r_id[1])
      hits[[length(hits) + 1]] <- tibble(breakpoint = i, left_parent = l_id[1], right_parent = r_id[1])
  }
  bind_rows(hits)
}
suspect_ids <- union(names(asv_char)[grepl(FWD, asv_char, fixed = TRUE)],
                     union(chim_ids, best$qseqid[best$ov3 >= 10]))
pool_ids  <- setdiff(names(asv_char), suspect_ids)
pool_seqs <- c(asv_char[pool_ids],
               setNames(all_12S$sequence, paste0("REF:", all_12S$sequence_id)))
pool_ids  <- c(pool_ids, paste0("REF:", all_12S$sequence_id))
bp <- map_dfr(suspect_ids, function(id)
  chimera_scan(asv_char[id], pool_seqs, pool_ids) |> mutate(asv_id = id))
bp_pairs <- bp |> count(asv_id, left_parent, right_parent)

rescue <- read.delim(list.files("exports/blastn_pass2", full.names = TRUE)[1],
                     header = FALSE)  # adjust col.names to rb_assign_edna outfmt
rescue_best <- rescue |> group_by(V1) |> slice_max(V5, n = 1, with_ties = FALSE)  # V1=qseqid, V5=bitscore
rescue_typ <- rescue_best |>
  left_join(rb_get_sequences(con, qc_flag = NULL) |> select(sequence_id, seq_type),
            by = c("V2" = "sequence_id")) |>
  select(asv_id = V1, rescue_via = seq_type)

asv_audit <- tibble(asv_id = names(asv_char)) |>
  left_join(abund, by = "asv_id") |>
  left_join(assign |> select(asv_id, assignment_status, species, pident, db_source), by = "asv_id") |>
  left_join(rescue_typ, by = "asv_id") |>
  mutate(len = nchar(asv_char),
         motif_fwd = grepl(FWD, asv_char, fixed = TRUE),
         ov3 = best$ov3[match(asv_id, best$qseqid)],
         is_bimera = asv_id %in% chim_ids,
         bp_join = asv_id %in% bp_pairs$asv_id,
         artifact_class = case_when(
           bp_join ~ "chimera_confirmed",
           motif_fwd ~ "primer_artifact",
           !is.na(ov3) & ov3 >= 10 & assignment_status == "no_match" ~ "non12S_tail",
           !is.na(ov3) & ov3 >= 10 & assignment_status != "no_match" & db_source == "full_reference" ~ "coverage_gap_rescued",
           is_bimera & !bp_join ~ "bimera_flag_unconfirmed",
           TRUE ~ "ok"))
write_csv(asv_audit, file.path(audit_dir, "asv_audit.csv"))
write_csv(bp_pairs, file.path(audit_dir, "breakpoint_joins.csv"))

# ---- 2. Reference-level audit ----
flag_incoherent <- function(refs, thresh = 0.25) {
  out <- list()
  for (sp in unique(refs$species)) {
    sub <- refs[refs$species == sp, , drop = FALSE]
    if (nrow(sub) < 2) next
    d <- pairdist(sub$sequence, sub$sequence_id)
    minD <- vapply(seq_len(nrow(d)), function(i) min(d[i, -i], na.rm = TRUE), numeric(1))
    if (any(minD > thresh))
      out[[sp]] <- tibble(species = sp, sequence_id = rownames(d)[minD > thresh],
                          min_intra_dist = minD[minD > thresh])
  }
  bind_rows(out)
}
ref_incoherent <- flag_incoherent(all_12S)

span_census <- all_12S |>
  mutate(spans = grepl(FWD, sequence, fixed = TRUE) & grepl(RCREV, sequence, fixed = TRUE)) |>
  group_by(species) |>
  summarise(n_12S_pass = sum(qc_flag == "pass"),
            n_12S_span = sum(qc_flag == "pass" & spans), .groups = "drop")
mcov  <- rb_marker_coverage(con, rank = "species")   # per-species, per-marker counts
occ_tbl   <- rb_get_sequences(con, qc_flag = NULL) |> distinct(species, occurrence)
alien_spp <- rb_get_sequences(con, occurrence = "alien", qc_flag = NULL) |>
  distinct(species) |> pull(species)
tax_cov <- rb_list_taxa(con, rank = "species")
gap12 <- rb_barcode_gap(con) |> filter(marker == "12S") |>
  select(species, taxonomic_resolution, assignment_risk)
qcsum <- rb_qc_summary(con) # DB-wide flag census (methods text)

# Inspect once and align names to what the attribution rules expect:
print(names(mcov)); print(head(mcov))
ref_audit <- mcov |>
  filter(seq_type %in% c("12S", "genome")) |>
  tidyr::pivot_wider(names_from = seq_type, values_from = n_sequences, values_fill = 0) |>
  dplyr::rename(n_12S_pass = `12S`, n_gen_pass = genome) |>
  left_join(span_census |> select(species, n_12S_span), by = "species") |>
  left_join(occ_tbl, by = "species") |>
  left_join(gap12, by = "species") |>
  mutate(n_12S_span = replace_na(n_12S_span, 0L))
write_csv(ref_incoherent, file.path(audit_dir, "reference_incoherent.csv"))
write_csv(ref_audit, file.path(audit_dir, "reference_coverage.csv"))

# ---- 3. Tie audit (diagnostic sites, distances, single-reference mislabel test) ----
tie_audit <- function(taxa, asv_ids) {
  refs <- all_12S |> filter(species %in% taxa)
  sel  <- names(sava_asvs) %in% asv_ids
  ids  <- c(paste0(refs$species, "|", refs$sequence_id), paste0("ASV|", names(sava_asvs)[sel]))
  seqs <- c(refs$sequence, asv_char[sel]); keep <- nchar(seqs) >= 100
  msa  <- DECIPHER::AlignSeqs(DNAStringSet(setNames(seqs[keep], ids[keep])), verbose = FALSE)
  mat  <- as.matrix(msa); is_asv <- startsWith(ids[keep], "ASV|")
  cols <- which(colSums(mat[is_asv, , drop = FALSE] != "-") > 0)
  frag <- Biostrings::subseq(msa, min(cols), max(cols))
  a <- as.matrix(frag); lab <- sub("\\|.*$", "", rownames(a)); foc <- lab %in% taxa
  sites <- which(vapply(seq_len(ncol(a)), function(j) {
    col <- a[foc, j]; l <- lab[foc]
    sa <- unique(col[l == taxa[1]]); sb <- unique(col[l == taxa[2]])
    length(sa) == 1 && length(sb) == 1 && sa != sb && !any(c(sa, sb) == "-")
  }, logical(1)))
  d <- as.matrix(ape::dist.dna(as.DNAbin(frag), model = "raw", pairwise.deletion = TRUE))
  n <- nrow(d); both <- matrix(foc, n, n) & matrix(foc, n, n, byrow = TRUE); low <- row(d) < col(d)
  intra <- both & low & outer(lab, lab, "=="); inter <- both & low & outer(lab, lab, "!=")
  tibble(pair = paste(taxa, collapse = " vs "), asv = paste(asv_ids, collapse = ","),
         frag_width = width(frag), n_refs = sum(foc),
         fixed_diagnostic_sites = length(sites),
         max_intra = if (any(intra)) max(d[intra], na.rm = TRUE) else NA_real_,
         min_inter = if (any(inter)) min(d[inter], na.rm = TRUE) else NA_real_)
}
single_ref_test <- function(taxa) {
  out <- list()
  for (sp in taxa) {
    refs <- all_12S |> filter(species == sp)
    if (nrow(refs) != 1) next
    g <- sub(" .*", "", sp)
    pool <- all_12S |> filter(species == sp | (sub(" .*", "", species) == g & species != sp) |
                                 species %in% setdiff(taxa, sp))
    if (nrow(pool) < 3) next
    d <- pairdist(pool$sequence, paste0(pool$species, "|", pool$sequence_id))
    l <- pool$species; r <- rownames(d)[l == sp][1]
    out[[sp]] <- tibble(species = sp, ref = refs$sequence_id,
                        d_to_congener  = min(d[r, l == g & l != sp], na.rm = TRUE),
                        d_to_other_taxon = min(d[r, l %in% setdiff(taxa, sp)], na.rm = TRUE))
  }
  bind_rows(out)
}
ties <- assign |> filter(assignment_status %in% c("check_full_tie", "ambiguous"))
tie_rows <- list(); mis_rows <- list()
if (nrow(ties)) {
  for (k in seq_len(nrow(ties))) {
    tx <- strsplit(ties$tie_taxa[k], "\\|")[[1]]
    tie_rows[[k]] <- tie_audit(tx, ties$asv_id[k])
    mis_rows[[k]] <- single_ref_test(tx)
  }
}
tie_audit_tbl  <- bind_rows(tie_rows); mislabel_tbl <- bind_rows(mis_rows)
write_csv(tie_audit_tbl, file.path(audit_dir, "tie_audit.csv"))
write_csv(mislabel_tbl, file.path(audit_dir, "single_reference_mislabel_test.csv"))
rb_disconnect(con)
message("Part 1 complete: asv_audit.csv, reference_*.csv, tie_audit.csv written.")

# ============================================================
#   Part 2: Direct Matrix-to-Matrix Comparison & Attribution
# ============================================================
# 1. Load the matrices
# our_mat is already in your environment as sava_spp_site
our_mat <- sava_spp_site 

# Load the article's matrix (ensure it has a 'species' column and site columns)
their_mat <- read_csv("exports/evaluation/article_species_matrix.csv", show_col_types = FALSE)

# 2. Pivot to long format for easy comparison
our_long <- our_mat |>
  pivot_longer(-species, names_to = "site", values_to = "abundance_ours") |>
  mutate(abundance_ours = replace_na(as.numeric(abundance_ours), 0))

their_long <- their_mat |>
  pivot_longer(-species, names_to = "site", values_to = "abundance_theirs") |>
  mutate(abundance_theirs = replace_na(as.numeric(abundance_theirs), 0))

# 3. Join the matrices
# We use a full join to keep species found by either study
comp_long <- full_join(our_long, their_long, by = c("species", "site"))

# 4. Summarize per species (Total reads and Site prevalence)
species_summary <- comp_long |>
  group_by(species) |>
  summarise(
    total_reads_ours = sum(abundance_ours, na.rm = TRUE),
    total_reads_theirs = sum(abundance_theirs, na.rm = TRUE),
    sites_detected_ours = sum(abundance_ours > 0, na.rm = TRUE),
    sites_detected_theirs = sum(abundance_theirs > 0, na.rm = TRUE),
    .groups = "drop"
  )

# 5. Load Reference Audit for Factor B attribution
# (We still need this to explain WHY a species might be missing in our matrix)
ref_audit <- read_csv(file.path(audit_dir, "reference_coverage.csv"), show_col_types = FALSE)

# 6. Define Overrides (Manual B/C factors)
overrides <- tibble(
  species = c("Chondrostoma soetta", 
              "Ballerus ballerus", "Ballerus sapa", 
              "Phoxinus phoxinus", "Phoxinus lumaireul", 
              "Clarias gariepinus"),
  factor  = c("B", "B", "B", "B", "B", "C"),
  note    = c("nasus/soetta complex; biogeography favours nasus",
              "B. ballerus/sapa: species-level swap within genus",
              "see ballerus row",
              "P. phoxinus/lumaireul: species-complex swap",
              "see phoxinus row",
              "genuine Clarias refs; alien/contaminant finding retained")
)

# 7. Attribution Logic (The 3 Factors)
# We apply this to the species_summary table
attribute_factor <- function(row) {
  # Check overrides first
  ov <- overrides |> filter(species == row$species)
  if (nrow(ov) > 0) return(tibble(factor = ov$factor[1], sub = paste0("audit override: ", ov$note[1])))
  
  # Concordant
  if (row$sites_detected_ours > 0 && row$sites_detected_theirs > 0) {
    return(tibble(factor = "-", sub = "concordant"))
  }
  
  # Study-Only (We found it, they didn't)
  if (row$sites_detected_ours > 0) {
    genus <- sub(" .*", "", row$species)
    g_in_paper <- genus %in% sub(" .*", "", unique(their_mat$species))
    
    if (row$species %in% c("Carassius auratus", "Salvelinus fontinalis")) {
      return(tibble(factor = "C", sub = "paper deleted alien/stocked taxon post hoc; study retains"))
    }
    if (row$species %in% c("Clarias gariepinus", "Neogobius fluviatilis", "Neogobius melanostomus", 
                           "Ctenopharyngodon idella", "Pseudorasbora parva", "Oncorhynchus mykiss")) {
      return(tibble(factor = "C", sub = "alien taxon: flag-vs-delete curation philosophy"))
    }
    if (g_in_paper) {
      return(tibble(factor = "B", sub = "congener detected where paper reports sibling species"))
    }
    # Process difference (Factor A)
    if (row$total_reads_ours >= 1000 && row$sites_detected_ours >= 5) {
      return(tibble(factor = "A", sub = "robust study detection absent from paper (ASV calling/assignment process)"))
    }
    return(tibble(factor = "A", sub = "low-abundance study detection (filter/depth sensitivity)"))
  }
  
  # Paper-Only (They found it, we didn't)
  if (row$sites_detected_theirs > 0) {
    genus <- sub(" .*", "", row$species)
    g_in_study <- genus %in% sub(" .*", "", unique(our_mat$species))
    
    if (row$species == "Barbus meridionalis") {
      return(tibble(factor = "C", sub = "paper retained range-implausible best hit; curated DB excludes"))
    }
    
    # Check reference coverage (Factor B)
    ref_info <- ref_audit |> filter(species == row$species)
    if (nrow(ref_info) > 0) {
      if (ref_info$n_12S_pass[1] == 0) {
        return(tibble(factor = "B", sub = "reference coverage gap: no passing 12S in study DB"))
      }
      if (ref_info$n_12S_span[1] == 0) {
        return(tibble(factor = "B", sub = "references present but not Teleo-amplicon compatible"))
      }
    }
    
    if (g_in_study) {
      return(tibble(factor = "B", sub = "congener swap: paper species vs study sibling"))
    }
    
    return(tibble(factor = "A", sub = "process difference (merging/filters/depth) - inferred"))
  }
  
  return(tibble(factor = "NA", sub = "not detected by either"))
}

# Apply the function row by row
attribution <- purrr::map_dfr(split(species_summary, seq_len(nrow(species_summary))), attribute_factor)
final_comparison <- bind_cols(species_summary, attribution)

# 8. Calculate Community Similarity Metrics
# Filter to shared species for correlation
shared_species <- final_comparison |> filter(factor == "-")
if (nrow(shared_species) > 5) {
  cor_test <- cor.test(shared_species$total_reads_ours, shared_species$total_reads_theirs, method = "spearman")
  message(sprintf("Spearman correlation of abundance for %d shared species: rho = %.3f, p = %.2e", 
                  nrow(shared_species), cor_test$estimate, cor_test$p.value))
}

# Jaccard Similarity (Presence/Absence)
species_ours <- unique(our_mat$species)
species_theirs <- unique(their_mat$species)
jaccard <- length(intersect(species_ours, species_theirs)) / length(union(species_ours, species_theirs))
message(sprintf("Species-level Jaccard Similarity: %.3f (%d shared / %d total unique)", 
                jaccard, length(intersect(species_ours, species_theirs)), length(union(species_ours, species_theirs))))

# 9. Save Outputs
write_csv(final_comparison, file.path(audit_dir, "matrix_comparison_attribution.csv"))

# Print summary table for the manuscript
message("\n--- Discordance Summary ---")
print(final_comparison |> 
        filter(factor != "-") |> 
        group_by(factor, sub) |> 
        summarise(species = paste(species, collapse = ", "), .groups = "drop"))

message("\nPart 2 complete. Comparison saved to matrix_comparison_attribution.csv")
