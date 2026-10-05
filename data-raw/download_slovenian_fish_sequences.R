# ============================================================
# data-raw/download_slovenian_fish_sequences.R
#
# Comprehensive, resumable, reproducible NCBI download for the
# Slovenian fish species list. Designed for the barcurateR
# curation workflow.
#
# Key design changes from previous versions:
#  - Query is SHORT and uses NCBI's native mitochondrial
#    filters (mitochondrion[Filter] and gene_in_mitochondrion[PROP])
#    instead of enumerating every possible title phrase.
#    This avoids E-utilities query truncation.
#  - Marker classification is performed POST-FETCH using
#    barcurateR's own marker patterns, not a fragile header regex.
#  - All mitochondrial records are retained; classification only
#    labels them for downstream use.
# ============================================================

# ------------------------------------------------------------
# 1. Dependencies
# ------------------------------------------------------------
required_pkgs <- c("rentrez", "readr", "dplyr", "tibble", "purrr", "tidyr", "stringr")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
})

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# ------------------------------------------------------------
# 2. Configuration
# ------------------------------------------------------------
species_csv    <- "data-raw/slovenian_fish_species.csv"
out_root       <- "data-raw/slovenian_fish"
cache_dir      <- file.path(out_root, "cache")
manifest_dir   <- file.path(out_root, "manifest")
raw_dir        <- file.path(out_root, "raw")
processed_dir  <- file.path(out_root, "processed")
for (d in c(cache_dir, manifest_dir, raw_dir, processed_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

overwrite     <- TRUE      # Set TRUE for a fresh comprehensive run
batch_size    <- 200
request_delay <- 0.4       # Anonymous limit is 3 req/s

# ------------------------------------------------------------
# 3. Read species table
# ------------------------------------------------------------
if (!file.exists(species_csv)) stop("Species table not found: ", species_csv)

species_table <- readr::read_csv(species_csv, show_col_types = FALSE)
if (!"species" %in% names(species_table)) stop("Need a 'species' column.")

species_table <- species_table %>% filter(!grepl(" sp\\.$", species))

alien_spp <- c(
  "Salvelinus fontinalis", "Salvelinus umbla", "Oncorhynchus mykiss",
  "Coregonus lavaretus", "Micropterus salmoides", "Lepomis gibbosus",
  "Ameiurus nebulosus", "Gambusia holbrooki", "Hypophthalmichthys molitrix",
  "Hypophthalmichthys nobilis", "Ctenopharyngodon idella",
  "Pseudorasbora parva", "Carassius gibelio", "Carassius auratus",
  "Clarias gariepinus", "Ameiurus melas", "Oreochromis niloticus",
  "Mylopharyngodon piceus", "Polyodon spathula", "Acipenser baerii",
  "Ponticola kessleri", "Oncorhynchus kisutch", "Babka gymnotrachelus",
  "Neogobius fluviatilis", "Neogobius melanostomus"
)

species_table <- species_table %>%
  mutate(
    genus = if ("genus" %in% names(.)) {
      ifelse(is.na(genus) | !nzchar(genus), sub(" .*", "", species), genus)
    } else {
      sub(" .*", "", species)
    },
    kingdom    = if ("kingdom" %in% names(.)) kingdom else "Animalia",
    phylum     = if ("phylum" %in% names(.)) phylum else "Chordata",
    occurrence = ifelse(species %in% alien_spp, "alien", "native")
  )

taxonomy_table <- species_table %>%
  select(dplyr::any_of(c(
    "species", "genus", "family", "order", "class",
    "phylum", "kingdom", "occurrence"
  )))

readr::write_csv(taxonomy_table, file.path(processed_dir, "slovenian_fish_taxonomy.csv"))

all_spp <- unique(species_table$species)
message("Loaded ", length(all_spp), " valid species.")

# ------------------------------------------------------------
# 4. Helpers
# ------------------------------------------------------------
normalize_sums <- function(sums, uids = NULL) {
  if (is.null(sums) || !is.list(sums) || length(sums) == 0) return(NULL)
  if (!is.list(sums[[1]])) sums <- list(sums)
  sums <- sums[vapply(sums, is.list, logical(1))]
  if (length(sums) == 0) return(NULL)
  if (is.null(names(sums)) || any(!nzchar(names(sums)))) {
    if (!is.null(uids) && length(uids) == length(sums)) names(sums) <- as.character(uids)
  }
  sums
}

dedup_rows <- function(d) {
  if (is.null(d) || nrow(d) == 0) return(d)
  d$.key <- ifelse(
    !is.na(d$accession), d$accession,
    ifelse(!is.na(d$uid), paste0("uid:", d$uid), paste0("seq:", d$sequence))
  )
  d <- d %>% arrange(desc(!is.na(uid))) %>% distinct(.key, .keep_all = TRUE)
  d$.key <- NULL
  d
}

with_retry <- function(expr, max_retries = 5, delay = 3) {
  msg <- NULL
  for (i in seq_len(max_retries)) {
    res <- tryCatch(expr(), error = function(e) e)
    if (!inherits(res, "error")) return(res)
    msg <- conditionMessage(res)
    if (grepl("^NON-TRANSIENT:", msg)) stop(msg, call. = FALSE)
    message("    Attempt ", i, " failed: ", msg)
    Sys.sleep(delay)
  }
  stop("Max retries reached: ", msg, call. = FALSE)
}

is_transient <- function(msg) {
  grepl("timeout|timed out|HTTP status 5|50[0-9]|temporarily unavailable|rate limit|too many requests|connection reset|could not resolve",
        msg, ignore.case = TRUE)
}

safe_search <- function(term, retmax) {
  res <- tryCatch(
    rentrez::entrez_search(db = "nuccore", term = term, retmax = retmax),
    error = function(e) e
  )
  if (!inherits(res, "error")) return(res)
  if (is_transient(conditionMessage(res))) stop(conditionMessage(res), call. = FALSE)
  stop("NON-TRANSIENT: ", conditionMessage(res), call. = FALSE)
}

# ------------------------------------------------------------
# SIMPLIFIED QUERY: short, robust, avoids truncation
# ------------------------------------------------------------
build_query <- function(sp) {
  paste0(
    '"', sp, '"[Organism] AND (',
    'mitochondrion[Filter] OR gene_in_mitochondrion[PROP]',
    ')'
  )
}

parse_fasta_text <- function(fasta_text) {
  if (!nzchar(fasta_text)) return(tibble())
  records <- strsplit(fasta_text, "\n>", fixed = TRUE)[[1]]
  records[1] <- sub("^>", "", records[1])
  purrr::map(records, function(rec) {
    lines  <- strsplit(rec, "\n", fixed = TRUE)[[1]]
    header <- lines[1]
    seq    <- paste(lines[-1], collapse = "")
    tibble(accession = sub("\\s.*$", "", header), header = header, sequence = seq)
  }) |> purrr::list_rbind()
}

# ------------------------------------------------------------
# POST-FETCH MARKER CLASSIFICATION using barcurateR patterns
# ------------------------------------------------------------
classify_marker <- function(header) {
  # Pattern set — any of these appearing in the header tags the record.
  # Order in case_when() determines priority (genome > multi > single).
  p_12s    <- "12s|mt-rnr1|12s ribosomal rna|rrns|small subunit ribosomal"
  p_16s    <- "16s|mt-rnr2|16s ribosomal rna|rrnl|large subunit ribosomal"
  p_coi    <- "\\bcoi\\b|\\bco1\\b|\\bcox1\\b|\\bcoxi\\b|cytochrome c oxidase subunit (1|i)|mt-co1"
  p_cytb   <- "\\bcytb\\b|cytochrome b|mt-cyb"
  p_nd     <- "\\bnd[1-6]\\b|nadh dehydrogenase"
  p_atp    <- "\\batp[68]\\b|atp synthase"
  p_genome <- "complete genome|complete mitochondrial|complete mtdna|mitogenome"
  
  # Vectorised hit counts (each returns an integer vector)
  n_12s  <- as.integer(grepl(p_12s,  header, ignore.case = TRUE))
  n_16s  <- as.integer(grepl(p_16s,  header, ignore.case = TRUE))
  n_coi  <- as.integer(grepl(p_coi,  header, ignore.case = TRUE))
  n_cytb <- as.integer(grepl(p_cytb, header, ignore.case = TRUE))
  n_nd   <- as.integer(grepl(p_nd,   header, ignore.case = TRUE))
  n_atp  <- as.integer(grepl(p_atp,  header, ignore.case = TRUE))
  is_genome <- grepl(p_genome, header, ignore.case = TRUE)
  
  total <- n_12s + n_16s + n_coi + n_cytb + n_nd + n_atp
  
  dplyr::case_when(
    is_genome          ~ "complete_genome",
    total >= 2         ~ "multi_marker",
    n_12s  == 1        ~ "12S",
    n_16s  == 1        ~ "16S",
    n_coi  == 1        ~ "COI",
    n_cytb == 1        ~ "CYTB",
    n_nd   == 1        ~ "ND",
    n_atp  == 1        ~ "ATP",
    TRUE               ~ "other_mito"
  )
}

parse_fasta_batch <- function(fasta_text, sums = NULL) {
  pf <- parse_fasta_text(fasta_text)
  if (nrow(pf) == 0) return(pf)
  
  if (!is.null(sums) && length(sums) > 0) {
    meta <- data.frame(
      uid        = names(sums) %||% rep(NA_character_, length(sums)),
      accession  = vapply(sums, function(s) s$caption %||% NA_character_, character(1)),
      ncbi_taxid = vapply(sums, function(s) as.integer(s$taxid %||% NA_integer_), integer(1)),
      stringsAsFactors = FALSE
    )
    pf <- dplyr::left_join(pf, meta, by = "accession")
  } else {
    pf$uid <- NA_character_
    pf$ncbi_taxid <- NA_integer_
  }
  pf$marker <- classify_marker(pf$header)
  pf
}

# ------------------------------------------------------------
# Species resolution (unchanged)
# ------------------------------------------------------------
resolve_species <- function(taxids) {
  map_file <- file.path(cache_dir, "taxid_species_map.rds")
  map <- if (file.exists(map_file)) readRDS(map_file) else list()
  new_ids <- setdiff(unique(as.character(taxids[!is.na(taxids)])), names(map))
  
  if (length(new_ids) > 0) {
    for (ch in split(new_ids, ceiling(seq_along(new_ids) / 100))) {
      sums <- with_retry(\() rentrez::entrez_summary(db = "taxonomy", ids = ch))
      sums  <- normalize_sums(sums, ch)
      Sys.sleep(request_delay)
      for (s in sums) {
        sci   <- s$scientificname %||% NA_character_
        lin   <- s$lineage %||% ""
        parts <- trimws(strsplit(lin, ";")[[1]])
        cand  <- tail(intersect(parts, all_spp), 1)
        if (!length(cand)) cand <- if (!is.na(sci) && sci %in% all_spp) sci else NA_character_
        if (!is.na(cand) && grepl("\\bx\\b|hybrid", sci, ignore.case = TRUE)) cand <- NA_character_
        map[[as.character(s$taxid %||% NA)]] <- cand
      }
    }
    saveRDS(map, map_file)
  }
  vapply(as.character(taxids),
         function(t) if (!is.null(map[[t]])) map[[t]] else NA_character_,
         character(1), USE.NAMES = FALSE)
}

# ------------------------------------------------------------
# 5. UID snapshot and resumable fetch
# ------------------------------------------------------------
capture_uids <- function(sp) {
  mf_path <- file.path(manifest_dir, paste0(gsub(" ", "_", sp), "_manifest.rds"))
  mf_old  <- if (!overwrite && file.exists(mf_path)) readRDS(mf_path) else NULL
  term    <- build_query(sp)
  
  if (!is.null(mf_old) && identical(mf_old$query, term)) return(mf_old)
  
  message("  Searching: ", sp)
  srch <- with_retry(\() safe_search(term, 0))
  count <- as.integer(srch$count %||% 0L)
  uids  <- integer(0)
  
  if (count > 0) {
    Sys.sleep(request_delay)
    # Use rentrez's automatic POST handling for large ID sets
    uids <- sort(as.integer(
      with_retry(\() safe_search(term, min(count, 100000L)))$ids))
  }
  
  mf <- list(
    species = sp, query = term,
    date = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    n_expected = count, uids = uids,
    failed_uids = character(0),
    n_parsed = mf_old$n_parsed, n_kept = mf_old$n_kept,
    date_retrieved = mf_old$date_retrieved
  )
  saveRDS(mf, mf_path)
  mf
}

placeholder_row <- function(sp) {
  tibble(
    species_query = sp, uid = NA_character_, accession = NA_character_,
    ncbi_taxid = NA_integer_, sequence = NA_character_,
    header = NA_character_, marker = NA_character_,
    species_resolved = NA_character_, species_check = NA_character_
  )
}

fetch_species <- function(sp) {
  cache_file <- file.path(cache_dir, paste0("ncbi_", gsub(" ", "_", sp), ".rds"))
  cached <- if (!overwrite && file.exists(cache_file)) readRDS(cache_file) else NULL
  
  mf <- capture_uids(sp)
  if (length(mf$uids) == 0) return(placeholder_row(sp))
  
  done <- if (!is.null(cached)) as.character(cached$uid) else character(0)
  todo <- setdiff(as.character(mf$uids), done)
  
  parsed <- list()
  failed <- character(0)
  
  if (length(todo) > 0) {
    for (ch in split(todo, ceiling(seq_along(todo) / batch_size))) {
      res <- tryCatch({
        fasta <- with_retry(\()
                            rentrez::entrez_fetch(db = "nuccore", id = ch, rettype = "fasta"))
        Sys.sleep(request_delay)
        sums <- tryCatch(
          with_retry(\() rentrez::entrez_summary(db = "nuccore", id = ch)),
          error = function(e) NULL
        )
        sums <- normalize_sums(sums, ch)
        if (!is.null(sums)) Sys.sleep(request_delay)
        parse_fasta_batch(fasta, sums)
      }, error = function(e) {
        message("    batch failed (", length(ch), " uids): ", conditionMessage(e))
        NULL
      })
      if (is.null(res)) failed <- c(failed, ch) else parsed[[length(parsed) + 1]] <- res
    }
  }
  
  if (length(parsed) > 0) {
    new_rows <- purrr::list_rbind(parsed) %>% mutate(species_query = sp)
    out <- dplyr::bind_rows(cached, new_rows)
  } else {
    out <- cached
  }
  
  if (is.null(out) || nrow(out) == 0) return(placeholder_row(sp))
  
  out <- dedup_rows(out)
  
  # Resolve species
  out$species_resolved <- resolve_species(out$ncbi_taxid)
  out$species_check <- dplyr::case_when(
    is.na(out$sequence)                          ~ "placeholder",
    is.na(out$species_resolved)                  ~ "off_list_or_hybrid",
    out$species_resolved == out$species_query    ~ "consistent",
    TRUE                                         ~ "subspecies_or_synonym"
  )
  
  n_parsed <- nrow(out)
  
  # KEEP ALL mitochondrial records — classification is for labelling only
  # No marker-based filtering here; barcurateR handles that downstream
  saveRDS(out, cache_file)
  
  mf$failed_uids    <- failed
  mf$n_parsed       <- n_parsed
  mf$n_kept         <- nrow(out)
  mf$date_retrieved <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  saveRDS(mf, file.path(manifest_dir, paste0(gsub(" ", "_", sp), "_manifest.rds")))
  
  if (length(failed) > 0) {
    message("    WARNING: ", length(failed), " uids failed for ", sp,
            " (recorded in manifest; rerun to top up).")
  }
  out
}

# ------------------------------------------------------------
# 6. Run
# ------------------------------------------------------------
failed_species <- character(0)

species_results <- purrr::map(all_spp, function(sp) {
  tryCatch(
    fetch_species(sp),
    error = function(e) {
      message("  SPECIES FAILED: ", sp, " — ", conditionMessage(e))
      failed_species <<- c(failed_species, sp)
      placeholder_row(sp)
    }
  )
})

if (length(failed_species) > 0) {
  readr::write_lines(failed_species, file.path(processed_dir, "failed_species.txt"))
  message("Species that failed entirely: ", paste(failed_species, collapse = ", "))
}

all_ncbi <- purrr::list_rbind(species_results) |>
  dplyr::rename(species = species_query) |>
  dedup_rows()

# ------------------------------------------------------------
# 7. Summary and outputs — shaped for barcurateR
# ------------------------------------------------------------
kept <- all_ncbi %>% filter(!is.na(sequence))

message("\n=== Download summary ===")
message("Total records fetched: ", nrow(all_ncbi))
message("Records with sequences: ", nrow(kept))
message("Species represented: ", length(unique(kept$species)))

# Marker distribution
marker_dist <- kept %>% count(marker, sort = TRUE)
print(marker_dist)

# Add a source column for barcurateR's rb_parse_source_table()
all_ncbi <- all_ncbi %>% mutate(source = "ncbi")

# barcurateR expects at minimum: species, sequence, seq_type, source
# We provide marker as a description column for rb_parse_source_table()
all_ncbi <- all_ncbi %>%
  mutate(
    seq_type = marker,  # will be refined by rb_parse_source_table if desired
    length = nchar(sequence)
  )

today <- format(Sys.Date(), "%Y%m%d")
stamped <- file.path(raw_dir, paste0("ncbi_sequence_raw_", today, ".csv"))
working <- file.path(raw_dir, "ncbi_sequence_raw.csv")

readr::write_csv(all_ncbi, stamped)
file.copy(stamped, working, overwrite = TRUE)

message("\nSaved snapshot: ", stamped)
message("Saved working copy: ", working)
message("This CSV is ready for the barcurateR workflow script.")
message("Finished script successfully!")
