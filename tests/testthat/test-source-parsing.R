# ============================================================
# tests/testthat/test-source-parsing.R
#
# Intention-based tests for source parsing and standardization
# in barcurateR.
#
# These tests are designed to check whether the current functions
# match the intended behavior discussed for:
#
#   - marker pattern defaults
#   - marker/completeness extraction from descriptions
#   - source table standardization
#   - sequence cleaning
#   - primer removal
#   - source combining
#   - ambiguous sequence resolution
#
# Some tests may fail until the implementation is updated.
# ============================================================

# Skip the whole file if required packages are unavailable.
# These should eventually be declared appropriately in DESCRIPTION.
testthat::skip_if_not_installed("stringr")
testthat::skip_if_not_installed("stringi")
testthat::skip_if_not_installed("dplyr")

# ------------------------------------------------------------
# rb_standardize_columns()
# ------------------------------------------------------------

test_that("rb_standardize_columns works directly", {
  
  raw <- data.frame(
    accession = "X1",
    scientific_name = "Danio rerio",
    sequence = "ACGT",
    stringsAsFactors = FALSE
  )
  
  out <- rb_standardize_columns(
    raw,
    c(
      sequence_id = "accession",
      species = "scientific_name",
      sequence = "sequence"
    )
  )
  
  expect_true(all(c("sequence_id", "species", "sequence") %in% names(out)))
  
  # Original raw columns should be removed after renaming.
  expect_false("accession" %in% names(out))
  expect_false("scientific_name" %in% names(out))
  
  expect_equal(out$sequence_id, "X1")
  expect_equal(out$species, "Danio rerio")
  expect_equal(out$sequence, "ACGT")
})

# ------------------------------------------------------------
# rb_default_marker_patterns()
# ------------------------------------------------------------

test_that("rb_default_marker_patterns includes expected vertebrate markers", {
  
  patterns <- rb_default_marker_patterns()
  
  expect_type(patterns, "list")
  
  # Expected standardized marker labels.
  # If you prefer "cytb" instead of "CYTB", update this and the code.
  expect_true(
    all(
      c("12S", "16S", "COI", "CYTB", "18S") %in% names(patterns)
    )
  )
  
  # The patterns should recognize common terms.
  expect_match(patterns[["12S"]], "12s", ignore.case = TRUE)
  expect_match(patterns[["16S"]], "16s", ignore.case = TRUE)
  expect_match(patterns[["COI"]], "coi", ignore.case = TRUE)
  expect_match(patterns[["CYTB"]], "cytochrome b", ignore.case = TRUE)
  expect_match(patterns[["18S"]], "18s", ignore.case = TRUE)
})


# ------------------------------------------------------------
# rb_extract_markers()
# ------------------------------------------------------------

test_that("rb_extract_markers identifies complete genome records first", {
  
  expect_equal(
    rb_extract_markers("mitochondrion, complete genome"),
    "genome_complete"
  )
  
  expect_equal(
    rb_extract_markers("mitogenome"),
    "genome_complete"
  )
  
  # Genome label should take precedence even if markers are mentioned.
  expect_equal(
    rb_extract_markers(
      "mitochondrion, complete genome, contains COI and 12S"
    ),
    "genome_complete"
  )
})


test_that("rb_extract_markers identifies partial marker records", {
  
  expect_equal(
    rb_extract_markers("Danio rerio 12S ribosomal RNA gene, partial sequence"),
    "12S_partial"
  )
  
  expect_equal(
    rb_extract_markers("16S ribosomal RNA gene, partial sequence"),
    "16S_partial"
  )
  
  expect_equal(
    rb_extract_markers("COI gene, partial sequence"),
    "COI_partial"
  )
  
  expect_equal(
    rb_extract_markers("cytochrome b gene, partial cds"),
    "CYTB_partial"
  )
  
  expect_equal(
    rb_extract_markers("cytb gene, partial cds"),
    "CYTB_partial"
  )
  
  expect_equal(
    rb_extract_markers("18S ribosomal RNA gene, partial sequence"),
    "18S_partial"
  )
})


test_that("rb_extract_markers identifies complete marker records", {
  
  # These tests are important for checking the regex grouping behavior.
  # If they fail, the completeness regex likely needs to be fixed.
  
  expect_equal(
    rb_extract_markers("Danio rerio 12S ribosomal RNA gene, complete sequence"),
    "12S_complete"
  )
  
  expect_equal(
    rb_extract_markers("16S ribosomal RNA gene, complete sequence"),
    "16S_complete"
  )
  
  expect_equal(
    rb_extract_markers("COI gene, complete sequence"),
    "COI_complete"
  )
  
  expect_equal(
    rb_extract_markers("cytochrome c oxidase subunit I gene, complete sequence"),
    "COI_complete"
  )
  
  expect_equal(
    rb_extract_markers("cytb gene, complete cds"),
    "CYTB_complete"
  )
  
  expect_equal(
    rb_extract_markers("cytochrome b gene, complete cds"),
    "CYTB_complete"
  )
  
  expect_equal(
    rb_extract_markers("18S ribosomal RNA gene, complete sequence"),
    "18S_complete"
  )
})


test_that("rb_extract_markers returns unknown completeness when not stated", {
  
  expect_equal(
    rb_extract_markers("COI gene"),
    "COI_unknown"
  )
  
  expect_equal(
    rb_extract_markers("12S ribosomal RNA gene"),
    "12S_unknown"
  )
  
  expect_equal(
    rb_extract_markers("16S ribosomal RNA gene"),
    "16S_unknown"
  )
  
  expect_equal(
    rb_extract_markers("cytochrome b gene"),
    "CYTB_unknown"
  )
  
  expect_equal(
    rb_extract_markers("18S ribosomal RNA gene"),
    "18S_unknown"
  )
})


test_that("rb_extract_markers uses explicit wording even for marker fragments", {
  
  # This documents the current intended behavior:
  # if the source description explicitly says complete, it is labelled complete,
  # even for fragments such as COI-5P.
  #
  # If you later decide that COI-5P should always be partial, this test
  # should be changed.
  
  expect_equal(
    rb_extract_markers("COI-5P gene, complete sequence"),
    "COI_complete"
  )
  
  expect_equal(
    rb_extract_markers("COI-5P gene, partial sequence"),
    "COI_partial"
  )
})


test_that("rb_extract_markers handles multi-gene records", {
  
  expect_equal(
    rb_extract_markers(
      "Species tRNA-Phe gene, 12S ribosomal RNA gene, and tRNA-Val gene, complete sequence"
    ),
    "12S_complete;other_complete"
  )
  
  expect_equal(
    rb_extract_markers(
      "genes for 12S rRNA and tRNA-Phe, partial sequence"
    ),
    "12S_partial;other_partial"
  )
})


test_that("rb_extract_markers supports custom marker patterns", {
  
  expect_equal(
    rb_extract_markers(
      "Some plant ITS region, partial sequence",
      marker_patterns = list(ITS = "its")
    ),
    "ITS_partial"
  )
  
  expect_equal(
    rb_extract_markers(
      "Some plant ITS region, complete sequence",
      marker_patterns = list(ITS = "its")
    ),
    "ITS_complete"
  )
})


test_that("rb_extract_markers handles custom patterns containing alternation", {
  
  # This is a specific regression test for regex grouping.
  # If marker_pattern contains "|", completeness detection should still work.
  
  custom_patterns <- list(TEST = "test1|test2")
  
  expect_equal(
    rb_extract_markers(
      "test2 gene, partial sequence",
      marker_patterns = custom_patterns
    ),
    "TEST_partial"
  )
  
  expect_equal(
    rb_extract_markers(
      "test2 gene, complete sequence",
      marker_patterns = custom_patterns
    ),
    "TEST_complete"
  )
})


test_that("rb_extract_markers returns other_unknown when nothing is recognized", {
  
  expect_equal(
    rb_extract_markers("some unknown sequence"),
    "other_unknown"
  )
  
  expect_equal(
    rb_extract_markers("environmental sample DNA"),
    "other_unknown"
  )
})


# ------------------------------------------------------------
# rb_clean_sequence()
# ------------------------------------------------------------

test_that("rb_clean_sequence uppercases sequences", {
  
  expect_equal(rb_clean_sequence("acgt"), "ACGT")
  expect_equal(rb_clean_sequence("ACGT"), "ACGT")
  expect_equal(rb_clean_sequence("aCgT"), "ACGT")
})


test_that("rb_clean_sequence removes whitespace but retains alignment gaps", {
  
  # Whitespace should be removed.
  expect_equal(rb_clean_sequence(" acgt "), "ACGT")
  expect_equal(rb_clean_sequence("AC GT"), "ACGT")
  expect_equal(rb_clean_sequence("AC\nGT"), "ACGT")
  expect_equal(rb_clean_sequence("AC\tGT"), "ACGT")
  expect_equal(rb_clean_sequence(" AC GT "), "ACGT")
  
  # Real alignment gaps should be retained.
  expect_equal(rb_clean_sequence("AC-GT"), "AC-GT")
  expect_equal(rb_clean_sequence(" AC-GT "), "AC-GT")
  expect_equal(rb_clean_sequence("AC - GT"), "AC-GT")
  expect_equal(rb_clean_sequence("AC--GT"), "AC--GT")
})


# ------------------------------------------------------------
# rb_parse_source_table()
# ------------------------------------------------------------

test_that("rb_parse_source_table standardizes a simple NCBI-like table", {
  
  raw <- data.frame(
    accession = c("X1", "X2", "X3"),
    scientific_name = c("Danio rerio", "Cyprinus carpio", "Carassius auratus"),
    sequence = c("acgtacgt", "ACGTACGTAA", "AC-GT"),
    definition = c(
      "12S ribosomal RNA gene, partial sequence",
      "mitochondrion, complete genome",
      "cytochrome b gene, partial cds"
    ),
    stringsAsFactors = FALSE
  )
  
  parsed <- rb_parse_source_table(
    raw,
    source_name = "ncbi_test",
    column_map = c(
      sequence_id = "accession",
      species = "scientific_name",
      sequence = "sequence"
    ),
    description_col = "definition"
  )
  
  expect_s3_class(parsed, "data.frame")
  
  # Standardized columns should be present.
  expect_true(
    all(
      c(
        "sequence_id",
        "species",
        "sequence",
        "source",
        "length",
        "seq_type"
      ) %in% names(parsed)
    )
  )
  
  # Source identifiers should become sequence_id.
  expect_equal(parsed$sequence_id, c("X1", "X2", "X3"))
  
  # Source label should be added.
  expect_equal(parsed$source, rep("ncbi_test", 3))
  
  # Sequences should be cleaned but gaps retained.
  expect_equal(
    parsed$sequence,
    c("ACGTACGT", "ACGTACGTAA", "AC-GT")
  )
  
  # Length should count retained gap characters.
  expect_equal(parsed$length, c(8, 10, 5))
  
  # seq_type should be inferred from the description column.
  expect_equal(
    parsed$seq_type,
    c("12S_partial", "genome_complete", "CYTB_partial")
  )
})


test_that("rb_parse_source_table works when seq_type is already supplied", {
  
  raw <- data.frame(
    processid = c("BOLD1", "BOLD2"),
    species_name = c("Danio rerio", "Cyprinus carpio"),
    nucleotides = c("ACGTACGT", "ACGTACGTAA"),
    marker = c("COI_partial", "12S_partial"),
    stringsAsFactors = FALSE
  )
  
  parsed <- rb_parse_source_table(
    raw,
    source_name = "bold_test",
    column_map = c(
      sequence_id = "processid",
      species = "species_name",
      sequence = "nucleotides",
      seq_type = "marker"
    ),
    description_col = NULL
  )
  
  expect_true("seq_type" %in% names(parsed))
  
  expect_equal(
    parsed$seq_type,
    c("COI_partial", "12S_partial")
  )
})


test_that("rb_parse_source_table assigns other_unknown when no marker info exists", {
  
  raw <- data.frame(
    accession = "X1",
    scientific_name = "Danio rerio",
    sequence = "ACGTACGT",
    stringsAsFactors = FALSE
  )
  
  parsed <- rb_parse_source_table(
    raw,
    source_name = "local_test",
    column_map = c(
      sequence_id = "accession",
      species = "scientific_name",
      sequence = "sequence"
    ),
    description_col = NULL
  )
  
  expect_equal(parsed$seq_type, "other_unknown")
})


test_that("rb_parse_source_table removes whitespace but retains gaps", {
  
  raw <- data.frame(
    accession = c("X1", "X2"),
    scientific_name = c("Danio rerio", "Cyprinus carpio"),
    sequence = c(" acgt acgt ", "AC-GT"),
    definition = c(
      "12S ribosomal RNA gene, partial sequence",
      "COI gene, partial sequence"
    ),
    stringsAsFactors = FALSE
  )
  
  parsed <- rb_parse_source_table(
    raw,
    source_name = "ncbi_test",
    column_map = c(
      sequence_id = "accession",
      species = "scientific_name",
      sequence = "sequence"
    ),
    description_col = "definition"
  )
  
  expect_equal(
    parsed$sequence,
    c("ACGTACGT", "AC-GT")
  )
  
  expect_equal(
    parsed$length,
    c(8, 5)
  )
})


# ------------------------------------------------------------
# rb_reverse_complement()
# ------------------------------------------------------------

test_that("rb_reverse_complement returns correct reverse complements", {
  
  expect_equal(rb_reverse_complement("ATCG"), "CGAT")
  expect_equal(rb_reverse_complement("AAAC"), "GTTT")
  expect_equal(rb_reverse_complement("ACGT"), "ACGT")
})


# ------------------------------------------------------------
# rb_remove_primers()
# ------------------------------------------------------------

test_that("rb_remove_primers removes forward primers", {
  
  expect_equal(
    rb_remove_primers("AAACGTACGT", f_primer_seqs = "AAA"),
    "CGTACGT"
  )
})


test_that("rb_remove_primers removes reverse primers", {
  
  # Reverse primer TTT has reverse complement AAA.
  expect_equal(
    rb_remove_primers("ACGTACGTAAA", r_primer_seqs = "TTT"),
    "ACGTACGT"
  )
})


test_that("rb_remove_primers removes both forward and reverse primers", {
  
  # Case 1:
  # Forward primer: AAA
  # Reverse primer: TTT
  # Reverse complement of TTT is AAA.
  #
  # Full amplicon:
  # AAA + CGTACGT + AAA
  expect_equal(
    rb_remove_primers(
      "AAACGTACGTAAA",
      f_primer_seqs = "AAA",
      r_primer_seqs = "TTT"
    ),
    "CGTACGT"
  )
  
  # Case 2:
  # Forward primer: AAA
  # Reverse primer: AAA
  # Reverse complement of AAA is TTT.
  #
  # Full amplicon:
  # AAA + CGTACGT + TTT
  expect_equal(
    rb_remove_primers(
      "AAACGTACGTTTT",
      f_primer_seqs = "AAA",
      r_primer_seqs = "AAA"
    ),
    "CGTACGT"
  )
})


test_that("rb_remove_primers returns sequence unchanged if no primer matches", {
  
  expect_equal(
    rb_remove_primers("ACGTACGT", f_primer_seqs = "AAA"),
    "ACGTACGT"
  )
  
  expect_equal(
    rb_remove_primers("ACGTACGT", f_primer_seqs = NA, r_primer_seqs = NA),
    "ACGTACGT"
  )
})


test_that("rb_remove_primers checks multiple possible forward primers", {
  
  expect_equal(
    rb_remove_primers("CCCGTACGT", f_primer_seqs = c("AAA", "CCC")),
    "GTACGT"
  )
})


test_that("rb_remove_primers retains alignment gaps after trimming", {
  
  expect_equal(
    rb_remove_primers("AAAAC-GT", f_primer_seqs = "AAA"),
    "AC-GT"
  )
  
  expect_equal(
    rb_remove_primers("AC-GTAAA", r_primer_seqs = "TTT"),
    "AC-GT"
  )
})


# ------------------------------------------------------------
# rb_combine_sources()
# ------------------------------------------------------------

test_that("rb_combine_sources combines and deduplicates records", {
  
  a <- data.frame(
    species = "sp1",
    sequence = "AC-GT",
    stringsAsFactors = FALSE
  )
  
  b <- data.frame(
    species = "sp1",
    sequence = "ac-gt",
    stringsAsFactors = FALSE
  )
  
  c <- data.frame(
    species = "sp2",
    sequence = "AC-GT",
    stringsAsFactors = FALSE
  )
  
  combined <- rb_combine_sources(list(a, b, c))
  
  expect_s3_class(combined, "data.frame")
  
  # sp1 / AC-GT duplicate should be collapsed.
  # sp2 / AC-GT should remain because the species differs.
  expect_equal(nrow(combined), 2)
  expect_true(all(c("sp1", "sp2") %in% combined$species))
  
  sp1_seq <- combined$sequence[combined$species == "sp1"]
  expect_equal(sp1_seq, "AC-GT")
})


test_that("rb_combine_sources cleans sequences before deduplication", {
  
  a <- data.frame(
    species = "sp1",
    sequence = " ACGT ",
    stringsAsFactors = FALSE
  )
  
  b <- data.frame(
    species = "sp1",
    sequence = "acgt",
    stringsAsFactors = FALSE
  )
  
  combined <- rb_combine_sources(list(a, b))
  
  expect_equal(nrow(combined), 1)
  expect_equal(combined$sequence, "ACGT")
})


# ------------------------------------------------------------
# rb_resolve_ambiguous()
# ------------------------------------------------------------

test_that("rb_resolve_ambiguous returns data unchanged when no ambiguity exists", {
  
  d <- data.frame(
    species = c("A", "B"),
    sequence = c("AAA", "CCC"),
    stringsAsFactors = FALSE
  )
  
  out <- rb_resolve_ambiguous(d)
  
  expect_equal(nrow(out), 2)
  expect_true(all(c("A", "B") %in% out$species))
})


test_that("rb_resolve_ambiguous warns when ambiguity exists and on_unresolved = 'warn'", {
  
  dup <- data.frame(
    species = c("A", "B"),
    sequence = c("ACGT", "ACGT"),
    stringsAsFactors = FALSE
  )
  
  expect_warning(
    out <- rb_resolve_ambiguous(dup, on_unresolved = "warn")
  )
  
  # With warn, data should be returned unchanged.
  expect_equal(nrow(out), 2)
})


test_that("rb_resolve_ambiguous errors when ambiguity exists and on_unresolved = 'stop'", {
  
  dup <- data.frame(
    species = c("A", "B"),
    sequence = c("ACGT", "ACGT"),
    stringsAsFactors = FALSE
  )
  
  expect_error(
    rb_resolve_ambiguous(dup, on_unresolved = "stop")
  )
})


test_that("rb_resolve_ambiguous drops ambiguous sequences when on_unresolved = 'drop'", {
  
  dup <- data.frame(
    species = c("A", "B"),
    sequence = c("ACGT", "ACGT"),
    stringsAsFactors = FALSE
  )
  
  expect_warning(
    out <- rb_resolve_ambiguous(dup, on_unresolved = "drop")
  )
  
  expect_equal(nrow(out), 0)
})


test_that("rb_resolve_ambiguous can write ambiguous records to a CSV file", {
  
  dup <- data.frame(
    species = c("A", "B"),
    sequence = c("ACGT", "ACGT"),
    stringsAsFactors = FALSE
  )
  
  tmp <- tempfile(fileext = ".csv")
  
  expect_warning(
    rb_resolve_ambiguous(
      dup,
      on_unresolved = "warn",
      flag_output_path = tmp
    )
  )
  
  expect_true(file.exists(tmp))
  
  unlink(tmp)
})


test_that("rb_resolve_ambiguous applies keep actions from a resolution table", {
  
  d <- data.frame(
    species = c("A", "B", "C"),
    sequence = c("ACGT", "ACGT", "TTTT"),
    source = "test",
    stringsAsFactors = FALSE
  )
  
  resolution <- data.frame(
    species = "A",
    sequence = "ACGT",
    source = "test",
    action = "keep",
    stringsAsFactors = FALSE
  )
  
  out <- rb_resolve_ambiguous(d, resolution_table = resolution)
  
  expect_s3_class(out, "data.frame")
  
  # Original ambiguous A/B rows are removed, then A is kept back.
  # C is not ambiguous and should remain.
  expect_equal(nrow(out), 2)
  expect_true("A" %in% out$species)
  expect_true("C" %in% out$species)
  expect_false("B" %in% out$species)
})


test_that("rb_resolve_ambiguous applies combine actions from a resolution table", {
  
  d <- data.frame(
    species = c("A", "B", "C"),
    sequence = c("ACGT", "ACGT", "TTTT"),
    source = "test",
    stringsAsFactors = FALSE
  )
  
  resolution <- data.frame(
    species = c("A", "B"),
    sequence = c("ACGT", "ACGT"),
    source = c("test", "test"),
    action = "combine",
    stringsAsFactors = FALSE
  )
  
  out <- rb_resolve_ambiguous(d, resolution_table = resolution)
  
  expect_s3_class(out, "data.frame")
  
  # ACGT should be combined into one row with species A|B.
  # TTTT should remain unchanged.
  expect_equal(nrow(out), 2)
  expect_true("A|B" %in% out$species)
  expect_true("C" %in% out$species)
})


test_that("rb_resolve_ambiguous treats blank, NA, and explicit drop actions as dropped", {
  
  d <- data.frame(
    species = c("A", "B", "C", "D"),
    sequence = c("ACGT", "ACGT", "TTTT", "GGGG"),
    stringsAsFactors = FALSE
  )
  
  resolution <- data.frame(
    species = c("A", "B", "C", "D"),
    sequence = c("ACGT", "ACGT", "TTTT", "GGGG"),
    action = c("keep", "", NA, "drop"),
    stringsAsFactors = FALSE
  )
  
  # Intended behavior: the function should warn that rows without
  # a usable action were dropped.
  expect_warning(
    out <- rb_resolve_ambiguous(d, resolution_table = resolution),
    "drop|missing|blank|action"
  )
  
  # A is kept.
  expect_true("A" %in% out$species)
  
  # B has blank action and should be dropped.
  expect_false("B" %in% out$species)
  
  # C has NA action and should be dropped.
  expect_false("C" %in% out$species)
  
  # D has explicit drop action and should be dropped.
  expect_false("D" %in% out$species)
})


test_that("rb_resolve_ambiguous reports the number of dropped blank/drop actions", {
  
  d <- data.frame(
    species = c("A", "B", "C"),
    sequence = c("ACGT", "ACGT", "TTTT"),
    stringsAsFactors = FALSE
  )
  
  resolution <- data.frame(
    species = c("A", "B", "C"),
    sequence = c("ACGT", "ACGT", "TTTT"),
    action = c("keep", NA, ""),
    stringsAsFactors = FALSE
  )
  
  # Intended behavior: warning/message should include a count.
  # The exact wording can be adjusted later; this checks that some
  # numeric count is emitted.
  expect_warning(
    out <- rb_resolve_ambiguous(d, resolution_table = resolution),
    "2|two"
  )
  
  # Only A should remain.
  expect_equal(nrow(out), 1)
  expect_true("A" %in% out$species)
})