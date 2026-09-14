#' barcurateR: Standardize, Curate, Audit, and Deploy DNA Reference Databases
#'
#' @description
#' Provides tools to standardize, curate, query, audit, and export DNA
#' reference databases, and assign amplicon sequence variants for
#' ecological metabarcoding in a reproducible manner.
#'
#' @section Main workflow:
#' 1. Parse and standardize source sequences ([rb_parse_source_table()])
#' 2. Combine sources and resolve ambiguities ([rb_combine_sources()], [rb_resolve_ambiguous()])
#' 3. Run the full curation pipeline ([rb_curate_reference()])
#' 4. Query and export the curated database ([rb_get_sequences()], [rb_export_fasta()])
#' 5. Assign eDNA sequences ([rb_assign_edna()])
#'
#' @docType package
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @keywords internal
#' @name barcurateR-package
#' @aliases barcurateR
"_PACKAGE"