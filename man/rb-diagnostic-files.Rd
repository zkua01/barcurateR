\name{rb_read_barcode_gap}
\alias{rb_read_barcode_gap}
\alias{rb_read_ambiguity}
\title{Read barcode gap and ambiguity diagnostics from CSV files}
\description{
Load barcode gap metrics or ambiguous sequence tables from standalone CSV files.
These functions are the file-based counterparts to \code{\link{rb_barcode_gap}} and
\code{\link{rb_ambiguity}}, which query a connected SQLite database.

Use these functions when working with intermediate pipeline outputs that have
not yet been loaded into a SQLite database (e.g., the CSV output of
\code{\link{rb_run_barcode_gap}} or a manually curated ambiguity resolution table).
}
\usage{
rb_read_barcode_gap(path = NULL, species = NULL, marker = NULL)
rb_read_ambiguity(path = NULL, action = NULL)
}
\arguments{
\item{path}{Path to the CSV file. If \code{NULL}, the function searches for the
default file (\code{barcode_gap_metrics_enhanced.csv} or
\code{ambiguous_sequences_tbl_cleaned.csv}) in the following locations, in order:
(1) the directory specified by the \env{RB_YZFISHDB_DATA_DIR} environment
variable, (2) the package's bundled \code{inst/extdata} directory, and
(3) \code{inst/extdata} relative to the current working directory.}
\item{species}{Optional character vector of species names to retain.}
\item{marker}{Optional character vector of marker names to retain.}
\item{action}{Optional character vector of resolution actions to retain
(e.g., \code{"keep"}, \code{"combine"}).}
}
\value{
A data frame containing the filtered contents of the CSV file. Returns all
rows if no filters are applied.
}
\seealso{
\code{\link{rb_barcode_gap}} and \code{\link{rb_ambiguity}} for the
SQLite-based equivalents.

\code{\link{rb_run_barcode_gap}} for generating barcode gap metrics from
scratch.

\code{\link{rb_resolve_ambiguous}} for resolving ambiguous sequences during
the curation pipeline.
}
\examples{
\dontrun{
# Read from a specific CSV
gap <- rb_read_barcode_gap("path/to/barcode_gap_metrics_enhanced.csv")

# Filter to COI results for two species
gap_subset <- rb_read_barcode_gap(
  species = c("Danio rerio", "Cyprinus carpio"),
  marker = "COI"
)

# Read ambiguity table, keeping only rows marked for manual review
amb <- rb_read_ambiguity(action = "keep")
}
}