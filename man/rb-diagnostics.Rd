\name{rb_marker_coverage}
\alias{rb_marker_coverage}
\alias{rb_source_coverage}
\alias{rb_qc_summary}
\alias{rb_barcode_gap}
\alias{rb_ambiguity}
\title{Summarize reference coverage, QC, barcode gaps, and ambiguity from SQLite}
\description{
Generate diagnostic summaries from a connected SQLite reference database.
These functions query named tables inside the database and are intended for
use after the database has been fully built and populated.
}
\usage{
rb_marker_coverage(con, rank = "species", table_name = "yzfishdb_final")
rb_source_coverage(con, table_name = "yzfishdb_final")
rb_qc_summary(con, table_name = "yzfishdb_final", qc_table_name = "qc_reference_p2")
rb_barcode_gap(con, species = NULL, marker = NULL, table_name = "barcode_gap_metrics")
rb_ambiguity(con, table_name = "ambiguous_sequences")
}
\arguments{
\item{con}{A DBI SQLite connection returned by \code{\link{rb_connect}}.}
\item{rank}{Taxonomic rank for marker coverage summaries. Supported values are
\code{"species"}, \code{"genus"}, \code{"family"}, and \code{"order"}.}
\item{species}{Optional character vector of species names to filter barcode gap results.}
\item{marker}{Optional character vector of marker names to filter barcode gap results.}
\item{table_name}{Name of the main reference table to query. Defaults to
\code{"yzfishdb_final"}. Users building custom databases should set this to
the name of their own curated sequence table.}
\item{qc_table_name}{Name of the QC reference table. Defaults to
\code{"qc_reference_p2"}. If this table does not exist in the database,
the function falls back to querying \code{table_name}.}
}
\details{
All five functions query tables inside a connected SQLite database. If a
required table does not exist (e.g., \code{barcode_gap_metrics} has not yet
been populated), the function returns an empty data frame rather than
raising an error.

For users who have barcode gap or ambiguity results stored as standalone CSV
files (e.g., from an intermediate analysis step before database assembly),
use \code{\link{rb_read_barcode_gap}} and \code{\link{rb_read_ambiguity}} instead.
}
\value{
Data frames summarizing marker coverage, source contribution, quality-control
states, barcode-gap metrics, or curated ambiguous records.
}
\examples{
con <- rb_connect()
head(rb_marker_coverage(con))
head(rb_source_coverage(con))
rb_qc_summary(con)
rb_barcode_gap(con, marker = "COI")
rb_ambiguity(con)
rb_disconnect(con)
}