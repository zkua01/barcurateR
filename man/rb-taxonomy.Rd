\name{rb_taxonomy}
\alias{rb_taxonomy}
\alias{rb_build_taxonomy_string}
\title{Read taxonomy and build pipeline taxonomy strings}
\description{
Inspect distinct taxonomy records and convert seven-rank taxonomy columns into
plain or rank-prefixed taxonomy strings for downstream export formats.
}
\usage{
rb_taxonomy(con, species = NULL, qc_flag = "pass", exact = FALSE,
            table_name = "yzfishdb_final")
rb_build_taxonomy_string(data, style = c("plain", "rank_prefix"))
}
\arguments{
\item{con}{A DBI SQLite connection returned by \code{\link{rb_connect}}.}
\item{species}{Optional character vector of species names to filter.}
\item{qc_flag}{Character vector of quality-control state(s) to retain. Defaults
to \code{"pass"}.}
\item{exact}{Logical. If \code{TRUE}, enforces strict equality for \code{qc_flag}
matching. If \code{FALSE} (default), allows partial matching.}
\item{table_name}{Name of the main reference table to query. Defaults to
\code{"yzfishdb_final"}. Users building custom databases should set this to
the name of their own curated sequence table.}
\item{data}{A data frame containing \code{kingdom}, \code{phylum}, \code{class},
\code{order}, \code{family}, \code{genus}, and \code{species} columns.}
\item{style}{Taxonomy string style. \code{"plain"} returns semicolon-separated
rank values; \code{"rank_prefix"} returns strings such as \code{k__Animalia}.}
}
\value{
\code{rb_taxonomy()} returns a data frame of distinct taxonomy records
containing the standard seven taxonomic ranks, along with \code{occurrence}
and \code{habitat} columns.

\code{rb_build_taxonomy_string()} returns a character vector of formatted
taxonomy strings, with one string per row in the input \code{data}.
}
\examples{
con <- rb_connect()

tax <- rb_taxonomy(con)
head(tax)

# Query a custom table
custom_tax <- rb_taxonomy(con, table_name = "my_curated_table")

# Build taxonomy strings
plain_strings <- rb_build_taxonomy_string(head(tax), style = "plain")
prefixed_strings <- rb_build_taxonomy_string(head(tax), style = "rank_prefix")

rb_disconnect(con)
}