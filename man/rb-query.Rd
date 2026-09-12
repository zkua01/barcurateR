\name{rb_get_sequences}
\alias{rb_get_sequences}
\alias{rb_list_taxa}
\alias{rb_filter_reference}
\title{Query and filter reference sequences}
\description{
Retrieve sequence records from the curated reference table and apply
reproducible filters for taxonomy, marker, source, occurrence status,
quality-control state, and sequence length.
}
\usage{
rb_get_sequences(con, species = NULL, genus = NULL, family = NULL,
                 marker = NULL, source = NULL, occurrence = NULL,
                 qc_flag = "pass", exact = FALSE, table_name = "yzfishdb_final")
rb_list_taxa(con, rank = "species", occurrence = NULL,
             marker = NULL, qc_flag = "pass", exact = FALSE,
             table_name = "yzfishdb_final")
rb_filter_reference(data, marker = NULL, source = NULL, occurrence = NULL,
                    min_length = NULL, max_length = NULL, qc_flag = "pass",
                    exact = FALSE)
}
\arguments{
\item{con}{A DBI SQLite connection returned by \code{\link{rb_connect}}.}
\item{species, genus, family}{Optional taxonomic filters.}
\item{marker}{Optional marker filter, for example \code{"12S"}, \code{"16S"},
or \code{"COI"}.}
\item{source}{Optional source database or collection filter.}
\item{occurrence}{Optional occurrence filter, for example \code{"native"} or
\code{"introduced"}.}
\item{qc_flag}{Character vector of quality-control state(s) to retain. Defaults
to \code{"pass"}.}
\item{exact}{Logical. If \code{TRUE}, enforces strict equality (\code{==}) for
\code{qc_flag} matching. If \code{FALSE} (default), allows partial matching for
composite flags.}
\item{rank}{Taxonomic rank to list. Supported ranks are kingdom, phylum, class,
order, family, genus, and species.}
\item{data}{A reference sequence data frame, usually returned by
\code{rb_get_sequences()}.}
\item{min_length, max_length}{Optional sequence length limits in base pairs.}
\item{table_name}{Name of the main reference table to query. Defaults to
\code{"yzfishdb_final"}. Users building custom databases should set this to
the name of their own curated sequence table.}
}
\details{
\code{rb_get_sequences} and \code{rb_list_taxa} execute SQL queries against the
SQLite database. \code{rb_filter_reference} operates on an in-memory data frame
and does not require a database connection (note: \code{rb_filter_reference} does
not accept a \code{table_name} argument because it operates on data already
loaded into R).

When \code{exact = FALSE} (the default), composite QC flags such as
\code{"tree_high_length_variation|manual_review_needed"} will be matched by
either component. When \code{exact = TRUE}, only rows with an exactly matching
\code{qc_flag} value are retained.
}
\value{
\itemize{
\item \code{rb_get_sequences}: A data frame of sequence records matching the
  specified filters.
\item \code{rb_list_taxa}: A data frame with two columns: the specified taxonomic
  \code{rank} and \code{n_sequences}.
\item \code{rb_filter_reference}: A filtered data frame retaining the same column
  structure as the input \code{data}.
}
}
\examples{
con <- rb_connect()

# Query with default settings
native_12s <- rb_get_sequences(con, marker = "12S", occurrence = "native")
head(native_12s)

# Partial QC flag matching (default)
multi_qc <- rb_get_sequences(con, qc_flag = c("pass", "manual_review_needed"))

# Exact QC flag matching
strict_qc <- rb_get_sequences(con, qc_flag = c("pass", "manual_review_needed"),
                               exact = TRUE)

# Query a custom table
custom_data <- rb_get_sequences(con, marker = "COI",
                                 table_name = "my_custom_reference")

# List taxa and filter in-memory
species_counts <- rb_list_taxa(con, rank = "species", marker = "12S")
filtered_refs <- rb_filter_reference(native_12s, min_length = 80, max_length = 250)

rb_disconnect(con)
}