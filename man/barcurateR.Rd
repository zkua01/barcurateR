\name{barcurateR}
\alias{barcurateR}
\docType{package}
\title{Standardize, Curate, Audit, and Deploy DNA reference databases}
\description{
\pkg{barcurateR} provides tools to standardize, curate, query, audit, and export DNA reference databases, and assign 
}
\details{
The package is organized into function families:
\itemize{
  \item database installation and cache helpers: \code{\link{rb_install_db}};
  \item SQLite connection and table inspection: \code{\link{rb_connect}};
  \item reference querying and filtering: \code{\link{rb_get_sequences}};
  \item taxonomy strings: \code{\link{rb_build_taxonomy_string}};
  \item pipeline-ready FASTA exports: \code{\link{rb_export_fasta}};
  \item database diagnostics: \code{\link{rb_marker_coverage}};
  \item ASV assignment and community matrices: \code{\link{rb_assign_edna}}.
}
}
\examples{
\dontrun{
library(barcurateR)
rb_install_db(timeout = 7200, retries = 5)
con <- rb_connect()
rb_tables(con)
rb_disconnect(con)
}
}
