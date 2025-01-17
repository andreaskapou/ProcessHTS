% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_enh_beatson.R
\name{read_enh_beatson}
\alias{read_enh_beatson}
\title{Read Beatson bed formatted enhancer file}
\usage{
read_enh_beatson(enh_file, gene_file, chr_discarded = NULL,
  is_GRanges = TRUE)
}
\arguments{
\item{enh_file}{File containing enhancer annotation data}

\item{gene_file}{File containing gene annotation data}

\item{chr_discarded}{A vector with chromosome names to be discarded.}

\item{is_GRanges}{Logical: if TRUE a GRanges object is returned, otherwise a
data.frame object is returned.}
}
\value{
A \code{\link[GenomicRanges]{GRanges}} object if \code{is_GRanges} is
  TRUE, otherwise a \code{\link[data.table]{data.table}} object.

  The GRanges object contains three additional metadata columns: \itemize{
  \item \code{ensembl_id}: Ensembl IDs of each gene promoter. \item
  \code{gene_name}: Gene name.} These columns can be accessed as follows:
  \code{granges_object$ensembl_id}
}
\description{
\code{read_enh_beatson} reads a file containing enhancer annotation data
using the \code{\link[data.table]{fread}} function.
}
\examples{
# Get the location of the RNA-Seq file
enh_beatson_file <- system.file("extdata", "enhancers_beatson.bed", package = "processHTS")
rna_beatson_file <- system.file("extdata", "rna_beatson2.bed", package = "processHTS")

enh_beatson_data <- read_enh_beatson(enh_beatson_file, rna_beatson_file, chr_discarded = "chrX", is_GRanges=TRUE)

}
\author{
C.A.Kapourani \email{C.A.Kapourani@ed.ac.uk}
}
\seealso{
\code{\link{read_annot_beatson}}, \code{\link{read_bs_bismark_cov}}
}

