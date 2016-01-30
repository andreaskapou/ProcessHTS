#' Read file containing Bismark Cov formatted BS-Seq data
#'
#' \code{read_bs_bismark_cov} reads a file containing methylation data from
#'  BS-Seq experiments using the \code{\link{scan}} function. The BS-Seq file
#' should be in Bismark Cov format.
#'
#' @inheritParams read_bs_encode_haib
#'
#' @return a \code{\link[GenomicRanges]{GRanges}} object if \code{is_GRanges}
#'  is TRUE, otherwise a data.frame object.
#'
#' @seealso \code{\link{pool_bs_bismark_cov_rep}},
#'  \code{\link{preprocess_bs_bismark_cov}}
#'
#' @references
#'   \url{http://rnbeads.mpi-inf.mpg.de/data/RnBeads.pdf}
#'
#' @examples
#' # Get the location of the bismark file
#' bismark_file <- system.file("extdata", "bism_rep1.bed", package = "processHTS")
#' bs_data <- read_bs_bismark_cov(file=bismark_file, is_GRanges=TRUE)
#'
#' @export
read_bs_bismark_cov <- function(file, chr_discarded = NULL, is_GRanges = TRUE){
  message("Reading file ", file, " ...")
  bs_data <- read.table(file = file,
                        header = FALSE,
                        sep = "\t",
                        col.names = c("chr", "start",
                                      "meth_reads", "unmeth_reads"),
                        comment.char = "")


  # Remove selected chromosomes  ------------------------------
  bs_data <- discard_chr(x = bs_data, chr_discarded = chr_discarded)


  # Sorting data -----------------------------------------------
  # With order priority: 1. chr, 2. start
  message("Sorting BS-Seq data ...")
  bs_data <- bs_data[with(bs_data, order(bs_data$chr,
                                         bs_data$start)), ]

  # Get sequential row numbers
  row.names(bs_data) <- NULL


  if (is_GRanges){
    # Create a GRanges object ---------------------------------
    message("Creating GRanges object ...")
    bs_data <- GenomicRanges::GRanges(seqnames = bs_data$chr,
                      ranges = IRanges::IRanges(start=bs_data$start, width=1),
                      meth_reads   = bs_data$meth_reads,
                      unmeth_reads = bs_data$unmeth_reads)
  }
  message("Done!\n")
  return(bs_data)
}
