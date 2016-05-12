#' Read Bismark Cov formatted BS-Seq file
#'
#' \code{read_bs_bismark_cov} reads a file containing methylation data from
#' BS-Seq experiments using the \code{\link[data.table]{fread}} function. The
#' BS-Seq file should be in Bismark Cov format. Read the Important section below
#' on when using this function.
#'
#' @inheritParams read_bs_encode_haib
#'
#' @return A \code{\link[GenomicRanges]{GRanges}} object if \code{is_GRanges} is
#'   TRUE, otherwise a \code{\link[data.table]{data.table}} object.
#'
#'   The GRanges object contains two additional metadata columns: \itemize{
#'   \item \code{total_reads}: total reads mapped to each genomic location.
#'   \item \code{meth_reads}: methylated reads mapped to each genomic location.
#'   } These columns can be accessed as follows:
#'   \code{granges_object$total_reads}
#'
#' @section Important: Unless you want to create a different workflow when
#'   processing the BS-Seq data, you should NOT call this function, since this
#'   is a helper function. Instead you should call the
#'   \code{\link{preprocess_bs_seq}} function.
#'
#' @seealso \code{\link{pool_bs_seq_rep}}, \code{\link{preprocess_bs_seq}}
#'
#' @references \url{http://rnbeads.mpi-inf.mpg.de/data/RnBeads.pdf}
#'
#' @examples
#' # Get the location of the bismark file
#' bism_file <- system.file("extdata", "bism_rep1.bed", package = "processHTS")
#' bs_data <- read_bs_bismark_cov(file = bism_file, is_GRanges = TRUE)
#'
#' @export
read_bs_bismark_cov <- function(file, chr_discarded = NULL, is_GRanges = TRUE){
  message("Reading file ", file, " ...")
  bs_data <- data.table::fread(input = file,
                               sep = "\t",
                               header = FALSE,
                               col.names = c("chr", "start", "meth_reads",
                                             "unmeth_reads"))


  # Remove selected chromosomes  -------------------------------
  bs_data <- discard_chr(x = bs_data, chr_discarded = chr_discarded)


  # Sorting data -----------------------------------------------
  # With order priority: 1. chr, 2. start
  message("Sorting BS-Seq data ...")
  bs_data <- bs_data[order(bs_data$chr, bs_data$start)]


  if (is_GRanges){
    # Create a GRanges object ---------------------------------
    message("Creating GRanges object ...")
    bs_data <- GenomicRanges::GRanges(seqnames = bs_data$chr,
                      ranges = IRanges::IRanges(start=bs_data$start, width=1),
                      total_reads = bs_data$meth_reads + bs_data$unmeth_reads,
                      meth_reads  = bs_data$meth_reads)
  }
  message("Finished reading BS-Seq file!\n")
  return(bs_data)
}
