#' Read file containing genome chromosome sizes.
#'
#' \code{read_chrom_size} reads a file containing genome chromosome sizes
#' using the \code{\link{scan}} function.
#'
#' @param file The name of the file to read data values from.
#'
#' @return data frame.
#'
#' @seealso \code{\link{read_rna_encode_caltech}}, \code{\link{read_bs_encode_haib}}
#'
#' @examples
#' # Get the location of the hg19 file
#' hg19_file <- system.file("extdata", "hg19", package = "processHTS")
#' data <- read_chrom_size(file=hg19_file)
#'
#' \dontrun{
#' read.chrom_size(20)
#' }
#'
#' @export
read_chrom_size <- function(file){
  message("Reading file ", file, " ...")
  data_raw <- scan(file=file,
                   sep="\t",
                   what=list("character",  # Reference chromosome
                             integer()     # Chromosome size
                   ))

  # Create a data frame with the chromosome sizes
  data <- as.data.frame(data_raw[[2]], row.names=data_raw[[1]])
  colnames(data) <- "size"

  message("Done!\n")
  return(data)
}
