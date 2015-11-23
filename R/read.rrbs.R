#' Read file containing \code{bed} formatted RRBS data
#'
#' \code{read.rrbs} reads a file containing methylation data from RRBS
#' experiments using the \code{\link{scan}} function. The RRBS file should be
#' in \code{bed} format.
#'
#' @param rrbs_cov Discard an optional integer to be used to disacrd low
#'   coverage reads.
#' @inheritParams read.rnaseq
#'
#' @return a list if \code{is_list} is TRUE, otherwise a
#'  \code{\link[GenomicRanges]{GRanges}} object.
#'
#' @seealso \code{\link{read.chrom_size}}, \code{\link{read.rnaseq}}
#'
#' @references
#'   \url{http://genome.ucsc.edu/cgi-bin/hgTables?db=hg19&hgta_group=regulation&hgta_track=wgEncodeHaibMethylRrbs&hgta_doSchema=describe+table+schema}
#'
#' @examples
#' # Get the location of the RRBS file
#' rrbs_file <- system.file("extdata", "rrbs.bed", package = "processHTS")
#' data <- read.rrbs(file=rrbs_file, is_list=FALSE)
#'
#' @export
read.rrbs <- function(file, is_del_chrom = FALSE, rrbs_cov = 0, is_list = TRUE){
  message("Reading file ", file, " ...")
  data_raw <- scan(file=file,
                   skip=1,
                   sep="\t",
                   what=list("character",  # Reference chromosome or scaffold
                             integer(),    # Start position in chromosome
                             NULL,         # End position in chromosome
                             NULL,         # Name of item
                             integer(),    # Score from 0-1000. Capped number of reads
                             "character",  # Strand : + or - or . for unknown
                             NULL,         # Start position
                             NULL,         # End position
                             NULL,         # Color value R,G,B
                             NULL,         # Number of reads or coverage
                             integer()     # Methylation percentage at this position
                   ))

  # Sorting data -----------------------------------------------
  # According to the following order priority:
  #     1: First 'chromosome',
  #     2: then 'position in chromosome' and
  #     3: finally 'strand'
  message("Sorting data...")
  entries <- c(1,2,5,6,11)  # Keep only data that are of interest

  Order <- with(data_raw, order(data_raw[[1]], data_raw[[2]], data_raw[[6]]))
  for (j in 1:length(entries)){
    entry <- entries[j]
    data_raw[[entry]] <- data_raw[[entry]][Order]
  }

  # Delete chromosomes X, Y and M -----------------------------
  if (is_del_chrom){
    message("Removing X, Y and M chromosomes...")
    chrom_ind <- which(data_raw[[1]] == "chrM")
    chrom_ind <- c(chrom_ind, which(data_raw[[1]] == "chrX"))
    chrom_ind <- c(chrom_ind, which(data_raw[[1]] == "chrY"))
    for (j in 1:length(entries)){
      entry <- entries[j]
      data_raw[[entry]] <- data_raw[[entry]][-chrom_ind]
    }
  }

  # Discard low coverage reads --------------------------------
  #   i.e CpGs with less than n total reads)
  if (rrbs_cov > 0){
    message("Discarding low coverage reads...")
    low_cov <- which(data_raw[[5]] < rrbs_cov)
    for (j in 1:length(entries)){
      entry <- entries[j]
      data_raw[[entry]] <- data_raw[[entry]][-low_cov]
    }
  }

  if (is_list){
    message("Done!\n")
    return(data_raw)
  }

  # Create a GRanges object -----------------------------------
  message("Creating GRanges object...")
  meth_data <- GenomicRanges::GRanges(seqnames = data_raw[[1]],
                      strand = data_raw[[6]],
                      ranges = IRanges::IRanges(start=data_raw[[2]], width=1),
                      total_reads = data_raw[[5]],
                      meth_reads = as.integer(round(0.01 * data_raw[[5]] * data_raw[[11]]))
               )

  message("Done!\n")
  return(meth_data)
}
