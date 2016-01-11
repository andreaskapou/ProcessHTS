#' Read file containing \code{bed} formatted RNA-Seq data
#'
#' \code{read.rnaseq} reads a file containing promoter annotation data
#' from RNA-Seq experiments using the \code{\link{scan}} function. The RNA-Seq
#' file should be in \code{bed} format, e.g. use \code{gtf2bed} tool if your
#' initial file is in \code{gtf} format.
#'
#' @param file The name of the file to read data values from.
#' @param is_del_chrom Logical: if TRUE, delete chromosomes X, Y and M.
#' @param is_list Logical: if TRUE a list is returned, otherwise
#'  a GRanges object.
#'
#' @return a list if \code{is_list} is TRUE, otherwise a
#'  \code{\link[GenomicRanges]{GRanges}} object.
#'
#' @seealso \code{\link{read.chrom_size}}, \code{\link{read.rrbs}}
#'
#' @examples
#' # Get the location of the RNA-Seq file
#' rnaseq_file <- system.file("extdata", "rnaseq.bed", package = "processHTS")
#' data <- read.rnaseq(file=rnaseq_file, is_list=FALSE)
#'
#' @export
read.rnaseq <- function(file, is_del_chrom = FALSE, is_list = TRUE){
  message("Reading file ", file, " ...")
  data_raw <- scan(file=file,
                   sep="\t",
                   what=list("character",  # Reference chromosome
                             integer(),    # Start position in chromosome
                             integer(),    # End position in chromosome
                             "character",  # Gene ENSEMBL id
                             numeric(),    # Expression level
                             "character",  # Strand : + or - or . for unknown
                             NULL,         # Source, e.g. HAVANA
                             NULL,         # Type of feature, e.g. gene
                             NULL,         # No information
                             "character"   # Metadata
                   ))

  # Extract FPKM from each gene --------------------------------
  message("Extracting FPKM...")
  data_raw[[11]] <- vector(mode = "numeric")
  for (i in 1:length(data_raw[[10]])){
    data_raw[[10]][i] <- extract_fpkm(data_raw[[10]][i])
  }

  # Sorting data -----------------------------------------------
  # According to the following order priority:
  #     1: First 'chromosome',
  #     2: then 'position in chromosome' and
  #     3: finally 'strand'
  message("Sorting data...")
  entries <- c(1,2,3,4,5,6,10)  # Keep only data that are of interest

  Order <- with(data_raw, order(data_raw[[1]], data_raw[[2]], data_raw[[6]]))
  for (j in 1:length(entries)){
    entry <- entries[j]
    data_raw[[entry]] <- data_raw[[entry]][Order]
  }

  # Delete chromosomes X, Y and M -----------------------------------
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

  if (is_list){
    message("Done!\n")
    return(data_raw)
  }

  # Create a GRanges object -----------------------------------
  message("Creating GRanges object...")
  expr_data <- GenomicRanges::GRanges(seqnames = data_raw[[1]],
                      strand = data_raw[[6]],
                      ranges = IRanges::IRanges(start=data_raw[[2]], end=data_raw[[3]]),
                      gene_id = data_raw[[4]],
                      gene_expr = data_raw[[5]],
                      gene_fpkm = data_raw[[10]]
                )

  message("Done!\n")
  return(expr_data)
}
