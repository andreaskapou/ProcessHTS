#' Wrapper method for processing RRBS and RNA-Seq HTS data.
#'
#' \code{meth_expr_data} is a wrapper method for processing HTS data and
#' returning the methylation promoter regions and the corresponding gene
#' expression levels for those promoter regions.
#'
#' @param rrbs_file The name of the RRBS '.bed' formatted file to read
#'  data values from.
#' @param rnaseq_file The name of the RNA-Seq '.bed' formatted file to
#'  read data values from.
#' @param chrom_size_file Optional name of the file containing genome
#'  chromosome sizes
#' @param is_del_chrom Optional logical: if TRUE, delete chromosomes X,
#'  Y and M.
#' @param rrbs_cov Optional integer to disacrd low coverage reads.
#' @param tss_data Logical, if the data read are TSS or gene
#' @inheritParams create_meth_reg
#'
#' @return A \code{methExpr} object which contains, among others the following
#'  information:
#'  \itemize{
#'    \item{ \code{meth_data}: A list containing the methylation data, where each
#'      each entry in the list consists of an L X 3 dimensional matrix, where:
#'      \enumerate{
#'        \item{ 1st column: Contains the locations of the CpGs relative to TSS,
#'          where the range (min, max) of possible values is given, by the
#'          inputs fmin and fmax.
#'        }
#'        \item{ 2nd column: The total reads of the CpG in the corresponding
#'          location.}
#'        \item{ 3rd column: The methylated reads of the CpG in the corresponding
#'          location.}
#'      }
#'    }
#'    \item{ \code{expr_data}: A \code{\link[GenomicRanges]{GRanges}} object containing the
#'    corresponding gene expression data for each entry of the \code{meth_data} list.}
#'  }
#'
#' @examples
#' # Get the location of the files
#' #rrbs_file <- system.file("extdata", "rrbsH1hESC.bed", package = "processHTS")
#' #rnaseq_file <- system.file("extdata", "rnaseqH1hESC.bed", package = "processHTS")
#' #data <- meth_expr_data(rrbs_file, rnaseq_file)
#'
#' @export
meth_expr_data <- function(bs_file, rna_file, chrom_size_file = NULL,
                           upstream = -100, downstream = 100,
                           num_CpG = 1, sd_thresh = 0, chr_discarded = NULL,
                           rrbs_cov = 0, ignore_strand = FALSE,
                           tss_data = FALSE, fmin = -1, fmax = 1){


  # Process BS-Seq file and return data in the required format
  bs_data <- read_bs_encode_haib(file          = bs_file,
                                 chr_discarded = chr_discarded,
                                 is_GRanges    = TRUE)
  # Read RNA-Seq BED file
  rna_data <- read_rna_encode_caltech(file          = rna_file,
                                      chr_discarded = chr_discarded,
                                      is_GRanges    = TRUE)

  # Read the chromosome size file, if it is supplied
  if (!is.null(chrom_size_file)){
    chrom_size <- read_chrom_size(file = chrom_size_file)
  }else{
    chrom_size <- NULL
  }

  # Read TSS BED file and return a promoter region n upstream
  # and m downstream of Transcription Start Site (TSS)
  if (tss_data){
    prom_data <- create_tss_reg(data       = rnaseq_data,
                                    chrom_size = chrom_size,
                                    upstream   = upstream,
                                    downstream = downstream)
  }else{
    prom_data <- create_prom_regions(data       = rnaseq_data,
                                     chrom_size = chrom_size,
                                     upstream   = upstream,
                                     downstream = downstream)
  }
  # Remove objects from workspace
  rm(chrom_size)
  rm(rnaseq_data)

  # Create methylation regions data
  meth_regions <- create_meth_reg(rrbs_data     = rrbs_data,
                                      promoter_data = prom_data,
                                      upstream      = upstream,
                                      downstream    = downstream,
                                      num_CpG       = num_CpG,
                                      sd_thresh     = sd_thresh,
                                      ignore_strand = ignore_strand,
                                      fmin          = -1,
                                      fmax          = 1)

  # Keep only the corresponding gene expression data
  expr_data <- prom_data[meth_regions$prom_ind]

  # Create object
  meth_expr <- list(meth_data  = meth_regions$meth_data,
                    expr_data  = expr_data,
                    upstream   = upstream,
                    downstream = downstream,
                    num_CpG    = num_CpG,
                    sd_thresh  = sd_thresh,
                    rrbs_cov   = rrbs_cov,
                    fmin       = fmin,
                    fmax       = fmax)

  class(meth_expr) <- "methExpr"
  return(meth_expr)
}
