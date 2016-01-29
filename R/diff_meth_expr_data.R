#' Wrapper method for processing RRBS and RNA-Seq HTS data.
#'
#' \code{diff_meth_expr_data} is a wrapper method for processing HTS data and
#' returning the methylation promoter regions and the corresponding gene
#' expression levels for those promoter regions.
#'
#' @param rrbs_files The names of the RRBS '.bed' formatted file to read
#'  data values from.
#' @param rnaseq_files The names of the RNA-Seq '.bed' formatted file to
#'  read data values from.
#' @param chrom_size_file Optional name of the file containing genome
#'  chromosome sizes
#' @param is_del_chrom Optional logical: if TRUE, delete chromosomes X,
#'  Y and M.
#' @param rrbs_cov Optional integer to disacrd low coverage reads.
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
#' @export
diff_meth_expr_data <- function(rrbs_files, rnaseq_files, chrom_size_file = NULL,
                           upstream = -100, downstream = 100, num_CpG = 1,
                           sd_thresh = 0, is_del_chrom = FALSE, rrbs_cov = 0,
                           fmin = -1, fmax = 1){

  hts_data <- list(contr = NULL, treat = NULL)
  message("Control ...")
  hts_data$contr <- meth_expr_data(rrbs_file       = rrbs_files[1],
                                   rnaseq_file     = rnaseq_files[1],
                                   chrom_size_file = chrom_size_file,
                                   upstream        = upstream,
                                   downstream      = downstream,
                                   num_CpG         = num_CpG,
                                   sd_thresh       = sd_thresh,
                                   is_del_chrom    = is_del_chrom,
                                   rrbs_cov        = rrbs_cov,
                                   fmin            = fmin,
                                   fmax            = fmax)
  message("Treatment ...")
  hts_data$treat <- meth_expr_data(rrbs_file       = rrbs_files[2],
                                   rnaseq_file     = rnaseq_files[2],
                                   chrom_size_file = chrom_size_file,
                                   upstream        = upstream,
                                   downstream      = downstream,
                                   num_CpG         = num_CpG,
                                   sd_thresh       = sd_thresh,
                                   is_del_chrom    = is_del_chrom,
                                   rrbs_cov        = rrbs_cov,
                                   fmin            = fmin,
                                   fmax            = fmax)

  # Find overlaps between promoter regions and RRBS data
  overlaps <- GenomicRanges::findOverlaps(query   = hts_data$contr$expr_data,
                                          subject = hts_data$treat$expr_data,
                                          ignore.strand = FALSE)
  hits <- as.data.frame(overlaps)   # Create a HITS matrix for faster computations
  hts_data$contr$expr_data <- hts_data$contr$expr_data[hits[ ,1]]
  hts_data$treat$expr_data <- hts_data$treat$expr_data[hits[ ,2]]

  contr_meth <- list()
  treat_meth <- list()
  for (i in 1:NROW(hits)){
    contr_meth[[i]] <- hts_data$contr$meth_data[[hits[i,1]]]
    treat_meth[[i]] <- hts_data$treat$meth_data[[hits[i,2]]]
  }
  hts_data$contr$meth_data <- contr_meth
  hts_data$treat$meth_data <- treat_meth

  class(hts_data) <- "diff_meth_expr"
  return(hts_data)
}
