#' Wrapper method for processing BS-Seq and RNA-Seq Beatson HTS data
#'
#' \code{process_beatson} is a wrapper method for processing HTS data and
#' returning the methylation promoter regions and the corresponding gene
#' expression data for those promoter regions.
#'
#' @inheritParams process_haib_caltech
#'
#' @return A \code{processHTS} object which contains among others the following
#'  information:
#'  \itemize{
#'    \item{ \code{methyl_region}: A list containing the methylation regions,
#'      where each each entry in the list consists of an L X 3 dimensional
#'      matrix, where:
#'      \enumerate{
#'        \item{ 1st col: Contains the locations of the CpGs relative to TSS,
#'          where the range (min, max) of possible values is given, by the
#'          inputs fmin and fmax.
#'        }
#'        \item{ 2nd col: The total reads of the CpG in the corresponding
#'          location.}
#'        \item{ 3rd col: The methylated reads of the CpG in the corresponding
#'          location.}
#'      }
#'    }
#'    \item{ \code{prom_region}: A \code{\link[GenomicRanges]{GRanges}} object
#'      containing corresponding annotated promoter regions for each entry of
#'      the \code{methyl_region} list..
#'
#'    }
#'    \item{ \code{rna_data}: A \code{\link[GenomicRanges]{GRanges}} object
#'      containing the corresponding RNA-Seq data for each entry of the
#'      \code{methyl_region} list.}
#'  }
#'
#' @examples
#' # Get the location of the files
#' rrbs_file <- system.file("extdata", "rrbsH1hESC.bed", package = "processHTS")
#' rnaseq_file <- system.file("extdata", "rnaseqH1hESC.bed", package = "processHTS")
#' data <- process_haib_caltech(rrbs_file, rnaseq_file)
#'
#' @export
process_beatson <- function(bs_files, rna_files, chrom_size_file = NULL,
                            chr_discarded = NULL, upstream = -100,
                            downstream = 100, min_bs_cov = 0,
                            max_bs_cov = 1000, cpg_density = 1,
                            sd_thresh = 0, ignore_strand = FALSE,
                                                    fmin = -1, fmax = 1){

  # Process BS-Seq file and return data in the required format
  bs_data <- preprocess_bs_bismark_cov(files         = bs_files,
                                       chr_discarded = chr_discarded,
                                       min_bs_cov    = min_bs_cov,
                                       max_bs_cov    = max_bs_cov)

  # Read the chromosome size file, if it is supplied
  if (!is.null(chrom_size_file)){
    chrom_size <- read_chrom_size(file = chrom_size_file)
  }else{
    chrom_size <- NULL
  }

  # Read RNA-Seq BED file
  rna_data <- read_rna_beatson(file          = rna_files,
                               chr_discarded = chr_discarded,
                               is_GRanges    = TRUE)

  # Create promoter regions
  prom_reg <- create_prom_region(annot_data = rna_data,
                                 chrom_size = chrom_size,
                                 upstream   = upstream,
                                 downstream = downstream)

  # Create methylation regions data
  methyl_reg <- create_methyl_region(bs_data       = bs_data,
                                     prom_region   = prom_reg,
                                     cpg_density   = cpg_density,
                                     sd_thresh     = sd_thresh,
                                     ignore_strand = ignore_strand,
                                     fmin          = -1,
                                     fmax          = 1)

  # Keep only the corresponding gene expression data
  rna_data <- rna_data[methyl_reg$prom_ind]
  # Keep only the corresponding gene annotation data
  prom_reg <- prom_reg[methyl_reg$prom_ind]

  # Create object
  obj <- structure(list(methyl_region = methyl_reg$meth_data,
                        prom_region   = prom_reg,
                        rna_data      = rna_data,
                        upstream      = upstream,
                        downstream    = downstream,
                        cpg_density   = cpg_density,
                        sd_thresh     = sd_thresh,
                        fmin          = fmin,
                        fmax          = fmax),
                   class = "processHTS")
  return(obj)
}
