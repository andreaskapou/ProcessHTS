#' Process differential BS-Seq and RNA-Seq Beatson HTS data
#'
#' \code{process_diff_beatson} is a wrapper method ...
#'
#' @param bs_contr_files The name of the control BS-Seq '.bed' formatted files
#'  to read data values from.
#' @param bs_treat_files The name of the treatment BS-Seq '.bed' formatted
#'  files to read data values from.
#' @param rna_contr_files The name of the control RNA-Seq '.bed' formatted
#'  files to read data values from.
#' @param rna_treat_files The name of the treatment RNA-Seq '.bed' formatted
#'  files to read data values from.
#' @param chrom_size_file Optional name of the file containing genome
#'  chromosome sizes
#' @param chr_discarded A vector with chromosome names to be discarded.
#' @param upstream Integer defining the length of bp upstream of TSS.
#' @param downstream Integer defining the length of bp downstream of TSS.
#' @param min_bs_cov The minimum number of reads mapping to each CpG site.
#' @param max_bs_cov The maximum number of reads mapping to each CpG site.
#' @inheritParams create_methyl_region
#'
#' @return A \code{diff_processHTS} object which contains among others the
#'  following information:
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
#' bs_file <- system.file("extdata", "bism_rep1.bed", package = "processHTS")
#' rna_file <- system.file("extdata", "rna_beatson.bed", package = "processHTS")
#' #data <- process_beatson(bs_file, rna_file)
#'
#' @export
process_diff_beatson <- function(bs_contr_files, bs_treat_files,
                                 rna_contr_files, rna_treat_files,
                                 chrom_size_file = NULL, chr_discarded = NULL,
                                 upstream = -100, downstream = 100,
                                 min_bs_cov = 0, max_bs_cov = 1000,
                                 cpg_density = 1, sd_thresh = 0,
                                 ignore_strand = FALSE, fmin = -1, fmax = 1){

  # Process control dataset
  contr_data <- process_beatson(bs_files        = bs_contr_files,
                                rna_files       = rna_contr_files,
                                chrom_size_file = chrom_size_file,
                                chr_discarded   = chr_discarded,
                                upstream        = upstream,
                                downstream      = downstream,
                                min_bs_cov      = min_bs_cov,
                                max_bs_cov      = max_bs_cov,
                                cpg_density     = cpg_density,
                                sd_thresh       = sd_thresh,
                                ignore_strand   = ignore_strand,
                                fmin            = fmin,
                                fmax            = fmax)

  # Process treatment dataset
  treat_data <- process_beatson(bs_files        = bs_treat_files,
                                rna_files       = rna_treat_files,
                                chrom_size_file = chrom_size_file,
                                chr_discarded   = chr_discarded,
                                upstream        = upstream,
                                downstream      = downstream,
                                min_bs_cov      = min_bs_cov,
                                max_bs_cov      = max_bs_cov,
                                cpg_density     = cpg_density,
                                sd_thresh       = sd_thresh,
                                ignore_strand   = ignore_strand,
                                fmin            = fmin,
                                fmax            = fmax)

  # Find overlaps between control and treatment data
  overlaps <- GenomicRanges::findOverlaps(query   = contr_data$rna_data,
                                          subject = treat_data$rna_data,
                                          select  = "first")

  # Get only the subset of overlapping control data
  keep_non_na <- which(!is.na(overlaps))

  # Create final object
  obj <- structure(list(contr       = list(),
                        treat       = list(),
                        upstream    = treat_data$upstream,
                        downstream  = treat_data$downstream,
                        cpg_density = treat_data$cpg_density,
                        sd_thresh   = treat_data$sd_thresh,
                        fmin        = treat_data$fmin,
                        fmax        = treat_data$fmax),
                   class = "diff_processHTS")

  # Keep only methylation regions that overlap on both control and treatment cases
  obj$contr$methyl_region <- contr_data$methyl_region[keep_non_na]
  obj$treat$methyl_region <- treat_data$methyl_region[overlaps[keep_non_na]]

  # Keep only promoter regions that overlap on both control and treatment cases
  obj$contr$prom_region <- contr_data$prom_region[keep_non_na]
  obj$treat$prom_region <- treat_data$prom_region[overlaps[keep_non_na]]

  # Keep only RNA-Seq data that overlap on both control and treatment cases
  obj$contr$rna_data <- contr_data$rna_data[keep_non_na]
  obj$treat$rna_data <- treat_data$rna_data[overlaps[keep_non_na]]

  return(obj)
}
