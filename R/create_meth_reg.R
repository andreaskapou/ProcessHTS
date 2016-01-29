#' Create methylation regions for each gene promoter.
#'
#' \code{create_meth_reg} creates methylation regions using a combination
#'  of BS-Seq and RNA-Seq data. The BS-Seq data give information for the
#'  methylated CpGs individually, and the RNA-Seq data are used to locate the
#'  TSS of each gene and its promoter region.
#'
#' @param bs_data \code{\link[GenomicRanges]{GRanges}} object containing the
#'  BS-Seq data.
#' @param prom_data \code{\link[GenomicRanges]{GRanges}} object containing
#'  the processed RNA-Seq data, i.e. keeping the promoter regions around TSS.
#' @param upstream Integer defining the length of bps upstream of TSS.
#' @param downstream Integer defining the length of bps downstream of TSS.
#' @param num_CpG Optional integer defining the minimum number of CpGs that
#'  have to be in a methylated region. Regions with less than \code{n} CpGs
#'  are discarded.
#' @param sd_thresh Optional numeric defining the minimum standard deviation
#'  of the methylation change in a region.
#' @param ignore_strand Logical, whether or not to ignore strand information.
#' @param fmin Optional minimum range value for region location scaling.
#' @param fmax Optional maximum range value for region location scaling.
#'
#' @return A \code{meth_regions} object containing the following information:
#'  \itemize{
#'    \item{ \code{meth_data}: A list containing methylation data, where each
#'      each entry in the list consists of an L X 3 dimensional matrix, where:
#'      \enumerate{
#'        \item{ 1st column: Contains the locations of CpGs relative to TSS,
#'          where the range (min, max) of possible values is given, by the
#'          inputs fmin and fmax.
#'        }
#'        \item{ 2nd column: Contains the total reads of the CpG in the
#'          corresponding location.}
#'        \item{ 3rd column: Contains the methylated reads of the CpG in the
#'          corresponding location.}
#'      }
#'    }
#'    \item{ \code{prom_ind}: A vector storing the corresponding promoter
#'      indices.}
#'  }
#'
#' @seealso \code{\link{read_rna_encode_caltech}},
#'  \code{\link{read_bs_encode_haib}},
#'  \code{\link{create_encode_caltech_prom_reg}}
#'
#' @examples
#' # Load the RNA-Seq example dataset
#' rna_data <- rnaseq_data
#' prom_data <- create_encode_caltech_prom_reg(rna_data, upstream=-2000, downstream=2000)
#' bs_data <- rrbs_data
#'
#' meth_regions <- create_meth_reg(bs_data, prom_data)
#'
#' @export
create_meth_reg <- function(bs_data, prom_data, upstream = -100,
                                downstream = 100, num_CpG = 1, sd_thresh = 0,
                                  ignore_strand = FALSE, fmin = -1, fmax = 1){

  message("Creating methylation regions ...")
  assertthat::assert_that(is(bs_data, "GRanges"))
  assertthat::assert_that(is(prom_data, "GRanges"))
  if (upstream > 0 ){
    upstream <- -upstream
  }

  # Find overlaps between promoter regions and BS-Seq data -------
  overlaps <- GenomicRanges::findOverlaps(query   = prom_data,
                                          subject = bs_data,
                                          ignore.strand = ignore_strand)

  if (length(overlaps) < 2){
    stop("Not enough matches between the BS-Seq and RNA-Seq data.")
  }


  # Convert data in vector format for faster lookup --------------
  query_hits <- S4Vectors::queryHits(overlaps)
  subj_hits  <- S4Vectors::subjectHits(overlaps)

  prom_loc   <- unique(query_hits)     # Indices of promoter locations
  tss_loc    <- prom_data$tss          # TSS locations
  cpg_loc    <- GenomicRanges::ranges(bs_data)@start  # CpG locations
  tot_reads  <- bs_data$total_reads    # Total reads
  meth_reads <- bs_data$meth_reads     # Methylated reads


  # Initialize variables -----------------------------------------
  n            <- 1                         # Data points counter
  LABEL        <- FALSE                     # Flag variable
  meth_data    <- list()                    # List where data will be stored
  prom_counter <- 0                         # Promoter counter
  prom_ind     <- vector(mode = "integer")  # Vector of promoter indices
  cpg_ind      <- vector(mode = "integer")  # Vector of CpG indices
  cpg_ind      <- c(cpg_ind, subj_hits[1])  # Add the first subject hit

  for (i in 2:NROW(query_hits)){
    # If query hits is the same as the previous one
    if (query_hits[i] == query_hits[i - 1]){
      cpg_ind <- c(cpg_ind, subj_hits[i])  # Add subject hit
    }else{
      prom_counter <- prom_counter + 1  # Increase promoter counter
      LABEL <- TRUE
    }

    if(LABEL){

      # Only keep regions that have more than 'n' CpGs
      if (length(cpg_ind) > num_CpG){

        # If standard deviation of the methylation level is above threshold
        if (sd(meth_reads[cpg_ind] / tot_reads[cpg_ind]) > sd_thresh){

          # Promoter indices
          prom_ind <- c(prom_ind, prom_loc[prom_counter])

          # Locations of CpGs in the genome
          region <- cpg_loc[cpg_ind]

          # TSS location for promoter 'promCount'
          tss <- tss_loc[prom_loc[prom_counter]]

          # Extract strand information, i.e. direction
          strand_direct <- as.character(GenomicRanges::strand(
                                  prom_data[prom_loc[prom_counter]]))

          # Shift CpG locations relative to TSS
          center_data  <- center_loc(region = region,
                                     tss = tss,
                                     strand_direction = strand_direct)

          # In the "-" strand the order of the locations should change
          Order <- order(center_data)

          meth_data[[n]] <- matrix(0, nrow = length(cpg_ind), ncol = 3)

          # Store normalized locations of methylated CpGs in (fmin,fmax) region
          meth_data[[n]][ ,1] <- minmax_scaling(data = center_data[Order],
                                                xmin = upstream,
                                                xmax = downstream,
                                                fmin = fmin,
                                                fmax = fmax)

          # Store total reads in the corresponding locations
          meth_data[[n]][ ,2] <- tot_reads[cpg_ind][Order]

          # Store methylated reads in the corresponding locations
          meth_data[[n]][ ,3] <- meth_reads[cpg_ind][Order]

          # Increase data points counter
          n <- n + 1
        }
      }
      LABEL   <- FALSE
      cpg_ind <- vector(mode="integer")
      cpg_ind <- c(cpg_ind, subj_hits[i])
    }
  }
  meth_regions <- structure(list(meth_data = meth_data,
                                 prom_ind = prom_ind),
                            class = "meth_regions")
  message("Done!\n")
  return(meth_regions)
}


#' Center CpG locations relative to TSS
#'
#' \code{center_loc} centera CpG locations relative to TSS
#'
#' @param region CpG locations
#' @param tss TSS location
#' @param strand_direction Strand direction
#'
#' @return Centered location data relative to TSS
#'
center_loc <- function(region, tss, strand_direction){
  assertthat::assert_that(is.character(strand_direction))

  center <- region - tss
  if (identical(strand_direction, "-")){
    center  <- (-center)  # If '-' strand, swap CpG locations
  }
  return(center)
}
