#' Create methylation regions for each gene promoter.
#'
#' \code{create_meth_regions} creates methylation regions data using a
#' combination of RRBS and RNA-Seq data. The RRBS data give information for
#' the methylated CpGs individually, and the RNA-Seq data are used to locate
#' the Transcription Start Sites (TSS) of each gene and its promoter region.
#' The RRBS and RNA-Seq data should be a \code{\link[GenomicRanges]{GRanges}}
#' object.
#'
#' @param rrbs_data \code{\link[GenomicRanges]{GRanges}} object containing the
#'  RRBS data.
#' @param promoter_data \code{\link[GenomicRanges]{GRanges}} object containing
#'  the processed RNA-Seq data, i.e. keeping the promoter regions around the TSS.
#' @param upstream Integer defining the length of bps upstream of TSS.
#' @param downstream Integer defining the length of bps downstream of TSS.
#' @param num_CpG Optional integer defining the minimum number of CpGs that
#'  have to be in a methylated region. Regions with less than \code{n} CpGs
#'  are discarded.
#' @param sd_thresh Opotional numeric defining the minimum standard deviation
#'  of the methylation change in a region.
#' @param fmin Optional minimum range value for region location scaling.
#' @param fmax Optional maximum range value for region location scaling.
#'
#' @return A \code{methRegions} object which contains the following information:
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
#'    \item{ \code{prom_ind}: A vector storing the corresponding promoter indices.}
#'  }
#'
#' @seealso \code{\link{read.rnaseq}}, \code{\link{read.rrbs}}, \code{\link{create_prom_regions}}
#'
#' @examples
#' # Load the RNA-Seq example dataset
#' data <- rnaseq_data
#' promoter_data <- create_prom_regions(data, upstream=-2000, downstream=2000)
#' rrbs_data <- rrbs_data
#'
#' meth_regions <- create_meth_regions(rrbs_data, promoter_data)
#'
#' @export
create_meth_regions <- function(rrbs_data, promoter_data, upstream = -100, downstream = 100,
                                num_CpG = 1, sd_thresh = 0, fmin = -1, fmax = 1){
  assertthat::assert_that(is(rrbs_data, "GRanges"))
  assertthat::assert_that(is(promoter_data, "GRanges"))
  if (upstream > 0 ){
    upstream <- -upstream
  }

  # Find overlaps between promoter regions and RRBS data
  overlaps <- GenomicRanges::findOverlaps(query   = promoter_data,
                                          subject = rrbs_data,
                                          ignore.strand = FALSE)

  if (length(overlaps) < 2){
    stop("Not enough matches between the RRBS data and RNA-Seq data. Check again the file names provided.")
  }

  message("Converting Granges objects in matrices...")
  hits <- as.data.frame(overlaps)   # Create a HITS matrix for faster computations
  prom_loc <- unique(hits[,1])  # Keep promoter locations

  # Convert data in matrix format for faster lookup
  tss_loc    <- as.vector(GenomicRanges::elementMetadata(promoter_data)$tss)
  cpg_loc    <- as.vector(GenomicRanges::ranges(rrbs_data)@start)
  tot_reads  <- as.vector(GenomicRanges::elementMetadata(rrbs_data)$total_reads)
  meth_reads <- as.vector(GenomicRanges::elementMetadata(rrbs_data)$meth_reads)

  n            <- 1                         # Data points counter
  LABEL        <- FALSE                     # Flag variable
  meth_data    <- list()                    # List where data will be stored
  prom_counter <- 0                         # Promoter counter
  prom_ind     <- vector(mode = "integer")  # Vector of promoter indices
  cpg_ind      <- vector(mode = "integer")  # Vector of CpG indices
  cpg_ind      <- c(cpg_ind, hits[1,2])     # Add the first subject hit


  message("Generating promoter regions data...")
  for (i in 2:NROW(hits)){
    # If query hits is the same as the previous one
    if (hits[i,1] == hits[(i - 1),1]){
      cpg_ind      <- c(cpg_ind, hits[i,2])  # Add subject hit in the vector
    }else{
      prom_counter <- prom_counter + 1  # Increase promoter counter
      LABEL        <- TRUE
    }

    if(LABEL){
      # Only keep regions that have more than 'n' CpGs
      if (length(cpg_ind) > num_CpG){
        # If the standard deviation of the methylation level is above threshold
        if (sd(meth_reads[cpg_ind] / tot_reads[cpg_ind]) > sd_thresh){
          prom_ind      <- c(prom_ind, prom_loc[prom_counter])
          # Locations of CpGs in the genome
          region        <- cpg_loc[cpg_ind]
          # TSS location for promoter 'promCount'
          tss           <- tss_loc[prom_loc[prom_counter]]
          # Extract strand information, i.e. direction
          strand_direct <- as.matrix(GenomicRanges::strand(promoter_data[prom_loc[prom_counter]])@values)[1,1]
          # Shift CpG locations relative to TSS
          centerd_data  <- center_loc(region = region,
                                      tss = tss,
                                      strand_direction = strand_direct)
          # In the "-" strand the order of the locations should change
          Order <- order(centerd_data)

          meth_data[[n]] <- matrix(0, nrow = length(cpg_ind), ncol = 3)

          # Store normalized locations of methylated CpGs in (fmin, fmax) region
          meth_data[[n]][ ,1] <- minmax_scaling(data = centerd_data[Order],
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
      cpg_ind <- c(cpg_ind, hits[i,2])
    }
  }
  message("Done!")
  meth_regions <- list(meth_data = meth_data, prom_ind = prom_ind)
  class(meth_regions) <- "methRegions"
  return(meth_regions)
}


# Center CpG locations relative to TSS
center_loc <- function(region, tss, strand_direction){
  assertthat::assert_that(is.character(strand_direction))

  center <- region - tss
  if (identical(strand_direction,"-")){
    center  <- (-center)  # If '-' strand, swap CpG locations
  }
  return(center)
}
