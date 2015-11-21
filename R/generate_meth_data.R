#' Generate methylation data

meth_data <- function(rrbs_data, promoter_data, upstream, downstream, num_CpG=0, sd_thresh=0, fmin=-1, fmax=1){
  assert_that(is(rrbs_data, "GRanges"))
  assert_that(is(promoter_data, "GRanges"))
  if (upstream > 0 ){
    upstream < -upstream
  }

  # Find overlaps between promoter regions and RRBS data
  overlaps <- GenomicRanges::findOverlaps(query   = promoter_data,
                                          subject = rrbs_data,
                                          ignore.strand = FALSE)
  message("Converting Granges objects in matrices...")
  hits <- as.matrix(overlaps)   # Create a HITS matrix for faster computations
  prom_loc <- unique(hits[,1])  # Keep promoter locations

  # Convert data in matrix format for faster lookup
  tss_loc    <- as.vector(elementMetadata(promoter_data)$tss)
  cpg_loc    <- as.vector(ranges(rrbs_data)@start)
  tot_reads  <- as.vector(elementMetadata(rrbs_data)$total_reads)
  meth_reads <- as.vector(elementMetadata(rrbs_data)$meth_reads)

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
          strand_direct <- as.matrix(strand(promoter_data[prom_loc[prom_counter]])@values)[1,1]
          # Shift CpG locations relative to TSS
          centerd_data  <- center_loc(region = region,
                                      tss = tss,
                                      strand_direction = strand_direct)
          # In the "-" strand the order of the locations should change
          Order <- order(centerd_data)

          meth_data[[n]] <- matrix(0, nrow = 3, ncol = length(cpg_ind))

          # Store normalized locations of methylated CpGs in (fmin, fmax) region
          meth_data[[n]][1, ] <- minmax_scaling(data = centerd_data[Order],
                                                xmin = upstream,
                                                xmax = downstream,
                                                fmin = fmin,
                                                fmax = fmax)

          # Store total reads in the corresponding locations
          meth_data[[n]][2, ] <- tot_reads[cpg_ind][Order]
          # Store methylated reads in the corresponding locations
          meth_data[[n]][3, ] <- meth_reads[cpg_ind][Order]

          # Increase data points counter
          n <- n + 1
        }
      }
      LABEL   <- FALSE
      cpg_ind <- vector(mode="integer")
      cpg_ind <- c(cpg_ind, hits[i,2])
    }
  }
}


# Center CpG locations relative to TSS
center_loc <- function(region, tss, strand_direction){
  assert_that(is.character(strand_direction))

  center <- region - tss
  if (identical(strand_direction,"-")){
    center  <- (-center)  # If '-' strand, swap CpG locations
  }

  return(center)
}
