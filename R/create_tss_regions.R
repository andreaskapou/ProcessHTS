#' Create TSS regions using RNA-Seq data
#'
#' \code{create_tss_regions} creates promoter regions data using Transcription
#' Start Sites (TSS) of RNA-Seq data as ground truth labels. Using the TSS
#' data we create promoter regions \code{n} bps upstream and \code{m} bps
#' downstream of TSS. The RNA-Seq data should be a \code{\link{list}} object.
#'
#' @inheritParams create_prom_regions
#'
#' @return The created promoters stored in a \code{\link[GenomicRanges]{GRanges}} object.
#'
#' @seealso \code{\link{read.rnaseq}}
#'
#' @examples
#' # Load the RNA-Seq example dataset
#' data <- rnaseq_data
#' promoters <- create_tss_regions(data)
#'
#' @export
create_tss_regions <- function(data, chrom_size = NULL, upstream = -100,
                               downstream = 100){
  assertthat::assert_that(is.list(data))
  N <- length(data[[1]])  # Length of the dataset
  if (upstream > 0 ){
    upstream <- -upstream
  }

  # Create list object for keeping unique TSS regions
  tss <- list()
  tss[[1]] <- vector(mode = "character") #data[[1]]  # Reference chromosome
  tss[[2]] <- vector(mode = "integer") #, N)  # Start position in chromosome
  tss[[3]] <- vector(mode = "integer") # N)  # End position in chromosome
  tss[[4]] <- vector(mode = "character") # data[[4]]  # Gene ENSEMBL id
  tss[[5]] <- vector(mode = "integer") # N)  # TSS position
  tss[[6]] <- vector(mode = "character") # data[[6]]  # Strand: + or - or . for unknown
  tss[[7]] <- vector(mode = "numeric") # data[[5]]  # Expression level
  tss[[8]] <- vector(mode = "numeric") # data[[10]]  # Expression level in FPKM

  n            <- 1                         # Data points counter
  LABEL        <- FALSE                     # Flag variable
  tss_counter  <- 0                         # TSS counter
  tss_ind      <- vector(mode = "integer")  # Vector of tss indices
  tss_ind      <- c(tss_ind, 1)             # Add the first index

  for (i in 2:N){
    if (data[[2]][i] > (data[[2]][i - 1] + downstream)){
      tss_counter  <- tss_counter + 1
      LABEL        <- TRUE
    }else{
      tss_ind     <- c(tss_ind, i)
    }

    if (LABEL){
      expr <- vector(mode = "numeric")
      for (j in 1:length(tss_ind)){
        expr[j] <- data[[5]][tss_ind[j]]
      }
      position <- which.max(expr)
      position <- tss_ind[position]
      #max_expr <- max(expr)

      tss[[1]][tss_counter] <- data[[1]][position]

      #if (identical(data[[6]][i], "+")){
      #  tss[[2]][tss_counter] <- max(0, data[[2]][position] + upstream)
      #  if (is.null(chrom_size)){
      #    tss[[3]][tss_counter] <- data[[2]][position] + downstream
      #  }else{
      #    tss[[3]][tss_counter] <- min(chrom_size[data[[1]][i], ], data[[2]][position] + downstream)
      #  }
      #}else if (identical(data[[6]][i], "-")){
      #  tss[[2]][tss_counter] <- max(0, data[[2]][position] - downstream)
      #  if (is.null(chrom_size)){
      #    tss[[3]][tss_counter] <- data[[2]][position] - upstream
      #  }else
      #    tss[[3]][tss_counter] <- min(chrom_size[data[[1]][i], ], data[[2]][position] - upstream)
      #}

      tss[[2]][tss_counter] <- data[[2]][position] + upstream
      tss[[3]][tss_counter] <- data[[2]][position] + downstream
      tss[[4]][tss_counter] <- data[[4]][position]
      tss[[5]][tss_counter] <- data[[2]][position]
      tss[[6]][tss_counter] <- data[[6]][position]
      tss[[7]][tss_counter] <- data[[5]][position]
      tss[[8]][tss_counter] <- data[[10]][position]

      LABEL   <- FALSE
      tss_ind <- vector(mode = "integer")  # Vector of tss indices
      tss_ind <- c(tss_ind, i)             # Add the first index
    }
  }

  rm(data) # Remove objects from workspace

  # Create a GRanges object
  promoter_data <- GenomicRanges::GRanges(seqnames = tss[[1]],
                                          ranges    = IRanges::IRanges(start = tss[[2]], end = tss[[3]]),
                                          strand    = tss[[6]],
                                          gene_id   = tss[[4]],
                                          tss       = tss[[5]],
                                          gene_expr = tss[[7]],
                                          gene_fpkm = tss[[8]]
  )
  message("Created GRanges object for promoter regions...")
  return(promoter_data)
}
