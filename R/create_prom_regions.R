#' Create promoter regions using RNA-Seq data
#'
#' \code{create_prom_regions} creates promoter regions data using Transcription
#' Start Sites (TSS) of RNA-Seq data as ground truth labels. Using the TSS
#' data we create promoter regions \code{n} bps upstream and \code{m} bps
#' downstream of TSS. The RNA-Seq data should be a \code{\link{list}} object.
#'
#' @param data List object containing the RNA-Seq data.
#' @param chrom_size Optional data frame containing the chromosome sizes.
#' @param upstream Integer defining the length of bps upstream of TSS.
#' @param downstream Integer defining the length of bps downstream of TSS.
#'
#' @return The created promoters stored in a \code{\link[GenomicRanges]{GRanges}} object.
#'
#' @seealso \code{\link{read.rnaseq}}
#'
#' @examples
#' # Load the RNA-Seq example dataset
#' data <- rnaseq_data
#' promoters <- create_prom_regions(data)
#'
#' @export
create_prom_regions <- function(data, chrom_size = NULL, upstream = -100,
                                                       downstream = 100){
  assertthat::assert_that(is.list(data))
  N <- length(data[[1]])  # Length of the dataset
  if (upstream > 0 ){
    upstream <- -upstream
  }

  # Create list object for keeping unique promoter regions
  prom <- list()
  prom[[1]] <- data[[1]]  # Reference chromosome
  prom[[2]] <- vector(mode = "integer", N)  # Start position in chromosome
  prom[[3]] <- vector(mode = "integer", N)  # End position in chromosome
  prom[[4]] <- data[[4]]  # Gene ENSEMBL id
  prom[[5]] <- vector(mode = "integer", N)  # TSS position
  prom[[6]] <- data[[6]]  # Strand: + or - or . for unknown
  prom[[7]] <- data[[5]]  # Expression level
  prom[[8]] <- data[[10]]  # Expression level in FPKM

  for (i in 1:N){
    # Depending on strand we change the regions upstream and downstream of TSS
    if (identical(data[[6]][i], "+")){
      prom[[2]][i] <- max(0, data[[2]][i] + upstream)
      if (is.null(chrom_size)){
        prom[[3]][i] <- data[[2]][i] + downstream
      }else{
        prom[[3]][i] <- min(chrom_size[data[[1]][i], ], data[[2]][i] + downstream)
      }
      prom[[5]][i] <- data[[2]][i]
    }else if (identical(data[[6]][i], "-")){
      prom[[2]][i] <- max(0, data[[3]][i] - downstream)
      if (is.null(chrom_size)){
        prom[[3]][i] <- data[[3]][i] - upstream
      }else{
        prom[[3]][i] <- min(chrom_size[data[[1]][i], ], data[[3]][i] - upstream)
      }
      prom[[5]][i] <- data[[3]][i]
    }
  }
  rm(data) # Remove objects from workspace

  # Create a GRanges object
  promoter_data <- GenomicRanges::GRanges(seqnames = prom[[1]],
                      ranges    = IRanges::IRanges(start = prom[[2]], end = prom[[3]]),
                      strand    = prom[[6]],
                      gene_id   = prom[[4]],
                      tss       = prom[[5]],
                      gene_expr = prom[[7]],
                      gene_fpkm = prom[[8]]
                )
  message("Created GRanges object for promoter regions...")
  return(promoter_data)
}
