#' Create promoter regions from ENCODE Caltech RNA-Seq data
#'
#' \code{create_encode_caltech_prom_reg} creates promoter regions from ENCODE
#'  Caltech RNA-Seq data. Using the TSS of RNA-Seq data as ground truth labels
#'  we create promoter regions \code{N} bp upstream and \code{M} bp
#'  downstream of TSS.
#'
#' @param rna_data data.frame object containing the RNA-Seq data.
#' @param chrom_size Optional data.frame containing the chromosome sizes.
#' @param upstream Integer defining the length of bp upstream of TSS.
#' @param downstream Integer defining the length of bp downstream of TSS.
#'
#' @return The created promoter regions stored in a
#'  \code{\link[GenomicRanges]{GRanges}} object.
#'
#' @seealso \code{\link{read_rna_encode_caltech}}
#'
#' @examples
#' # Load the RNA-Seq example dataset
#' rna_data <- rnaseq_data
#' promoters <- create_encode_caltech_prom_reg(rna_data)
#'
#' @export
create_encode_caltech_prom_reg <- function(rna_data, chrom_size = NULL,
                                           upstream = -100, downstream = 100){

  message("Creating promoter regions ...")
  assertthat::assert_that(is.data.frame(rna_data))
  N <- NROW(rna_data)  # Length of the dataset
  if (upstream > 0 ){
    upstream <- -upstream
  }

  # Create empty vectors
  up_prom   <- vector(mode = "integer", N)  # Start position in chromosome
  down_prom <- vector(mode = "integer", N)  # End position in chromosome
  tss       <- vector(mode = "integer", N)  # TSS position

  for (i in 1:N){
    # Depending on the strand we change regions upstream and downstream of TSS
    if (identical(rna_data$strand[i], "+")){
      up_prom[i] <- max(0, rna_data$start[i] + upstream)
      if (is.null(chrom_size)){
        down_prom[i] <- rna_data$start[i] + downstream
      }else{
        down_prom[i] <- min(chrom_size[rna_data$chr[i], ],
                            rna_data$start[i] + downstream)
      }
      tss[i] <- rna_data$start[i]
    }else if (identical(rna_data$strand[i], "-")){
      up_prom[i] <- max(0, rna_data$end[i] - downstream)
      if (is.null(chrom_size)){
        down_prom[i] <- rna_data$end[i] - upstream
      }else{
        down_prom[i] <- min(chrom_size[rna_data$chr[i], ],
                            rna_data$end[i] - upstream)
      }
      tss[i] <- rna_data$end[i]
    }
  }

  # Create GRanges object
  message("Creating GRanges object for promoter regions ...")
  rna_data <- GenomicRanges::GRanges(seqnames = rna_data$chr,
                     ranges     = IRanges::IRanges(start = up_prom,
                                                   end = down_prom),
                     strand     = rna_data$strand,
                     tss        = tss,
                     ensembl_id = rna_data$ensembl_id,
                     gene_name  = rna_data$gene_name,
                     gene_expr  = rna_data$gene_expr,
                     gene_fpkm  = rna_data$fpkm)
  message("Done!\n")
  return(rna_data)
}
