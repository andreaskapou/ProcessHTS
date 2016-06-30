#' Create enhancer regions from gene annotation data.
#'
#' \code{create_enh_region} creates enhancer region from gene annotation data.
#'
#' @param annot_data A \code{\link[GenomicRanges]{GRanges}} object containing
#'   the enhancer annotation data.
#' @param chrom_size Optional \code{\link[data.table]{data.table}} containing
#'   chromosome sizes, e.g. using the \code{\link{read_chrom_size}} function.
#' @param upstream Integer defining the length of bp upstream of enhancer
#'   centre.
#' @param downstream Integer defining the length of bp downstream of enhancer
#'   cetre.
#'
#' @return A \code{\link[GenomicRanges]{GRanges}} object containing the enhancer
#'   regions data.
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @examples
#' enh_beatson_file <- system.file("extdata", "enhancers_beatson.bed", package = "processHTS")
#' gene_beatson_file <- system.file("extdata", "rna_beatson2.bed", package = "processHTS")
#' enh_beatson_data <- read_enh_beatson(enh_beatson_file, gene_beatson_file)
#' enh_beatson_data <- enh_beatson_data$enhancer_data
#' enh_regions <- create_enh_region(enh_beatson_data)
#'
#' @importFrom methods is
#' @export
create_enh_region <- function(annot_data, chrom_size = NULL, upstream = -15000,
                               downstream = 15000){

  message("Creating promoter regions ...")
  assertthat::assert_that(methods::is(annot_data, "GRanges"))
  N <- NROW(annot_data)  # Number of enhancers
  if (upstream > 0 ){
    upstream <- -upstream
  }

  # Create empty vectors
  enh       <- vector(mode = "integer", N)  # Enhancer centre location
  up_enh    <- vector(mode = "integer", N)  # Start location in chromosome
  down_enh  <- vector(mode = "integer", N)  # End location in chromosome

  # Extract chromosome information
  annot_chr    <- as.character(annot_data@seqnames)
  # Extract strand information
  annot_strand <- as.character(GenomicRanges::strand(annot_data))
  # Extract start information
  annot_start  <- GenomicRanges::ranges(annot_data)@start
  # Extract end information
  annot_end    <- annot_start + GenomicRanges::ranges(annot_data)@width - 1
  # Centre of enhancer region
  centre <- as.integer( (annot_end - annot_start) / 2) + 1

  for (i in 1:N){
    # Set enhancer centre location
    enh[i] <- annot_start[i] + centre[i]

    # Depending on the strand we change regions
    if (identical(annot_strand[i], "+")){
      # Set upstream bp promoter region
      up_enh[i] <- max(0, enh[i] + upstream)
      # Set downstream bp promoter region
      if (is.null(chrom_size)){
        down_enh[i] <- enh[i] + downstream
      }else{
        down_enh[i] <- min(chrom_size[chrom_size$chr == annot_chr[i]]$size,
                           enh[i] + downstream)
      }
    }else if (identical(annot_strand[i], "-")){
      up_enh[i] <- max(0, enh[i] - downstream)
      # Set upstream bp promoter region
      if (is.null(chrom_size)){
        down_enh[i] <- enh[i] - upstream
      }else{
        down_enh[i] <- min(chrom_size[chrom_size$chr == annot_chr[i]]$size,
                            enh[i] - upstream)
      }
    }
  }

  # Create GRanges object
  message("Creating GRanges object for enhancer regions ...")
  enh_region <- GenomicRanges::GRanges(seqnames = annot_chr,
                              ranges  = IRanges::IRanges(start = up_enh,
                                                         end = down_enh),
                              strand  = annot_strand,
                              tss     = enh)
  message("Finished creating enhancer regions!\n")
  return(enh_region)
}
