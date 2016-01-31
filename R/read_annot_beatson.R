#' Read file containing gene annotation data from Beatson
#'
#' \code{read_annot_beatson} reads a file containing gene annotation data from
#' Beatson, using the \code{\link{scan}} function.
#'
#' @inheritParams read_bs_encode_haib
#'
#' @return a \code{\link[GenomicRanges]{GRanges}} object if \code{is_GRanges}
#'  is TRUE, otherwise a data.frame object.
#'
#' @seealso \code{\link{read_chrom_size}}, \code{\link{read_bs_encode_haib}}
#'
#' @examples
#' # Get the location of the annotated Beatson file
#' annot_beatson_file <- system.file("extdata", "annot_beatson.bed", package = "processHTS")
#' annot_beatson_data <- read_annot_beatson(file=annot_beatson_file, is_GRanges=TRUE)
#'
#' @export
read_annot_beatson <- function(file, chr_discarded = NULL,
                                        is_GRanges = TRUE){
  message("Reading file ", file, " ...")
  data_raw <- scan(file=file,
                   sep="\t",
                   what=list("character",  # Reference chromosome
                             NULL,         # M kb upstream
                             NULL,         # N kb downstream
                             "character",  # Gene name
                             numeric(),    # Strand direction, 1 or -1
                             "character",  # Ensembl ID
                             NULL,         # Reference chromosome (again)
                             integer(),    # Start position in chromosome
                             integer(),    # End position in chromosome
                             NULL,         # Strand direction (again)
                             NULL,         # Gene name (again)
                             NULL,         # Gene type
                             NULL          # Gene status
                   ))

  # Convert strand direction to '+' and '-'
  strand_dir <- vector(mode = "character", length(data_raw[[5]]))
  strand_dir[which(data_raw[[5]] == 1)] <- "+"
  strand_dir[which(data_raw[[5]] == -1)] <- "-"

  # Store only required fields
  annot_data <- data.frame(chr = data_raw[[1]], start = data_raw[[8]],
                           end = data_raw[[9]], strand = strand_dir,
                           ensembl_id = data_raw[[6]],
                           gene_name = data_raw[[4]],
                           stringsAsFactors = FALSE)
  rm(data_raw)


  # Remove selected chromosomes  -------------------------------
  annot_data <- discard_chr(x = annot_data, chr_discarded = chr_discarded)


  # Sorting data -----------------------------------------------
  # With order priority: 1. chr, 2. start, 3. strand
  message("Sorting annotation data ...")
  annot_data <- annot_data[with(annot_data, order(annot_data$chr,
                                                  annot_data$start,
                                                  annot_data$strand)), ]
  # Get sequential row numbers
  row.names(annot_data) <- NULL


  if (is_GRanges){
    # Create a GRanges object -----------------------------------
    message("Creating GRanges object ...")
    annot_data <- GenomicRanges::GRanges(seqnames = annot_data$chr,
                         strand = annot_data$strand,
                         ranges = IRanges::IRanges(start = annot_data$start,
                                                   end = annot_data$end),
                         ensembl_id = annot_data$ensembl_id,
                         gene_name  = annot_data$gene_name)
  }
  message("Done!\n")
  return(annot_data)
}
