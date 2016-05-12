#' Read gene annotation file from Beatson
#'
#' \code{read_annot_beatson} reads a file containing gene annotation data from
#' Beatson, using the \code{\link{scan}} function.
#'
#' @inheritParams read_bs_encode_haib
#'
#' @return A \code{\link[GenomicRanges]{GRanges}} object if \code{is_GRanges} is
#'   TRUE, otherwise a \code{\link[data.table]{data.table}} object.
#'
#'   The GRanges object contains two additional metadata columns: \itemize{
#'   \item \code{ensembl_id}: Ensembl IDs of each gene promoter. \item
#'   \code{gene_name}: Gene name. } These columns can be accessed as follows:
#'   \code{granges_object$ensembl_id}
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @seealso \code{\link{read_chrom_size}}, \code{\link{read_bs_encode_haib}}
#'
#' @examples
#' # Get the location of the annotated Beatson file
#' annot_beatson_file <- system.file("extdata", "annot_beatson.bed", package = "processHTS")
#' annot_beatson_data <- read_annot_beatson(file=annot_beatson_file, is_GRanges=TRUE)
#'
#' @export
read_annot_beatson <- function(file, chr_discarded = NULL, is_GRanges = TRUE){
  message("Reading file ", file, " ...")
  annot_data <- data.table::fread(input = file,
                                  sep = "\t",
                                  header = FALSE,
                                  col.names = c("chr", "gene_name", "strand",
                                                "ensembl_id", "start", "end"),
                                  drop = c(2:3, 7, 10:13))

  # Convert strand direction to '+' and '-'
  annot_data$strand[annot_data[, annot_data$strand == 1]] <- "+"
  annot_data$strand[annot_data[, annot_data$strand == -1]] <- "-"


  # Remove selected chromosomes  -------------------------------
  annot_data <- discard_chr(x = annot_data, chr_discarded = chr_discarded)


  # Sorting data -----------------------------------------------
  # With order priority: 1. chr, 2. start, 3. strand
  message("Sorting annotation data ...")
  annot_data <- annot_data[order(annot_data$chr,
                                 annot_data$start,
                                 annot_data$strand)]


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
  message("Finished reading annotation file!\n")
  return(annot_data)
}
