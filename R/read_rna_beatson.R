#' Read file containing Beatson \code{bed} formatted RNA-Seq data
#'
#' \code{read_rna_beatson} reads a file containing promoter annotation
#' data together with gene expression levels from RNA-Seq experiments using the
#' \code{\link{scan}} function.
#'
#' @inheritParams read_bs_encode_haib
#'
#' @return a \code{\link[GenomicRanges]{GRanges}} object if \code{is_GRanges}
#'  is TRUE, otherwise a data.frame object.
#'
#' @seealso \code{\link{read_annot_beatson}}, \code{\link{read_bs_bismark_cov}}
#'
#' @examples
#' # Get the location of the RNA-Seq file
#' rna_beatson_file <- system.file("extdata", "rna_beatson.bed", package = "processHTS")
#' rna_beatson_data <- read_rna_beatson(file=rna_beatson_file, is_GRanges=TRUE)
#'
#' @export
read_rna_beatson <- function(file, chr_discarded = NULL,
                                    is_GRanges = TRUE){
  message("Reading file ", file, " ...")
  rna_data <- read.table(file = file,
                         header = TRUE,
                         sep = "\t",
                         col.names = c("chr", "start", "end", "strand",
                                       "gene_name", "ensembl_id", "gene_fpkm"),
                         comment.char = "",
                         stringsAsFactors = FALSE)


  # Remove selected chromosomes  -------------------------------
  rna_data <- discard_chr(x = rna_data, chr_discarded = chr_discarded)


  # Sorting data -----------------------------------------------
  # With order priority: 1. chr, 2. start, 3. strand
  message("Sorting RNA-Seq data ...")
  rna_data <- rna_data[with(rna_data, order(rna_data$chr,
                                            rna_data$start,
                                            rna_data$strand)), ]

  # Get sequential row numbers
  row.names(rna_data) <- NULL


  if (is_GRanges){
    # Create a GRanges object -----------------------------------
    message("Creating GRanges object ...")
    rna_data <- GenomicRanges::GRanges(seqnames = rna_data$chr,
                                       strand = rna_data$strand,
                                       ranges = IRanges::IRanges(start = rna_data$start,
                                                                 end = rna_data$end),
                                       ensembl_id = rna_data$ensembl_id,
                                       gene_name  = rna_data$gene_name,
                                       gene_fpkm  = rna_data$gene_fpkm)
  }
  message("Done!\n")
  return(rna_data)
}
