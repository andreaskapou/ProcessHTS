#' Read Beatson bed formatted RNA-Seq file
#'
#' \code{read_rna_beatson} reads a file containing promoter annotation data
#' together with gene expression levels from RNA-Seq experiments using the
#' \code{\link[data.table]{fread}} function.
#'
#' @inheritParams read_bs_encode_haib
#'
#' @return A \code{\link[GenomicRanges]{GRanges}} object if \code{is_GRanges} is
#'   TRUE, otherwise a \code{\link[data.table]{data.table}} object.
#'
#'   The GRanges object contains three additional metadata columns: \itemize{
#'   \item \code{ensembl_id}: Ensembl IDs of each gene promoter. \item
#'   \code{gene_name}: Gene name. \item \code{gene_fpkm}: Expression level in
#'   FPKM. } These columns can be accessed as follows:
#'   \code{granges_object$ensembl_id}
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @seealso \code{\link{read_annot_beatson}}, \code{\link{read_bs_bismark_cov}}
#'
#' @examples
#' # Get the location of the RNA-Seq file
#' rna_beatson_file <- system.file("extdata", "rna_beatson.bed", package = "processHTS")
#' rna_beatson_data <- read_rna_beatson(file=rna_beatson_file, chr_discarded = "chrX", is_GRanges=TRUE)
#'
#' @export
read_rna_beatson <- function(file, chr_discarded = NULL, is_GRanges = TRUE){
  message("Reading file ", file, " ...")
  rna_data <- data.table::fread(input = file,
                                sep = "\t",
                                header = TRUE,
                                col.names = c("chr", "start", "end", "strand",
                                              "gene_name", "ensembl_id",
                                              "gene_fpkm"))


  # Remove selected chromosomes  -------------------------------
  rna_data <- discard_chr(x = rna_data, chr_discarded = chr_discarded)


  # Sorting data -----------------------------------------------
  # With order priority: 1. chr, 2. start, 3. strand
  message("Sorting RNA-Seq data ...")
  rna_data <- rna_data[order(rna_data$chr, rna_data$start, rna_data$strand)]


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
  message("Finished reading RNA-Seq file!\n")
  return(rna_data)
}
