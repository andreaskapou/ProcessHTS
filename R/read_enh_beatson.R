#' Read Beatson bed formatted enhancer file
#'
#' \code{read_enh_beatson} reads a file containing enhancer annotation data
#' using the \code{\link[data.table]{fread}} function.
#'
#' @param enh_file File containing enhancer annotation data
#' @param gene_file File containing gene annotation data
#' @inheritParams read_bs_encode_haib
#'
#' @return A \code{\link[GenomicRanges]{GRanges}} object if \code{is_GRanges} is
#'   TRUE, otherwise a \code{\link[data.table]{data.table}} object.
#'
#'   The GRanges object contains three additional metadata columns: \itemize{
#'   \item \code{ensembl_id}: Ensembl IDs of each gene promoter. \item
#'   \code{gene_name}: Gene name.} These columns can be accessed as follows:
#'   \code{granges_object$ensembl_id}
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @seealso \code{\link{read_annot_beatson}}, \code{\link{read_bs_bismark_cov}}
#'
#' @examples
#' # Get the location of the RNA-Seq file
#' enh_beatson_file <- system.file("extdata", "enhancers_beatson.bed", package = "processHTS")
#' rna_beatson_file <- system.file("extdata", "rna_beatson2.bed", package = "processHTS")
#'
#' enh_beatson_data <- read_enh_beatson(enh_beatson_file, rna_beatson_file, chr_discarded = "chrX", is_GRanges=TRUE)
#'
#' @export
read_enh_beatson <- function(enh_file, gene_file, chr_discarded = NULL, is_GRanges = TRUE){
  message("Reading file ", enh_file, " ...")
  enhancer_data <- data.table::fread(input = enh_file,
                               sep = "\t",
                               header = TRUE,
                               col.names = c("chr", "start", "end", "strand",
                                             "gene_name", "ensembl_id"))


  # Remove duplicate gene names
  dupl <- which(duplicated(enhancer_data$gene_name))
  enhancer_data <- enhancer_data[-dupl, ]
  # Remove duplicate enhancer regions
  dupl <- which(duplicated(enhancer_data$start))
  enhancer_data <- enhancer_data[-dupl, ]
  rownames(enhancer_data) <- NULL


  # Read RNA-Seq BED file and return data.table format
  gene_data <- read_rna_beatson(file          = gene_file,
                                chr_discarded = chr_discarded,
                                is_GRanges    = FALSE)

  # Remove duplicate gene names
  gene_data <- gene_data[-which(duplicated(gene_data$gene_name)),]
  # Keep only genes that overlap with enhancers
  gene_data <- gene_data[gene_data$gene_name %in% enhancer_data$gene_name]
  rownames(gene_data) <- NULL

  # Get strand information in order to create the appropriate methylation profiles
  enhancer_data <- enhancer_data[enhancer_data$gene_name %in% gene_data$gene_name]
  enhancer_data$strand <- gene_data$strand
  rownames(enhancer_data) <- NULL

  # Remove selected chromosomes  -------------------------------
  enhancer_data <- discard_chr(x = enhancer_data, chr_discarded = chr_discarded)


  # Sorting data -----------------------------------------------
  # With order priority: 1. chr, 2. start, 3. strand
  message("Sorting Enhancer data ...")
  enhancer_data <- enhancer_data[order(enhancer_data$chr,
                                       enhancer_data$start,
                                       enhancer_data$strand)]

  if (is_GRanges){
    # Create GRanges object
    message("Creating GRanges object for enhancer data ...")
    enhancer_data <- GenomicRanges::GRanges(seqnames = enhancer_data$chr,
                 ranges     = IRanges::IRanges(start = enhancer_data$start,
                                               end   = enhancer_data$end),
                 strand     = enhancer_data$strand,
                 ensembl_id = enhancer_data$ensembl_id,
                 gene_name  = enhancer_data$gene_name)

    # Create GRanges object also for RNA-Seq data
    gene_data <- GenomicRanges::GRanges(seqnames = gene_data$chr,
                                        strand = gene_data$strand,
                                        ranges = IRanges::IRanges(start = gene_data$start,
                                                                  end = gene_data$end),
                                        ensembl_id = gene_data$ensembl_id,
                                        gene_name  = gene_data$gene_name,
                                        gene_fpkm  = gene_data$gene_fpkm)
  }
  message("Finished reading Enhancer file!\n")
  return(list(enhancer_data = enhancer_data, gene_data = gene_data))
}
