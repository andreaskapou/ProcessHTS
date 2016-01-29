#' Compute the min-max scaling
#'
#' \code{minmax_scaling} normalizes a given vector using the the min-max
#' scaling method. More formally:
#' \deqn{scaled = \frac{data -x_{min}}{x_{max} - x_{min}} \times (f_{max} - f_{min}) + f_{min}}
#'
#' @param data Vector with numeric data to be scaled.
#' @param xmin Optional minimum value, otherwise \code{min(data)} will be used.
#' @param xmax Optional maximum value, otherwise \code{max(data)} will be used.
#' @param fmin Optional minimum range value, default is 0.
#' @param fmax Optional maximum range value, default is 1.
#'
#' @return The scaled data in the given range, default is between (0, 1). If
#'  xmin = xmax the input vector \code{data} is returned.
#'
#' @examples
#' data <- c(-20, 0, 15, 20)
#' scaled <- minmax_scaling(data)
#'
#' @export
minmax_scaling <- function(data, xmin = NULL, xmax = NULL, fmin = 0, fmax = 1){
  if (is.null(xmin)){
    xmin <- min(data)
  }
  if (is.null(xmax)){
    xmax <- max(data)
  }
  if ( (xmin - xmax) == 0){
    return(data)
  }
  minmax <- (data - xmin) / (xmax - xmin)
  minmax_scaled <- minmax * (fmax - fmin) + fmin

  return(minmax_scaled)
}


#' Extract FPKM from string
#'
#' \code{extract_fpkm} Extracts FPKM value from a string
#'
#' @param x a string containing FPKM information
#'
#' @return The FPKM numeric value
#'
#' @examples
#' data <- 'gene_id "72"; transcr "ENST00000456328"; FPKM "0.0736851531";'
#' scaled <- extract_fpkm(data)
#'
#' @export
extract_fpkm <- function(x){
  # TODO test when no FPKM is available
  fpkm <- gsub(".* FPKM ([^;]+);.*", "\\1", x)
  return(as.numeric(fpkm))
}


#' Extract gene name from string
#'
#' \code{extract_gene_name} Extracts gene name from a string
#'
#' @param x a string containing gene name information
#'
#' @return The gene name as a string
#'
#' @examples
#' data <- 'gene_name "Bnt1.1"; transcr "ENST00000456328"; FPKM "0.0736831";'
#' scaled <- extract_gene_name(data)
#'
#' @export
extract_gene_name <- function(x){
  # TODO test when no gene name is available
  gene_name <- gsub(".* gene_name ([^;]+);.*", "\\1", x)
  return(gene_name)
}


#' Discard selected chromosomes
#'
#' \code{discard_chr} Discards selected chromosomes
#'
#' @param x The HTS data stored in a data frame object
#' @param chr_discarded A vector with chromosome names to be discarded.
#'
#' @return The reduced HTS data.
#'
#' @export
discard_chr <- function(x, chr_discarded = NULL){
  if (!is.null(chr_discarded)){
    message("Removing selected chromosomes ...")
    for (i in 1:length(chr_discarded)){
      x <- subset(x, x$chr != chr_discarded[i])
    }
  }
  return(x)
}


#' Discard BS-Seq noisy reads
#'
#' \code{discard_bs_noise_reads} discards low coverage and (really) high reads
#'  from BS-Seq experiments. These reads can be thought as noise of the
#'  experiment.
#'
#' @param bs_data A Granges object containing the BS-Seq data.
#' @param min_bs_cov The minimum number of reads mapping to each CpG site.
#' @param max_bs_cov The maximum number of reads mapping to each CpG site.
#'
#' @return The reduced Granges object without noisy observations
#'
#' @export
discard_bs_noise_reads <- function(bs_data, min_bs_cov = 2, max_bs_cov = 1000){

  message("Discarding noisy reads ...")
  bs_data <- subset(bs_data,
                    bs_data$total_reads > min_bs_cov &
                    bs_data$total_reads < max_bs_cov)
  return(bs_data)
}
