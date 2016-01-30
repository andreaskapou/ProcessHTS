#' Plot promoter CpG density
#'
#' Plot a density of number of CpGs that are found in each promoter region.
#'
#' @param obj Object of type 'methExpr'
#'
#' @export
plot_prom_cpg_density <- function(obj){
  prom <- vector(mode = "numeric")
  for (i in 1:length(obj$methyl_region)){
    prom[i] <- length(obj$methyl_region[[i]])
  }
  plot(density(prom),
       col = "darkblue",
       main = "Distribution of # CpGs in promoter regions")
}


#' Plot promoter CpG histogram
#'
#' Plot a histogram of number of CpGs that are found in each promoter region.
#'
#' @param obj Object of type 'methExpr'
#'
#' @export
plot_prom_cpg_hist <- function(obj){
  prom <- vector(mode = "numeric")
  for (i in 1:length(obj$methyl_region)){
    prom[i] <- length(obj$methyl_region[[i]])
  }
  hist(prom,
       border = "darkblue",
       freq = TRUE,
       col = "#CCCCFF",
       xlab = "# of CpGs",
       main = "Histogram of # CpGs in promoter regions")
}
