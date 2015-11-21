
# Create RRBS data ------------------------------
rrbs_file <- system.file("extdata", "rrbs.bed", package = "processHTS")
rrbs_data <- read.rrbs(file=rrbs_file)
devtools::use_data(rrbs_data)


# Create RNA-Seq data ------------------------------
rnaseq_file <- system.file("extdata", "rnaseq.bed", package = "processHTS")
rnaseq_data <- read.rnaseq(file=rnaseq_file)
devtools::use_data(rnaseq_data)


# Create hg19 data ------------------------------
hg19_file <- system.file("extdata", "hg19", package = "processHTS")
hg19_data <- read.chrom_size(file=hg19_file)
devtools::use_data(hg19_data)
