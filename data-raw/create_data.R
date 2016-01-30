
# Create RRBS data --------------------------------
rrbs_file <- system.file("extdata", "rrbs.bed", package = "processHTS")
rrbs_data <- read_bs_encode_haib(file=rrbs_file, is_GRanges = TRUE)
devtools::use_data(rrbs_data, overwrite = TRUE)


# Create RNA-Seq data -----------------------------
rnaseq_file <- system.file("extdata", "rnaseq.bed", package = "processHTS")
rnaseq_data <- read_rna_encode_caltech(file=rnaseq_file, is_GRanges = FALSE)
devtools::use_data(rnaseq_data, overwrite = TRUE)

# Create annot data in GRanges format -----------
rnaseq_file <- system.file("extdata", "rnaseq.bed", package = "processHTS")
annot_data <- read_rna_encode_caltech(file=rnaseq_file, is_GRanges = TRUE)
devtools::use_data(annot_data, overwrite = TRUE)


# Create hg19 data --------------------------------
hg19_file <- system.file("extdata", "hg19", package = "processHTS")
hg19_data <- read_chrom_size(file=hg19_file)
devtools::use_data(hg19_data, overwrite = TRUE)

# Create bismark data -----------------------------
bismark_file <- system.file("extdata", "bism_rep1.bed", package = "processHTS")
bs_data      <- read_bs_bismark_cov(file=bismark_file, is_GRanges=TRUE)
devtools::use_data(bs_data, overwrite = TRUE)
