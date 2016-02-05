context("Reading files")

rnaseq_file <- "rnaseq.bed"

test_that("erros are thrown", {
  expect_error(read_chrom_size("wrong_file_name"))
})

test_that("returns required output", {
   out <- read_rna_encode_caltech(rnaseq_file, is_GRanges=TRUE)
   expect_is(out, "GRanges")

   out <- read_rna_encode_caltech(rnaseq_file, is_GRanges=FALSE)
   expect_is(out, "data.table")
})
