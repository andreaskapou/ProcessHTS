context("Reading files")

rnaseq_file <- "rnaseq.bed"

test_that("erros are thrown", {
  expect_error(read.chrom_size("wrong_file_name"))
})

test_that("returns required output", {
  out <- read.rnaseq(rnaseq_file, is_list=TRUE)
  expect_is(out, "list")

  out <- read.rnaseq(rnaseq_file, is_list=FALSE)
  expect_is(out, "GRanges")
})
