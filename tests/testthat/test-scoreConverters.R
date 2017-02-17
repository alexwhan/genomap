library(dplyr)
library(tidyr)
library(purrr)
context("converting genotyping scores")

load("../test_data/genotype_converter_data.rda")
test_that("isHet works", {
  expect_false(isHet("AA"))
  expect_true(isHet("AG"))
  expect_error(isHet('u'))
  expect_error(isHet("GT", "[AB]", "error"))
  expect_error(isHet("GT", "[AB]", "other"))
  expect_warning(isHet(12))
})

test_that("isHomo works", {
  expect_false(isHomo("AG"))
  expect_true(isHomo("AA"))
  expect_error(isHomo('u'))
  expect_error(isHomo("TT", "[AB]", "error"))
  expect_error(isHomo("GT", "[AB]", "other"))
  expect_warning(isHomo(11))
})

test_that("convert_rel converts correctly", {
  expect_error(convert_rel(genotype_raw_df))
  expect_error(convert_rel(genotype_raw_df, markerName, parent1, parent2, contains("prog"), missingString = "--"))
})

test_that("convertScore works manually", {
  expect_equal_to_reference({
    suppressWarnings(manual_convert <- genotype_raw_df %>%
      convert_rel(markerName, parent1, parent2, dplyr::contains("prog")) %>%
      dplyr::select(markerName, converted) %>%
      tidyr::unnest())
  }, "../test_data/genotype_rel_df.rds")
})

test_that("check_parents_ works", {
  expect_equal_to_reference({
    genotype_raw_df %>%
      check_parents_('parent1', 'parent2') %>%
      ungroup()
  }, "../test_data/genotype_parent_status_df.rds")
})
