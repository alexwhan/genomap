library(dplyr)
library(tidyr)
library(purrr)
context("converting genotyping scores")

load("../test_data/genotype_raw_df.rda")
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
})

test_that("convertScore works manually", {
  expect_equal_to_reference({
    manual_convert <- genotype_raw_df %>%
      dplyr::rowwise() %>%
      mutate(new = map(parent1, genomap:::convertScore, paternal = parent2, progeny = c(prog1, prog2))) %>%
      select(markerName, new) %>%
      unnest()
  }, "convertScore_manual.rds")
})
