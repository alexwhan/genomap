context("converting genotyping scores")

load("../tests/test-data/genotype_raw_df.rda")
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

temp <- genotype_raw_df %>% convert_rel(markerName, parent1, parent2)

test_that("convert_rel converts correctly") {
  expect_error(convert_rel(genotype_raw_df))
}
