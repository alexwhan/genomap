context("Testing holmans functions")

load("../test_data/holmans_data.rda")
test_that("trixy converts accurately", {
  expect_error(trixy(mtcars))
  expect_error(trixy(mean))
  expect_warning(trixy(tri_df * 3))
  expect_identical(genomap:::trixy(tri_df), tri_trixy_df)
})

test_that("ggholmans producing plots correctly", {
  expect_error(ggholmans(mtcars))
  expect_warning(ggholmans(tri_df * 3))
  expect_equal(ggholmans(tri_df), holmans_tri)
})
