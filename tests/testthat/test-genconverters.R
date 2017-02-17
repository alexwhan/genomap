context("Testing genoconverter functions")

load("../test_data/map_and_cross.rda")

test_that("map2df converts correctly", {
  expect_equal(map2df(test_map),
               test_mapdf)
  expect_equal(map2df(test_cross),
               test_crossdf)
  expect_error(map2df(mtcars))
  expect_silent(map2df(test_map))
  expect_silent(map2df(test_cross))
})
