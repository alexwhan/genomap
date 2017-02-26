library(genomap)
context("get_map_coords")

starts <- c(0, 61, 173, 283, 394, 527)
ends <- c(61, 173, 283, 394, 527, 674)

test_that("lg starts are right from map, cross, mpcross and tidy_gen_map", {
  expect_equal(get_map_coords(bp_map)$lg_start, starts)
  expect_equal(get_map_coords(bp_cross)$lg_start, starts)
  expect_equal(get_map_coords(m4_cross)$lg_start, starts)
  expect_equal(get_map_coords(map2df(bp_map))$lg_start, starts)
})

test_that("lg ends are right from map, cross, mpcross and tidy_gen_map", {
  expect_equal(get_map_coords(bp_map)$lg_end, ends)
  expect_equal(get_map_coords(bp_cross)$lg_end, ends)
  expect_equal(get_map_coords(m4_cross)$lg_end, ends)
  expect_equal(get_map_coords(map2df(bp_map))$lg_end, ends)
})
