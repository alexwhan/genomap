library(genomap)
context("get_lg_offsets")

offsets <- c(0, 61, 173, 283, 394, 527)

test_that("lg offsets are right from map, cross, mpcross and tidy_gen_map", {
  expect_equal(get_lg_offsets(bp_map)$lg_offset, offsets)
  expect_equal(get_lg_offsets(bp_cross)$lg_offset, offsets)
  expect_equal(get_lg_offsets(m4_cross)$lg_offset, offsets)
  expect_equal(get_lg_offsets(map2df(bp_map))$lg_offset, offsets)
})

