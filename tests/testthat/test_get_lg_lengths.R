library(genomap)
context("get_lg_lengths")

lengths <- c(61, 112, 110, 111, 133, 147)

test_that("lg lengths are right from map, cross, mpcross and tidy_gen_map", {
  expect_equal(get_lg_lengths(bp_map)$lg_length, lengths)
  expect_equal(get_lg_lengths(bp_cross)$lg_length, lengths)
  expect_equal(get_lg_lengths(m4_cross)$lg_length, lengths)
  expect_equal(get_lg_lengths(map2df(bp_map))$lg_length, lengths)
})

