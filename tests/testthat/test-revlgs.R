context("reversing linkage groups")

test_that("revvec reverses vectors correctly", {
  expect_equal(revvec(c(a = 1, b = 2, c = 5)), c(c = 1, b = 4, a = 5))
})
