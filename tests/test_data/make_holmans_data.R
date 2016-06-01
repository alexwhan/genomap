library(dplyr)

set.seed(1981)
tri_df <- data_frame(AA = abs(rnorm(10, mean = 0, sd = 0.15)),
                     BB = abs(rnorm(10, mean = 0, sd = 0.15))) %>%
  rowwise() %>%
  mutate(AB = 1 - sum(AA, BB))

tri_trixy_df <- genomap:::trixy(tri_df)

holmans_tri <- genomap::ggholmans(tri_df)

save(tri_df, tri_trixy_df, holmans_tri, file = "tests/test_data/holmans_data.rda")
