library(dplyr)

set.seed(1981)
tri_df <- data_frame(AA = abs(rnorm(10, mean = 0, sd = 0.15)),
                     BB = abs(rnorm(10, mean = 0, sd = 0.15))) %>%
  rowwise() %>%
  mutate(AB = 1 - sum(AA, BB))

holmans_tri <- genomap::ggholmans(tri_df)
devtools::use_data(tri_df)
devtools::use_data(holmans_tri, internal = TRUE)
