library(dplyr)
library(readr)

genotype_raw_df <- data_frame(markerName = paste0("M", 1:6),
                              parent1 = c("AA", "TT", "GG", "CC", "AG", "TC"),
                              parent2 = c("GG", "CC", "AA", "CC", "AA", "TC"),
                              prog1 = c("AA", "TT", "AA", "CC", "AA", "TC"),
                              prog2 = c("AA", "CC", "AA", "CC", "AG", "TC"),
                              prog3 = c("GG", "CC", "GG", "CC", "AG", "TC"),
                              prog4 = c("AG", "TT", "AG", "CC", "AG", "TC"))

saveRDS(genotype_raw_df, file = "tests/test_data/genotype_raw_df.rds")

genotype_rel_df <- data_frame(markerName = paste0("M", 1:6),
                              prog1 = c("AA", "AA", "BB", rep(NA, 3)),
                              prog2 = c("AA", "BB", "BB", rep(NA, 3)),
                              prog3 = c("BB", "BB", "AA", rep(NA, 3)),
                              prog4 = c("AB", "AA", "AB", rep(NA, 3)))

saveRDS(genotype_rel_df, file = "tests/test_data/genotype_rel_df.rds")

genotype_parent_status_df <- genotype_raw_df %>%
  mutate(parent_status = c("informative", "informative", "informative",
                            "nonseg", "heterozygous", "heterozygous")) %>%
  select(parent1, parent2, parent_status, markerName, contains("prog")) %>%
  rowwise() %>%
  ungroup()

saveRDS(genotype_parent_status_df, file = "tests/test_data/genotype_parent_status_df.rds")
