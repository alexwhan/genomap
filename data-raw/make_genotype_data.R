library(dplyr)
library(readr)

genotype_raw_df <- data_frame(markerName = paste0("M", 1:6),
                           parent1 = c("AA", "TT", "GG", "CC", "AG", "TC"),
                           parent2 = c("GG", "CC", "AA", "CC", "AA", "TC"),
                           prog1 = c("AA", "TT", "AA", "CC", "AA", "TC"),
                           prog2 = c("AA", "CC", "AA", "CC", "AG", "TC"),
                           prog3 = c("GG", "CC", "GG", "CC", "AG", "TC"),
                           prog4 = c("AG", "TT", "AG", "CC", "AG", "TC"))

genotype_rel_df <- data_frame(markerName = paste0("M", 1:3),
                              prog1 = c("AA", "AA", "BB"),
                              prog2 = c("AA", "BB", "BB"),
                              prog3 = c("BB", "BB", "AA"),
                              prog4 = c("AB", "AA", "AB"))

write_csv(genotype_raw_df, "data-raw/genotype_raw_df.csv")
write_csv(genotype_rel_df, "data-raw/genotype_rel_df.csv")
