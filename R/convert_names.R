library(tidyr)
library(ASMap)
library(dplyr)
data(mapDH)
mapDH$pheno <-
old_pheno <- mapDH$pheno
old_pheno$obs <- rnorm(nrow(old_pheno))
gc <- genClones(mapDH)
mapDHf <- fixClones(mapDH, gc$cgd, consensus = TRUE)

rename_pheno <- function(pheno, cgd, .id = "Genotype") {
  new_names <- cgd %>%
    select(G1, G2) %>%
    mutate(new_name = paste0(G1, "_", G2)) %>%
    gather(pos, old_name, G1, G2) %>%
    select(-pos)

  names(new_names)[names(new_names) == "old_name"] <- .id

  if(is.factor(pheno[[.id]])) {
    pheno[[.id]] <- as.character(pheno[[.id]])
  }

  stopifnot(all(new_names$old_name %in% pheno[[.id]]))

  new_pheno <- old_pheno %>%
    left_join(new_names) %>%
    mutate(new_name = ifelse(is.na(new_name), Genotype, new_name))
}

test <- rename_pheno(old_pheno, gc$cgd)
