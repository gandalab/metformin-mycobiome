### fungal and bacterial alpha div
# EVS 2/2023

library(phyloseq)
library(microViz)
library(tidyverse)
library(rstatix)
x
### ---- fungi ----

# get data
load("data/fungi-phyloseq-filtered.RData")
psf <- tax_glom(psf, "Genus") %>% 
  ps_mutate(phenmet = paste(phenotype, metformin, sep = "_")) %>% 
  ps_mutate(Treatment = case_when(
    phenmet %in% "T2D_metformin" ~ "T2D-MET",
    phenmet %in% "T2D_no.metformin" ~ "T2D-NOMET",
    phenmet %in% "healthy_metformin" ~ "NORM-MET",
    phenmet %in% "healthy_no.metformin" ~ "NORM"
  )) %>% 
  ps_mutate(Treatment =  replace_na(Treatment, "Other/Unknown"))
taxa_names(psf) <- psf@tax_table[,'Genus']

# get data
alpha <- estimate_richness(psf, measures = "Shannon") %>% 
  rownames_to_column(var = "ID") %>% 
  left_join(sample_data(psf) %>% data.frame() %>% 
              rownames_to_column(var = "ID") %>% 
              dplyr::select(ID, metformin, phenmet, Treatment, phenotype, paper, Sample))

# get summary stats
alpha %>% get_summary_stats(Shannon, type = "median_iqr")

## save
save(alpha, file = "results/fungi-shannons.RData")

### ---- bacteria -----

load("data/bac-phyloseq-filtered.RData")

psb <- tax_glom(psb, "Genus") %>% 
  ps_mutate(phenmet = paste(phenotype, metformin, sep = "_")) %>% 
  ps_mutate(Treatment = case_when(
    phenmet %in% "T2D_metformin" ~ "T2D-MET",
    phenmet %in% "T2D_no.metformin" ~ "T2D-NOMET",
    phenmet %in% "healthy_metformin" ~ "NORM-MET",
    phenmet %in% "healthy_no.metformin" ~ "NORM"
  )) 
taxa_names(psb) <- psb@tax_table[,'Genus']

# get data
alphab <- estimate_richness(psb, measures = "Shannon") %>% 
  rownames_to_column(var = "ID") %>% 
  left_join(sample_data(psf) %>% data.frame() %>% 
              rownames_to_column(var = "ID") %>% 
              dplyr::select(ID, Treatment, metformin, phenotype, paper, Sample))

## summary stats
alphab %>% get_summary_stats(Shannon, type = "median_iqr")

## save
save(alphab, file = "results/bacteria-shannons.RData")

### ---- correlate ----


both <- alpha %>% 
  mutate(fun.rich = Shannon) %>% 
  select(-Shannon) %>% 
  full_join(alphab %>% mutate(bac.rich = Shannon) %>% select(-Shannon)) %>%
  drop_na(bac.rich)

cor.test(both$fun.rich, both$bac.rich, method = "spearman", exact = FALSE)   # 0.3

ggscatter(both, x = "bac.rich", y = "fun.rich",
          add = "reg.line")
