## sample preparation - mouse
# EVS 7/2022

require(tidyverse)
require(rstatix)
require(phyloseq)
require(microViz)

## ---- fungi filter ----

## define new filtering parameters from filter-kraken.R
myNewKraken <- function(biom, rank) {
  
  # read biom table
  ps <- biom
  # fix column names
  colnames(tax_table(ps)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  # tax fix
  s1 <- tax_fix(ps)
  # filter out NA phylum
  s2 <- subset_taxa(s1, !is.na(Phylum) & !Phylum %in% c("", "NA"))
  # remove taxa not present in any samples
  s3 <- prune_taxa(taxa_sums(s2) > 0, s2)
  # tax glom at specified rank
  glom <- tax_glom(s3, taxrank = rank)
  
  ## filter by prevalence 
  filt <- glom %>% tax_filter(
    min_prevalence = 0.1,
    min_total_abundance = 0.001)
  
  return(filt)
}

## ---- bacteria filter ----

## make filtering function
myFilter <- function(biom, rank, abun_filt = 1e-5) {
  
  # read biom table
  ps <- biom
  
  # fix column names
  colnames(tax_table(ps)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  
  # tax fix
  s1 <- tax_fix(ps)
  
  # table and number of features for each phyla and kingdom
  rank_names(s1)
  table(tax_table(s1)[, "Phylum"], exclude = NULL) 
  table(tax_table(s1)[, "Kingdom"], exclude = NULL)
  
  # filter out NA phylum
  s2 <- subset_taxa(s1, !is.na(Phylum) & !Phylum %in% c("", "NA"))
  
  # remove taxa not present in any samples
  s3 <- prune_taxa(taxa_sums(s2) > 0, s2)
  
  ## if rank is not "species", tax_glom
  if(rank != "Species") {
    
    cat("\n glomming... this may take a while... \n")
    # tax glom at specified rank
    glom <- tax_glom(s3, taxrank = rank)
    cat("\n done glomming \n")
    
    ## filter by prevalence 
    filt <- glom %>% tax_filter(
      min_prevalence = 0.1,
      min_total_abundance = abun_filt)
    
  } else {
    
    # no glom
    
    # just filter
    filt <- s3 %>% tax_filter(
      min_prevalence = 0.1,
      min_total_abundance = abun_filt)
    
  }
  
  return(filt)
  
}

### ---- FUNGI ----

# kraken
psM <- import_biom("~/OneDrive - The Pennsylvania State University/Shared-Trello-Projects/Fungi Meta-sync issue/metformin-mice/Bioinformatics/bioms/fung_mouse.biom",
                     parseFunction = parse_taxonomy_greengenes)
# fix sample names
sample_names(psM) <- str_remove(sample_names(psM), "_krakenfung")

# make sample data
sample_data(psM) <- data.frame(row.names = sample_names(psM),
                               ID = sapply(str_split(sample_names(psM), "_"), `[`, 1),
                               fakedata = "fakedata")

## add together duplicate samples by making a column in the sample data table to merge by
psmerge <- merge_samples(psM, group = "ID", fun = sum)

# remove fake sample data and add treatment assignment
sample_data(psmerge) <- data.frame(row.names = sample_names(psmerge),
                                   ID = str_remove_all(sample_names(psmerge), "M"))

psmerge <- psmerge %>% 
  ps_mutate(treatment = case_when(
    ID %in% c("1", "2", "3", "4", "5", "6") ~ "placebo",
    ID %in% c("7", "8", "9", "10", "11", "12") ~ "metformin",
    ID %in% "13" ~ "PC_Ext",
    ID %in% "14" ~ "NC_Ext"
  ))
                                   

## filter at Genus level
psf<- myNewKraken(biom = psmerge,
                      rank = "Genus") # 43 taxa

# save
save(psf, file = "R/mouse/mouse-fungi-ps-filtered.RData")

## ---- BACTERIA ----

# kraken
psb <- import_biom("~/OneDrive - The Pennsylvania State University/Shared-Trello-Projects/Fungi Meta-sync issue/metformin-mice/Bioinformatics/bioms/std_mouse.biom",
                   parseFunction = parse_taxonomy_greengenes)
# fix sample names
sample_names(psb) <- str_remove(sample_names(psb), "_krakenst")

# make sample data
sample_data(psb) <- data.frame(row.names = sample_names(psb),
                               ID = sapply(str_split(sample_names(psb), "_"), `[`, 1),
                               fakedata = "fakedata")

## add together duplicate samples by making a column in the sample data table to merge by
psmerge <- merge_samples(psb, group = "ID", fun = sum)

# remove fake sample data and add treatment assignment
sample_data(psmerge) <- data.frame(row.names = sample_names(psmerge),
                                   ID = str_remove_all(sample_names(psmerge), "M"))

psmerge <- psmerge %>% 
  ps_mutate(treatment = case_when(
    ID %in% c("1", "2", "3", "4", "5", "6") ~ "placebo",
    ID %in% c("7", "8", "9", "10", "11", "12") ~ "metformin",
    ID %in% "13" ~ "PC_Ext",
    ID %in% "14" ~ "NC_Ext"
  ))

## raw - 7727 taxa
hist(taxa_sums(psmerge))
hist(sample_sums(psmerge))
## filter at Genus level
psfb<- myFilter(biom = psmerge,
                  rank = "Genus",
                ## abundnace filter is way too low
                abun_filt = 1e-4) 

# save
save(psfb, file = "R/mouse/mouse-bac-ps-filtered.RData")

