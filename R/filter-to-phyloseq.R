## filtering kraken fungi tables

# EVS 11/2021, updated 1/2022

# get data
require(dada2)
require(tidyverse)
require(rstatix)
require(phyloseq)
require(microViz)


## ----- define function -----

## define filtering parameters
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


# get metadata
load("all-metadata-objects.RData")


## ---- sun2018 ----

# kraken
psSun <- import_biom("../Bioinformatics/kraken-bioms/fung_sun2018.biom",
                     parseFunction = parse_taxonomy_greengenes)
# fix sample names
sample_names(psSun) <- str_remove(sample_names(psSun), "_krakenfung")

# add metadata
rownames(sun) <- sun$Run
sample_data(psSun) <- sun

# filter
psSunSp <- myNewKraken(biom = psSun, rank = "Species")

## ---- karlsson2013 ----

# kraken biom
psKarl <- import_biom("../Bioinformatics/kraken-bioms/fung_karlsson2013.biom",
                      parseFunction = parse_taxonomy_greengenes)
# fix sample names
sample_names(psKarl) <- str_remove(sample_names(psKarl), "_krakenfung")

# add metadata
rownames(karl) <- karl$Run
sample_data(psKarl) <- karl

### filter
psKarlSp <- myNewKraken(biom = psKarl, rank = "Species")

## ---- elbere2020 ----

# kraken biom
psElb <- import_biom("../Bioinformatics/kraken-bioms/fung_elbere2020.biom",
                     parseFunction = parse_taxonomy_greengenes)
# fix sample names
sample_names(psElb) <- str_remove(sample_names(psElb), "_krakenfung")

# add metadata
elb <- elb %>% column_to_rownames(var = "Run")
sample_data(psElb) <- elb

## filter
psElbSp <- myNewKraken(biom = psElb, rank = "Species")

## ---- wu2017 ----

# kraken biom
psWu <- import_biom("../Bioinformatics/kraken-bioms/fung_wu2017.biom",
                    parseFunction = parse_taxonomy_greengenes)
# fix sample names
sample_names(psWu) <- str_remove(sample_names(psWu), "_krakenfung")

# add metadata
wus <- wu %>% column_to_rownames(var = "Run")
sample_data(psWu) <- wus

## merge by sample name (131 unique samples)
psWum <- merge_samples(psWu, "Sample")

# fix sample data
psWum <- phyloseq(tax_table(psWum),
                  otu_table(psWum, taxa_are_rows = FALSE))
wuc <- wu %>% 
  select(-Run) %>% 
  distinct() %>% 
  column_to_rownames(var = "Sample")
sample_data(psWum) <- wuc

## filter
psWuSp <- myNewKraken(biom = psWum, rank = "Species")


## ---- forslund2015 ----

# get Kraken biom
psFors <- import_biom("../Bioinformatics/kraken-bioms/fung_forslund2015.biom",
                      parseFunction = parse_taxonomy_greengenes)
sample_names(psFors) <- str_remove_all(sample_names(psFors), "_krakenfung")

# get metadata
forsm <- fors %>% column_to_rownames(var = "Run")
sample_data(psFors) <- forsm

# merge by samples (124 unique samples)
psFors <- merge_samples(psFors, 'Sample')
psForsm <- phyloseq(tax_table(psFors),
                    otu_table(psFors, taxa_are_rows = FALSE))

# get new metadata
forsn <- fors %>% 
  select(-Run) %>% 
  distinct() %>% 
  column_to_rownames(var = "Sample")
sample_data(psForsm) <- forsn

## filter
psForSp <- myNewKraken(biom = psForsm, rank = "Species")

## ---- li2014 ----

# kraken
psLi <- import_biom("../Bioinformatics/kraken-bioms/fung_li2014.biom",
                    parseFunction = parse_taxonomy_greengenes)
sample_names(psLi) <- str_remove_all(sample_names(psLi), "_krakenfung")

# get single reads
psLiS <- import_biom("../Bioinformatics/kraken-bioms/fung_SINGLE_li2014.biom",
                     parseFunction = parse_taxonomy_greengenes)
sample_names(psLiS) <- str_remove_all(sample_names(psLiS), "_krakenfung")

# get metadata
lim <- li[li$Run %in% sample_names(psLi),] 
rownames(lim) <- lim$Run
sample_data(psLi) <- lim

# merge by sample (106 unique samples)
psLi <- merge_samples(psLi, "Sample")
psLim <- phyloseq(tax_table(psLi),
                  otu_table(psLi, taxa_are_rows = FALSE))

# merge single by sample
liS <- li[li$Run %in% sample_names(psLiS),]
rownames(liS) <- liS$Run
sample_data(psLiS) <- liS
psLiS <- merge_samples(psLiS, "Sample")

## merge both
psLiB <- merge_phyloseq(psLim, psLiS)

# get new metadata
lin <- li %>% 
  select(-Run) %>% 
  distinct() %>% 
  column_to_rownames(var = "Sample")
sample_data(psLiB) <- lin

## filter
psLiSp <- myNewKraken(biom = psLiB, rank = "Species")

## ---- qin2012 ----

# get biom 
psQin <- import_biom("../Bioinformatics/kraken-bioms/fung_qin2012.biom",
                     parseFunction = parse_taxonomy_greengenes)
sample_names(psQin) <- str_remove_all(sample_names(psQin), "_krakenfung")

# get metadata
qinm <- qin %>% column_to_rownames(var = "Run")
sample_data(psQin) <- qinm

## filter
psQinSp <- myNewKraken(biom = psQin, rank = "Species")


## ---- nielsen2014 ----

# kraken biom
psNiel <- import_biom("../Bioinformatics/kraken-bioms/fung_nielsen2014.biom",
                      parseFunction = parse_taxonomy_greengenes) 
sample_names(psNiel) <- str_remove_all(sample_names(psNiel), "_krakenfung")

# get single reads
psNielS <- import_biom("../Bioinformatics/kraken-bioms/fung_SINGLE_nielsen2014.biom",
                       parseFunction = parse_taxonomy_greengenes)
sample_names(psNielS) <- str_remove_all(sample_names(psNielS), "_krakenfung")

# get metadata for PE
nielP <- niel[niel$Run %in% sample_names(psNiel),]
rownames(nielP) <- nielP$Run
sample_data(psNiel) <- nielP

# merge PE
psNiel <- merge_samples(psNiel, "Sample")
psNielm <- phyloseq(tax_table(psNiel),
                    otu_table(psNiel, taxa_are_rows = FALSE))

# get metdata for SE
nielS <- niel[niel$Run %in% sample_names(psNielS),]
rownames(nielS) <- nielS$Run
sample_data(psNielS) <- nielS

# merge Se
psNielS <- merge_samples(psNielS, "Sample")
psNielSm <- phyloseq(tax_table(psNielS),
                     otu_table(psNielS, taxa_are_rows = FALSE))

# merge PE and SE
psNielB <- merge_phyloseq(psNielm, psNielSm)

# get new metadata
nieln <- niel %>% 
  select(-Run) %>% 
  distinct() %>% 
  column_to_rownames(var = "Sample")
sample_data(psNielB) <- nieln

## filter
psNielSp <- myNewKraken(biom = psNielB, rank = "Species")

## ---- lechatlier2013 ----

# kraken
psLec <- import_biom("../Bioinformatics/kraken-bioms/fung_lechatlier2013.biom",
                     parseFunction = parse_taxonomy_greengenes)
sample_names(psLec) <- str_remove_all(sample_names(psLec), "_krakenfung")

# metadata
lechm <- lech %>% column_to_rownames(var = "Run")
sample_data(psLec) <- lechm

# merge by sample (292 unique samples)
psLec <- merge_samples(psLec, "Sample")
psLecm <- phyloseq(tax_table(psLec),
                   otu_table(psLec, taxa_are_rows = FALSE))

# get new metadata
lechn <- lech %>% 
  select(-Run) %>% 
  distinct() %>% 
  column_to_rownames(var = "Sample")
sample_data(psLecm) <- lechn

## filter
psLecSp <- myNewKraken(biom = psLecm, rank = "Species")


## ---- merge all phyloseqs ----

#### 2/2023; remove matching sample names before merging
library(microViz)
psLecREM <- psLecSp %>% 
  ps_filter(!sample_names(psLecSp) %in% sample_names(psForSp))
psNielREM <- psNielSp %>% 
  ps_filter(!sample_names(psNielSp) %in% sample_names(psForSp))
psNielREM <- psNielREM %>% 
  ps_filter(!sample_names(psNielREM) %in% sample_names(psLecREM))

## 2/2023; update all NA metaHIT ("MH") samples to healthy, no metformin (checked methods)
psLecREM <- psLecREM %>% 
  ps_mutate(metformin = "no.metformin", phenotype = "healthy")
psNielREM <- psNielREM %>% 
  ps_mutate(metformin = if_else(phenotype == "healthy", "no.metformin", metformin))


## remove Forslund - now, all NA samples (unknown phenotypes)
## species level
psf <- merge_phyloseq(psSunSp, psKarlSp, psElbSp, psForSp,
                      psLecREM, psLiSp, psNielREM, psQinSp, psWuSp)
## keep only known phenotype and T2D
psf.filt <- psf %>% 
  ps_filter(metformin %in% c("metformin", "no.metformin"),
            phenotype %in% c("healthy", "T2D", "IGT"))

## which are empty?
empt <- prune_samples(sample_sums(psf) == 0, psf)
sd <- samdat_tbl(empt)
sd %>% group_by(paper) %>% count()
### forslund samples were empty in SRA

# remove empty samples (samples with 0 reads)
psf <- prune_samples(sample_sums(psf.filt) > 0, psf.filt)
psf


## fix duplicate species names
psf@tax_table["460519","Species"] <- "Komagataella.phaffii"
psf@tax_table["113608", "Species"] <- "Tetrapisispora.phaffii"
psf@tax_table["5062", "Species"] <- "Aspergillus.oryzae"
psf@tax_table["318829", "Species"] <- "Pyricularia.oryzae"


# save
save(psf, file = "fungi-phyloseq-filtered.RData")
