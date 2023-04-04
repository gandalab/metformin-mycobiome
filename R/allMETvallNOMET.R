### all MET vs all NOMET 

# EVS 2/2023


library(tidyverse)
library(fido)
library(microViz)
library(phyloseq)
library(ggpubr)

set.seed(123)

# get data
load("data/fungi-phyloseq-filtered.RData")

psf <- psf %>% tax_glom("Genus")
taxa_names(psf) <- tax_table(psf)[,"Genus"]

### get only taxa that are shared between the three papers
shared <- psf %>% 
  ps_filter(!paper %in% c("lechatlier2013", "forslund2015", "nielsen2014")) %>% 
  merge_samples(group = "paper") %>% 
  tax_transform("pa") %>% 
  tax_filter(min_prevalence = 6)

# filter original ps by shared taxa
psh <- prune_taxa(taxa_names(shared), psf)

### get data
psh <- psh %>% 
  ps_filter(!paper %in% c("lechatlier2013", "forslund2015", "nielsen2014")) %>% 
  ps_mutate(phenmet = paste(phenotype, metformin, sep = "_")) %>% 
  # 1/23; dummy variables
  ps_mutate(metformin = if_else(metformin == "metformin", 1, 0)) 
# 845 samples, 36 taxa

#### ----- get paired samples in their own column ----

# first note which samples are paired
## have to do this by paper since each has their own naming conventions

metadata <- sample_data(psh)
metadata <- data.frame(metadata)

# keep only columns that we need
meta.paired <- metadata %>% 
  dplyr::select(Sample, metformin, paper) %>% 
  # create new column of unique names for each paper
  mutate(newID = case_when(
    paper %in% "sun2018" ~ paste("sun2018", str_remove_all(Sample, "\\w+M"), sep = "_"),
    paper %in% "elbere2020" ~ paste("elbere2020", str_extract(Sample, "S(\\d){1,3}"), sep = "_"),
    paper %in% "wu2017" ~ paste("wu2017", str_extract(row.names(metadata), "(\\d){1,3}"), sep = "_"),
    # li2014 has sample names in rownames due to merging phyloseqs earlier
    paper %in% "li2014" ~ paste("li2014", row.names(metadata), sep = "_"),
    # same with forslund2015 and lechatlier 2013
    #paper %in% "forslund2015" ~ paste("forslund2015", row.names(metadata), sep = "_"),
    # paper %in% "lechatlier2013" ~ paste("lechatlier2013", row.names(metadata), sep = "_"),
    # for papers without paired, add the paper name and sample 
    paper %in% c("karlsson2013", "qin2012") ~ paste(paper, Sample, sep = "_")
  ))

# (lazy) do in for loop
samps <- unique(meta.paired$newID)
samps <- samps[!is.na(samps)]
for(i in 1:length(samps)) {
  
  # add column
  new = if_else(meta.paired$newID %in% samps[i], 1, 0)
  
  # add to dataframe as a column
  meta.paired[, ncol(meta.paired) + 1] <- new
  
  # change column name
  colnames(meta.paired)[ncol(meta.paired)] <- samps[i]
  
}

### keep only the columns we want in the model matrix
meta <- meta.paired %>% 
  dplyr::select(metformin, paper, 5:ncol(meta.paired)) 

#### ----- get priors and posterior ----


### get OTU table
otu <- otu_table(psh)
X <- t(model.matrix(~ ., data = meta))
Y <- otu[,colnames(otu) %in% colnames(X)]
# Specify priors
upsilon <- ntaxa(psh)+3 
Omega <- diag(ntaxa(psh))
G <- cbind(diag(ntaxa(psh)-1), -1)
Xi <- (upsilon-ntaxa(psh))*G%*%Omega%*%t(G)
Theta <- matrix(0, ntaxa(psh)-1, nrow(X))
Gamma <- diag(nrow(X))

## prior checking
priors <- pibble(NULL, X, upsilon, Theta, Gamma, Xi)  
#print(priors)

## convert to CLR
priors <- to_clr(priors)  
names_covariates(priors) <- rownames(X)
priors$Y <- Y

## save to run on ROAR
save.image(file = "R/T2DMETvNORM/allMETvallNOMET-priors.RData")

## read in posterior from R
load("R/T2DMETvNORM/fido-posterior-allMETvNOMET.RData")
posterior <- mod

# check posterior
ppc_summary(posterior) # 97%

# get params
posterior <- to_clr(posterior)
names_categories(posterior) <- taxa_names(psh)
params = summary(posterior, pars = "Lambda")

# get significant params
sig.params = params$Lambda %>% 
  filter(covariate == "metformin") %>%
  filter(sign(p2.5) == sign(p97.5))
sig.params # 6

# get r2
source("R/r2_emily.R")
rm <- r2_emily(posterior, covariates.include = 2, covariates.exclude = c(3, 4, 5, 6, 7))
quantile(rm, probs = c(0.025, 0.5, 0.975))

## save parameters
pars <- params$Lambda
save(pars, file = "results/allMETvallNOMET-params.RData")

### ----- per - paper -----

source("R/myFIDO.R")

### SUN 2018
ps <- psh %>% 
  ps_filter(paper == "sun2018")

df <- meta %>% 
  filter(paper == "sun2018") %>% 
  dplyr::select(metformin, starts_with("sun2018"))

sun <- myFIDO(ps, df, 2)

### WU 2017
wu <- myFIDO(psh %>% 
               ps_filter(paper == "wu2017"),
             meta %>% 
               filter(paper == "wu2017") %>% 
               dplyr::select(metformin, starts_with("wu2017")),
             2)

### ELBERE 2020
elb <- myFIDO(psh %>% 
                ps_filter(paper == "elbere2020"),
              meta %>% 
                filter(paper == "elbere2020") %>% 
                dplyr::select(metformin, starts_with("elbere2020")),
              2)

### QIN 2012
qin <- myFIDO(psh %>% 
                ps_filter(paper == "qin2012"),
              meta %>% 
                filter(paper == "qin2012") %>% 
                dplyr::select(metformin, starts_with("qin2012")),
              2)

### LI 2014
li <- myFIDO(psh %>% 
                ps_filter(paper == "li2014"),
              meta %>% 
                filter(paper == "li2014") %>% 
                dplyr::select(metformin, starts_with("li2014")),
              2)

### KARLSSON 2013
karl <- myFIDO(psh %>% 
                ps_filter(paper == "karlsson2013"),
              meta %>% 
                filter(paper == "karlsson2013") %>% 
                dplyr::select(metformin, starts_with("karlsson2013")),
              2)

### combine params and save
all_par <- rbind(sun %>% mutate(Study = "Sun 2018"), 
                 wu %>% mutate(Study = "Wu 2017"), 
                 elb %>% mutate(Study = "Elbere 2020"),
                 qin %>% mutate(Study = "Qin 2012"),
                 li %>% mutate(Study = "Li 2014"),
                 karl %>% mutate(Study = "Karlsson 2013"))

save(all_par, file = "results/alLMETvallNOMET-bypaper-allparams.RData")

## ----- filter candidate taxa ----

## the significant coords from global comparison
sigs <- sig.params$coord

## filter the full dataset to only significant taxa
bypaper <- all_par %>% 
  filter(coord %in% sigs) %>% 
  # get only metformin covariate
  filter(covariate == "metformin") %>% 
  mutate(Genus = str_remove(coord, "clr_"))

### THRESHOLD: there were 6 studies in the comparison,
# so consider significant if 4 have same direction and have magnitude over 0.1


taxon <- unique(bypaper$Genus)
outdf <- data.frame()
for(i in 1:length(taxon)) {
  
  # subset
  sub <- bypaper %>% filter(Genus == taxon[i])
  outdf[i, "Genus"] <- taxon[i]
  
  # summarize effect size direction (sign of mean)
  
  if(length(which(sign(sub$mean) == 1)) >= 4 |
     length(which(sign(sub$mean) == -1)) >= 4) {
    
    # create new column in outdf
    outdf[i, "signs"] <- TRUE
  } else { outdf[i, "signs"] <- FALSE}
  
  
  # summarize strength of effect size
  if(length(which(abs(sub$mean) > 0.1)) >= 4) {
    
    outdf[i, "size"] <- TRUE
  } else { outdf[i, "size"] <- FALSE}
  
}

outdf

## save taxa for new plot
mod1 <- outdf %>% 
  filter(signs == TRUE & size == TRUE) 
