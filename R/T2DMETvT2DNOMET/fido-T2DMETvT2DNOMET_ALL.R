### FIDO; all data T2DMET v T2DNOMET
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

### get only taxa that are shared between the 6 papers
shared <- psf %>% 
  ps_filter(!paper %in% c("lechatlier2013", "nielsen2014", "forslund2015")) %>% 
  merge_samples(group = "paper") %>% 
  tax_transform("pa") %>% 
  tax_filter(min_prevalence = 6)

# filter original ps by shared taxa
psh <- prune_taxa(taxa_names(shared), psf)

### get only RCTs
psh <- psh %>% 
  ps_mutate(phenmet = paste(phenotype, metformin, sep = "_")) %>% 
  ps_filter(phenmet %in% c("T2D_metformin", "T2D_no.metformin")) %>% 
  # 1/23; dummy variables
  ps_mutate(metformin = if_else(phenmet == "T2D_metformin", 1, 0))
# 466 samples, 36 taxa


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
    paper %in% "wu2017" ~ paste("wu2017", str_extract(row.names(metadata), "(\\d){1,3}V"), sep = "_"),
    # li2014 has sample names in rownames due to merging phyloseqs earlier
    paper %in% "li2014" ~ paste("li2014", row.names(metadata), sep = "_"),
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

## get posterior
posterior <- refit(priors, optim_method="lbfgs", n_samples = 5000)

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
sig.params # 12

# get r2
source("R/r2_emily.R")
rm <- r2_emily(posterior, covariates.include = 2, covariates.exclude = c(3, 4))
quantile(rm, probs = c(0.025, 0.5, 0.975))

## save parameters
pars <- params$Lambda
save(pars, file = "results/T2DMETvNOMET_ALL-params.RData")

### ----- per - paper -----

### This has already been ran in the two separate analyses - just combine
load("results/T2DMETvT2DNOMET_RCT-bypaper-allparams.RData")
rct <- all_par
load("results/T2DMETvT2DNOMET_CS-bypaper-allparams.RData")

#### EXCEPT WU 2017; additional sampling times were added to this global analysis 
#(only wanted short-term for RCT)
# and placebo samples
source("R/myFIDO.R")

# run
wu <- myFIDO(
  phyloseq = psh %>% 
    ps_filter(paper == "wu2017"),
  metadata = meta %>% 
    filter(paper == "wu2017") %>% 
    dplyr::select(metformin, starts_with("wu2017")),
  2
)

## save
all_par <- rbind(all_par, rct %>% filter(!Study == "Wu 2017"), wu %>% mutate(Study = "Wu 2017"))
save(all_par, file = "results/T2DMETvT2DNOMET_ALL-bypaper.RData")

## ----- filter candidate taxa ----

## the significant coords from global comparison
sigs <- sig.params$coord # 3

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

# all three remain