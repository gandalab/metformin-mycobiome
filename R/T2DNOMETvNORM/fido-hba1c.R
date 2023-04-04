### FIDO - covariates - by paper

## HbA1c

# EVS  2/2023

library(fido)
library(tidyverse)
library(phyloseq)
library(ggpubr)
library(microViz)
library(rstatix)

# get myFIDO function
source("R/myFIDO.R")

## get OTU table
load("data/fungi-phyloseq-filtered.RData")

psf <- psf %>% tax_glom("Genus")
taxa_names(psf) <- tax_table(psf)[,"Genus"]

# how many papers have insulin metrics?
psf %>% samdat_tbl() %>% 
  group_by(paper) %>% 
  get_summary_stats(hba1c.perc, type = "mean_sd")
# forslund, karlsson, lechatlier, li, qin, sun, elbere2020

### get only taxa that are shared between all papers
shared <- psf %>% ps_filter(!paper %in% c("wu2017", "nielsen2014")) %>% merge_samples(group = "paper") %>% 
  tax_transform("pa") %>% 
  tax_filter(min_prevalence = 7)

# filter original ps by shared taxa
psh <- prune_taxa(taxa_names(shared), psf)

## ---- global test ----

# NO PAIRED SAMPLES - use only baseline values from RCTs (only provided one set of clinical vars - not long enough to see differences in BMI)
ps <- psh

metadata <- sample_data(psh)
df <- data.frame(metadata) %>% 
  drop_na(hba1c.perc) %>% 
  # remove post-metformin samples from RCTs (sun and elbere)
  mutate(papermet = paste0(paper, metformin)) %>% 
  filter(!papermet %in% c("sun2018metformin")) %>% 
  filter(!papermet %in% c("elbere2020metformin")) %>% 
  # standardize
  
  mutate(hba.stand = (hba1c.perc - mean(hba1c.perc) ) / sd(hba1c.perc)) %>% 
  mutate(metformin = if_else(metformin == "metformin", 1, 0))

# get summary stats
df %>% get_summary_stats(hba1c.perc, type = "common")

# get otu tabl
otu <- otu_table(ps)

## make model matrix
X <- t(model.matrix(~ hba.stand + metformin + paper, data = df))

# set Y
Y <- otu[,colnames(otu) %in% colnames(X)]

# Specify priors
upsilon <- ntaxa(ps)+3 
Omega <- diag(ntaxa(ps))
G <- cbind(diag(ntaxa(ps)-1), -1)
Xi <- (upsilon-ntaxa(ps))*G%*%Omega%*%t(G)
Theta <- matrix(0, ntaxa(ps)-1, nrow(X))
Gamma <- diag(nrow(X))

## prior checking
priors <- pibble(NULL, X, upsilon, Theta, Gamma, Xi)  
#print(priors)

## convert to CLR
priors <- to_clr(priors)  

# set names
names_covariates(priors) <- rownames(X)

priors$Y <- Y

# refit with priors to get posterior
posterior <- refit(priors, optim_method="lbfgs", n_samples = 5000)

# convert to CLR
posterior <- to_clr(posterior)

# name the posterior categories (taxa)
names_categories(posterior) <- taxa_names(ps)
params = summary(posterior, pars = "Lambda")

# get sigs
sig.params <- params$Lambda %>% 
  filter(covariate == "hba.stand") %>% 
  filter(sign(p2.5) == sign(p97.5)) # 9

# save results
hbapar <- params
save(hbapar, file = "results/fido-hba1c-params.RData")

# get r2
quantile(r2(posterior, covariate = 2),
         probs = c(0.025, 0.5, 0.975)
)


## ----- BY PAPER ----

### FORSLUND 2015
### NO PAIRED SAMPLES

## subset metadata and OTU table
ps <- psf %>% 
  ps_filter(paper == "forslund2015") 

df <- data.frame(metadata) %>% 
  filter(paper == "forslund2015") %>% 
  # drop NAs
  drop_na(hba1c.perc) %>% 
  # standardize fbg
  mutate(hba.stand = (hba1c.perc - mean(hba1c.perc) ) / sd(hba1c.perc)) %>% 
  dplyr::select(hba.stand)

## run 
fors <- myFIDO(ps, df, covariates.include = 2)

### KARLSSON

## subset metadata and OTU table
ps <- psf %>% 
  ps_filter(paper == "karlsson2013") 

df <- data.frame(metadata) %>% 
  filter(paper == "karlsson2013") %>% 
  drop_na(hba1c.perc) %>% 
  # standardize fbg
  mutate(hba.stand = (hba1c.perc - mean(hba1c.perc) ) / sd(hba1c.perc)) %>% 
  dplyr::select(metformin, hba.stand)

karl <- myFIDO(ps, df, get_r2 = FALSE)

### LECHATLIER 

ps <- psf %>% 
  ps_filter(paper == "lechatlier2013") 

df <- data.frame(metadata) %>% 
  filter(paper == "lechatlier2013") %>% 
  drop_na(hba1c.perc) %>% 
  # standardize fbg
  mutate(hba.stand = (hba1c.perc - mean(hba1c.perc) ) / sd(hba1c.perc)) %>% 
  dplyr::select(hba.stand)

lech <- myFIDO(ps, df, get_r2 = FALSE)

### LI 2014

ps <- psf %>% 
  ps_filter(paper == "li2014") 

df <- data.frame(metadata) %>% 
  filter(paper == "li2014") %>% 
  drop_na(hba1c.perc) %>% 
  # standardize fbg
  mutate(hba.stand = (hba1c.perc - mean(hba1c.perc) ) / sd(hba1c.perc)) %>% 
  dplyr::select(hba.stand, metformin)

li <- myFIDO(ps, df, get_r2 = FALSE)


### QIN 2012
ps <- psf %>% 
  ps_filter(paper == "qin2012") 

df <- data.frame(metadata) %>% 
  filter(paper == "qin2012") %>% 
  drop_na(hba1c.perc) %>% 
  # standardize fbg
  mutate(hba.stand = (hba1c.perc - mean(hba1c.perc) ) / sd(hba1c.perc)) %>% 
  dplyr::select(hba.stand, metformin)

qin <- myFIDO(ps, df, get_r2 = FALSE)

## SUN 2018 

ps <- psf %>% 
  ps_filter(paper == "sun2018") 

df <- data.frame(metadata) %>% 
  filter(paper == "sun2018") %>% 
  ### GET ONLY BEFORE METFORMIN 
  filter(metformin == "no.metformin") %>% 
  drop_na(hba1c.perc) %>% 
  # standardize fbg
  mutate(hba.stand = (hba1c.perc - mean(hba1c.perc) ) / sd(hba1c.perc)) %>% 
  dplyr::select(hba.stand)

sun <- myFIDO(ps, df, get_r2 = FALSE)

### ELBERE 2020
ps <- psf %>% 
  ps_filter(paper == "elbere2020") 

df <- data.frame(metadata) %>% 
  filter(paper == "elbere2020") %>% 
  ### GET ONLY BEFORE METFORMIN 
  filter(metformin == "no.metformin") %>% 
  drop_na(hba1c.perc) %>% 
  # standardize fbg
  mutate(hba.stand = (hba1c.perc - mean(hba1c.perc) ) / sd(hba1c.perc)) %>% 
  dplyr::select(hba.stand)

elb <- myFIDO(ps, df, get_r2 = FALSE)


## ---- save all ----

## add together
all_pars <- rbind(fors %>% mutate(Study = "Forslund 2014"),
                  karl %>% mutate(Study = "Karlsson 2013"),
                  lech %>% mutate(Study = "LeChatlier 2013"),
                  li %>% mutate(Study = "Li 2014"),
                  qin %>% mutate(Study = "Qin 2012"),
                  sun %>% mutate(Study = "Sun 2018"),
                  elb %>% mutate(Study = "Elbere 2020"))

save(all_pars, file = "results/fido-hba1c-bypaper-params.RData")

## ---- FILTER ----


## the significant coords from global comparison
sigs <- sig.params$coord # 9

## filter the full dataset to only significant taxa
bypaper <- all_pars %>% 
  filter(coord %in% sigs) %>% 
  # get only insulin
  filter(covariate == "hba.stand") %>% 
  mutate(Genus = str_remove(coord, "clr_"))

### THRESHOLD: there were 7 studies in the comparison,
# so consider significant if 5 have same direction and have magnitude over 0.1


taxon <- unique(bypaper$Genus)
outdf <- data.frame()
for(i in 1:length(taxon)) {
  
  # subset
  sub <- bypaper %>% filter(Genus == taxon[i])
  outdf[i, "Genus"] <- taxon[i]
  
  # summarize effect size direction (sign of mean)
  
  if(length(which(sign(sub$mean) == 1)) >= 5 |
     length(which(sign(sub$mean) == -1)) >= 5) {
    
    # create new column in outdf
    outdf[i, "signs"] <- TRUE
  } else { outdf[i, "signs"] <- FALSE}
  
  
  # summarize strength of effect size
  if(length(which(abs(sub$mean) > 0.1)) >= 5) {
    
    outdf[i, "size"] <- TRUE
  } else { outdf[i, "size"] <- FALSE}
  
}

outdf

## save taxa for new plot
mod1 <- outdf %>% 
  filter(signs == TRUE & size == TRUE) # 0


