### FIDO - covariates - by paper

## fasting insulin

# EVS 2/2023

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
  get_summary_stats(fasting.insulin.pmol.L, type = "mean_sd")
# forslund, karlsson, lechatlier, li, qin, sun

### get only taxa that are shared between all papers
shared <- psf %>% ps_filter(!paper %in% c("wu2017", "nielsen2014", "elbere2020")) %>% merge_samples(group = "paper") %>% 
  tax_transform("pa") %>% 
  tax_filter(min_prevalence = 6)

# filter original ps by shared taxa
psh <- prune_taxa(taxa_names(shared), psf)

## ---- global test ----

# NO PAIRED SAMPLES - use only baseline values from RCTs (only provided one set of clinical vars - not long enough to see differences in BMI)
ps <- psh

metadata <- sample_data(psh)
df <- data.frame(metadata) %>% 
  drop_na(fasting.insulin.pmol.L) %>% 
  # remove post-metformin samples from RCTs (sun and elbere)
  mutate(papermet = paste0(paper, metformin)) %>% 
  filter(!papermet %in% c("sun2018metformin")) %>% 
  # standardize
  
  mutate(fpi.stand = (fasting.insulin.pmol.L - mean(fasting.insulin.pmol.L) ) / sd(fasting.insulin.pmol.L)) %>% 
  mutate(metformin = if_else(metformin == "metformin", 1, 0))

# get summary stats
df %>% get_summary_stats(fasting.insulin.pmol.L, type = "common")

# get otu tabl
otu <- otu_table(ps)

## make model matrix
X <- t(model.matrix(~ fpi.stand + metformin + paper, data = df))

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
  filter(covariate == "fpi.stand") %>% 
  filter(sign(p2.5) == sign(p97.5)) # 10

# save results
inspar <- params
save(inspar, file = "results/fido-insulin-params.RData")

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
  drop_na(fasting.insulin.pmol.L) %>% 
  # standardize fbg
  mutate(insulin.stand = (fasting.insulin.pmol.L - mean(fasting.insulin.pmol.L) ) / sd(fasting.insulin.pmol.L)) %>% 
  dplyr::select(insulin.stand)

## run 
fors <- myFIDO(ps, df, covariates.include = 2)

### KARLSSON

## subset metadata and OTU table
ps <- psf %>% 
  ps_filter(paper == "karlsson2013") 

df <- data.frame(metadata) %>% 
  filter(paper == "karlsson2013") %>% 
  drop_na(fasting.insulin.pmol.L) %>% 
  # standardize fbg
  mutate(insulin.stand = (fasting.insulin.pmol.L - mean(fasting.insulin.pmol.L) ) / sd(fasting.insulin.pmol.L)) %>% 
  dplyr::select(metformin, insulin.stand)

karl <- myFIDO(ps, df, get_r2 = FALSE)

### LECHATLIER 

ps <- psf %>% 
  ps_filter(paper == "lechatlier2013") 

df <- data.frame(metadata) %>% 
  filter(paper == "lechatlier2013") %>% 
  drop_na(fasting.insulin.pmol.L) %>% 
  # standardize fbg
  mutate(insulin.stand = (fasting.insulin.pmol.L - mean(fasting.insulin.pmol.L) ) / sd(fasting.insulin.pmol.L)) %>% 
  dplyr::select(insulin.stand)

lech <- myFIDO(ps, df, get_r2 = FALSE)

### LI 2014

ps <- psf %>% 
  ps_filter(paper == "li2014") 

df <- data.frame(metadata) %>% 
  filter(paper == "li2014") %>% 
  drop_na(fasting.insulin.pmol.L) %>% 
  # standardize fbg
  mutate(insulin.stand = (fasting.insulin.pmol.L - mean(fasting.insulin.pmol.L) ) / sd(fasting.insulin.pmol.L)) %>% 
  dplyr::select(insulin.stand, metformin)

li <- myFIDO(ps, df, get_r2 = FALSE)


### QIN 2012
ps <- psf %>% 
  ps_filter(paper == "qin2012") 

df <- data.frame(metadata) %>% 
  filter(paper == "qin2012") %>% 
  drop_na(fasting.insulin.pmol.L) %>% 
  # standardize fbg
  mutate(insulin.stand = (fasting.insulin.pmol.L - mean(fasting.insulin.pmol.L) ) / sd(fasting.insulin.pmol.L)) %>% 
  dplyr::select(insulin.stand, metformin)

qin <- myFIDO(ps, df, get_r2 = FALSE)

## SUN 2018 

ps <- psf %>% 
  ps_filter(paper == "sun2018") 

df <- data.frame(metadata) %>% 
  filter(paper == "sun2018") %>% 
  ### GET ONLY BEFORE METFORMIN 
  filter(metformin == "no.metformin") %>% 
  drop_na(fasting.insulin.pmol.L) %>% 
  # standardize fbg
  mutate(insulin.stand = (fasting.insulin.pmol.L - mean(fasting.insulin.pmol.L) ) / sd(fasting.insulin.pmol.L)) %>% 
  dplyr::select(insulin.stand)

sun <- myFIDO(ps, df, get_r2 = FALSE)

## ---- save all ----

## add together
all_pars <- rbind(fors %>% mutate(Study = "Forslund 2014"),
                  karl %>% mutate(Study = "Karlsson 2013"),
                  lech %>% mutate(Study = "LeChatlier 2013"),
                  li %>% mutate(Study = "Li 2014"),
                  qin %>% mutate(Study = "Qin 2012"),
                  sun %>% mutate(Study = "Sun 2018"))

save(all_pars, file = "results/fido-insulin-bypaper-params.RData")

## ---- FILTER ----


## the significant coords from global comparison
sigs <- sig.params$coord # 10

## filter the full dataset to only significant taxa
bypaper <- all_pars %>% 
  filter(coord %in% sigs) %>% 
  # get only insulin
  filter(covariate == "insulin.stand") %>% 
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
  filter(signs == TRUE & size == TRUE) # 4
