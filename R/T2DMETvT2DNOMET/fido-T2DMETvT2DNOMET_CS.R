### FIDO; cross-sectional T2D-MET vs T2D-NOMET
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
  ps_filter(paper %in% c("qin2012", "li2014", "karlsson2013")) %>% 
  merge_samples(group = "paper") %>% 
  tax_transform("pa") %>% 
  tax_filter(min_prevalence = 3)

# filter original ps by shared taxa
psh <- prune_taxa(taxa_names(shared), psf)

### get only RCTs
psh <- psh %>% 
  ps_mutate(phenmet = paste(phenotype, metformin, sep = "_")) %>% 
  ps_filter(phenmet %in% c("T2D_metformin", "T2D_no.metformin")) %>% 
  # 1/23; dummy variables
  ps_mutate(metformin = if_else(phenmet == "T2D_metformin", 1, 0)) %>% 
  # get only cross-sectioal studies
  ps_filter(paper %in% c("qin2012", "li2014", "karlsson2013")) 
# 191 samples, 40 taxa

#### NO PAIRED SAMPLES
metadata <- sample_data(psh)
meta <- data.frame(metadata) %>% 
  dplyr::select(metformin, paper)

### what is the distribution of clinical markers in this population?
metadata %>% group_by(metformin) %>% 
  get_summary_stats(fasting.glucose.mg.dL, type = "median_iqr")
metadata %>% group_by(metformin) %>% 
  get_summary_stats(hba1c.perc, type = "mean_sd")

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
save(pars, file = "results/T2DMETvNOMET_CS-params.RData")

### ----- per - paper -----

source("R/myFIDO.R")

### QIN 2012
qin <- myFIDO(psh %>% 
               ps_filter(paper == "qin2012"),
             meta %>% 
               filter(paper == "qin2012") %>% 
               dplyr::select(metformin, starts_with("qin2012")),
             2)

### KARL 2013
karl <- myFIDO(psh %>% 
                ps_filter(paper == "karlsson2013"),
              meta %>% 
                filter(paper == "karlsson2013") %>% 
                dplyr::select(metformin, starts_with("karlsson2013")),
              2)

### LI 2014
li <- myFIDO(psh %>% 
                 ps_filter(paper == "li2014"),
               meta %>% 
                 filter(paper == "li2014") %>% 
                 dplyr::select(metformin, starts_with("li2014")),
               2)

### combine params and save
all_par <- rbind(qin %>% mutate(Study = "Qin 2012"),
                 karl %>% mutate(Study = "Karlsson 2013"),
                 li %>% mutate(Study = "Li 2014"))

save(all_par, file = "results/T2DMETvT2DNOMET_CS-bypaper-allparams.RData")

## ----- filter candidate taxa ----

## the significant coords from global comparison
sigs <- sig.params$coord

## filter the full dataset to only significant taxa
bypaper <- all_par %>% 
  filter(coord %in% sigs) %>% 
  # get only metformin covariate
  filter(covariate == "metformin") %>% 
  mutate(Genus = str_remove(coord, "clr_"))

### THRESHOLD: there were 3 studies in the comparison,
# so consider significant if 2 have same direction and have magnitude over 0.1


taxon <- unique(bypaper$Genus)
outdf <- data.frame()
for(i in 1:length(taxon)) {
  
  # subset
  sub <- bypaper %>% filter(Genus == taxon[i])
  outdf[i, "Genus"] <- taxon[i]
  
  # summarize effect size direction (sign of mean)
  
  if(length(which(sign(sub$mean) == 1)) >= 2 |
     length(which(sign(sub$mean) == -1)) >= 2) {
    
    # create new column in outdf
    outdf[i, "signs"] <- TRUE
  } else { outdf[i, "signs"] <- FALSE}
  
  
  # summarize strength of effect size
  if(length(which(abs(sub$mean) > 0.1)) >= 2) {
    
    outdf[i, "size"] <- TRUE
  } else { outdf[i, "size"] <- FALSE}
  
}

outdf

## save taxa for new plot
mod1 <- outdf %>% 
  filter(signs == TRUE & size == TRUE) 
