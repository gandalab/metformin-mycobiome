### elbere: NORM+MET v NORM

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

### get elbere NORM v NORM+MET
psh <- psf %>% 
  ps_filter(paper == "elbere2020") %>% 
  ps_mutate(phenmet = paste(phenotype, metformin, sep = "_")) %>% 
  ps_filter(phenmet %in% c("healthy_metformin", "healthy_no.metformin")) %>% 
  # 1/23; dummy variables
  ps_mutate(metformin = if_else(phenmet == "healthy_metformin", 1, 0)) %>% 
  # remove 24H sampling
  ps_filter(!treatment == "metformin.24H")

# 69 samples, 41 taxa

## ---- get paired samples ----

metadata <- sample_data(psh)
metadata <- data.frame(metadata) 

# keep only columns that we need
meta.paired <- metadata %>% 
  dplyr::select(Sample, metformin) %>% 
  # create new column of unique names for each paper
  mutate(newID = paste("elbere2020", str_extract(Sample, "S(\\d){1,3}"), sep = "_")
  )

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
  dplyr::select(metformin, 4:ncol(meta.paired)) 


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
sig.params # 3

# get r2
quantile(r2(posterior, covariate = 2), probs = c(0.025, 0.5, 0.975))

## save parameters
elbpars <- params$Lambda
save(elbpars, file = "results/NORMMETvNORM-params.RData")
