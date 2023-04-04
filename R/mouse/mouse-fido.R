## FIDO on mouse data for metformin effect in T2D


library(fido) 
library(phyloseq)
library(plyr)
library(dplyr)
library(microViz)

## set seed
set.seed(123)

# get data
load("data/mouse-phyloseq-filtered.RData")
psf <- pssave

## ---- filter and glom to Genus ----

## filter to only T2D and metformin/no metformin (remove NA's)
psfilt <- psf %>% 
  # glom to Genus
  tax_glom("Genus") 

# for ease downstream, rename taxa names to genus
taxa_names(psfilt) <- psfilt@tax_table[,"Genus"]

# 1/23; to numeric
psfilt <- psfilt %>% 
  ps_mutate(metformin = if_else(treatment == "metformin", 1, 0)) %>% 
  ps_select(metformin)


## ---- set priors ----

##Extracting the OTU data
otu <- t(otu_table(psfilt))
dim(otu)


##Extracting the metadata
metadata <- sample_data(psfilt)
metadata <- data.frame(metadata) #

X <- t(model.matrix(~ ., data = metadata))
Y <- otu[,colnames(otu) %in% colnames(X)]

##Specifying priors
upsilon <- ntaxa(psfilt)+3 
Omega <- diag(ntaxa(psfilt))
G <- cbind(diag(ntaxa(psfilt)-1), -1)
Xi <- (upsilon-ntaxa(psfilt))*G%*%Omega%*%t(G)
Theta <- matrix(0, ntaxa(psfilt)-1, nrow(X))
Gamma <- diag(nrow(X))

# get priors
priors <- pibble(NULL, X, upsilon, Theta, Gamma, Xi)  
#print(priors)

priors <- to_clr(priors)  
names_covariates(priors) <- rownames(X)

priors$Y <- Y

## get posterior
posterior <- refit(priors, optim_method="lbfgs",
                   n_samples = 5000)

ppc_summary(posterior) # 99% 

posterior <- to_clr(posterior)

# name the posterior categories (taxa)
names_categories(posterior) <- taxa_names(psfilt)
params = summary(posterior, pars = "Lambda")

# any significant?
sig.params = params$Lambda %>% 
  filter(covariate == "metformin") %>%
  filter(sign(p2.5) == sign(p97.5))
sig.params # 1

# save params
mouse_pars <- params$Lambda
save(mouse_pars, file = "results/mouse-fido-params.RData")

## ---- r2 ----

# global
quantile(r2(posterior), probs = c(0.025, 0.5, 0.975))

# partial metformin
quantile(r2(posterior, covariate = 2), probs = c(0.025, 0.5, 0.975))
