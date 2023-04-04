### Figure: T2DNOMET v NORM and clinical indicators

# EVS 2/2023

## 
library(ggpubr)
library(tidyverse)
library(microViz)
library(phyloseq)

# get colors
source("R/colors.R")

# get taxonomy table
load("data/fungi-phyloseq-filtered.RData")
psf <- tax_glom(psf, "Genus")
tax <- psf %>% tt_get() %>% as.data.frame()

## get T2DNOMET v NORM sigs
tsigs <- c("Nakaseomyces", "Saccharomyces", "Zygosaccharomyces")
# get params
load("results/T2DNOMETvNORM-params.RData")
tpars <- pars %>% 
  mutate(Genus = str_remove(coord, "clr_")) %>% 
  filter(covariate == "phenotype") %>% 
  filter(Genus %in% tsigs) %>% 
  left_join(tax)

## get insulin sigs
isigs <- c("Tetrapisispora", "Ustilago",
           "Zygosaccharomyces", "Zygotorulaspora")
load("results/fido-insulin-params.RData")
inspar <- inspar$Lambda %>% 
  mutate(Genus = str_remove(coord, "clr_")) %>% 
  filter(covariate == "fpi.stand") %>% 
  filter(Genus %in% isigs) %>% 
  mutate(Variable = "Insulin")

## get glucose sigs
gsigs <- c("Brettanomyces", "Saccharomyces")
load("results/fido-glucose-params.RData")
glucpar <- glucpar$Lambda %>% 
  mutate(Genus = str_remove(coord, "clr_")) %>% 
  filter(covariate == "fbg.stand") %>% 
  filter(Genus %in% gsigs) %>% 
  mutate(Variable = "Glucose")

## get BMI sigs
bsigs <- c("Ustilago", "Zygotorulaspora")
load("results/fido-bmi-params.RData")
bmipar <- bmipar$Lambda %>% 
  mutate(Genus = str_remove(coord, "clr_")) %>% 
  filter(covariate == "bmi.stand") %>% 
  filter(Genus %in% bsigs) %>% 
  mutate(Variable = "BMI")


## ---- panel A: T2DNOMET v NORM ----

pa <- ggplot(data = tpars,
       aes(y = reorder(Genus, mean), # changes order of y axis
           x = mean, xmin = p2.5, xmax = p97.5)) +
  geom_linerange(linewidth = 1.3) +
  geom_point(aes(color = Class), size = 5) +
  theme_pubr() +
  labs(y = NULL, x = "Centered Log-Ratio Value") +
  # dashed line at 0
  geom_vline(xintercept = 0, linetype = "dashed", color = "darkgrey") +
  # various theme stuff
  theme(axis.text.y = element_text(face = "italic"),
        strip.background = element_rect(fill = "white"),
        text = element_text(size = 14),
        legend.position = "right") +
  scale_color_manual(values = c_cols) +
  xlim(-1, 1)

## ---- panel B: all clinical covs ----

## add together
all <- bmipar %>% rbind(glucpar) %>% rbind(inspar) %>% 
  mutate(Variable = factor(Variable, ordered = TRUE,
                           levels = c("Glucose", "Insulin", "BMI")))

pb <- ggplot(data = all,
       aes(y = reorder(Genus, mean),
           x = mean, xmin = p2.5, xmax = p97.5)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "darkgrey") +
  geom_linerange(aes(group = Variable), 
                 position = position_dodge(width = 0.7)) +
  
  geom_point(aes(color = Variable), size = 5, shape = 15,
             position = position_dodge(width = 0.7)) +
  theme_pubr() +
  labs(y = NULL, x = "Centered Log-Ratio Value") +
  # various theme stuff
  theme(axis.text.y = element_text(face = "italic"),
        text = element_text(size = 14),
        legend.position = "right") +
  scale_color_manual(values = accent_cols[6:9], name = "") +
  scale_shape_discrete(name = "")

## --- add and save ----

ggarrange(pa, pb, labels = c("A.", "B."), 
          font.label = list(size = 18, face = "bold"),
          legend = "right")

ggsave(filename = "figures/fig-T2DNOMET_clinical.png", dpi = 600,
       height = 5, width = 12, units = "in")
