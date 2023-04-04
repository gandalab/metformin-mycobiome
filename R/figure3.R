### FIG: T2D-MET v NORM

library(ggpubr)
library(tidyverse)

# get colors
source("R/colors.R")

# get taxonomy table
load("data/fungi-phyloseq-filtered.RData")
psf <- tax_glom(psf, "Genus")
tax <- psf %>% tt_get() %>% as.data.frame()

# get sigs for T2DMET v NORM
sigs <- c("Cryptococcus", "Fusarium",
          "Saccharomycetales Order", "Thermothielavioides")

# get mouse results
load("results/mouse-fido-params.RData")
mouse <- mouse_pars %>% 
  mutate(Genus = str_remove(coord, "clr_")) %>% 
  filter(covariate == "metformin") %>% 
  filter(Genus %in% sigs) %>% 
  left_join(tax) %>% 
  mutate(Species = "Mouse")

# get T2DMET v NORM results
load("results/T2DMETvNORM-params.RData")
tpars <- pars %>% 
  mutate(Genus = str_remove(coord, "clr_")) %>% 
  filter(covariate == "metformin") %>% 
  filter(Genus %in% sigs) %>% 
  left_join(tax) %>% 
  mutate(Species = "Human")

# get sigs for all MET v all NOMET
asigs <- c("Fusarium", "Nakaseomyces", "Scheffersomyces",
           "Tetrapisispora")

# All MET results 
load("results/allMETvallNOMET-params.RData")
apars <- pars %>% 
  mutate(Genus = str_remove(coord, "clr_")) %>% 
  filter(covariate == "metformin") %>% 
  filter(Genus %in% asigs) %>% 
  left_join(tax) %>% 
  mutate(Species = "Human")

## ---- T2DMET v NORM ----

both <- mouse %>% rbind(tpars) %>% 
  ## arrange by hand
  mutate(Genus = factor(Genus, ordered = TRUE,
                        levels = c("Saccharomycetales Order", "Cryptococcus","Thermothielavioides","Fusarium"
                                    )))

pa <- ggplot(data = both,
       aes(y = Genus, # changes order of y axis
           x = mean, xmin = p2.5, xmax = p97.5)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "darkgrey") +
  geom_linerange(aes(group = Species), position = position_dodge(width = 0.4)) +

  geom_point(aes(color = Class, shape = Species), size = 6,
             position = position_dodge(width = 0.4)) +
  theme_pubr() +
  labs(y = NULL, x = "Centered Log-Ratio Value") +
  # various theme stuff
  theme(axis.text.y = element_text(face = "italic"),
        text = element_text(size = 18),
        legend.position = "right") +
  scale_color_manual(values = c_cols, name = "") +
  scale_shape_discrete(name = "")


## save
ggsave(plot = pa, filename = "figures/fig-T2DMETvNORM.png", dpi = 600,
       height = 5, width = 9, units = "in")

## ---- all MET v all NORM ----


# add mouse
both1 <- apars %>% 
  rbind(mouse_pars %>% 
          mutate(Genus = str_remove(coord, "clr_")) %>% 
          filter(covariate == "metformin") %>% 
          filter(Genus %in% asigs) %>% 
          left_join(tax) %>% 
          mutate(Species = "Mouse"))


## plot
pb <- ggplot(data = both1,
       aes(y = reorder(Genus, mean), # changes order of y axis
           x = mean, xmin = p2.5, xmax = p97.5)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "darkgrey") +
  geom_linerange(aes(group = Species), position = position_dodge(width = 0.4)) +
  
  geom_point(aes(color = Class, shape = Species), size = 6,
             position = position_dodge(width = 0.4)) +
  theme_pubr() +
  labs(y = NULL, x = "Centered Log-Ratio Value") +
  # various theme stuff
  theme(axis.text.y = element_text(face = "italic"),
        text = element_text(size = 18),
        legend.position = "right") +
  scale_color_manual(values = accent_cols, name = "") +
  scale_shape_discrete(name = "")


## ---- add together and save

ggarrange(pa, pb, common.legend = TRUE, legend = "right",
          labels = c("A.", "B."), font.label = list(size = 20, face = "bold"))

## save
ggsave(filename = "figures/fig-T2DMETvNORM.png", dpi = 600,
       height = 6, width = 14, units = "in")
