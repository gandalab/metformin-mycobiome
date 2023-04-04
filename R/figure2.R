### Figure; T2DMET v T2DNOMET cross-sectional

# EVS 2/2023


library(ggpubr)
library(tidyverse)

# get colors
source("R/colors.R")

# get taxonomy table
load("data/fungi-phyloseq-filtered.RData")
psf <- tax_glom(psf, "Genus")
tax <- psf %>% tt_get() %>% as.data.frame()

# CS sigs
sigcs <- c("Aspergillus", "Brettanomyces", "Fusarium",
           "Saccharomycetaceae Family", "Saccharomycetales Order", "Scheffersomyces",
           "Tetrapisispora", "Thermothelomyces", "Zygotorulaspora")
# read in CS pars
load("results/T2DMETvNOMET_CS-params.RData")
cspars <- pars %>% 
  mutate(Genus = str_remove(coord, "clr_")) %>% 
  filter(covariate == "metformin") %>% 
  filter(Genus %in% sigcs) %>% 
  left_join(tax)

globsig <- c("Fusarium", "Nakaseomyces", "Tetrapisispora")
# read in global pars
load("results/T2DMETvNOMET_ALL-params.RData")
globpars <- pars %>% 
  mutate(Genus = str_remove(coord, "clr_")) %>% 
  filter(covariate == "metformin") %>% 
  filter(Genus %in% globsig) %>% 
  mutate(panel = "Global", Study = "Global") %>% 
  left_join(tax)

# read in global by-paper
load("results/T2DMETvT2DNOMET_ALL-bypaper.RData")
bypaper <- all_par %>% 
  mutate(Genus = str_remove(coord, "clr_")) %>% 
  filter(covariate == "metformin") %>% 
  filter(Genus %in% globsig) %>% 
  mutate(panel = "bypaper") %>% 
  left_join(tax)

## ---- panel A; CS ----

# plot
pa <- ggplot(data = cspars,
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
  xlim(-0.8, 0.8)

## ---- panel B; all with paper ----

# combine
both <- rbind(globpars, bypaper) %>% 
  mutate(Study = factor(Study, ordered = TRUE,
                        levels = c("Elbere 2020", "Karlsson 2013",
                                   "Li 2014", "Qin 2012", "Sun 2018", "Wu 2017", "Global")))
## change order of genera to match panel A
both <- both %>% 
  mutate(Genus = factor(Genus, ordered = TRUE,
                        levels = c("Fusarium", "Tetrapisispora",
                                   "Nakaseomyces")))
  
pb <- ggplot(data = both,
       aes(y = reorder(Genus, mean), # changes order of y axis
           x = mean)) +
  # dashed line at 0
  geom_vline(xintercept = 0, linetype = "dashed", color = "darkgrey") +
  
  geom_linerange(data = both %>% filter(panel == "Global"),
                 aes(y = reorder(Genus, mean), x = mean,
                     xmin = p2.5, xmax = p97.5), linewidth = 1.5,
                 color = "black"
  ) +
  geom_point(aes(color = Study, shape = panel, 
                 fill = Study, size = panel, alpha = panel)) +
  theme_pubr() +
  labs(y = NULL, x = "Centered Log-Ratio Value") +
  # various theme stuff
  theme(axis.text.y = element_text(face = "italic"),
        #strip.background = element_rect(fill = "white"),
        text = element_text(size = 14),
        legend.position = "right") +
  scale_color_manual(values = c(pap_cols, Global = "black"), name = "Study") +
  scale_fill_manual(values = c(pap_cols, Global = "black"), name = "Study") +
  scale_shape_manual(values = c(23, 24), guide = "none") +
  scale_size_manual(values = c(5, 7), guide = "none") +
  scale_alpha_manual(values = c(Global = 0.9, bypaper = 0.6), guide = "none")


## ---- combine ----

ggarrange(pa, pb, ncol = 2, legend = "bottom" #labels = c("A.", "B."), 
          #font.label = list(size = 18, face = "bold")
          ) +
  theme(plot.background = element_rect(fill = "white", color = "white"))

# save
ggsave(filename = "figures/fig-T2DMETvT2DNOMET_CS.png", dpi = 600,
       height = 6, width = 14, units = "in")
