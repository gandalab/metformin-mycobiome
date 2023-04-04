### figure showing low and 'normal' Shannons

## EVS 2/2023

library(phyloseq)
library(tidyverse)
library(microViz)
library(ggpubr)
library(cowplot)

source("R/colors.R")

# get taxonomy table
load("data/fungi-phyloseq-filtered.RData")
psf <- tax_glom(psf, "Genus")


# get alpha diversity
load("results/fungi-shannons.RData")

# get alpha under 1
alpha <- alpha %>% 
  mutate(is.low = if_else(Shannon < 1, "low", "normal")) %>% 
  column_to_rownames(var = "ID")

sample_data(psf) <- alpha

## add Study
psf <- psf %>% 
  ps_mutate(Study = str_to_title(str_replace(paper, "20", " 20"))) %>% 
  ps_mutate(Study = if_else(str_detect(Study, "Chat"), "LeChatlier 2013", Study))



## ----- 3/7; diff colors and sorting ----
source("figures/functions-palettes.R")

# sort 
#l <- psf %>% 
 # ps_filter(is.low == "low") 
lp <- myShades(psf, nHue = 5, nShade = 3, rank1 = "Class",
               rank1pl = "Class",
               rank2 = "Family")
pal <- myPal(psf, nHue = 5, nShade = 3, rank1 = "Class",
             rank1pl = "Class",
             rank2 = "Family")

pa <- lp %>%
  ps_filter(is.low == "low", .keep_all_taxa = TRUE) %>% 
  #tax_agg("Family:Genus") %>% 
  tax_transform("compositional") %>%
  ps_arrange(desc(Saccharomycetaceae), .target = "otu_table") %>%
  comp_barplot(tax_level = "Class:Family", 
               sample_order = "asis", 
               tax_order = "asis",
               label = NULL,
               n_taxa = length(pal),
               palette = pal, bar_width =1,
               bar_outline_color = NA) +
  
  theme_pubr(base_size = 14, legend = "top") +
  theme(axis.ticks.y = element_blank()) +
  coord_flip() +
  ylab("Relative Abundance") 

## ----- 3/7; panel B ----

# sort 
#n <- psf %>% 
 # ps_filter(is.low == "normal") 
#np <- myShades(psf, nHue = 5, nShade = 3, rank1 = "Class",
 #              rank1pl = "Class",
  #             rank2 = "Family")
#pal <- myPal(n, nHue = 5, nShade = 3, rank1 = "Class",
 #            rank1pl = "Class",
  #           rank2 = "Family")

pb <- lp %>%
  ps_filter(is.low == "normal", .keep_all_taxa = TRUE) %>% 
  #tax_agg("Family:Genus") %>% 
  tax_transform("compositional") %>%
  ps_arrange(desc(Saccharomycetaceae), .target = "otu_table") %>%
  comp_barplot(tax_level = "Class:Family", 
               sample_order = "asis", 
               tax_order = "asis",
               label = NULL,
               n_taxa = length(pal),
               palette = pal, bar_width = 1,
               bar_outline_colour = NA) +
  
  theme_pubr(base_size = 14, legend = "top") +
  theme(axis.ticks.y = element_blank()) +
  coord_flip() +
  ylab("Relative Abundance") 

## ---- combine and save ----

ggarrange(pa, pb, ncol = 2, 
          #labels = c("A.", "B."), font.label = list(size = 20, face = "bold"),
          #widths = c(1, 3),
          common.legend = TRUE, legend = "bottom") 

## save
ggsave(filename = "figures/fig-adiv-relabun.png", dpi = 300,
       height = 8, width = 18, units = "in") 


## ----- OLD -----

## panel
pb <- psf %>% 
  ps_filter(is.low == "low") %>% 
  tax_fix() %>% 
  comp_barplot("Genus",
               n_taxa = 8,
               label = NULL,
               bar_outline_width = 0.005)+
  theme_pubr(base_size = 14, legend = "top") +
  theme(axis.ticks.x = element_blank())

pb <- psf %>% 
  ps_filter(is.low == "normal") %>% 
  tax_fix() %>% 
  comp_barplot("Genus",
               n_taxa = 8,
               label = NULL,
               bar_outline_width = 0.005) +
  theme_pubr(base_size = 18, legend = "top") +
  theme(axis.ticks.x = element_blank())

## arrange
ggarrange(pa, pb, ncol = 2, 
          #labels = c("A.", "B."), font.label = list(size = 20, face = "bold"),
          widths = c(1, 3),
          common.legend = TRUE, legend = "right") 

## save
ggsave(filename = "figures/fig-adiv-relabun.png", dpi = 300,
       height = 8, width = 18, units = "in") +
  theme(plot.background = element_rect(fill = "white", color = "white"))


## ---- ordination ---

### plot ordination
# load residuals
load("data/fungi-resids-avgs-phyloseq.RData")
sample_data(psresid.avgs) <- alpha
# plot
psresid.avgs %>% 
  tax_fix() %>% 
  #tax_transform("clr") %>% 
  ord_calc() %>% ord_plot(color = "is.low",
                          alpha = 0.5, auto_caption = NA,
                          plot_taxa = 1:5,
                          tax_lab_style = 
                            tax_lab_style(size = 4,
                                         position = position_jitter(width = 0.3,
                                                                   height = 1)),
                          tax_vec_length = 0.6) +
  stat_ellipse(aes(color = is.low), size = 1) +
  scale_color_manual(values = accent_cols[4:5], name = "Shannon's H",
                     labels = c("H < 1", "H > 1")) +
  bgcolor("white") + theme_pubr(legend = "right")



# save
ggsave(filename = "figures/fig-adiv-ordination.png", dpi = 600,
       height = 5, width = 8, units = "in")
