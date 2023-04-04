### FIGURE: T2DMET v T2DNOMET - RESULTS

# EVS 2/2023

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

## get RCT params
load("results/T2DMETvNOMET_RCT-params.RData")
rctpars <- pars

## get CS pars
load("results/T2DMETvNOMET_CS-params.RData")
cspars <- pars

## get "all" pars
load("results/T2DMETvNOMET_ALL-params.RData")
globpars <- pars

## get residuals for ordination
load("data/fungi-resids-onedraw-phyloseq.RData")
load("data/fungi-resids-avgs-phyloseq.RData")


# get just RCT and shared taxa
taxa_names(psf) <- tax_table(psf)[,"Genus"]

### get only taxa that are shared between the three papers
shared <- psf %>% 
  ps_filter(paper %in% c("wu2017", "elbere2020", "sun2018")) %>% 
  merge_samples(group = "paper") %>% 
  tax_transform("pa") %>% 
  tax_filter(min_prevalence = 3)

# filter original ps by shared taxa
psh <- prune_taxa(taxa_names(shared), psf)

### get only RCTs
forplotc <- psh %>% 
  ps_mutate(phenmet = paste(phenotype, metformin, sep = "_")) %>% 
  ps_filter(phenmet %in% c("T2D_metformin", "T2D_no.metformin")) %>% 
  # 1/23; dummy variables
  ps_mutate(metformin = if_else(phenmet == "T2D_metformin", 1, 0)) %>% 
  # get only RCTs and only one time point from Wu 2017
  ps_filter(paper %in% c("wu2017", "elbere2020", "sun2018")) %>% 
  ps_filter(!treatment == "metformin.24H") %>% 
  ps_filter(!treatment %in% c("metformin.T4mo", "placebo.to.metformin.6mo")) %>% 
  # remove wu placebo (want only paired before and after metformin)
  ps_filter(!treatment %in% c("placebo.T0", "placebo.T2mo", "placebo.T4mo")) %>% 
  # re-name
  ps_mutate(Timing = case_when(
    treatment %in% "before.metformin" ~ "Baseline",
    treatment %in% "after.3d.metformin" ~ "MET 3d",
    treatment %in% "metformin.7d" ~ "MET 7d",
    treatment %in% "metformin.baseline" ~ "Baseline",
    treatment %in% "metformin.T0" ~ "Baseline",
    treatment %in% "metformin.T2mo" ~ "MET 2mo"
  )) %>% 
  ps_mutate(Timing = factor(Timing, ordered = TRUE, 
                            levels = c("Baseline", "MET 3d", "MET 7d", "MET 2mo"))) %>% 
  ps_arrange(paper, Timing)


## ----- panel A: RCT ----

# filter to significants coords and wrangle
rctpars <- rctpars %>% 
  filter(covariate == "metformin") %>% 
  mutate(Genus = str_remove(coord, "clr_")) %>% 
  filter(Genus %in% c("Botrytis", "Tetrapisispora", "Zygotorulaspora")) %>% 
  # add taxonomy
  left_join(tax)

## plot
pa <- ggplot(data = rctpars,
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
  scale_color_manual(values = c_cols)



## ---- panel B: PCA of RCT ----

## get only shared taxxa
psb <- prune_taxa(taxa_names(shared), psresid)

# get correct sample data and RCT data
sample_data(psb) <- sample_data(forplotc)

# get unique IDs for arrows

sam <- psb %>% samdat_tbl() 
sam <- sam %>% 
  mutate(newID = case_when(
    paper %in% "sun2018" ~ paste("sun", str_remove_all(Sample, "\\w+M"), sep = "-"),
    paper %in% "elbere2020" ~ paste("elb", str_extract(Sample, "S(\\d){1,3}"), sep = "-"),
    paper %in% "wu2017" ~ paste("wu", str_extract(.sample_name, "(\\d){1,3}"), sep = "-"))) %>% 
  column_to_rownames(var = ".sample_name")
sample_data(psb) <- sam

# ordinate and beautify
pb <- psb %>% 
  dist_calc("euclidean") %>% 
  ps_mutate(Treatment = factor(if_else(metformin == 1, "MET", "Baseline"))) %>% 
  ps_mutate(Study = str_to_title(str_replace(paper, "20", " 20"))) %>% 
  ord_calc() %>% 
  ord_plot(shape = "Treatment", alpha = 0.6, color = "Study", size = 2.5,
           auto_caption = NA) %>% 
  add_paths(
    id_var = "newID",
    id_values = unique(sam$newID),
    alpha = 0.5,
    arrow = grid::arrow(length = grid::unit(1.5, units = "mm"))
  ) +
  scale_color_manual(values = pap_cols) + theme_pubr(base_size = 14, legend = "right") + bgcolor("white")


## ---- panel C: relative abundance barplot RCT ----

### 3/7/2023 update; facet by paper
source("figures/functions-palettes.R")

### NEW PLOTS
## ----- panel C plot A ----
ps1 <- forplotc %>% 
  ps_filter(paper == "sun2018") %>% 
  ps_mutate(ID = as.numeric(
    str_extract(Sample, "(\\d){1,2}"))) %>% 
  ps_arrange(desc(ID), .target = "sample_data") %>%
  tax_fix() 

plota <- myShades(phylo = ps1)
pal <- myPal(phylo = ps1)

pca <- plota %>% 
  comp_barplot(
    tax_level = "Class:Family",
    n_taxa = length(hierarchicalPal),
    tax_order = "asis", 
    palette = pal, bar_width = 0.975,
    sample_order = "asis",
    label = NULL) +
  facet_wrap(~Timing, scales = "free_y") +
  theme_transparent(base_size = 16) +
  theme(axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle = 80,
                                   vjust = 0.6),
        strip.text = element_text(size = 16),
        strip.background = element_blank()) +
  ylab("Relative Abundance") +
  coord_flip() 

## ---- panel C plot B ----

ps1 <- forplotc %>% 
  ps_filter(paper == "elbere2020") %>% 
  ps_mutate(ID = as.numeric(
    str_extract(Sample, "(\\d){1,2}"))) %>% 
  ps_arrange(desc(ID), .target = "sample_data") %>%
  tax_fix() 

plotb <- myShades(phylo = ps1)
pal <- myPal(phylo = ps1)

pcb <- plotb %>% 
  comp_barplot(
    tax_level = "Class:Family",
    n_taxa = length(hierarchicalPal),
    tax_order = "asis", 
    palette = pal, bar_width = 0.975,
    sample_order = "asis",
    label = NULL) +
  facet_wrap(~Timing, scales = "free_y") +
  theme_transparent(base_size = 16) +
  theme(axis.ticks.y = element_blank(),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 80,
                                   vjust = 0.6),
        strip.text = element_text(size = 16)) +
  ylab("Relative Abundance") +
  coord_flip() 

## ---- panel C plot C ----


ps1 <- forplotc %>% 
  ps_filter(paper == "wu2017") 
ps1 <- ps1 %>% 
  ps_mutate(ID = as.numeric(
    str_extract(sample_names(ps1), 
                "(\\d){1,2}"))) %>% 
  ps_arrange(desc(ID), .target = "sample_data") %>%
  tax_fix() 

plotc <- myShades(ps1)
pal <- myPal(ps1)

pcc <- plotc %>% 
  comp_barplot(
    tax_level = "Class:Family",
    n_taxa = length(hierarchicalPal),
    tax_order = "asis", 
    palette = pal, bar_width = 0.975,
    sample_order = "asis",
    label = NULL) +
  facet_wrap(~Timing, scales = "free_y") +
  theme_transparent(base_size = 16) +
  theme(axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle = 80,
                                   vjust = 0.6),
        strip.background = element_blank(),
        strip.text = element_text(size = 16)) +
  ylab("Relative Abundance") +
  coord_flip() 


## ---- OLD  plot ----
forplotc %>% 
  ps_filter(paper == "sun2018") %>% 
  tax_fix() %>% 
  comp_barplot("Class",
               label = NULL,
               sample_order = "asis") +
  facet_wrap(~Timing, scales = "free_x") +
  theme_transparent(base_size = 14, legend = "right") +
  theme(axis.ticks.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 16)) +
  ylab("Relative Abundance")

### ---- add together ----

pc <- ggarrange(pca, pcb, pcc, common.legend = TRUE,
          nrow = 1, legend = "bottom")
pc1 <- ggpar(pc, legend.title = list(fill = ""))
## save
ggsave(plot = pc, filename = "figures/fig-T2DMETvT2DNOMET_RCT-C.png", dpi = 600,
       height = 8, width = 18, units = "in")


ggarrange(pa, pb#, #labels = c("A.", "B.", "C."), 
         # font.label = list(size = 16, face = "bold"),
          #widths = c(1, 1.5)#, heights = c(0.7, 1)#,
          #vjust = c(0.3, 0.3, 1.5)
         ) +
  #theme(plot.margin = margin(0.5, 0, 0, 0, unit = "in")) +
  theme(plot.background = element_rect(fill = "white", color = "white"))

# save as-is for now (manipulate in Inkscape/Illustrator later)
ggsave("figures/fig-T2DMETvT2DNOMET_RCT.png", dpi = 600,
       height = 4, width = 12, units = "in")

### save just panel C to add together in Illustrator
ggsave(plot = pc, filename = "figures/fig-T2DMETvT2DNOMET_RCT-C.png", dpi = 600,
       height = 5, width = 11, units = "in")
