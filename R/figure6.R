### FUNGI HEATMAP
library(phyloseq)
library(microViz)
library(ggpubr)
library(tidyverse)
library(grid)
library(ComplexHeatmap)

source("R/colors.R")


# get data
load("data/fungi-phyloseq-filtered.RData")
psf <- tax_glom(psf, "Genus")
taxa_names(psf) <- psf@tax_table[,"Genus"]


### get data for heatmap
mps <- psf %>% 
  # make phenmet column and make pretty for plot
  ps_mutate(phenmet = paste(phenotype, metformin, sep = "_")) %>% 
  # ps_filter(phenmet %in% c("T2D_metformin", "T2D_no.metformin", "healthy_no.metformin"),
  #  .keep_all_taxa = TRUE) %>% 
  ps_mutate(forplot = case_when(
    phenmet %in% "T2D_metformin" ~ "T2D-MET",
    phenmet %in% "T2D_no.metformin" ~ "T2D-NOMET",
    phenmet %in% "healthy_no.metformin" ~ "NORM"
  )) %>% 
  ps_filter(!is.na(forplot)) %>% 
  # beautify Study names
  ps_mutate(
    Study = str_to_title(paper),
    Study = str_replace(Study, "2", " 2")) %>% 
  ps_mutate(Study = if_else(str_detect(Study, "atlier"), "LeChatlier 2013", Study)) %>% 
  ps_mutate(both = paste(Study, forplot, sep = "_")) %>% 
  # show at Class level
  tax_glom("Class") %>% 
  merge_samples(group = "both", fun = sum)

taxa_names(mps) <- tax_table(mps)[,"Class"]

## add sample data back
sample_data(mps) <- data.frame(
  sample = sample_names(mps),
  row.names = sample_names(mps)
)

hm <- mps %>%
  ps_mutate(
    Treatment = sapply(str_split(sample, "_"), `[`, 2),
    Study = sapply(str_split(sample, "_"), `[`, 1)) %>% 
  tax_transform("clr") %>% 
  comp_heatmap(
    tax_anno = taxAnnotation(
      Prevalence = anno_tax_prev(ylim = 0:1),
      `RA` = anno_tax_box(trans = "compositional")
    ),
    sample_anno = sampleAnnotation(
      Treatment = anno_sample_cat(var = "Treatment",
                                  col = treat_cols,
                                  legend_title = "Treatment"),
      Study = anno_sample_cat(var = "Study",
                              col = pap_cols,
                              legend_title = "Study")
    ),
    colors = heat_palette(palette = "Blue-Red 3"),
    # change legend title
    heatmap_legend_param = list(title = "CLR"),
    # make fungi names italic
    row_names_gp = gpar(fontface = "italic")
  )

### SAVE
png(filename = "figures/fungi-heatmap.png", width = 9, height = 6, units = "in", res = 600)
draw(hm)
dev.off()


### BACTERIA HEATMAP
library(phyloseq)
library(microViz)
library(ggpubr)
library(tidyverse)
library(grid)

source("R/colors.R")
treat_cols <- treat_cols[1:3]

# get data
load("data/bac-phyloseq-filtered.RData")
psb <- tax_glom(psb, "Genus")
taxa_names(psb) <- psb@tax_table[,"Genus"]


### get data for heatmap
mps <- psb %>% 
  # make phenmet column and make pretty for plot
  ps_mutate(phenmet = paste(phenotype, metformin, sep = "_")) %>% 
  # ps_filter(phenmet %in% c("T2D_metformin", "T2D_no.metformin", "healthy_no.metformin"),
  #  .keep_all_taxa = TRUE) %>% 
  ps_mutate(forplot = case_when(
    phenmet %in% "T2D_metformin" ~ "T2D-MET",
    phenmet %in% "T2D_no.metformin" ~ "T2D-NOMET",
    phenmet %in% "healthy_no.metformin" ~ "NORM"
  )) %>% 
  ps_filter(!is.na(forplot)) %>% 
  # beautify Study names
  ps_mutate(
    Study = str_to_title(paper),
    Study = str_replace(Study, "2", " 2")) %>% 
  ps_mutate(Study = if_else(str_detect(Study, "atlier"), "LeChatlier 2013", Study)) %>% 
  ps_mutate(both = paste(Study, forplot, sep = "_")) %>% 
  # show at Class level
  tax_glom("Class") %>% 
  merge_samples(group = "both", fun = sum)

taxa_names(mps) <- tax_table(mps)[,"Class"]

## add sample data back
sample_data(mps) <- data.frame(
  sample = sample_names(mps),
  row.names = sample_names(mps)
)

hmb <- mps %>%
  ps_mutate(
    Treatment = sapply(str_split(sample, "_"), `[`, 2),
    Study = sapply(str_split(sample, "_"), `[`, 1)) %>% 
  tax_transform("clr") %>% 
  comp_heatmap(
    tax_anno = taxAnnotation(
      Prevalence = anno_tax_prev(ylim = 0:1)#,
      #`RA` = anno_tax_box(trans = "compositional")
    ),
    sample_anno = sampleAnnotation(
      Treatment = anno_sample_cat(var = "Treatment",
                                  col = treat_cols,
                                  legend_title = "Treatment"),
      Study = anno_sample_cat(var = "Study",
                              col = pap_cols,
                              legend_title = "Study")
    ),
    colors = heat_palette(palette = "Blue-Red 3"),
    # change legend title
    heatmap_legend_param = list(title = "CLR"),
    # make fungi names italic
    row_names_gp = gpar(fontface = "italic", fontsize = 9)
  )

### SAVE
png(filename = "figures/sfig-full-bac-heatmap.png", 
    width = 7, height = 9, units = "in", res = 600)
hmb %>% ComplexHeatmap::draw(annotation_legend_list = attr(hmb, "AnnoLegends"))
dev.off()



## ---- subset to top 10 classes to match fungi ----

## subset to top 10


sub <- prune_taxa(tax_top(mps, n = 10, rank = "Class"), mps)

## plot
hmsub <- sub %>%
  ps_mutate(
    Treatment = sapply(str_split(sample, "_"), `[`, 2),
    Study = sapply(str_split(sample, "_"), `[`, 1)) %>% 
  tax_transform("clr") %>% 
  comp_heatmap(
    tax_anno = taxAnnotation(
      Prevalence = anno_tax_prev(ylim = 0:1),
      `RA` = anno_tax_box(trans = "compositional")
    ),
    sample_anno = sampleAnnotation(
      Treatment = anno_sample_cat(var = "Treatment",
                                  col = treat_cols,
                                  legend_title = "Treatment"),
      Study = anno_sample_cat(var = "Study",
                              col = pap_cols,
                              legend_title = "Study")
    ),
    colors = heat_palette(palette = "Blue-Red 3"),
    # change legend title
    heatmap_legend_param = list(title = "CLR"),
    # make fungi names italic
    row_names_gp = gpar(fontface = "italic")
  )

### SAVE
png(filename = "figures/fig-subset-bacteria-heatmap.png",
    res = 600, height = 6, width = 8, units = "in")
hmsub %>% ComplexHeatmap::draw(annotation_legend_list = attr(hmsub, "AnnoLegends"))
dev.off()
