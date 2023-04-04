## mouse bacteria controls 

# EVS 8/2022

# 
library(phyloseq)
library(microViz)
library(tidyverse)
source("R/colors.R")
## ---- Bacteria ----

# get data
load("R/mouse/mouse-bac-ps-filtered.RData")
psf <- psfb

# what is the composition of the negative control?
psf <- tax_fix(psf)

psf %>% 
  comp_barplot("Family")

# get just controls
psf %>% 
  ps_filter(str_detect(treatment, "Ext")) %>% 
  ps_mutate(Control = if_else(treatment == "PC_Ext", "Mock", "Blank")) %>% 
  comp_barplot("Genus",
               label = "Control")

# save
ggsave(filename = "R/updated-silvermanlab-analyses/working-plots/mouse-controls-bac.png", dpi = 600)

## plot absolute abundance
df <- psf %>% 
  ps_filter(str_detect(treatment, "Ext")) %>% 
  ps_mutate(Control = if_else(treatment == "PC_Ext", "Mock", "Blank")) %>% 
  ps_melt() %>% as_tibble() %>% 
  mutate(Abundance = as.numeric(Abundance)) %>% 
  arrange(desc(Abundance)) %>% 
  mutate(Control = factor(Control, ordered = TRUE, levels = c("Mock", "Blank")))

p <- ggbarplot(df %>% slice_head(n = 10), x = "Genus", y = "Abundance", 
          fill = "Control",
          #sort.val = "desc",
          #label = TRUE,
          position = position_dodge(width = 0.8),
          x.text.angle = 45,
          #ylim = c(0, 3500),
          xlab = "") +
  scale_fill_manual(values = accent_cols)
ggpar(p, yscale = "log10", format.scale = TRUE)

#ggsave(filename = "R/updated-silvermanlab-analyses/working-plots/mouse-control-bac-counts.png", dpi = 600, plot = last_plot())

# decontam based on negative control
library(decontam)

# ---- filter for Bacteria & remove NA phyla ----

# get only bacteria
psf <- subset_taxa(psf, Kingdom == "Bacteria")

# remove NA phylum
psf1 <- subset_taxa(psf, !is.na(Phylum) & !Phylum %in% c("", "NA"))

# validate
psf2 <- tax_fix(psf1)
psf3 <- phyloseq_validate(psf2, remove_undetected = TRUE)

## ---- decontam ----

## workflow from: https://benjjneb.github.io/decontam/vignettes/decontam_intro.html

# add sampling variables
psd <- psf3 %>% 
  ps_mutate(
    SampleBinary = if_else(str_detect(treatment,"NC_Ext"), true = "Control", false = "Sample")
  ) 

## inspect library sizes 
sampdf <- data.frame(sample_data(psd))
sampdf$LibrarySize <- sample_sums(psd)
sampdf <- sampdf[order(sampdf$LibrarySize), ]
sampdf$Index <- seq(nrow(sampdf))
ggplot(data = sampdf, aes(x = Index, y = LibrarySize, color = SampleBinary)) + geom_point()


## use negative controls to identify contaminants
pds <- psd %>% ps_mutate(is.neg = if_else(SampleBinary == "Control", true = TRUE, false = FALSE))
contamdf.prev <- isContaminant(pds, method="prevalence", neg="is.neg", threshold = 0.5) # more aggressive threshold is 0.5
table(contamdf.prev$contaminant) 

## identify the contaminant taxa and remove them from the phyloseq object
psdecon <- prune_taxa(!contamdf.prev$contaminant, pds)

## remove the positive controls and save the decontaminated phyloseq object
pssave <- ps_filter(psdecon, SampleBinary == "Sample") %>% 
  # remove positive controls
  ps_filter(!str_detect(treatment, "PC")) %>% 
  ps_select(-c(is.neg, SampleBinary))
psb <- pssave

# save
taxa_names(psb) <- psb@tax_table[,"Genus"]
save(psb, file = "data/mouse-bac-phyloseq-filtered.RData")
