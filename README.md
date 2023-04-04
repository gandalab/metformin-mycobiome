# metformin-mycobiome

Repo for data and analysis for the manuscript Van Syoc et al 2023.

*Some data used in these analyses cannot be publicly shared since they were obtained from study authors or institutional Data Use Agreements.* Bash scripts to obtain raw metagenomics files and obtain fungal taxonomy are fully reproducible. R scripts for downstream analyses are not reproducible due to data sharing limitations; all files start with a phyloseq object containing fungal taxonomy and metadata but the phyloseq objects are not provided.

The following analyses are provided:
1. Bash scripts to [obtain raw metagenomics files](bash/fastq-dump.sh) from Sequence Read Archive, [remove host reads](bash/prepare-reads.sh) and obtain bacterial and [fungal taxonomy.](bash/kraken.sh). 

2. [Filtering fungal taxonomy for relative abundance.](R/filter-to-phyloseq.R). 

3. Bayesian multinomial logistic-normal linear models (MLN) for effects of [metformin treatment in T2D](R/T2DMETvT2DNOMET/). 

3.1 Building [figure 1](R/figure1.R). 

3.2 Building [figure 2](R/figure2.R). 

4. MLN for comparison of [T2D-MET to NORM](R/T2DMETvNORM/). 

4.1 [Mouse](R/mouse/) study fungal taxonomy, filtering/decontamination, and MLN comparison. 

4.2 Building [Figure 3](R/figure3.R). 

5. MLN comparison of [T2D-NOMET to NORM](R/T2DNOMETvNORM/). 

5.1 MLN comparison of [obesity and glycemic metrics](R/T2DNOMETvNORM/). 

5.2 Building [Figure 4](R/figure4.R). 

6. [Fungal alpha diversity](R/alpha-div.R). 

6.1 Building [Figure 5](R/figure5.R). 

6.2 Building [Figure 6](R/figure6.R). 

