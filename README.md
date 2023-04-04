# metformin-mycobiome

Repo for data and analysis for the manuscript Van Syoc et al 2023.

*Some data used in these analyses cannot be publicly shared since they were obtained from study authors or institutional Data Use Agreements.* Bash scripts to obtain raw metagenomics files and obtain fungal taxonomy are fully reproducible. R scripts for downstream analyses are not reproducible due to data sharing limitations; all files start with a phyloseq object containing fungal taxonomy and metadata but the phyloseq objects are not provided. 

The following analyses are provided:
1. Bash scripts to obtain raw metagenomics files from Sequence Read Archive and obtain bacterial and fungal taxonomy.
2. Filtering bacterial and fungal taxonomy for relative abundance.
3. Bayesian multinomial logistic-normal linear models (MLN) for effects of metformin treatment in T2D
3.1 Building figure 1
3.2 Building figure 2
4. MLN for comparison of T2D-MET to NORM 
4.1 Mouse study fungal taxonomy, filtering/decontamination, and MLN comparison
4.2 Building Figure 3
5. MLN comparison of T2D-NOMET to NORM
5.1 MLN comparison of obesity and glycemic metrics
5.2 Building Figure 4
6. Fungal alpha diversity 
6.1 Building Figure 5
6.2 Building Figure 6
