#!/bin/bash


## kraken 2 on standard db and then kraken 2 fungi on the unclassified reads

## ---- VARIABLES ----
conda activate bioinfo

set -uex

# output name
OUTNAME=studyID

# working directory with list of run IDs
DIR=/path

# input directory (fastq files)
INDIR=~/scratch/lengthtrim/$OUTNAME/

# kraken2 database - standard build
KSTAN=/refs/krakdb

# kraken2 database - fungi only
KFUNG=/refs/krakfun

# output directory for kraken reports - standard
OUTST=/path
# make directory if it doesn't exist
mkdir -p $OUTST

# output directory for kraken unclassified
UNCLASS=~/path
mkdir -p $UNCLASS

# output directory for kraken reports - fungi
REFUN=/path
# make directory if it doesn't exist
mkdir -p $REFUN

# output directory for kraken fungal output
OUTFUN=~/scratch/kraken-output/jan22/$OUTNAME
mkdir -p $OUTFUN

## ---- kraken standard database ----

## write unclassified reads to file for next step
cat $DIR/runids_$OUTNAME.txt | parallel -j10 /storage/work/epb5360/miniconda3/envs/bioinfo/bin/kraken2 --paired --db $KSTAN --threads 10 --unclassified-out $UNCLASS/{}#.unclas.fastq --output "-" --report $OUTST/{}_krakenst.txt $INDIR/{}_trim_1.fastq $INDIR/{}_trim_2.fastq

# print finished message
echo "kraken standard complete"

## ---- kraken fungal database on unclassified reads ----

cat $DIR/runids_$OUTNAME.txt | parallel -j10 /storage/work/epb5360/miniconda3/envs/bioinfo/bin/kraken2 --paired --db $KFUNG --threads 10 --output $OUTFUN/{} --report $REFUN/{}_krakenfung.txt $UNCLASS/{}_1.unclas.fastq $UNCLASS/{}_2.unclas.fastq

# print finished message
echo "kraken fungi complete"
