#!/bin/bash

## script to get runinfo and test accession number

# activate conda
conda activate bioinfo

# catch errors
set -uex

## ---- CHANGE THESE VARIABLES ----

# accession number
ACC=PRJEB39500

# input directory
DIR=/path

# general output name
OUTNAME=studyID

# set prefix (European is ERR, American is SRR)
PREF=ERR

## ---- workflow: get runinfo and test accession ----
# from Biostars Handbook

# get runinfo
esearch -db sra -query $ACC | efetch -format runinfo > $DIR/runinfo_$OUTNAME.csv

# get runids
cat $DIR/runinfo_$OUTNAME.csv | cut -f 1 -d ',' | grep $PREF > $DIR/runids_$OUTNAME.txt

# get one sample with 10,000 reads to test connection
cat $DIR/runids_$OUTNAME.txt | head -1 | parallel /storage/work/epb5360/miniconda3/envs/bioinfo/bin/fastq-dump -X 10000 --split-files --outdir $DIR {}

# print stats for the reads
seqkit stat $DIR/*.fastq

# print finished message
echo "done!"

