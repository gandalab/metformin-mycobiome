#!/bin/bash

## --- get reads from SRA

conda activate bioinfo

set -uex

## ---- CHANGE THESE VARIABLES ----

# general output name
OUTNAME=studyID

# input scratch directory
DIR=/path/

# output directory for raw reads
RAWDIR=~/scratch/allreads/

# output directory for fastqc and multiqc quality
QUALOUT=~/scratch/quality/

## ---- get reads from SRA with fastq-dump ---

# make output directory if it doesn't exist
mkdir -p $RAWDIR

# fastq-dump
cat $DIR/runids_$OUTNAME.txt | parallel -j10 /storage/work/epb5360/miniconda3/envs/bioinfo/bin/fastq-dump {} --split-files --outdir $RAWDIR

echo "read dump complete"

## ---- run fastqc and multiqc for quality ----

# make output directory if it doesn't exist
mkdir -p $QUALOUT/$OUTNAME

# run fastqc on forward reads 
cat $DIR/runids_$OUTNAME.txt | parallel -j10 /storage/work/epb5360/miniconda3/envs/bioinfo/bin/fastqc $RAWDIR/{}_1.fastq -o $QUALOUT/$OUTNAME

# run fastqc on reverse reads 
cat $DIR/runids_$OUTNAME.txt | parallel -j10 /storage/work/epb5360/miniconda3/envs/bioinfo/bin/fastqc $RAWDIR/{}_1.fastq -o $QUALOUT/$OUTNAME

# run multiqc on all reads
multiqc $QUALOUT/$OUTNAME/* -o $QUALOUT -n $OUTNAME.multiqc

echo "quality check complete"
