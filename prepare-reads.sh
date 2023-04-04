#!/bin/bash


## full workflow 

conda activate bioinfo
set -uex

## ----- set variables ----

# name
OUTNAME=studyID

# minimum length to trim (2/3 of average read length)
MINLEN=44

# working directory
DIR=/path/

# raw read directory
#RAWDIR=~/scratch/allreads

# adapter directory
ADDIR=~/scratch/adapter/$OUTNAME
mkdir -p $ADDIR

# human read removal directory
BOWDIR=~/scratch/hostremoval/$OUTNAME
mkdir -p $BOWDIR

# directory to bowtie database
BTDB=/path/to/hg37dec_v0.1

# quality trim directory
TRIMDIR=~/scratch/lengthtrim/$OUTNAME
mkdir -p $TRIMDIR

# final fastq and multiqc directory
QUALDIR=~/scratch/fastpquality/$OUTNAME
mkdir -p $QUALDIR

## ---- detect and remove adapters with fastp ----

cat $DIR/runids_$OUTNAME.txt | parallel -j10 /storage/work/epb5360/miniconda3/envs/bioinfo/bin/fastp -i $RAWDIR/{}_1.fastq -I $RAWDIR/{}_2.fastq -o $ADDIR/{}_1_adap.fastq -O $ADDIR/{}_2_adap.fastq --detect_adapter_for_pe
echo "adapter removal complete"

## ---- remove human reads with bowtie2 ----

cat $DIR/runids_$OUTNAME.txt | parallel -j10 /storage/work/epb5360/miniconda3/envs/bioinfo/bin/bowtie2 --quiet -x $BTDB -1 $ADDIR/{}_1_adap.fastq -2 $ADDIR/{}_2_adap.fastq --un-conc $BOWDIR/{}_hostremoved_%.fastq

echo "host removal complete"

## ---- quality and length trim with fastp ----

cat $DIR/runids_$OUTNAME.txt | parallel -j10 /storage/work/epb5360/miniconda3/envs/bioinfo/bin/fastp -i $BOWDIR/{}_hostremoved_1.fastq -I $BOWDIR/{}_hostremoved_2.fastq -o $TRIMDIR/{}_trim_1.fastq -O $TRIMDIR/{}_trim_2.fastq --cut_tail --cut_mean_quality 20 --length_required $MINLEN

## ---- quality with fastqc and multiqc ----

# fastqc on forward reads
cat $DIR/runids_$OUTNAME.txt | parallel -j10 /storage/work/epb5360/miniconda3/envs/bioinfo/bin/fastqc $TRIMDIR/{}_trim_1.fastq -o $QUALDIR

# multiqc on fastqc reports
multiqc $QUALDIR/* -o $QUALDIR -n $OUTNAME.multiqc

echo "quality check complete"

