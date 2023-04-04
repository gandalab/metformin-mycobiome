#!/bin/bash
conda activate base

## script to revert new report (with minimizer) back to old report and build biom tables
 
# name in an outside variable
OUTNAME=$1

# working directory
DIR=/path

# new directory - standard 
NEWDIR=/path
mkdir -p $NEWDIR
 
# new directory - fungi
NEWDIRF=/path
mkdir -p $NEWDIRF

# old directory - standard
OLDDIR=/path

# old directory - fungi
OLDDIRF=/path

# biom directory - standard
BIOM=/path

## ---- SCRIPT ----

# revert standard kraken run
cat $DIR/runids_$OUTNAME.txt | parallel cut -f1-3,6-8 $OLDDIR/{}_krakenst.txt \> $NEWDIR/{}_stdreport.txt

echo "standard done"

# revert fungi kraken run
cat $DIR/runids_$OUTNAME.txt | parallel cut -f1-3,6-8 $OLDDIRF/{}_krakenfung.txt \> $NEWDIRF/{}_fungreport.txt

echo "fungi done"

## make bioms out of old reports
kraken-biom -o $BIOM/std_$OUTNAME.biom --fmt json $NEWDIR/*.txt

# fungi
kraken-biom -o $BIOM/fung_$OUTNAME.biom --fmt json $NEWDIRF/*.txt

echo "bioms done"
