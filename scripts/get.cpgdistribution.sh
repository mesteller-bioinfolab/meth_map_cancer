#!/bin/bash

# get input file
infile=$1
outdir=$2

# get sample name
name=`basename $infile .smooth.bed.gz`

zcat $infile  | awk '{u=$5; m=$6; t=u+m; if(t >= 10){ a[int(m/t*50+.5)/50] += 1; b[int($7*50+.5)/50] += 1}}END{for(i=0; i<=50; i++) print i/50, a[i/50], b[i/50]}' > $outdir/$name.cpgdensity.txt

