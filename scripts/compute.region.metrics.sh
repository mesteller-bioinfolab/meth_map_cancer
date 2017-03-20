#!/bin/bash

if [ $# -lt 4 ]
then
        echo
        echo "Takes two tabix smoothed methylation files and a bed file of regions (dmr)"
        echo
        echo "... and returns a bed file like the input one plus some area metrics"
        echo
        echo "    output: region_files + total_area max_difference N_CpG average_difference"
        echo
        echo "usage:"
        echo
        echo "     compute.real.dmr.metrics.sh <sample1_methylation.gz> <sample2_methylation.gz> <region.bed> <summary.script>"
        echo 
        exit 1
fi


# get arguments
smp1=$1
smp2=$2
region=$3
script=$4

# get number of fiels of region file
Nfields=`awk -F "\t" '{print NF; exit}' $region`

# do the monkey
paste <(tabix -B $smp1 $region) <(tabix -B $smp2 $region) | bedtools intersect -wao -a $region -b - | $script -v n=$Nfields
