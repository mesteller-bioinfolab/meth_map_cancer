#!/bin/bash

# get arguments
smp=$1
chr=$2

# set sample name and location
name=`basename $smp`
loc=`dirname $smp`

# do!
hmr <(tabix $smp $chr | awk 'BEGIN{OFS="\t"}{t = $5 + $6; if(t > 0) print $1, $2, $2 + 1, "CpG:"t, $6/t, "+"}') >\
 $loc/$name.hmr.chr$chr
