#!/bin/bash

# get arguments
genename=$1
outPref=$2
geneannotation=$3
exprDir=$4

##  Get methylation

# get promoter coordinates
tmpProm=`mktemp`
grep -i $genename\
 $geneannotation |\
 bedtools groupby -g 1,2,3,6 -c 5 -o distinct | tr "\t" " " >\
 $tmpProm

# for each transcript, get methylation
while read line
do
# get info
lineARR=($line)
region=${lineARR[0]}:${lineARR[1]}-${lineARR[2]}
gene=${lineARR[3]}
transcript=${lineARR[4]}

# for each tcga tissue
for file in /home/user/bsdata/tcga_meth/*.bed.gz
do
# put sample id, average methyltion and transcript into a file
paste <(zcat $file | head -n 1 | cut -f 5- | tr "\t" "\n")\
 <(tabix $file $region | cut -f 5- | awk '{for(i=1; i<= NF; i ++){if($i != "NA"){m[i] += $i; n[i] += 1}}}END{for(i=1; i <= NF; i++) {if(n[i] == 0){ out="NA"}else{out=m[i]/n[i]}; print out}}') |\
 awk -v tt=$transcript '{OFS="\t"; print $0, tt}'
done
done < $tmpProm > $outPref.$genename.meth.tab

# remove tpm file
rm $tmpProm

##  Get expression

# for each tcga tissue
for file in $exprDir/*.rsem
do
# put sample id, expression value and gene ID
cat <(head -n 1 $file) <(grep -i $genename $file) |\
 awk '{OFS="\t"; if(gsub("_", "_", $1) != 0){for(i=2; i<= NF; i ++) sample[i] = $i}else{for(i=2; i<= NF; i ++) print sample[i], $i, $1}}'
done > $outPref.$genename.expr.tab






