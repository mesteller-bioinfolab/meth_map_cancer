#!/bin/bash

# get arguments
genename=$1
outPref=$2
tissuetable=$3
geneannotation=$4

##  Get methylation

# get promoter coordinates
tmpProm=`mktemp`

for file in /home/user/bsdata/tcga_expr/*.rsem
do
bssample=$(grep $(basename $file .rsem) $tissuetable | cut -d " " -f 2)
bedtools intersect -wo -b <(grep $genename $geneannotation) -a /home/user/bsdata/hmr/$bssample.common.hmr.bed |\
 bedtools groupby -g 1,2,3,17 -c 16 -o distinct > $tmpProm.$(basename $file .rsem)
done

# for each tcga tissue
for file in /home/user/bsdata/tcga_meth/*.bed.gz
do

# for each transcript, get methylation
while read line
do
# get info
lineARR=($line)
region=${lineARR[0]}:${lineARR[1]}-${lineARR[2]}
gene=${lineARR[3]}
transcript=${lineARR[4]}

# put sample id, average methyltion and transcript into a file
paste <(zcat $file | head -n 1 | cut -f 5- | tr "\t" "\n")\
 <(tabix $file $region | cut -f 5- | awk '{for(i=1; i<= NF; i ++){if($i != "NA"){m[i] += $i; n[i] += 1}}}END{for(i=1; i <= NF; i++) {if(n[i] == 0){ out="NA"}else{out=m[i]/n[i]}; print out}}') |\
 awk -v tt=$transcript '{OFS="\t"; print $0, tt}'
done < $tmpProm.$(basename $file .bed.gz | awk '{print toupper($1)}')
done > $outPref.$genename.meth.tab

# remove tpm file
rm $tmpProm*

##  Get expression

# for each tcga tissue
for file in /home/user/bsdata/tcga_expr/*.rsem
do
# put sample id, expression value and gene ID
grep -i -e gene_id -e $genename $file |\
 awk '{OFS="\t"; if($1 == "gene_id"){for(i=2; i<= NF; i ++) sample[i] = $i}else{for(i=2; i<= NF; i ++) print sample[i], $i, $1}}'
done > $outPref.$genename.expr.tab






