#!/bin/bash

# get samples
samples=/home/user/bsdata/*smooth.bed.gz

# make results directory
OUTDIR=/home/user/bsdata/hmr
mkdir -p $OUTDIR

# set wd
cd /home/user/bsdata/hmr/

#########################
######      Variant sites
#########################

##  Compute distribution of SD per CpG site
# normal samples
for region in `seq 1 22`
do
command=`for i in G145 W145 617N NY108 22A B11-27235 CP10 PL5 Laia
do
	echo "<(tabix CpG_context_$i.smooth.bed.gz $region | cut -f 7)"
done | tr "\n" " "`

echo "eval paste $command | awk -v chr=$region '{for(i=1; i<= NF; i++){if(\$i != \".\") {m += \$i; n += 1}}; for(i=1; i<=NF; i++){if(\$i != \".\"){a += (\$i - m/n)^2}}; if(n > 5) out[int(sqrt(a/n)*1000 +.5)/1000] += 1; a=0; m=0; n=0}END{for(j in out) print chr, j, out[j]}'"
done | parallel --gnu > /home/user/bsdata/cpg.sd.distribution.normals.txt

# tumor samples
for region in `seq 1 22`
do
command=`for i in SH5YSY U87MG H1437 H1672 H157 M673 M697 22B 22C 468PT 468LN PC3 22RV1
do
	echo "<(tabix CpG_context_$i.smooth.bed.gz $region | cut -f 7)"
done | tr "\n" " "`

echo "eval paste $command | awk -v chr=$region '{for(i=1; i<= NF; i++){if(\$i != \".\") {m += \$i; n += 1}}; for(i=1; i<=NF; i++){if(\$i != \".\"){a += (\$i - m/n)^2}}; if(n > 5) out[int(sqrt(a/n)*1000 +.5)/1000] += 1; a=0; m=0; n=0}END{for(j in out) print chr, j, out[j]}'"
done | parallel --gnu > /home/user/bsdata/cpg.sd.distribution.cancers.txt

##  Compute distribution of SD per CpG site
# normal samples
for region in `seq 1 22`
do
command=`for i in G145 W145 617N NY108 22A B11-27235 CP10 PL5 Laia
do
	echo "<(tabix CpG_context_$i.smooth.bed.gz $region | cut -f 7)"
done | tr "\n" " "`

echo "eval paste $command | awk -v chr=$region '{for(i=1; i<= NF; i++){if(\$i != \".\") {m += \$i; n += 1}}; for(i=1; i<=NF; i++){if(\$i != \".\"){a += (\$i - m/n)^2}}; if(n > 5) out[int(sqrt(a/n)*1000 +.5)/1000\" \"int((m/n)*1000 +.5)/1000] += 1; a=0; m=0; n=0}END{for(j in out) print chr, j, out[j]}'"
done | parallel --gnu > /home/user/bsdata/cpg.sdandmethylation.distribution.normals.txt

# tumor samples
for region in `seq 1 22`
do
command=`for i in SH5YSY U87MG H1437 H1672 H157 M673 M697 22B 22C 468PT 468LN PC3 22RV1
do
	echo "<(tabix CpG_context_$i.smooth.bed.gz $region | cut -f 7)"
done | tr "\n" " "`

echo "eval paste $command | awk -v chr=$region '{for(i=1; i<= NF; i++){if(\$i != \".\") {m += \$i; n += 1}}; for(i=1; i<=NF; i++){if(\$i != \".\"){a += (\$i - m/n)^2}}; if(n > 5) out[int(sqrt(a/n)*1000 +.5)/1000\" \"int((m/n)*1000 +.5)/1000] += 1; a=0; m=0; n=0}END{for(j in out) print chr, j, out[j]}'"
done | parallel --gnu > /home/user/bsdata/cpg.sdandmethylation.distribution.cancers.txt


##  Compute distribution of SD per CpG site, only inside CGI
region=/home/user/resources/human/hg19.ucsc.cgi.nochr.bed

# normal samples
command=`for i in G145 W145 617N NY108 22A B11-27235 CP10 PL5 Laia
do
	echo "<(tabix -B CpG_context_$i.smooth.bed.gz $region | cut -f 7)"
done | tr "\n" " "`

eval paste $command | awk -v chr=$region '{for(i=1; i<= NF; i++){if($i != ".") {m += $i; n += 1}}; for(i=1; i<=NF; i++){if($i != "."){a += ($i - m/n)^2}}; if(n > 5) out[int(sqrt(a/n)*1000 +.5)/1000" "int((m/n)*1000 +.5)/1000] += 1; a=0; m=0; n=0}END{for(j in out) print j, out[j]}' > /home/user/bsdata/cpg.sdandmethylation.distribution.normals.CGI.txt

# tumor samples
command=`for i in SH5YSY U87MG H1437 H1672 H157 M673 M697 22B 22C 468PT 468LN PC3 22RV1
do
	echo "<(tabix -B CpG_context_$i.smooth.bed.gz $region | cut -f 7)"
done | tr "\n" " "`

eval paste $command | awk -v chr=$region '{for(i=1; i<= NF; i++){if($i != ".") {m += $i; n += 1}}; for(i=1; i<=NF; i++){if($i != "."){a += ($i - m/n)^2}}; if(n > 5) out[int(sqrt(a/n)*1000 +.5)/1000" "int((m/n)*1000 +.5)/1000] += 1; a=0; m=0; n=0}END{for(j in out) print j, out[j]}' > /home/user/bsdata/cpg.sdandmethylation.distribution.cancers.CGI.txt



##  Get raw SD per position
# normal samples
for region in `seq 1 22`
do
command=`for i in G145 W145 617N NY108 22A B11-27235 CP10 PL5 Laia
do
	echo "<(tabix CpG_context_$i.smooth.bed.gz $region | cut -f 1,2,7)"
done | tr "\n" " "`

echo "eval paste $command | awk 'BEGIN{OFS=\"\t\"}{for(i=1; i<= NF/3; i++){if(\$(3*i) != \".\") {m += \$(3*i); n += 1}}; for(i=1; i<=NF/3; i++){if(\$(3*i) != \".\"){a += (\$(3*i) - m/n)^2}}; print \$1, \$2, \$2+1, a, m, n; a=0; m=0; n=0}'"
done | parallel --gnu > /home/user/bsdata/cpg.sd.perposition.normals.bed

# compress and index
bgzip /home/user/bsdata/cpg.sd.perposition.normals.bed
tabix -p bed /home/user/bsdata/cpg.sd.perposition.normals.bed.gz

# get GpG in CGI islands and outside
tabix -B cpg.sd.perposition.normals.bed.gz \
 <(bedtools complement -i \
 <(sed 's/^chr//1' /home/user/resources/human/hg19.ucsc.cgi.bed) -g /home/user/ref/human/hg19.sizes) >\
 cpg.sd.perposition.normals.nonCGI.bed

tabix -B cpg.sd.perposition.normals.bed.gz \
 <(sed 's/^chr//1' /home/user/resources/human/hg19.ucsc.cgi.bed) |\
 awk 'BEGIN{OFS="\t"}{if($4 == "") $4=0; if($5 == "") $5 = 0; if($6 == "") $6 = 0; print}' |\
 bedtools sort -i - > cpg.sd.perposition.normals.CGI.bed


# tumor samples
for region in `seq 1 22`
do
command=`for i in SH5YSY U87MG H1437 H1672 H157 M673 M697 22B 22C 468PT 468LN PC3 22RV1
do
	echo "<(tabix CpG_context_$i.smooth.bed.gz $region | cut -f 1,2,7)"
done | tr "\n" " "`

echo "eval paste $command | awk 'BEGIN{OFS=\"\t\"}{for(i=1; i<= NF/3; i++){if(\$(3*i) != \".\") {m += \$(3*i); n += 1}}; for(i=1; i<=NF/3; i++){if(\$(3*i) != \".\"){a += (\$(3*i) - m/n)^2}}; print \$1, \$2, \$2+1, a, m, n; a=0; m=0; n=0}'"
done | parallel --gnu > /home/user/bsdata/cpg.sd.perposition.cancers.bed

# compress and index
bgzip /home/user/bsdata/cpg.sd.perposition.cancers.bed
tabix -p bed /home/user/bsdata/cpg.sd.perposition.cancers.bed.gz

# get GpG in CGI islands and outside
tabix -B cpg.sd.perposition.cancers.bed.gz \
 <(bedtools complement -i \
 <(sed 's/^chr//1' /home/user/resources/human/hg19.ucsc.cgi.bed) -g /home/user/ref/human/hg19.sizes) >\
 cpg.sd.perposition.cancers.nonCGI.bed

tabix -B cpg.sd.perposition.cancers.bed.gz \
 <(sed 's/^chr//1' /home/user/resources/human/hg19.ucsc.cgi.bed) |\
 awk 'BEGIN{OFS="\t"}{if($4 == "") $4=0; if($5 == "") $5 = 0; if($6 == "") $6 = 0; print}' |\
 bedtools sort -i - > cpg.sd.perposition.cancers.CGI.bed

##  Get comparison per site
# SD
for i in `seq 1 22`
do
echo "eval paste <(tabix cpg.sd.perposition.normals.bed.gz $i) <(tabix cpg.sd.perposition.cancers.bed.gz $i) | awk '\$6 >= 5 && \$12 >= 5{sdN=int(sqrt(\$4/\$6)*1000 + .5)/1000; sdC=int(sqrt(\$10/\$12)*1000 + .5)/1000; out[sdN\" \"sdC] += 1}END{for(i in out) print i, out[i]}'"
done | parallel --gnu | awk '{out[$1" "$2] += $3}END{for(i in out) print i, out[i]}' > cpg.sd.perposition.txt

# methylation
for i in `seq 1 22`
do
echo "eval paste <(tabix cpg.sd.perposition.normals.bed.gz $i) <(tabix cpg.sd.perposition.cancers.bed.gz $i) | awk '\$6 >= 5 && \$12 >= 5{sdN=int((\$5/\$6)*1000 + .5)/1000; sdC=int((\$11/\$12)*1000 + .5)/1000; out[sdN\" \"sdC] += 1}END{for(i in out) print i, out[i]}'"
done | parallel --gnu | awk '{out[$1" "$2] += $3}END{for(i in out) print i, out[i]}' > cpg.methylation.perposition.txt

# only CGI sites
paste cpg.sd.perposition.normals.CGI.bed cpg.sd.perposition.cancers.CGI.bed |\
 awk '$6 >= 5 && $12 >= 5{sdN=int(sqrt($4/$6)*1000 + .5)/1000; sdC=int(sqrt($10/$12)*1000 + .5)/1000; out[sdN" "sdC] += 1}END{for(i in out) print i, out[i]}' >\
 cpg.sd.perposition.onlyCGI.txt

paste cpg.sd.perposition.normals.CGI.bed cpg.sd.perposition.cancers.CGI.bed |\
 awk '$6 >= 5 && $12 >= 5{sdN=int(($5/$6)*1000 + .5)/1000; sdC=int(($11/$12)*1000 + .5)/1000; out[sdN" "sdC] += 1}END{for(i in out) print i, out[i]}' >\
 cpg.methylation.perposition.onlyCGI.txt



##  snake plot
# normal samples
for region in `seq 1 22`
do
command=`for i in G145 W145 617N NY108 22A B11-27235 CP10 PL5 Laia
do
	echo "<(tabix CpG_context_$i.smooth.bed.gz $region | cut -f 7)"
done | tr "\n" " "`

echo "eval paste $command | gawk -v chr=$region '{for(i=1; i<= NF; i++){if(\$i != \".\") {m[i] = \$i}}; n=asort(m); if (n % 2) { median = m[(n + 1) / 2]} else { median = (m[(n / 2)] + m[(n / 2) + 1]) / 2.0};low=(m[int(5*n/100)]+m[int(5*n/100) + 1]) / 2.0; high = (m[int(95*n/100)] + m[int(95*n/100) + 1]) / 2.0; split(\"\", m); if(n > 5) {print low, median, high}}'"
done | parallel --gnu > /home/user/bsdata/cpg.snake90.normals.txt

for region in `seq 1 22`
do
command=`for i in SH5YSY U87MG H1437 H1672 H157 M673 M697 22B 22C 468PT 468LN PC3 22RV1
do
	echo "<(tabix CpG_context_$i.smooth.bed.gz $region | cut -f 7)"
done | tr "\n" " "`

echo "eval paste $command | gawk -v chr=$region '{for(i=1; i<= NF; i++){if(\$i != \".\") {m[i] = \$i}}; n=asort(m); if (n % 2) { median = m[(n + 1) / 2]} else { median = (m[(n / 2)] + m[(n / 2) + 1]) / 2.0};low=(m[int(5*n/100)]+m[int(5*n/100) + 1]) / 2.0; high = (m[int(95*n/100)] + m[int(95*n/100) + 1]) / 2.0; split(\"\", m); if( n > 5) print low, median, high}'"
done | parallel --gnu > /home/user/bsdata/cpg.snake90.cancers.txt

##################
###### Compute PCA
##################

# make results directory
mkdir -p /home/user/bsdata/pca

##  For each chormososme, get raw CpG count of all samples in a file
for chr in `seq 1 22` X
do
	# compose command
	command="paste `for i in ../CpG_context_*smooth.bed.gz; do echo "<(tabix $i $chr | cut -f 5,6 )"; done | tr "\n" " "`"

	# evaluate command & filter out by coverage & compute methylation & filter out by proportion of samples missing
	eval $command |\
     awk -v cove=4 -v pna=.8\
     '{p=0; t=$1 + $2; if(t > cove) {out=$2/t; p += 1}else{out="."}; for(i=2; i<= NF/2; i++) {t=$(i*2 - 1) + $(i*2); if(t > cove) {out=out"\t"$(i*2)/t; p += 1}else{out=out"\t."}}; if(p/(NF/2) > pna) print out}' >\
     ../pca/CpG.raw4cove.allsamples.chr$chr.tab
done

# get sample names in order
for i in ../CpG_context_*smooth.bed.gz
do
 basename $i .smooth.bed.gz | sed 's/CpG_context_//1'
done > ../pca/sample.names

# perform PCA
scripts/pca.r




##################
###### Compute HMR
##################

# compute HMR per sample
for sample in $samples
do
	# get sample name
	name=`basename $sample .smooth.bed.gz`

    # get chromosomes
	chromosomes=(`tabix -l $sample`)
	chromosomes=${chromosomes[@]:0:23}

	# compute HMR
	time echo $chromosomes | sed 's/ /\n/g' | parallel -j0 scripts/hmr.frommethylationbed.sh $sample {}

	# merge all chr in one file
	rm $OUTDIR/$name.hmr.bed

	for i in $chromosomes
	do
		cat $sample.hmr.chr$i >> $OUTDIR/$name.hmr.bed
	done
	
	rm $sample.hmr.chr*

	# compress & index
	bgzip -c $OUTDIR/$name.hmr.bed > $OUTDIR/$name.hmr.bed.gz
	tabix -p bed $OUTDIR/$name.hmr.bed.gz
done

##################
###### HMR metrics
##################


# compute HMR metrics (summary) per sample
for sample in $samples
do
	# get sample name
	name=`basename $sample .smooth.bed.gz`

	# compute hmr metrics (summary)
	time scripts/get.average.arrayregion.v2.sh $sample $OUTDIR/$name.hmr.bed scripts/hmr.summary.awk > $OUTDIR/$name.hmr.metrics.bed

	# compress & index
	bgzip -c $OUTDIR/$name.hmr.metrics.bed > $OUTDIR/$name.hmr.metrics.bed.gz
	tabix -p bed $OUTDIR/$name.hmr.metrics.bed.gz
done

# the correct way (removing CpG sites without coverage)
for sample in $samples
do
	# get sample name
	name=`basename $sample .smooth.bed.gz`

	# compute hmr metrics (summary)
	echo "tabix -B $sample $OUTDIR/$name.hmr.bed | awk '\$7 != \".\" {print}' | bedtools intersect -wo -a $OUTDIR/$name.hmr.bed -b - | bedtools groupby -g 1,2,3,4,5,6 -c 11,11,12,13,13 -o count,sum,sum,sum,mean > $OUTDIR/$name.hmr.metrics2.bed"
done | parallel --gnu -j 40

#######################
###### Genomic features
#######################

# split per genomic feature
mkdir $OUTDIR/genomicFeatures

# get samples again
samples=*.hmr.metrics.bed

# intersection with genomic features
echo $samples | sed 's/ /\n/g' | xargs -n 1 -P 40 -I {} /home/evidal/scripts/bsdata/split.hmr.genomicFeature.sh {} $OUTDIR/genomicFeatures


####################
###### HMR in Normal
####################

### consensus HMR

# get all bases present in any normal HMR
time for i in `awk '$3 ~ /ormal/ && $1 != "PL5" && $1 != "NB" && $1 != "Julia" && $1 != "103"  && $1 != "W145" {ORS=" "; print $1}'  /home/user/bsdata/samples.info.txt`
do
cat CpG_context_$i.hmr.bed
done | bedtools sort -i - | bedtools merge -i - > normal.all.hmr.bed


# get counts per sample of all normal hmr bases
awk '$3 ~ /ormal/ && $1 != "PL5" && $1 != "NB" && $1 != "Julia" && $1 != "103"  && $1 != "W145"  {print $1}'  /home/user/bsdata/samples.info.txt |
xargs -n 1 -P 10 -I {} sh -c\
 "tabix -B /home/user/bsdata/CpG_context_{}.gz normal.all.hmr.bed | awk 'BEGIN{OFS=\"\t\"}{print \$1, \$2 -1, \$2}' | bedtools intersect -wao -a - -b CpG_context_{}.hmr.bed | cut -f 1-3,10 > {}.allnormal.counts.bed"

# sum up all normal sample counts per base
paste *allnormal.counts.bed | awk 'BEGIN{OFS="\t"}{for(i=1; i <= NF/4; i++){a += $(i*4)} print $1, $2, $3, a; a = 0}' > normal.all.hmr.counts.bed

# get consensus hmr
scripts/consensus.normal.hmr.r

# merge consensus hmr closer than 100 bases
bedtools merge -d 100 -i normal.common.hmr.bed > normal.common.hmr.merged100b.bed
bedtools merge -d 100 -i normal.frequent.hmr.bed > normal.frequent.hmr.merged100b.bed
bedtools intersect -v -a normal.frequent.hmr.merged100b.bed -b normal.common.hmr.merged100b.bed > normal.frequent2.hmr.merged100b.bed
mv normal.frequent2.hmr.merged100b.bed normal.frequent.hmr.merged100b.bed

# make frequent HMRs presence-absence matrix
# get normal samples
normals=`awk '$3 ~ /ormal/ && $1 != "PL5" && $1 != "NB" && $1 != "Julia" && $1 != "103"  && $1 != "W145" {ORS=" "; print $1}'  /home/user/bsdata/samples.info.txt`

#overlap promoters with HMRs
for i in $normals
do
bedtools intersect -c -a normal.frequent.hmr.merged100b.bed -b CpG_context_$i.hmr.bed | cut -f 4 > $i.frequent.hmr.matrix.txt
done

paste *.frequent.hmr.matrix.txt > frequent.hmr.matrix
ls *.frequent.hmr.matrix.txt > frequent.hmr.matrix.names


# make frequent and tissue specific HMRs presence-absence matrix
# merge all tissue specific HMRs in a file
for i in $normals
do
awk -v t=$i '{OFS="\t"; $4 = t"_"NR; print}' $i.hmr.specific.bed
done | bedtools sort -i - > all.specific.bed
cp all.specific.bed normal.specific.hmr.bed

#overlap promoters with HMRs
for i in $normals
do
cat all.specific.bed normal.frequent.hmr.merged100b.bed | cut -f 1-3 | bedtools sort -i - | bedtools intersect -c -a - -b CpG_context_$i.hmr.bed | cut -f 4 > $i.frequentandspecific.hmr.matrix.txt
done

paste *.frequentandspecific.hmr.matrix.txt > frequentandspecific.hmr.matrix
ls *.frequentandspecific.hmr.matrix.txt > frequentandspecific.hmr.matrix.names


# make tissue specific HMRs presence-absence matrix
#overlap promoters with HMRs
for i in $normals
do
cut -f 1-3 all.specific.bed | bedtools intersect -c -a - -b CpG_context_$i.hmr.bed | cut -f 4 > $i.specific.hmr.matrix.txt
done

paste *.specific.hmr.matrix.txt > specific.hmr.matrix
ls *.specific.hmr.matrix.txt > specific.hmr.matrix.names

# make all HMRs presence-absence matrix
#overlap promoters with HMRs
for i in $normals
do
bedtools intersect -c -a normal.all.hmr.bed -b CpG_context_$i.hmr.bed | cut -f 4 > $i.all.hmr.matrix.txt
done

paste *.all.hmr.matrix.txt > all.hmr.matrix
ls *.all.hmr.matrix.txt > all.hmr.matrix.names

# make all but common HMRs presence-absence matrix
#overlap promoters with HMRs
for i in $normals
do
bedtools intersect -v -a normal.all.hmr.bed -b normal.common.hmr.merged100b.bed | bedtools intersect -c -a - -b CpG_context_$i.hmr.bed | cut -f 4 > $i.allbutcommon.hmr.matrix.txt
done

paste *.allbutcommon.hmr.matrix.txt > allbutcommon.hmr.matrix
ls *.allbutcommon.hmr.matrix.txt > allbutcommon.hmr.matrix.names

# remove sexual chromosomes
# make frequent HMRs presence-absence matrix
# get normal samples
normals=`awk '$3 ~ /ormal/ && $1 != "PL5" && $1 != "NB" && $1 != "Julia" && $1 != "103"  && $1 != "W145" {ORS=" "; print $1}'  /home/user/bsdata/samples.info.txt`

#overlap promoters with HMRs
for i in $normals
do
bedtools intersect -c -a normal.frequent.hmr.merged100b.bed -b CpG_context_$i.hmr.bed | grep -v "^X" | cut -f 4 > $i.frequent.hmr.matrixNOX.txt
done

paste *.frequent.hmr.matrixNOX.txt > frequent.hmr.matrixNOX
ls *.frequent.hmr.matrixNOX.txt > frequent.hmr.matrixNOX.names


# make frequent and tissue specific HMRs presence-absence matrix
# merge all tissue specific HMRs in a file
rm all.specific.bed
cat *specific.bed | bedtools sort -i - > all.specific.bed

#overlap promoters with HMRs
for i in $normals
do
cat all.specific.bed normal.frequent.hmr.merged100b.bed | cut -f 1-3 | bedtools sort -i - | bedtools intersect -c -a - -b CpG_context_$i.hmr.bed | grep -v "^X" | cut -f 4 > $i.frequentandspecific.hmr.matrixNOX.txt
done

paste *.frequentandspecific.hmr.matrixNOX.txt > frequentandspecific.hmr.matrixNOX
ls *.frequentandspecific.hmr.matrixNOX.txt > frequentandspecific.hmr.matrixNOX.names


# make tissue specific HMRs presence-absence matrix
#overlap promoters with HMRs
for i in $normals
do
cut -f 1-3 all.specific.bed | bedtools intersect -c -a - -b CpG_context_$i.hmr.bed | grep -v "^X" | cut -f 4 > $i.specific.hmr.matrixNOX.txt
done

paste *.specific.hmr.matrixNOX.txt > specific.hmr.matrixNOX
ls *.specific.hmr.matrixNOX.txt > specific.hmr.matrixNOX.names

# make all HMRs presence-absence matrix
#overlap promoters with HMRs
for i in $normals
do
bedtools intersect -c -a normal.all.hmr.bed -b CpG_context_$i.hmr.bed | grep -v "^X" | cut -f 4 > $i.all.hmr.matrixNOX.txt
done

paste *.all.hmr.matrixNOX.txt > all.hmr.matrixNOX
ls *.all.hmr.matrixNOX.txt > all.hmr.matrixNOX.names

# get genes with common HMR in their promoters
# get only coging promoters
coding=/home/user/resources/human/gencode.v16.alltranscripts.coding.prom.nochr.bed

bedtools intersect -u -a $coding -b normal.common.hmr.merged100b.bed > genomicFeatures/promoters.overlapping.common.hmr.bed
bedtools intersect -u -a $coding -b normal.frequent.hmr.merged100b.bed > genomicFeatures/promoters.overlapping.frequent.hmr.bed


# get consensus per sample
time for i in `awk '$3 ~ /ormal/ && $1 != "PL5" && $1 != "NB" && $1 != "Julia" && $1 != "103"  && $1 != "W145" {ORS=" "; print $1}'  /home/user/bsdata/samples.info.txt`
do
bedtools intersect -u -a CpG_context_$i.hmr.metrics.bed -b normal.common.hmr.merged100b.bed > $i.common.hmr.bed
bedtools intersect -u -a CpG_context_$i.hmr.metrics.bed -b normal.frequent.hmr.merged100b.bed > $i.frequent.hmr.bed
done

# intersection with genomic features                                                                                                                                                                                                        
/home/evidal/scripts/bsdata/split.hmr.genomicFeature.sh normal.common.hmr.merged100b.bed $OUTDIR/genomicFeatures
/home/evidal/scripts/bsdata/split.hmr.genomicFeature.sh normal.frequent.hmr.merged100b.bed $OUTDIR/genomicFeatures


# get samples again
samples=*.common.hmr.bed

# intersection with genomic features                                                                                                                                                                                                        
echo $samples | sed 's/ /\n/g' | xargs -n 1 -P 40 -I {} /home/evidal/scripts/bsdata/split.hmr.genomicFeature.sh {} $OUTDIR/genomicFeatures

# get samples again
samples=*.frequent.hmr.bed

# intersection with genomic features                                                                                                                                                                                                        
echo $samples | sed 's/ /\n/g' | xargs -n 1 -P 40 -I {} /home/evidal/scripts/bsdata/split.hmr.genomicFeature.sh {} $OUTDIR/genomicFeatures

###########################
####### tissue specific HMR
###########################

# get tissue specific HMR (not present in any other tissue)
scripts/tissue.specific.r

# compress and index for quick access
for i in *.hmr.specific.bed
do
bgzip -c $i > $i.gz
tabix -p bed $i.gz
done

# get numbers
allhmrs=/home/user/bsdata/hmr/normal.all.hmr.bed
commonhmrs=/home/user/bsdata/hmr/normal.common.hmr.merged100b.bed
frequenthmrs=/home/user/bsdata/hmr/normal.frequent.hmr.merged100b.bed
specifichmrs=/home/user/bsdata/hmr/normal.specific.hmr.bed
promoters=/home/user/resources/human/gencode.v16.alltranscripts.prom.nochr.bed
transcripts=/home/user/resources/human/gencode.v16.alltranscripts.nochr.bed
intergenic=/home/user/resources/human/gencode.v16.alltranscripts.intergenic.nochr.bed
array=/home/user/resources/human/450k.manifest.1to22andX.onlypos.bed


# number of HMR per type and coverage
echo hmr genome cover prom N > hmr.bsseqandarray.numbers
echo all all bsseq all `cat $allhmrs | wc -l` >> hmr.bsseqandarray.numbers
echo all all bsseq promoter `cat $allhmrs | bedtools intersect -u -a - -b $promoters | wc -l` >> hmr.bsseqandarray.numbers
echo all all array all `cat $allhmrs | bedtools intersect -u -a - -b $array | wc -l` >> hmr.bsseqandarray.numbers
echo all all array promoter `cat $allhmrs | bedtools intersect -u -a - -b $promoters | bedtools intersect -u -a - -b $array | wc -l` >> hmr.bsseqandarray.numbers

echo all autosomes bsseq all `grep -v X $allhmrs | wc -l` >> hmr.bsseqandarray.numbers
echo all autosomes bsseq promoter `grep -v X $allhmrs | bedtools intersect -u -a - -b $promoters | wc -l` >> hmr.bsseqandarray.numbers
echo all autosomes array all `grep -v X $allhmrs | bedtools intersect -u -a - -b $array | wc -l` >> hmr.bsseqandarray.numbers
echo all autosomes array promoter `grep -v X $allhmrs | bedtools intersect -u -a - -b $promoters | bedtools intersect -u -a - -b $array | wc -l` >> hmr.bsseqandarray.numbers

echo common all bsseq all `cat $allhmrs | bedtools intersect -u -a - -b $commonhmrs | wc -l` >> hmr.bsseqandarray.numbers
echo common all bsseq promoter `cat $allhmrs | bedtools intersect -u -a - -b $commonhmrs | bedtools intersect -u -a - -b $promoters | wc -l` >> hmr.bsseqandarray.numbers
echo common all array all `cat $allhmrs | bedtools intersect -u -a - -b $commonhmrs | bedtools intersect -u -a - -b $array | wc -l` >> hmr.bsseqandarray.numbers
echo common all array promoter `cat $allhmrs | bedtools intersect -u -a - -b $commonhmrs | bedtools intersect -u -a - -b $promoters | bedtools intersect -u -a - -b $array | wc -l` >> hmr.bsseqandarray.numbers

echo common autosomes bsseq all `grep -v X $allhmrs | bedtools intersect -u -a - -b $commonhmrs | wc -l` >> hmr.bsseqandarray.numbers
echo common autosomes bsseq promoter `grep -v X $allhmrs | bedtools intersect -u -a - -b $commonhmrs | bedtools intersect -u -a - -b $promoters | wc -l` >> hmr.bsseqandarray.numbers
echo common autosomes array all `grep -v X $allhmrs | bedtools intersect -u -a - -b $commonhmrs | bedtools intersect -u -a - -b $array | wc -l` >> hmr.bsseqandarray.numbers
echo common autosomes array promoter `grep -v X $allhmrs | bedtools intersect -u -a - -b $commonhmrs | bedtools intersect -u -a - -b $promoters | bedtools intersect -u -a - -b $array | wc -l` >> hmr.bsseqandarray.numbers

echo frequent all bsseq all `cat $allhmrs | bedtools intersect -u -a - -b $frequenthmrs | wc -l` >> hmr.bsseqandarray.numbers
echo frequent all bsseq promoter `cat $allhmrs | bedtools intersect -u -a - -b $frequenthmrs | bedtools intersect -u -a - -b $promoters | wc -l` >> hmr.bsseqandarray.numbers
echo frequent all array all `cat $allhmrs | bedtools intersect -u -a - -b $frequenthmrs | bedtools intersect -u -a - -b $array | wc -l` >> hmr.bsseqandarray.numbers
echo frequent all array promoter `cat $allhmrs | bedtools intersect -u -a - -b $frequenthmrs | bedtools intersect -u -a - -b $promoters | bedtools intersect -u -a - -b $array | wc -l` >> hmr.bsseqandarray.numbers

echo frequent autosomes bsseq all `grep -v X $allhmrs | bedtools intersect -u -a - -b $frequenthmrs | wc -l` >> hmr.bsseqandarray.numbers
echo frequent autosomes bsseq promoter `grep -v X $allhmrs | bedtools intersect -u -a - -b $frequenthmrs | bedtools intersect -u -a - -b $promoters | wc -l` >> hmr.bsseqandarray.numbers
echo frequent autosomes array all `grep -v X $allhmrs | bedtools intersect -u -a - -b $frequenthmrs | bedtools intersect -u -a - -b $array | wc -l` >> hmr.bsseqandarray.numbers
echo frequent autosomes array promoter `grep -v X $allhmrs | bedtools intersect -u -a - -b $frequenthmrs | bedtools intersect -u -a - -b $promoters | bedtools intersect -u -a - -b $array | wc -l` >> hmr.bsseqandarray.numbers

echo specific all bsseq all `cat $allhmrs | bedtools intersect -u -a - -b $specifichmrs | wc -l` >> hmr.bsseqandarray.numbers
echo specific all bsseq promoter `cat $allhmrs | bedtools intersect -u -a - -b $specifichmrs | bedtools intersect -u -a - -b $promoters | wc -l` >> hmr.bsseqandarray.numbers
echo specific all array all `cat $allhmrs | bedtools intersect -u -a - -b $specifichmrs | bedtools intersect -u -a - -b $array | wc -l` >> hmr.bsseqandarray.numbers
echo specific all array promoter `cat $allhmrs | bedtools intersect -u -a - -b $specifichmrs | bedtools intersect -u -a - -b $promoters | bedtools intersect -u -a - -b $array | wc -l` >> hmr.bsseqandarray.numbers

echo specific autosomes bsseq all `grep -v X $allhmrs | bedtools intersect -u -a - -b $specifichmrs | wc -l` >> hmr.bsseqandarray.numbers
echo specific autosomes bsseq promoter `grep -v X $allhmrs | bedtools intersect -u -a - -b $specifichmrs | bedtools intersect -u -a - -b $promoters | wc -l` >> hmr.bsseqandarray.numbers
echo specific autosomes array all `grep -v X $allhmrs | bedtools intersect -u -a - -b $specifichmrs | bedtools intersect -u -a - -b $array | wc -l` >> hmr.bsseqandarray.numbers
echo specific autosomes array promoter `grep -v X $allhmrs | bedtools intersect -u -a - -b $specifichmrs | bedtools intersect -u -a - -b $promoters | bedtools intersect -u -a - -b $array | wc -l` >> hmr.bsseqandarray.numbers

# number of HMR per feature
echo hmr genome feature N > hmr.features.numbers

echo all all all `cat $allhmrs | wc -l` >> hmr.features.numbers
echo all all promoters `cat $allhmrs | bedtools intersect -u -a - -b $promoters | wc -l` >> hmr.features.numbers
echo all all intragenic `cat $allhmrs | bedtools intersect -v -a - -b $promoters | bedtools intersect -u -a - -b $transcripts | wc -l` >> hmr.features.numbers
echo all all intergenic `cat $allhmrs | bedtools intersect -v -a - -b $promoters | bedtools intersect -v -a - -b $transcripts | wc -l` >> hmr.features.numbers

echo all autosomes all `grep -v X $allhmrs | wc -l` >> hmr.features.numbers
echo all autosomes promoters `grep -v X $allhmrs | bedtools intersect -u -a - -b $promoters | wc -l` >> hmr.features.numbers
echo all autosomes intragenic `grep -v X $allhmrs | bedtools intersect -v -a - -b $promoters | bedtools intersect -u -a - -b $transcripts | wc -l` >> hmr.features.numbers
echo all autosomes intergenic `grep -v X $allhmrs | bedtools intersect -v -a - -b $promoters | bedtools intersect -v -a - -b $transcripts | wc -l` >> hmr.features.numbers

echo common all all `cat $allhmrs | bedtools intersect -u -a - -b $commonhmrs | wc -l` >> hmr.features.numbers
echo common all promoters `cat $allhmrs | bedtools intersect -u -a - -b $commonhmrs | bedtools intersect -u -a - -b $promoters | wc -l` >> hmr.features.numbers
echo common all intragenic `cat $allhmrs | bedtools intersect -u -a - -b $commonhmrs | bedtools intersect -v -a - -b $promoters | bedtools intersect -u -a - -b $transcripts | wc -l` >> hmr.features.numbers
echo common all intergenic `cat $allhmrs | bedtools intersect -u -a - -b $commonhmrs | bedtools intersect -v -a - -b $promoters | bedtools intersect -v -a - -b $transcripts | wc -l` >> hmr.features.numbers

echo common autosomes all `grep -v X $allhmrs | bedtools intersect -u -a - -b $commonhmrs | wc -l` >> hmr.features.numbers
echo common autosomes promoters `grep -v X $allhmrs | bedtools intersect -u -a - -b $commonhmrs | bedtools intersect -u -a - -b $promoters | wc -l` >> hmr.features.numbers
echo common autosomes intragenic `grep -v X $allhmrs | bedtools intersect -u -a - -b $commonhmrs | bedtools intersect -v -a - -b $promoters | bedtools intersect -u -a - -b $transcripts | wc -l` >> hmr.features.numbers
echo common autosomes intergenic `grep -v X $allhmrs | bedtools intersect -u -a - -b $commonhmrs | bedtools intersect -v -a - -b $promoters | bedtools intersect -v -a - -b $transcripts | wc -l` >> hmr.features.numbers

echo frequent all all `cat $allhmrs | bedtools intersect -u -a - -b $frequenthmrs | wc -l` >> hmr.features.numbers
echo frequent all promoters `cat $allhmrs | bedtools intersect -u -a - -b $frequenthmrs | bedtools intersect -u -a - -b $promoters | wc -l` >> hmr.features.numbers
echo frequent all intragenic `cat $allhmrs | bedtools intersect -u -a - -b $frequenthmrs | bedtools intersect -v -a - -b $promoters | bedtools intersect -u -a - -b $transcripts | wc -l` >> hmr.features.numbers
echo frequent all intergenic `cat $allhmrs | bedtools intersect -u -a - -b $frequenthmrs | bedtools intersect -v -a - -b $promoters | bedtools intersect -v -a - -b $transcripts | wc -l` >> hmr.features.numbers

echo frequent autosomes all `grep -v X $allhmrs | bedtools intersect -u -a - -b $frequenthmrs | wc -l` >> hmr.features.numbers
echo frequent autosomes promoters `grep -v X $allhmrs | bedtools intersect -u -a - -b $frequenthmrs | bedtools intersect -u -a - -b $promoters | wc -l` >> hmr.features.numbers
echo frequent autosomes intragenic `grep -v X $allhmrs | bedtools intersect -u -a - -b $frequenthmrs | bedtools intersect -v -a - -b $promoters | bedtools intersect -u -a - -b $transcripts | wc -l` >> hmr.features.numbers
echo frequent autosomes intergenic `grep -v X $allhmrs | bedtools intersect -u -a - -b $frequenthmrs | bedtools intersect -v -a - -b $promoters | bedtools intersect -v -a - -b $transcripts | wc -l` >> hmr.features.numbers

echo specific all all `cat $allhmrs | bedtools intersect -u -a - -b $specifichmrs | wc -l` >> hmr.features.numbers
echo specific all promoters `cat $allhmrs | bedtools intersect -u -a - -b $specifichmrs | bedtools intersect -u -a - -b $promoters | wc -l` >> hmr.features.numbers
echo specific all intragenic `cat $allhmrs | bedtools intersect -u -a - -b $specifichmrs | bedtools intersect -v -a - -b $promoters | bedtools intersect -u -a - -b $transcripts | wc -l` >> hmr.features.numbers
echo specific all intergenic `cat $allhmrs | bedtools intersect -u -a - -b $specifichmrs | bedtools intersect -v -a - -b $promoters | bedtools intersect -v -a - -b $transcripts | wc -l` >> hmr.features.numbers

echo specific autosomes all `grep -v X $allhmrs | bedtools intersect -u -a - -b $specifichmrs | wc -l` >> hmr.features.numbers
echo specific autosomes promoters `grep -v X $allhmrs | bedtools intersect -u -a - -b $specifichmrs | bedtools intersect -u -a - -b $promoters | wc -l` >> hmr.features.numbers
echo specific autosomes intragenic `grep -v X $allhmrs | bedtools intersect -u -a - -b $specifichmrs | bedtools intersect -v -a - -b $promoters | bedtools intersect -u -a - -b $transcripts | wc -l` >> hmr.features.numbers
echo specific autosomes intergenic `grep -v X $allhmrs | bedtools intersect -u -a - -b $specifichmrs | bedtools intersect -v -a - -b $promoters | bedtools intersect -v -a - -b $transcripts | wc -l` >> hmr.features.numbers

# number of HMR @ promoters per gene type
codingproms=/home/user/resources/human/gencode.v16.alltranscripts.coding.prom.nochr.bed
noncodingproms=/home/user/resources/human/gencode.v16.alltranscripts.noncoding.prom.nochr.bed
pseudoproms=/home/user/resources/human/gencode.v16.alltranscripts.pseudo.prom.nochr.bed
otherproms=/home/user/resources/human/gencode.v16.alltranscripts.others.prom.nochr.bed

echo hmr genome type N > hmr.promoters.numbers

echo all all all `cat $allhmrs | bedtools intersect -u -a - -b $promoters | wc -l` >> hmr.promoters.numbers
echo all all coding `cat $allhmrs | bedtools intersect -u -a - -b $codingproms | wc -l` >> hmr.promoters.numbers
echo all all noncoding `cat $allhmrs | bedtools intersect -v -a - -b $codingproms | bedtools intersect -u -a - -b $noncodingproms | wc -l` >> hmr.promoters.numbers
echo all all pseudo `cat $allhmrs | bedtools intersect -v -a - -b $codingproms | bedtools intersect -v -a - -b $noncodingproms | bedtools intersect -u -a - -b $pseudoproms | wc -l` >> hmr.promoters.numbers
echo all all others `cat $allhmrs | bedtools intersect -v -a - -b $codingproms | bedtools intersect -v -a - -b $noncodingproms | bedtools intersect -v -a - -b $pseudoproms | bedtools intersect -u -a - -b $otherproms | wc -l` >> hmr.promoters.numbers

echo all autosomes all `grep -v X $allhmrs | bedtools intersect -u -a - -b $promoters | wc -l` >> hmr.promoters.numbers
echo all autosomes coding `grep -v X $allhmrs | bedtools intersect -u -a - -b $codingproms | wc -l` >> hmr.promoters.numbers
echo all autosomes noncoding `grep -v X $allhmrs | bedtools intersect -v -a - -b $codingproms | bedtools intersect -u -a - -b $noncodingproms | wc -l` >> hmr.promoters.numbers
echo all autosomes pseudo `grep -v X $allhmrs | bedtools intersect -v -a - -b $codingproms | bedtools intersect -v -a - -b $noncodingproms | bedtools intersect -u -a - -b $pseudoproms | wc -l` >> hmr.promoters.numbers
echo all autosomes others `grep -v X $allhmrs | bedtools intersect -v -a - -b $codingproms | bedtools intersect -v -a - -b $noncodingproms | bedtools intersect -v -a - -b $pseudoproms | bedtools intersect -u -a - -b $otherproms | wc -l` >> hmr.promoters.numbers

echo common all all `cat $allhmrs | bedtools intersect -u -a - -b $commonhmrs | bedtools intersect -u -a - -b $promoters | wc -l` >> hmr.promoters.numbers
echo common all coding `cat $allhmrs | bedtools intersect -u -a - -b $commonhmrs | bedtools intersect -u -a - -b $codingproms | wc -l` >> hmr.promoters.numbers
echo common all noncoding `cat $allhmrs | bedtools intersect -u -a - -b $commonhmrs | bedtools intersect -v -a - -b $codingproms | bedtools intersect -u -a - -b $noncodingproms | wc -l` >> hmr.promoters.numbers
echo common all pseudo `cat $allhmrs | bedtools intersect -u -a - -b $commonhmrs | bedtools intersect -v -a - -b $codingproms | bedtools intersect -v -a - -b $noncodingproms | bedtools intersect -u -a - -b $pseudoproms | wc -l` >> hmr.promoters.numbers
echo common all others `cat $allhmrs | bedtools intersect -u -a - -b $commonhmrs | bedtools intersect -v -a - -b $codingproms | bedtools intersect -v -a - -b $noncodingproms | bedtools intersect -v -a - -b $pseudoproms | bedtools intersect -u -a - -b $otherproms | wc -l` >> hmr.promoters.numbers

echo common autosomes all `grep -v X $allhmrs | bedtools intersect -u -a - -b $commonhmrs | bedtools intersect -u -a - -b $promoters | wc -l` >> hmr.promoters.numbers
echo common autosomes coding `grep -v X $allhmrs | bedtools intersect -u -a - -b $commonhmrs | bedtools intersect -u -a - -b $codingproms | wc -l` >> hmr.promoters.numbers
echo common autosomes noncoding `grep -v X $allhmrs | bedtools intersect -u -a - -b $commonhmrs | bedtools intersect -v -a - -b $codingproms | bedtools intersect -u -a - -b $noncodingproms | wc -l` >> hmr.promoters.numbers
echo common autosomes pseudo `grep -v X $allhmrs | bedtools intersect -u -a - -b $commonhmrs | bedtools intersect -v -a - -b $codingproms | bedtools intersect -v -a - -b $noncodingproms | bedtools intersect -u -a - -b $pseudoproms | wc -l` >> hmr.promoters.numbers
echo common autosomes others `grep -v X $allhmrs | bedtools intersect -u -a - -b $commonhmrs | bedtools intersect -v -a - -b $codingproms | bedtools intersect -v -a - -b $noncodingproms | bedtools intersect -v -a - -b $pseudoproms | bedtools intersect -u -a - -b $otherproms | wc -l` >> hmr.promoters.numbers

echo frequent all all `cat $allhmrs | bedtools intersect -u -a - -b $frequenthmrs | bedtools intersect -u -a - -b $promoters | wc -l` >> hmr.promoters.numbers
echo frequent all coding `cat $allhmrs | bedtools intersect -u -a - -b $frequenthmrs | bedtools intersect -u -a - -b $codingproms | wc -l` >> hmr.promoters.numbers
echo frequent all noncoding `cat $allhmrs | bedtools intersect -u -a - -b $frequenthmrs | bedtools intersect -v -a - -b $codingproms | bedtools intersect -u -a - -b $noncodingproms | wc -l` >> hmr.promoters.numbers
echo frequent all pseudo `cat $allhmrs | bedtools intersect -u -a - -b $frequenthmrs | bedtools intersect -v -a - -b $codingproms | bedtools intersect -v -a - -b $noncodingproms | bedtools intersect -u -a - -b $pseudoproms | wc -l` >> hmr.promoters.numbers
echo frequent all others `cat $allhmrs | bedtools intersect -u -a - -b $frequenthmrs | bedtools intersect -v -a - -b $codingproms | bedtools intersect -v -a - -b $noncodingproms | bedtools intersect -v -a - -b $pseudoproms | bedtools intersect -u -a - -b $otherproms | wc -l` >> hmr.promoters.numbers

echo frequent autosomes all `grep -v X $allhmrs | bedtools intersect -u -a - -b $frequenthmrs | bedtools intersect -u -a - -b $promoters | wc -l` >> hmr.promoters.numbers
echo frequent autosomes coding `grep -v X $allhmrs | bedtools intersect -u -a - -b $frequenthmrs | bedtools intersect -u -a - -b $codingproms | wc -l` >> hmr.promoters.numbers
echo frequent autosomes noncoding `grep -v X $allhmrs | bedtools intersect -u -a - -b $frequenthmrs | bedtools intersect -v -a - -b $codingproms | bedtools intersect -u -a - -b $noncodingproms | wc -l` >> hmr.promoters.numbers
echo frequent autosomes pseudo `grep -v X $allhmrs | bedtools intersect -u -a - -b $frequenthmrs | bedtools intersect -v -a - -b $codingproms | bedtools intersect -v -a - -b $noncodingproms | bedtools intersect -u -a - -b $pseudoproms | wc -l` >> hmr.promoters.numbers
echo frequent autosomes others `grep -v X $allhmrs | bedtools intersect -u -a - -b $frequenthmrs | bedtools intersect -v -a - -b $codingproms | bedtools intersect -v -a - -b $noncodingproms | bedtools intersect -v -a - -b $pseudoproms | bedtools intersect -u -a - -b $otherproms | wc -l` >> hmr.promoters.numbers

echo specific all all `cat $allhmrs | bedtools intersect -u -a - -b $specifichmrs | bedtools intersect -u -a - -b $promoters | wc -l` >> hmr.promoters.numbers
echo specific all coding `cat $allhmrs | bedtools intersect -u -a - -b $specifichmrs | bedtools intersect -u -a - -b $codingproms | wc -l` >> hmr.promoters.numbers
echo specific all noncoding `cat $allhmrs | bedtools intersect -u -a - -b $specifichmrs | bedtools intersect -v -a - -b $codingproms | bedtools intersect -u -a - -b $noncodingproms | wc -l` >> hmr.promoters.numbers
echo specific all pseudo `cat $allhmrs | bedtools intersect -u -a - -b $specifichmrs | bedtools intersect -v -a - -b $codingproms | bedtools intersect -v -a - -b $noncodingproms | bedtools intersect -u -a - -b $pseudoproms | wc -l` >> hmr.promoters.numbers
echo specific all others `cat $allhmrs | bedtools intersect -u -a - -b $specifichmrs | bedtools intersect -v -a - -b $codingproms | bedtools intersect -v -a - -b $noncodingproms | bedtools intersect -v -a - -b $pseudoproms | bedtools intersect -u -a - -b $otherproms | wc -l` >> hmr.promoters.numbers

echo specific autosomes all `grep -v X $allhmrs | bedtools intersect -u -a - -b $specifichmrs | bedtools intersect -u -a - -b $promoters | wc -l` >> hmr.promoters.numbers
echo specific autosomes coding `grep -v X $allhmrs | bedtools intersect -u -a - -b $specifichmrs | bedtools intersect -u -a - -b $codingproms | wc -l` >> hmr.promoters.numbers
echo specific autosomes noncoding `grep -v X $allhmrs | bedtools intersect -u -a - -b $specifichmrs | bedtools intersect -v -a - -b $codingproms | bedtools intersect -u -a - -b $noncodingproms | wc -l` >> hmr.promoters.numbers
echo specific autosomes pseudo `grep -v X $allhmrs | bedtools intersect -u -a - -b $specifichmrs | bedtools intersect -v -a - -b $codingproms | bedtools intersect -v -a - -b $noncodingproms | bedtools intersect -u -a - -b $pseudoproms | wc -l` >> hmr.promoters.numbers
echo specific autosomes others `grep -v X $allhmrs | bedtools intersect -u -a - -b $specifichmrs | bedtools intersect -v -a - -b $codingproms | bedtools intersect -v -a - -b $noncodingproms | bedtools intersect -v -a - -b $pseudoproms | bedtools intersect -u -a - -b $otherproms | wc -l` >> hmr.promoters.numbers



# get samples again
samples=*.specific.bed

# intersection with genomic features                                                                                                                                                                                                        
echo $samples | sed 's/ /\n/g' | xargs -n 1 -P 40 -I {} /home/evidal/scripts/bsdata/split.hmr.genomicFeature.sh {} $OUTDIR/genomicFeatures


# get only coging promoters
coding=/home/user/resources/human/gencode.v16.alltranscripts.coding.prom.nochr.bed

# get normal samples
normals=`awk '$3 ~ /ormal/ && $1 != "PL5" && $1 != "NB" && $1 != "Julia" && $1 != "103"  && $1 != "W145" {ORS=" "; print $1}'  /home/user/bsdata/samples.info.txt`

#overlap promoters with specific HMRs (average smoothed methylation below .5)
for i in $normals
do
awk '$10/$7 <= .5 {print}' $i.hmr.specific.bed | bedtools intersect -u -a $coding -b - > genomicFeatures/promoters.overlapping.$i.hmr.specific.bed
awk '$10/$7 <= .5 {print}' $i.hmr.specific.bed | bedtools intersect -wao -a $coding -b - | awk '{out[$6] += $NF}END{OFS="\t"; for(i in out) print i, out[i] != 0}' > genomicFeatures/genesymbol.promoters.overlapping.$i.hmr.specific.bed
done

# count the number of different t-HMR overlapping each gene promoter (gene symbol level)
paste genomicFeatures/genesymbol.promoters.overlapping.*.hmr.specific.bed | awk 'BEGIN{OFS="\t"}{for(i=1; i <= NF/2; i++){a += $(i*2)} print $1, a; a = 0}' > genomicFeatures/genesymbol.promoters.overlapping.all.hmr.counts.bed


#promoter enrichment
> normals.prom.enrichment.txt
wc -l genomicFeatures/*.hmr.specific.bed.prom.bed | grep -v total >> normals.prom.enrichment.txt
wc -l genomicFeatures/*.hmr.specific.bed.genomic.bed | grep -v total >> normals.prom.enrichment.txt
wc -l genomicFeatures/*.common.hmr.bed.genomic.bed | grep -v total >> normals.prom.enrichment.txt
wc -l genomicFeatures/*.common.hmr.bed.prom.bed | grep -v total >> normals.prom.enrichment.txt
wc -l genomicFeatures/*.frequent.hmr.bed.genomic.bed | grep -v total >> normals.prom.enrichment.txt
wc -l genomicFeatures/*.frequent.hmr.bed.prom.bed | grep -v total >> normals.prom.enrichment.txt


#########################################################
##  SD and methylation of c & t HMRs in Normal and Cancer
#########################################################

###   c-HMRs
# set region
region=/home/user/bsdata/hmr/normal.common.hmr.merged100b.bed

##  Compute distribution of SD per CpG site
# normal samples
command=`for i in G145 W145 617N NY108 22A B11-27235 CP10 PL5 Laia
do
	echo "<(tabix -B CpG_context_$i.smooth.bed.gz $region | cut -f 7)"
done | tr "\n" " "`

time eval paste $command |\
 awk '{for(i=1; i<= NF; i++){if($i != ".") {m += $i; n += 1}}; for(i=1; i<=NF; i++){if($i != "."){a += ($i - m/n)^2}}; if(n > 5) out[int(sqrt(a/n)*1000 +.5)/1000" "int((m/n)*1000 +.5)/1000] += 1; a=0; m=0; n=0}END{for(j in out) print j,  out[j]}'  > /home/user/bsdata/hmr/cpgincommonhmr.sdandmethylation.distribution.normals.txt


# tumor samples
# normal samples
command=`for i in SH5YSY U87MG H1437 H1672 H157 M673 M697 22B 22C 468PT 468LN PC3 22RV1
do
	echo "<(tabix -B CpG_context_$i.smooth.bed.gz $region | cut -f 7)"
done | tr "\n" " "`

time eval paste $command |\
 awk '{for(i=1; i<= NF; i++){if($i != ".") {m += $i; n += 1}}; for(i=1; i<=NF; i++){if($i != "."){a += ($i - m/n)^2}}; if(n > 5) out[int(sqrt(a/n)*1000 +.5)/1000" "int((m/n)*1000 +.5)/1000] += 1; a=0; m=0; n=0}END{for(j in out) print j,  out[j]}'  > /home/user/bsdata/hmr/cpgincommonhmr.sdandmethylation.distribution.cancers.txt


##  Compute distribution of SD per region site
# normal samples
command=`for i in G145 W145 617N NY108 22A B11-27235 CP10 PL5 Laia
do
	echo "<(tabix -B CpG_context_$i.smooth.bed.gz $region | bedtools intersect -wao -a $region -b - | bedtools groupby -g 1,2,3 -c 8,9 -o sum,sum | awk '{out=\".\"; t=\\$4 + \\$5; if(t>0) out=\\$5/t; print out}')"
done | tr "\n" " "`

time eval paste $command |\
 awk '{for(i=1; i<= NF; i++){if($i != ".") {m += $i; n += 1}}; for(i=1; i<=NF; i++){if($i != "."){a += ($i - m/n)^2}}; if(n > 5) out[int(sqrt(a/n)*1000 +.5)/1000" "int((m/n)*1000 +.5)/1000] += 1; a=0; m=0; n=0}END{for(j in out) print j,  out[j]}'  > /home/user/bsdata/hmr/commonhmr.sdandmethylation.distribution.normals.txt


# tumor samples
command=`for i in SH5YSY U87MG H1437 H1672 H157 M673 M697 22B 22C 468PT 468LN PC3 22RV1
do
	echo "<(tabix -B CpG_context_$i.smooth.bed.gz $region | bedtools intersect -wao -a $region -b - | bedtools groupby -g 1,2,3 -c 8,9 -o sum,sum | awk '{out=\".\"; t=\\$4 + \\$5; if(t>0) out=\\$5/t; print out}')"
done | tr "\n" " "`

time eval paste $command |\
 awk '{for(i=1; i<= NF; i++){if($i != ".") {m += $i; n += 1}}; for(i=1; i<=NF; i++){if($i != "."){a += ($i - m/n)^2}}; if(n > 5) out[int(sqrt(a/n)*1000 +.5)/1000" "int((m/n)*1000 +.5)/1000] += 1; a=0; m=0; n=0}END{for(j in out) print j,  out[j]}'  > /home/user/bsdata/hmr/commonhmr.sdandmethylation.distribution.cancers.txt


###   t-HMRs
# set region
region=/home/user/bsdata/hmr/normal.specific.hmr.onlycoords.bed

##  Compute distribution of SD per CpG site
# normal samples
command=`for i in G145 W145 617N NY108 22A B11-27235 CP10 PL5 Laia
do
	echo "<(tabix -B CpG_context_$i.smooth.bed.gz $region | cut -f 7)"
done | tr "\n" " "`

time eval paste $command |\
 awk '{for(i=1; i<= NF; i++){if($i != ".") {m += $i; n += 1}}; for(i=1; i<=NF; i++){if($i != "."){a += ($i - m/n)^2}}; if(n > 5) out[int(sqrt(a/n)*1000 +.5)/1000" "int((m/n)*1000 +.5)/1000] += 1; a=0; m=0; n=0}END{for(j in out) print j,  out[j]}'  > /home/user/bsdata/hmr/cpginspecifichmr.sdandmethylation.distribution.normals.txt


# tumor samples
command=`for i in SH5YSY U87MG H1437 H1672 H157 M673 M697 22B 22C 468PT 468LN PC3 22RV1
do
	echo "<(tabix -B CpG_context_$i.smooth.bed.gz $region | cut -f 7)"
done | tr "\n" " "`

time eval paste $command |\
 awk '{for(i=1; i<= NF; i++){if($i != ".") {m += $i; n += 1}}; for(i=1; i<=NF; i++){if($i != "."){a += ($i - m/n)^2}}; if(n > 5) out[int(sqrt(a/n)*1000 +.5)/1000" "int((m/n)*1000 +.5)/1000] += 1; a=0; m=0; n=0}END{for(j in out) print j,  out[j]}'  > /home/user/bsdata/hmr/cpginspecifichmr.sdandmethylation.distribution.cancers.txt

##  Exclude each tissue
region=/home/user/bsdata/hmr/normal.specific.hmr.bed
# normal samples
for i in G145 W145 617N NY108 22A B11-27235 CP10 PL5 Laia
do
	command=`for j in G145 W145 617N NY108 22A B11-27235 CP10 PL5 Laia
	do
		if [[ "$j" != "$i" ]]
		then
			echo "<(tabix -B CpG_context_$j.smooth.bed.gz <(grep $i $region) | cut -f 7)"
		fi
	done | tr "\n" " "`
	echo eval paste $command
done | parallel --gnu | \
 awk '{for(i=1; i<= NF; i++){if($i != ".") {m += $i; n += 1}}; for(i=1; i<=NF; i++){if($i != "."){a += ($i - m/n)^2}}; if(n > 5) out[int(sqrt(a/n)*1000 +.5)/1000" "int((m/n)*1000 +.5)/1000] += 1; a=0; m=0; n=0}END{for(j in out) print j,  out[j]}'  > /home/user/bsdata/hmr/cpginspecifichmr.sdandmethylation.distribution.normals.excludetissue.txt



# tumor samples
for i in G145 W145 617N NY108 22A B11-27235 CP10 PL5 Laia
do
	command=`for j in SH5YSY U87MG H1437 H1672 H157 M673 M697 22B 22C 468PT 468LN PC3 22RV1
	do
	    line=($(grep $j samples.info.newnames.control.csv | tr "," " "))
        jj=${line[5]}
		if [[ "$jj" != "$i" ]]
		then
			echo "<(tabix -B CpG_context_$j.smooth.bed.gz <(grep $i $region) | cut -f 7)"
		fi
	done | tr "\n" " "`
	echo eval paste $command
done | parallel --gnu | \
 awk '{for(i=1; i<= NF; i++){if($i != ".") {m += $i; n += 1}}; for(i=1; i<=NF; i++){if($i != "."){a += ($i - m/n)^2}}; if(n > 5) out[int(sqrt(a/n)*1000 +.5)/1000" "int((m/n)*1000 +.5)/1000] += 1; a=0; m=0; n=0}END{for(j in out) print j,  out[j]}'  > /home/user/bsdata/hmr/cpginspecifichmr.sdandmethylation.distribution.cancers.excludetissue.txt



######################


##  Compute distribution of SD per region site
# normal samples
command=`for i in G145 W145 617N NY108 22A B11-27235 CP10 PL5 Laia
do
	echo "<(tabix -B CpG_context_$i.smooth.bed.gz $region | bedtools intersect -wao -a $region -b - | bedtools groupby -g 1,2,3 -c 8,9 -o sum,sum | awk '{out=\".\"; t=\\$4 + \\$5; if(t>0) out=\\$5/t; print out}')"
done | tr "\n" " "`

time eval paste $command |\
 awk '{for(i=1; i<= NF; i++){if($i != ".") {m += $i; n += 1}}; for(i=1; i<=NF; i++){if($i != "."){a += ($i - m/n)^2}}; if(n > 5) out[int(sqrt(a/n)*1000 +.5)/1000" "int((m/n)*1000 +.5)/1000] += 1; a=0; m=0; n=0}END{for(j in out) print j,  out[j]}'  > /home/user/bsdata/hmr/specifichmr.sdandmethylation.distribution.normals.txt


# tumor samples
command=`for i in SH5YSY U87MG H1437 H1672 H157 M673 M697 22B 22C 468PT 468LN PC3 22RV1
do
	echo "<(tabix -B CpG_context_$i.smooth.bed.gz $region | bedtools intersect -wao -a $region -b - | bedtools groupby -g 1,2,3 -c 8,9 -o sum,sum | awk '{out=\".\"; t=\\$4 + \\$5; if(t>0) out=\\$5/t; print out}')"
done | tr "\n" " "`

time eval paste $command |\
 awk '{for(i=1; i<= NF; i++){if($i != ".") {m += $i; n += 1}}; for(i=1; i<=NF; i++){if($i != "."){a += ($i - m/n)^2}}; if(n > 5) out[int(sqrt(a/n)*1000 +.5)/1000" "int((m/n)*1000 +.5)/1000] += 1; a=0; m=0; n=0}END{for(j in out) print j,  out[j]}'  > /home/user/bsdata/hmr/specifichmr.sdandmethylation.distribution.cancers.txt


####################
####    Smiley plots
####################

###   c-HMRs

##  Compute distribution of methylation per relative position inside c-HMRs (c-HMRs scaled)
# set region
region=/home/user/bsdata/hmr/normal.common.hmr.merged100b.bed

for i in G145 W145 617N NY108 22A B11-27235 CP10 PL5 Laia SH5YSY U87MG H1437 H1672 H157 M673 M697 22B 22C 468PT 468LN PC3 22RV1
do

echo "tabix -B CpG_context_$i.smooth.bed.gz $region | bedtools intersect -wao -a $region -b - | awk '{start=\$2; end=\$3; pos=\$5; m=\$10; size=end - start; mid=(start + end)/2; i=int((pos - mid)/size*1000 + .5)/1000; if(m != \".\") {out[i] += m; n[i] += 1}}END{for(j in n) print j, out[j], n[j]}' > /home/user/bsdata/hmr/$i.common.hmr.smiley.txt"
done | parallel --gnu

##  Compute distribution of methylation per relative position inside c-HMRs (c-HMRs scaled + 50% flanking region)
# set region
awk '{OFS="\t"; start=$2; end=$3; size=end - start; flank=int(size/2 +.5); start=start - flank; if(start < 1) start=0; end=end + flank; print $1, start, end}'\
 /home/user/bsdata/hmr/normal.common.hmr.merged100b.bed > /home/user/bsdata/hmr/normal.common.hmr.merged100b.flank50.bed
region=/home/user/bsdata/hmr/normal.common.hmr.merged100b.flank50.bed

for i in G145 W145 617N NY108 22A B11-27235 CP10 PL5 Laia SH5YSY U87MG H1437 H1672 H157 M673 M697 22B 22C 468PT 468LN PC3 22RV1
do

echo "tabix -B CpG_context_$i.smooth.bed.gz $region | bedtools intersect -wao -a $region -b - | awk '{start=\$2; end=\$3; pos=\$5; m=\$10; size=end - start; mid=(start + end)/2; i=int((pos - mid)/size*1000 + .5)/1000; if(m != \".\") {out[i] += m; n[i] += 1}}END{for(j in n) print j, out[j], n[j]}' > /home/user/bsdata/hmr/$i.commonwithflanks.hmr.smiley.txt"
done | parallel --gnu -j 20


##  Compute distribution of methylation per relative position inside c-HMRs (2kb flanking region ti mid of c-HMRs)
# set region
awk -v flank=2000 '{OFS="\t"; start=$2; end=$3; mid=int((end + start)/2 +.5); start=mid - flank; if(start < 1) start=0; end=mid + flank; print $1, start, end}'\
 /home/user/bsdata/hmr/normal.common.hmr.merged100b.bed > /home/user/bsdata/hmr/normal.common.hmr.merged100b.mid2kb.bed
region=/home/user/bsdata/hmr/normal.common.hmr.merged100b.mid2kb.bed

for i in G145 W145 617N NY108 22A B11-27235 CP10 PL5 Laia SH5YSY U87MG H1437 H1672 H157 M673 M697 22B 22C 468PT 468LN PC3 22RV1
do

echo "tabix -B CpG_context_$i.smooth.bed.gz $region | bedtools intersect -wao -a $region -b - | awk '{start=\$2; end=\$3; pos=\$5; m=\$10; size=end - start; mid=(start + end)/2; i=int((pos - mid)/size*1000 + .5)/1000; if(m != \".\") {out[i] += m; n[i] += 1}}END{for(j in n) print j, out[j], n[j]}' > /home/user/bsdata/hmr/$i.common2kbfrommid.hmr.smiley.txt"
done | parallel --gnu -j 20


###   t-HMRs

##  Compute distribution of methylation per relative position inside t-HMRs (t-HMRs scaled)
# set region
region=/home/user/bsdata/hmr/normals.hmr.specific.bed

for i in G145 W145 617N NY108 22A B11-27235 CP10 Laia SH5YSY U87MG H1437 H1672 H157 M673 M697 22B 22C 468PT 468LN PC3 22RV1
do

control=`awk -F "," -v sam=$i '$1==sam {print $6}' /home/user/bsdata/samples.info.newnames.control.csv`

echo "tabix -B CpG_context_$i.smooth.bed.gz <(grep $control $region) | bedtools intersect -wao -a <(grep $control $region) -b - | awk '{start=\$2; end=\$3; pos=\$13; m=\$17; u=\$16; t=m+u; size=end - start; mid=(start + end)/2; i=int((pos - mid)/size*1000 + .5)/1000; if(t> 0) {out[i] += m/t; n[i] += 1}}END{for(j in n) print j, out[j], n[j]}' > /home/user/bsdata/hmr/$i.specific.hmr.smiley.txt"
done | parallel --gnu

##  Compute distribution of methylation per relative position inside t-HMRs (t-HMRs scaled + 50% flanking region)

# set region

for i in G145 W145 617N NY108 22A B11-27235 CP10 Laia SH5YSY U87MG H1437 H1672 H157 M673 M697 22B 22C 468PT 468LN PC3 22RV1
do

control=`awk -F "," -v sam=$i '$1==sam {print $6}' /home/user/bsdata/samples.info.newnames.control.csv`

grep $control $region | awk '{OFS="\t"; start=$2; end=$3; size=end - start; flank=int(size/2 +.5); start=start - flank; if(start < 1) start=0; end=end + flank; print $1, start, end}' >\
 /home/user/bsdata/hmr/$control.hmr.specific.flank50.bed

done


for i in G145 W145 617N NY108 22A B11-27235 CP10 Laia SH5YSY U87MG H1437 H1672 H157 M673 M697 22B 22C 468PT 468LN PC3 22RV1
do

control=`awk -F "," -v sam=$i '$1==sam {print $6}' /home/user/bsdata/samples.info.newnames.control.csv`

echo "tabix -B CpG_context_$i.smooth.bed.gz  /home/user/bsdata/hmr/$control.hmr.specific.flank50.bed | bedtools intersect -wao -a  /home/user/bsdata/hmr/$control.hmr.specific.flank50.bed -b - | awk '{start=\$2; end=\$3; pos=\$5; u=\$8; m=\$9; t=m+u; size=end - start; mid=(start + end)/2; i=int((pos - mid)/size*1000 + .5)/1000; if(t > 0) {out[i] += m/t; n[i] += 1}}END{for(j in n) print j, out[j], n[j]}' > /home/user/bsdata/hmr/$i.specificwithflanks.hmr.smiley.txt"
done | parallel --gnu -j 20




##  Compute distribution of methylation per relative position inside t-HMRs (2kb flanking region ti mid of c-HMRs)

# set region

for i in G145 W145 617N NY108 22A B11-27235 CP10 Laia SH5YSY U87MG H1437 H1672 H157 M673 M697 22B 22C 468PT 468LN PC3 22RV1
do

control=`awk -F "," -v sam=$i '$1==sam {print $6}' /home/user/bsdata/samples.info.newnames.control.csv`

grep $control $region |awk -v flank=2000 '{OFS="\t"; start=$2; end=$3; mid=int((end + start)/2 +.5); start=mid - flank; if(start < 1) start=0; end=mid + flank; print $1, start, end}' >\
 /home/user/bsdata/hmr/$control.hmr.specific.mid2kb.bed

done


for i in G145 W145 617N NY108 22A B11-27235 CP10 Laia SH5YSY U87MG H1437 H1672 H157 M673 M697 22B 22C 468PT 468LN PC3 22RV1
do

control=`awk -F "," -v sam=$i '$1==sam {print $6}' /home/user/bsdata/samples.info.newnames.control.csv`

echo "tabix -B CpG_context_$i.smooth.bed.gz /home/user/bsdata/hmr/$control.hmr.specific.mid2kb.bed | bedtools intersect -wao -a /home/user/bsdata/hmr/$control.hmr.specific.mid2kb.bed -b - | awk '{start=\$2; end=\$3; pos=\$5; u=\$8; m=\$9; t=u + m; size=end - start; mid=(start + end)/2; i=int((pos - mid)/size*1000 + .5)/1000; if(t > 0) {out[i] += m/t; n[i] += 1}}END{for(j in n) print j, out[j], n[j]}' > /home/user/bsdata/hmr/$i.specific2kbfrommid.hmr.smiley.txt"
done | parallel --gnu -j 20


####
##  ENCODE DNase1 overlap
####

# seet up codes
echo NY108 Hepg2 Liver >> chromatineState/codes.txt
echo Laia Gm12878 Blood >> chromatineState/codes.txt
echo B11-27235 Hmec Breast >> chromatineState/codes.txt

# arrange ENCODE data
while read line 
do
arr=($line)
code=${arr[1]}
tissue=${arr[2]}
sed 's/^chr//1' /home/user/resources/human/wgEncodeAwgDnaseUwduke${code}UniPk.narrowPeak > /home/user/resources/human/wgEncodeAwgDnaseUwduke${code}UniPk.narrowPeak.bed
awk '{OFS="\t"; if($1 != "Y" && $7 > 10) print }' /home/user/resources/human/wgEncodeAwgDnaseUwduke${code}UniPk.narrowPeak.bed > /home/user/resources/human/encode.dnasei.$tissue.noY.bed
done < chromatineState/codes.txt

# set up out file
echo HMR DNaseI specific all > encode.dnaseI.overlap.txt

# intersection with ENCODE data
while read line1
do
arr1=($line1)
sample1=${arr1[0]}
code1=${arr1[1]}
tissue1=${arr1[2]}

while read line2
do
arr2=($line2)
sample2=${arr2[0]}
code2=${arr2[1]}
tissue2=${arr2[2]}

overA=`bedtools intersect -u -a $sample1.hmr.specific.bed -b  /home/user/resources/human/encode.dnasei.$tissue2.noY.bed | wc -l`
overB=`bedtools intersect -u -a CpG_context_$sample1.hmr.bed -b  /home/user/resources/human/encode.dnasei.$tissue2.noY.bed | wc -l`
totA=`cat $sample1.hmr.specific.bed | wc -l`
totB=`cat CpG_context_$sample2.hmr.bed | wc -l`
rA=`echo "scale=3; $overA/$totA" | bc `
rB=`echo "scale=3; $overB/$totB" | bc `

echo $tissue1 $tissue2 $rA $rB

done < chromatineState/codes.txt

done < chromatineState/codes.txt >> encode.dnaseI.overlap.txt



####################################
###### tissue specific binding sites
####################################

# get only coging promoters
coding=/home/user/resources/human/gencode.v16.alltranscripts.coding.prom.nochr.bed

# get normal samples
normals=`awk '$3 ~ /ormal/ && $1 != "PL5" && $1 != "NB" && $1 != "Julia" && $1 != "103"  && $1 != "W145" {ORS=" "; print $1}'  /home/user/bsdata/samples.info.txt`

#overlap promoters with HMRs (average smoothed methylation below .1)
for i in $normals
do
awk '$10/$7 <= .1 {print}' CpG_context_$i.hmr.metrics.bed | bedtools intersect -wao -a $coding -b - | awk '{a[$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6] += $NF}END{for(i in a) print i "\t" 1-(a[i]==0)}' | bedtools sort -i > genomicFeatures/promoters.overlapping.$i.hmr.bed
done

# sum up all normal sample counts per promoter
paste genomicFeatures/promoters.overlapping.*.hmr.bed | awk 'BEGIN{OFS="\t"}{for(i=1; i <= NF/7; i++){a += $(i*7)} print $1, $2, $3, $4, $5, $6, a; a = 0}' > genomicFeatures/promoters.overlapping.all.hmr.counts.bed

# select promoters overlapping all normal samples and none of the common HMR
awk '$7>=7{print}' genomicFeatures/promoters.overlapping.all.hmr.counts.bed | bedtools intersect -v -a - -b normal.common.hmr.merged100b.bed > genomicFeatures/promoters.binding.sites.allsamples.nocommon.bed

# select promoters overlapping all normal samples and none of the frequent HMR
awk '$7>=7{print}' genomicFeatures/promoters.overlapping.all.hmr.counts.bed | bedtools intersect -v -a - -b normal.frequent.hmr.merged100b.bed > genomicFeatures/promoters.binding.sites.allsamples.nofrequent.bed


# select promoters overlapping at least ahlf of the normal samples and none of the common HMR
awk '$7>=4{print}' genomicFeatures/promoters.overlapping.all.hmr.counts.bed | bedtools intersect -v -a - -b normal.common.hmr.merged100b.bed > genomicFeatures/promoters.binding.sites.halfsamples.nocommon.bed

# select promoters overlapping at least ahlf of the normal samples and none of the frequent HMR
awk '$7>=4{print}' genomicFeatures/promoters.overlapping.all.hmr.counts.bed | bedtools intersect -v -a - -b normal.frequent.hmr.merged100b.bed > genomicFeatures/promoters.binding.sites.halfsamples.nofrequent.bed



##############
###### piRNA  
##############

# get only coging promoters
coding=/home/user/resources/human/gencode.v16.alltranscripts.coding.prom.nochr.bed

# get normal samples
normals=`awk '$3 ~ /ormal/ && $1 != "PL5" && $1 != "NB" && $1 != "Julia" && $1 != "103"  && $1 != "W145" {ORS=" "; print $1}'  /home/user/bsdata/samples.info.txt`

for i in $normals
do
o=`bedtools intersect -wao -a /home/user/resources/human/pirna.bed -b $i.hmr.specific.bed | awk '{a += $NF}END{print a}'`
u1=`bedtools intersect -u -b /home/user/resources/human/pirna.bed -a $i.hmr.specific.bed | wc -l`
u2=`bedtools intersect -u -a /home/user/resources/human/pirna.bed -b $i.hmr.specific.bed | wc -l`
r1=`cat $i.hmr.specific.bed | wc -l`
r2=`cat /home/user/resources/human/pirna.bed | wc -l`
b1=`awk '{a += $3 - $2}END{print a}' $i.hmr.specific.bed`
b2=`awk '{a += $3 - $2}END{print a}' /home/user/resources/human/pirna.bed `
echo -e "$i\t$o\t$u1\t$u2\t$r1\t$r2\t$b1\t$b2"
done > pirna.overlap.txt

#######################################
####    Distribution of CpG intensities
#######################################

echo /home/user/bsdata/CpG_context_*smooth.bed.gz | sed 's/ /\n/g' |  parallel -j10 scripts/get.cpgdistribution.sh {} ../results/

#### All CpG in genome
samples=`cat ../samples.info.txt | cut -d " " -f 1 | tr "\n" " "`
chromosomes=`seq -s " " 1 22`
chromosomes=`echo $chromosomes X`

# raw counts
# all
for i in $samples
do
time echo $chromosomes | tr " " "\n" | xargs -n 1 -P 40 -I {} sh -c "tabix /home/user/bsdata/CpG_context_$i.smooth.bed.gz {} | awk '{u = \$5; m = \$6; t = m + u; if( t > 0){ key = int(m/t*500); out[key] += 1}}END{for(i in out) print i, out[i]}'" | awk '{out[$1] += $2}END{for(i in out) print i, out[i]}' > ../results/$i.allcpg.rawmeth.distribution
done

# coverage above 10
for i in $samples
do
time echo $chromosomes | tr " " "\n" | xargs -n 1 -P 40 -I {} sh -c "tabix /home/user/bsdata/CpG_context_$i.smooth.bed.gz {} | awk '{u = \$5; m = \$6; t = m + u; if( t > 10){ key = int(m/t*500); out[key] += 1}}END{for(i in out) print i, out[i]}'" | awk '{out[$1] += $2}END{for(i in out) print i, out[i]}' > ../results/$i.allcpg.rawmeth.10.distribution
done

# coverage above 20
for i in $samples
do
time echo $chromosomes | tr " " "\n" | xargs -n 1 -P 40 -I {} sh -c "tabix /home/user/bsdata/CpG_context_$i.smooth.bed.gz {} | awk '{u = \$5; m = \$6; t = m + u; if( t > 20){ key = int(m/t*500); out[key] += 1}}END{for(i in out) print i, out[i]}'" | awk '{out[$1] += $2}END{for(i in out) print i, out[i]}' > ../results/$i.allcpg.rawmeth.20.distribution
done

# coverage above 100
for i in $samples
do
time echo $chromosomes | tr " " "\n" | xargs -n 1 -P 40 -I {} sh -c "tabix /home/user/bsdata/CpG_context_$i.smooth.bed.gz {} | awk '{u = \$5; m = \$6; t = m + u; if( t > 100){ key = int(m/t*500); out[key] += 1}}END{for(i in out) print i, out[i]}'" | awk '{out[$1] += $2}END{for(i in out) print i, out[i]}' > ../results/$i.allcpg.rawmeth.100.distribution
done

# smoothed methylation
# all
for i in $samples
do
time echo $chromosomes | tr " " "\n" | xargs -n 1 -P 40 -I {} sh -c "tabix /home/user/bsdata/CpG_context_$i.smooth.bed.gz {} | awk '{u = \$5; m = \$6; t = m + u; if( t > 0){ key = int(\$7*500); out[key] += 1}}END{for(i in out) print i, out[i]}'" | awk '{out[$1] += $2}END{for(i in out) print i, out[i]}' > ../results/$i.allcpg.smoothmeth.distribution
done

# coverage above 10
for i in $samples
do
time echo $chromosomes | tr " " "\n" | xargs -n 1 -P 40 -I {} sh -c "tabix /home/user/bsdata/CpG_context_$i.smooth.bed.gz {} | awk '{u = \$5; m = \$6; t = m + u; if( t > 10){ key = int(\$7*500); out[key] += 1}}END{for(i in out) print i, out[i]}'" | awk '{out[$1] += $2}END{for(i in out) print i, out[i]}' > ../results/$i.allcpg.smoothmeth.10.distribution
done

# coverage above 20
for i in $samples
do
time echo $chromosomes | tr " " "\n" | xargs -n 1 -P 40 -I {} sh -c "tabix /home/user/bsdata/CpG_context_$i.smooth.bed.gz {} | awk '{u = \$5; m = \$6; t = m + u; if( t > 20){ key = int(\$7*500); out[key] += 1}}END{for(i in out) print i, out[i]}'" | awk '{out[$1] += $2}END{for(i in out) print i, out[i]}' > ../results/$i.allcpg.smoothmeth.20.distribution
done

############
####    DMRs
############

#################
##  Colon triplet
#################

# make results directory
mkdir -p ../dmr

# compute dmr
scripts/smooth.dmr.r ../CpG_context_22B.smooth.bed.gz ../CpG_context_22A.smooth.bed.gz ../dmr/22B.vs.22A.dmr.bed $chromosomes
scripts/smooth.dmr.r ../CpG_context_22C.smooth.bed.gz ../CpG_context_22A.smooth.bed.gz ../dmr/22C.vs.22A.dmr.bed $chromosomes

# get consensus dmr (bases present in both dmrs)
# split by sign
awk '$5==1{print}' ../dmr/22B.vs.22A.dmr.bed > ../dmr/22B.vs.22A.dmr.bed.pos
awk '$5!=1{print}' ../dmr/22B.vs.22A.dmr.bed > ../dmr/22B.vs.22A.dmr.bed.neg
awk '$5==1{print}' ../dmr/22C.vs.22A.dmr.bed > ../dmr/22C.vs.22A.dmr.bed.pos
awk '$5!=1{print}' ../dmr/22C.vs.22A.dmr.bed > ../dmr/22C.vs.22A.dmr.bed.neg

# intersect
bedtools intersect -a ../dmr/22B.vs.22A.dmr.bed.pos -b ../dmr/22C.vs.22A.dmr.bed.pos > ../dmr/22A.consensus.dmr.bed.pos
bedtools intersect -a ../dmr/22B.vs.22A.dmr.bed.neg -b ../dmr/22C.vs.22A.dmr.bed.neg > ../dmr/22A.consensus.dmr.bed.neg

# merge and order
cat ../dmr/22A.consensus.dmr.bed.pos ../dmr/22A.consensus.dmr.bed.neg | bedtools sort -i - > ../dmr/22A.consensus.dmr.bed

# remove split files
rm ../dmr/22?.vs.22A.dmr.bed.??? ../dmr/22A.consensus.dmr.bed.???

# compute dmr real metrics
scripts/compute.region.metrics.sh ../CpG_context_22B.smooth.bed.gz ../CpG_context_22A.smooth.bed.gz ../dmr/22B.vs.22A.dmr.bed scripts/smoothed.real.dmr.metrics.awk > ../dmr/22B.vs.22A.dmr.real.metrics.bed
scripts/compute.region.metrics.sh ../CpG_context_22C.smooth.bed.gz ../CpG_context_22A.smooth.bed.gz ../dmr/22C.vs.22A.dmr.bed scripts/smoothed.real.dmr.metrics.awk > ../dmr/22C.vs.22A.dmr.real.metrics.bed
scripts/compute.region.metrics.sh ../CpG_context_22B.smooth.bed.gz ../CpG_context_22A.smooth.bed.gz ../dmr/22A.consensus.dmr.bed scripts/smoothed.real.dmr.metrics.awk > ../dmr/22A.consensus.dmr.real.metrics.bed

# compute dmr interval metrics
scripts/compute.region.metrics.sh ../CpG_context_22B.smooth.bed.gz ../CpG_context_22A.smooth.bed.gz ../dmr/22B.vs.22A.dmr.bed scripts/smoothed.interval.dmr.metrics.awk > ../dmr/22B.vs.22A.dmr.interval.metrics.bed
scripts/compute.region.metrics.sh ../CpG_context_22C.smooth.bed.gz ../CpG_context_22A.smooth.bed.gz ../dmr/22C.vs.22A.dmr.bed scripts/smoothed.interval.dmr.metrics.awk > ../dmr/22C.vs.22A.dmr.interval.metrics.bed
scripts/compute.region.metrics.sh ../CpG_context_22B.smooth.bed.gz ../CpG_context_22A.smooth.bed.gz ../dmr/22A.consensus.dmr.bed scripts/smoothed.interval.dmr.metrics.awk > ../dmr/22A.consensus.dmr.interval.metrics.bed

# get methylation level of all samples at B vs A dmrs
tabix -B ../CpG_context_22A.smooth.bed.gz ../dmr/22B.vs.22A.dmr.interval.metrics.bed | bedtools intersect -wao -a ../dmr/22B.vs.22A.dmr.interval.metrics.bed -b - | awk '{key=$1" "$2" "$3" "$6" "$7" "$8" "$9; u[key] += $14; m[key] += $15; s[key] += $16; n[key] += 1}END{for(i in n) print i, u[i], m[i], s[i], n[i]}' | sed 's/ /\t/g' | bedtools sort -i - > ../dmr/22A.at.22B.vs.22A.dmr.bed

tabix -B ../CpG_context_22B.smooth.bed.gz ../dmr/22B.vs.22A.dmr.interval.metrics.bed | bedtools intersect -wao -a ../dmr/22B.vs.22A.dmr.interval.metrics.bed -b - | awk '{key=$1" "$2" "$3" "$6" "$7" "$8" "$9; u[key] += $14; m[key] += $15; s[key] += $16; n[key] += 1}END{for(i in n) print i, u[i], m[i], s[i], n[i]}' | sed 's/ /\t/g' | bedtools sort -i - > ../dmr/22B.at.22B.vs.22A.dmr.bed

tabix -B ../CpG_context_22C.smooth.bed.gz ../dmr/22B.vs.22A.dmr.interval.metrics.bed | bedtools intersect -wao -a ../dmr/22B.vs.22A.dmr.interval.metrics.bed -b - | awk '{key=$1" "$2" "$3" "$6" "$7" "$8" "$9; u[key] += $14; m[key] += $15; s[key] += $16; n[key] += 1}END{for(i in n) print i, u[i], m[i], s[i], n[i]}' | sed 's/ /\t/g' | bedtools sort -i - > ../dmr/22C.at.22B.vs.22A.dmr.bed

# get methylation level of all samples at consensus dmrs
tabix -B ../CpG_context_22A.smooth.bed.gz ../dmr/22A.consensus.dmr.interval.metrics.bed | bedtools intersect -wao -a ../dmr/22A.consensus.dmr.interval.metrics.bed -b - | awk '{key=$1" "$2" "$3" "$6" "$7" "$8" "$9; u[key] += $14; m[key] += $15; s[key] += $16; n[key] += 1}END{for(i in n) print i, u[i], m[i], s[i], n[i]}' | sed 's/ /\t/g' | bedtools sort -i - > ../dmr/22A.at.consensus.dmr.interval.bed

tabix -B ../CpG_context_22B.smooth.bed.gz ../dmr/22A.consensus.dmr.interval.metrics.bed | bedtools intersect -wao -a ../dmr/22A.consensus.dmr.interval.metrics.bed -b - | awk '{key=$1" "$2" "$3" "$6" "$7" "$8" "$9; u[key] += $14; m[key] += $15; s[key] += $16; n[key] += 1}END{for(i in n) print i, u[i], m[i], s[i], n[i]}' | sed 's/ /\t/g' | bedtools sort -i - > ../dmr/22B.at.consensus.dmr.interval.bed

tabix -B ../CpG_context_22C.smooth.bed.gz ../dmr/22A.consensus.dmr.interval.metrics.bed | bedtools intersect -wao -a ../dmr/22A.consensus.dmr.interval.metrics.bed -b - | awk '{key=$1" "$2" "$3" "$6" "$7" "$8" "$9; u[key] += $14; m[key] += $15; s[key] += $16; n[key] += 1}END{for(i in n) print i, u[i], m[i], s[i], n[i]}' | sed 's/ /\t/g' | bedtools sort -i - > ../dmr/22C.at.consensus.dmr.interval.bed


#################
##  Breast triplet
#################

# make results directory
mkdir -p ../dmr

# compute dmr
scripts/smooth.dmr.r ../CpG_context_468PT.smooth.bed.gz ../CpG_context_B11-27235.smooth.bed.gz ../dmr/468PT.vs.B11-27235.dmr.bed $chromosomes
scripts/smooth.dmr.r ../CpG_context_468LN.smooth.bed.gz ../CpG_context_B11-27235.smooth.bed.gz ../dmr/468LN.vs.B11-27235.dmr.bed $chromosomes

# get consensus dmr (bases present in both dmrs)
# split by sign
awk '$5==1{print}' ../dmr/468PT.vs.B11-27235.dmr.bed > ../dmr/468PT.vs.B11-27235.dmr.bed.pos
awk '$5!=1{print}' ../dmr/468PT.vs.B11-27235.dmr.bed > ../dmr/468PT.vs.B11-27235.dmr.bed.neg
awk '$5==1{print}' ../dmr/468LN.vs.B11-27235.dmr.bed > ../dmr/468LN.vs.B11-27235.dmr.bed.pos
awk '$5!=1{print}' ../dmr/468LN.vs.B11-27235.dmr.bed > ../dmr/468LN.vs.B11-27235.dmr.bed.neg

# intersect
bedtools intersect -a ../dmr/468PT.vs.B11-27235.dmr.bed.pos -b ../dmr/468LN.vs.B11-27235.dmr.bed.pos > ../dmr/B11-27235.consensus.dmr.bed.pos
bedtools intersect -a ../dmr/468PT.vs.B11-27235.dmr.bed.neg -b ../dmr/468LN.vs.B11-27235.dmr.bed.neg > ../dmr/B11-27235.consensus.dmr.bed.neg

# merge and order
cat ../dmr/B11-27235.consensus.dmr.bed.pos ../dmr/B11-27235.consensus.dmr.bed.neg | bedtools sort -i - > ../dmr/B11-27235.consensus.dmr.bed

# remove split files
rm ../dmr/468??.vs.B11-27235.dmr.bed.??? ../dmr/B11-27235.consensus.dmr.bed.???

# compute dmr real metrics
scripts/compute.region.metrics.sh ../CpG_context_468PT.smooth.bed.gz ../CpG_context_B11-27235.smooth.bed.gz ../dmr/468PT.vs.B11-27235.dmr.bed scripts/smoothed.real.dmr.metrics.awk > ../dmr/468PT.vs.B11-27235.dmr.real.metrics.bed
scripts/compute.region.metrics.sh ../CpG_context_468LN.smooth.bed.gz ../CpG_context_B11-27235.smooth.bed.gz ../dmr/468LN.vs.B11-27235.dmr.bed scripts/smoothed.real.dmr.metrics.awk > ../dmr/468LN.vs.B11-27235.dmr.real.metrics.bed
scripts/compute.region.metrics.sh ../CpG_context_468PT.smooth.bed.gz ../CpG_context_B11-27235.smooth.bed.gz ../dmr/B11-27235.consensus.dmr.bed scripts/smoothed.real.dmr.metrics.awk > ../dmr/B11-27235.consensus.dmr.real.metrics.bed

# compute dmr interval metrics
scripts/compute.region.metrics.sh ../CpG_context_468PT.smooth.bed.gz ../CpG_context_B11-27235.smooth.bed.gz ../dmr/468PT.vs.B11-27235.dmr.bed scripts/smoothed.interval.dmr.metrics.awk > ../dmr/468PT.vs.B11-27235.dmr.interval.metrics.bed
scripts/compute.region.metrics.sh ../CpG_context_468LN.smooth.bed.gz ../CpG_context_B11-27235.smooth.bed.gz ../dmr/468LN.vs.B11-27235.dmr.bed scripts/smoothed.interval.dmr.metrics.awk > ../dmr/468LN.vs.B11-27235.dmr.interval.metrics.bed
scripts/compute.region.metrics.sh ../CpG_context_468PT.smooth.bed.gz ../CpG_context_B11-27235.smooth.bed.gz ../dmr/B11-27235.consensus.dmr.bed scripts/smoothed.interval.dmr.metrics.awk > ../dmr/B11-27235.consensus.dmr.interval.metrics.bed

# get methylation level of all samples at B vs A dmrs
tabix -B ../CpG_context_B11-27235.smooth.bed.gz ../dmr/468PT.vs.B11-27235.dmr.interval.metrics.bed | bedtools intersect -wao -a ../dmr/468PT.vs.B11-27235.dmr.interval.metrics.bed -b - | awk '{key=$1" "$2" "$3" "$6" "$7" "$8" "$9;         u[key] += $14; m[key] += $15; s[key] += $16; n[key] += 1}END{for(i in n) print i, u[i], m[i], s[i], n[i]}' | sed 's/ /\t/g' | bedtools sort -i - > ../dmr/B11-27235.at.468PT.vs.B11-27235.dmr.bed

tabix -B ../CpG_context_468PT.smooth.bed.gz ../dmr/468PT.vs.B11-27235.dmr.interval.metrics.bed | bedtools intersect -wao -a ../dmr/468PT.vs.B11-27235.dmr.interval.metrics.bed -b - | awk '{key=$1" "$2" "$3" "$6" "$7" "$8" "$9;             u[key] += $14; m[key] += $15; s[key] += $16; n[key] += 1}END{for(i in n) print i, u[i], m[i], s[i], n[i]}' | sed 's/ /\t/g' | bedtools sort -i - > ../dmr/468PT.at.468PT.vs.B11-27235.dmr.bed

tabix -B ../CpG_context_468LN.smooth.bed.gz ../dmr/468PT.vs.B11-27235.dmr.interval.metrics.bed | bedtools intersect -wao -a ../dmr/468PT.vs.B11-27235.dmr.interval.metrics.bed -b - | awk '{key=$1" "$2" "$3" "$6" "$7" "$8" "$9;             u[key] += $14; m[key] += $15; s[key] += $16; n[key] += 1}END{for(i in n) print i, u[i], m[i], s[i], n[i]}' | sed 's/ /\t/g' | bedtools sort -i - > ../dmr/468LN.at.468PT.vs.B11-27235.dmr.bed

# get methylation level of all samples at consensus dmrs
tabix -B ../CpG_context_B11-27235.smooth.bed.gz ../dmr/B11-27235.consensus.dmr.interval.metrics.bed | bedtools intersect -wao -a ../dmr/B11-27235.consensus.dmr.interval.metrics.bed -b - | awk '{key=$1" "$2" "$3" "$6" "$7" "$8" "$9;       u[key] += $14; m[key] += $15; s[key] += $16; n[key] += 1}END{for(i in n) print i, u[i], m[i], s[i], n[i]}' | sed 's/ /\t/g' | bedtools sort -i - > ../dmr/B11-27235.at.consensus.dmr.interval.bed

tabix -B ../CpG_context_468PT.smooth.bed.gz ../dmr/B11-27235.consensus.dmr.interval.metrics.bed | bedtools intersect -wao -a ../dmr/B11-27235.consensus.dmr.interval.metrics.bed -b - | awk '{key=$1" "$2" "$3" "$6" "$7" "$8" "$9;           u[key] += $14; m[key] += $15; s[key] += $16; n[key] += 1}END{for(i in n) print i, u[i], m[i], s[i], n[i]}' | sed 's/ /\t/g' | bedtools sort -i - > ../dmr/468PT.at.consensus.dmr.interval.bed

tabix -B ../CpG_context_468LN.smooth.bed.gz ../dmr/B11-27235.consensus.dmr.interval.metrics.bed | bedtools intersect -wao -a ../dmr/B11-27235.consensus.dmr.interval.metrics.bed -b - | awk '{key=$1" "$2" "$3" "$6" "$7" "$8" "$9;           u[key] += $14; m[key] += $15; s[key] += $16; n[key] += 1}END{for(i in n) print i, u[i], m[i], s[i], n[i]}' | sed 's/ /\t/g' | bedtools sort -i - > ../dmr/468LN.at.consensus.dmr.interval.bed


############################
####    Distance correlation
############################

# create results directory
mkdir -p /home/user/bsdata/cgdistance

# get chromosomes as regions                                                                                                                                                                                                                
chromosomes=`awk '{print $1":"1"-"$2}' /home/user/resources/human/hg19.chrsizes`

# get samples
samples=`grep -i normal ../samples.info.txt | cut -d " " -f 1`

### compute correlation
# for each sample
for smp in $samples
do
echo
echo "    Sample" $smp                                                                                                                                                                                                                       
# compute distance correlation for each chromosome
echo
echo "Computing unsigned distance correlation (w = 2000) ..."
time echo $chromosomes | tr " " "\n" | parallel -j 10\
 "tabix /home/user/bsdata/CpG_context_$smp.smooth.bed.gz {} | scripts/cg.distance.correlation.v3.raw.awk -v w=2000 -v cove=4" |\
 scripts/cg.distance.correlation.processV3output.awk > /home/user/bsdata/cgdistance/$smp.original.pairedend.2000.unsigned.txt
echo
echo "... Done"
don


############################
####    Closest promoter
############################

# set TSS
awk 'BEGIN{OFS="\t"}{tss=($2+$3)/2; $2=tss; $3=tss; print}' /home/user/resources/human/gencode.v16.alltranscripts.prom.nochr.bed > /home/user/resources/human/gencode.v16.alltranscripts.TSS.nochr.bed

# all HMRs (include overlap with promoters)

for i in /home/user/bsdata/hmr/CpG_context_*hmr.metrics.bed
do
name=`basename $i .bed | sed 's/CpG_context_//1'`
bedtools closest -D -a $i -b /home/user/resources/human/gencode.v16.alltranscripts.prom.nochr.bed > $name.closestsGene.bed
done

# only normal HMRs (distance to TSS, not to promoter ends)
normals=`awk '$3 ~ /ormal/ && $1 != "PL5" && $1 != "NB" && $1 != "Julia" && $1 != "103"  && $1 != "W145" {ORS=" "; print $1}'  /home/user/bsdata/samples.info.txt`

for i in $normals
do
bedtools closest -d -a /home/user/bsdata/hmr/CpG_context_$i.hmr.metrics.bed -b /home/user/resources/human/gencode.v16.alltranscripts.TSS.nochr.bed > $i.closestsGene.TSS.bed
bedtools shuffle -i /home/user/bsdata/hmr/CpG_context_$i.hmr.metrics.bed -g /home/user/resources/human/hg19.chrsizes |\
 bedtools closest -d -a - -b /home/user/resources/human/gencode.v16.alltranscripts.TSS.nochr.bed > $i.closestsGene.TSS.shuffle.bed
done


# tissue HMRs (outside of promoters)

for i in *hmr.specific.bed
do
name=`basename $i .bed`
bedtools closest -io -d -a $i -b /home/user/resources/human/gencode.v16.alltranscripts.prom.nochr.bed > $name.closestsGene.bed
bedtools closest -D "b" -a $i -b /home/user/resources/human/gencode.v16.alltranscripts.prom.nochr.bed > $name.closestGene.includeproms.bed
done

# Distal-proximal tissue specifc HMR proportion

prom=/home/user/resources/human/gencode.v16.alltranscripts.prom.nochr.bed
normals=`awk '$3 ~ /ormal/ && $1 != "PL5" && $1 != "NB" && $1 != "Julia" && $1 != "103"  && $1 != "W145" {ORS=" "; print $1}'  /home/user/bsdata/samples.info.txt`

## take into account all distal HMRs
# take into account only non-common proximal HMRs
for i in $normals
do
for j in $normals
do
bedtools intersect -v -a CpG_context_$j.hmr.bed -b normal.common.hmr.bed | bedtools intersect -u -a $prom -b - |\
 cut -f 5 | awk -v genes=$i.hmr.specific.closestsGene.bed -v i=$i -v j=$j 'BEGIN{while(( getline < genes ) > 0){if(!($16 in trans)) {trans[$16] = 1; n += 1 }}}{if($1 in trans) out += 1}END{OFS="\t"; print i, j, out, n, NR}' 
done
done > distal.noncommon_proximal.tissue.hmrs.txt

# take into account all proximal HMRs
for i in $normals
do
for j in $normals
do
bedtools intersect -u -a $prom -b CpG_context_$j.hmr.bed | cut -f 5 |\
 awk -v genes=$i.hmr.specific.closestsGene.bed -v i=$i -v j=$j 'BEGIN{while(( getline < genes ) > 0){if(!($16 in trans)) {trans[$16] = 1; n += 1 }}}{if($1 in trans) out += 1}END{OFS="\t"; print i, j, out, n, NR}' 
done
done > distal.all_proximal.tissue.hmrs.txt

## take into account only distal HMRs closer than 10Kb
# take into account only non-common proximal HMRs
for i in $normals
do
for j in $normals
do
bedtools intersect -v -a CpG_context_$j.hmr.bed -b normal.common.hmr.bed | bedtools intersect -u -a $prom -b - |\
 cut -f 5 | awk -v genes=$i.hmr.specific.closestsGene.bed -v i=$i -v j=$j\
 'BEGIN{while(( getline < genes ) > 0){if(!($16 in trans) && ($19 <= (10000 - 2000))) {trans[$16] = 1; n += 1 }}}{if($1 in trans) out += 1}END{OFS="\t"; print i, j, out, n, NR}' 
done
done > distal10kb.noncommon_proximal.tissue.hmrs.txt

# take into account all proximal HMRs
for i in $normals
do
for j in $normals
do
bedtools intersect -u -a $prom -b CpG_context_$j.hmr.bed | cut -f 5 |\
 awk -v genes=$i.hmr.specific.closestsGene.bed -v i=$i -v j=$j\
 'BEGIN{while(( getline < genes ) > 0){if(!($16 in trans) && ($19 <= (10000 - 2000))) {trans[$16] = 1; n += 1 }}}{if($1 in trans) out += 1}END{OFS="\t"; print i, j, out, n, NR}' 
done
done > distal10kb.all_proximal.tissue.hmrs.txt

##################
####    Super HMRs
##################

#set wd
cd /home/user/bsdata/hmr/superhmr

bedtools intersect -v -a 22B.superHMR.closestTSS.bed -b 22A.superHMR.closestTSS.bed > 22Bnotin22A.superHMR.closestTSS.bed
bedtools intersect -v -a 22C.superHMR.closestTSS.bed -b 22A.superHMR.closestTSS.bed > 22Cnotin22A.superHMR.closestTSS.bed
bedtools intersect -v -a 22C.superHMR.closestTSS.bed -b 22B.superHMR.closestTSS.bed > 22Cnotin22B.superHMR.closestTSS.bed


##  primary tumor superhmr occupancy of colon primary tumor superhmr
tumors="SH5YSY U87MG H1437 H1672 H157 M673 M697 468PT PC3 22RV1 "

tumorFiles=`for i in $tumors
do
	echo /home/user/bsdata/hmr/superhmr/$i.superHMR.closestTSS.bed
done | tr "\n" " "`

# counts and proportion
bedtools annotate\
 -i /home/user/bsdata/hmr/superhmr/22B.superHMR.closestTSS.bed\
 -files $tumorFiles -names $tumors>\
  /home/user/bsdata/hmr/superhmr/22B.superHMR.closestTSS.allprimarysuperhmrs.bed

# filter cancer super-HMRs by its corresponding normal HMRs
while read line
do
	line2=(`echo $line | tr "," " "`)
	sample=${line2[0]}
	control=${line2[5]}
	if [[ "$sample" != "$control" ]]
	then
		bedtools annotate -i /home/user/bsdata/hmr/superhmr/$sample.superHMR.closestTSS.bed -files /home/user/bsdata/hmr/CpG_context_${control}.hmr.bed |\
         awk '{if($NF < .25) print}' > /home/user/bsdata/hmr/superhmr/$sample.superHMR.closestTSS.normalfilt25.bed
		bedtools annotate -i /home/user/bsdata/hmr/superhmr/$sample.superHMR.closestTSS.bed -files /home/user/bsdata/hmr/CpG_context_${control}.hmr.bed |\
         awk '{if($NF < .5) print}' > /home/user/bsdata/hmr/superhmr/$sample.superHMR.closestTSS.normalfilt50.bed
	fi
done < /home/user/bsdata/samples.info.newnames.control.csv


#######################
####    Super Enhancers
#######################

# set wd
mkdir -p /home/user/bsdata/superenhancers
cd /home/user/bsdata/superenhancers

# Function to characterize superenhancers in terms of methylation and annotation
function superenhancers {

	# Get arguments
	allargs=($@)
	tissue=${allargs[0]}
	super=${allargs[1]}
	suf=${allargs[2]}
	samples=${allargs[@]:3:${#allargs[@]}}

	# set threshold for HMR methylation
	th=100

	for sample in $samples
	do
		
		# Count the number of bases covered by HMRs
		bedtools intersect -wo\
         -a $super\
         -b <(awk -v thres=$th '$10/$7 < thres/100 {print}' ../hmr/CpG_context_${sample}.hmr.metrics.bed) |\
         bedtools groupby -g 1,2,3,4,5,6 -c 18 -o sum >\
         $tissue.superenhancers.$sample.overlapbases.thres$th.${suf}
	done

	# Compute HMR overlap of superenhancers
	scripts/superenhancers.r\
     $tissue.superenhancers.methgain.hmroverlap.$th.${suf}\
     $tissue.superenhancers.*.overlapbases.thres$th.${suf}

	## annotate with the closest gene TSS
	# first print header
	fields=`awk -v file=$tissue.superenhancers.methgain.hmroverlap.$th.annotation.closestTSS.$suf '{OFS="\t"; print $0, "transcript", "gene", "type", "distance" > file; print NF; exit}' $tissue.superenhancers.methgain.hmroverlap.$th.$suf`

	# fields to group by
	grouping=`seq -s "," 1 $fields`
	cols=`seq -s "," \`expr $fields + 5\` \`expr $fields + 8\``

	# find closest TSS and arrange result
	bedtools closest -d\
     -a $tissue.superenhancers.methgain.hmroverlap.$th.$suf\
     -b $tss |\
     bedtools groupby -g $grouping -c $cols -o collapse,collapse,collapse,collapse >>\
     $tissue.superenhancers.methgain.hmroverlap.$th.annotation.closestTSS.$suf

	# annotate with all the TSS in a 1Mb flanking region to the mid of superenhancer
	# first print header
	fields=`awk -v file=$tissue.superenhancers.methgain.hmroverlap.$th.annotation.1MBflankingMid.$suf '{OFS="\t"; print $0, "gene" > file; print NF; exit}' $tissue.superenhancers.methgain.hmroverlap.$th.${suf}`

	# fields to group by
	grouping=`seq -s "," 1 $fields`
	cols=`expr $fields + 6`

	# find closest TSS and arrange result
	bedtools intersect -wao\
     -a <(awk 'BEGIN{OFS="\t"}{mid=int(($2 + $3)/2 +.5); start = mid - 1000000; end = mid + 1000000; if(start < 1) start = 1; $2 = start; $3 = end; print}' $tissue.superenhancers.methgain.hmroverlap.$th.$suf)\
     -b $tss |\
     bedtools groupby -g $grouping -c $cols -o distinct >>\
     $tissue.superenhancers.methgain.hmroverlap.$th.annotation.1MBflankingMid.${suf}

	# find closest TSS and arrange result
	bedtools intersect -wao\
     -a <(awk 'BEGIN{OFS="\t"}{start = $2- 100000; end = $3 + 100000; if(start < 1) start = 1; $2 = start; $3 = end; print}' $tissue.superenhancers.methgain.hmroverlap.$th.$suf)\
     -b $tss |\
     bedtools groupby -g $grouping -c $cols -o distinct >>\
     $tissue.superenhancers.methgain.hmroverlap.$th.annotation.100KBflanking.${suf}
}

## prepare SE files

awk 'BEGIN{FS=","; OFS="\t"}$7 == 1{print substr($2, 4, 100), $3, $4, $1, $8, $9}'\
 /home/user/resources/human/histones/breast/HMEC.csv |\
 bedtools sort -i - >\
 /home/user/resources/human/histones/breast/HMEC.bed

awk 'BEGIN{FS=","; OFS="\t"}$7 == 1{print substr($2, 4, 100), $3, $4, $1, $8, $9}'\
 /home/user/resources/human/histones/colon/UCSD_Sigmoid_Colon.csv |\
 bedtools sort -i - >\
 /home/user/resources/human/histones/colon/UCSD_Sigmoid_Colon.bed

awk 'BEGIN{FS=","; OFS="\t"}$7 == 1{print substr($2, 4, 100), $3, $4, $1, $8, $9}'\
 /home/user/resources/human/histones/lung/UCSD_Lung.csv |\
 bedtools sort -i - >\
 /home/user/resources/human/histones/lung/UCSD_Lung.bed

awk 'BEGIN{FS=","; OFS="\t"}$7 == 1{print substr($2, 4, 100), $3, $4, $1, $8, $9}'\
 /home/user/resources/human/histones/brain/BI_Brain_Inferior_Temporal_Lobe.csv |\
 bedtools sort -i - >\
 /home/user/resources/human/histones/brain/BI_Brain_Inferior_Temporal_Lobe.bed

awk 'BEGIN{FS=","; OFS="\t"}$7 == 1{print substr($2, 4, 100), $3, $4, $1, $8, $9}'\
 /home/user/resources/human/histones/brain/glioma/u87.csv |\
 bedtools sort -i - >\
 /home/user/resources/human/histones/brain/glioma/u87.bed

awk 'BEGIN{FS=","; OFS="\t"}$7 == 1{print substr($2, 4, 100), $3, $4, $1, $8, $9}'\
 /home/user/resources/human/histones/colon/tumor/HCT-116.csv |\
 bedtools sort -i - >\
 /home/user/resources/human/histones/colon/tumor/HCT-116.bed

## all genes
# TSS of genes file
tss=/home/user/resources/human/gencode.v16.alltranscripts.TSS.nochr.bed


#  Breast

super=/home/user/resources/human/histones/breast/HMEC.bed

time superenhancers breast $super "bed" B11-27235 468LN 468PT

#  Colon

super=/home/user/resources/human/histones/colon/UCSD_Sigmoid_Colon.bed

time superenhancers colon $super "bed" 22A 22B 22C

#  Lung

super=/home/user/resources/human/histones/lung/UCSD_Lung.bed

time superenhancers lung $super "bed" 617N H1437 H1672 H157

#  Brain

super=/home/user/resources/human/histones/brain/BI_Brain_Inferior_Temporal_Lobe.bed

time superenhancers brain $super "bed" G145 W145 SH5YSY U87MG

#  Glioblastoma

super=/home/user/resources/human/histones/brain/glioma/u87.bed

time superenhancers glioma $super "bed" W145 U87MG

#  Colon Cancer

super=/home/user/resources/human/histones/colon/tumor/HCT-116.bed

time superenhancers colontumor $super "bed" 22A 22B 22C

## all genes
# TSS of genes file
tss=/home/user/resources/human/gencode.v16.alltranscripts.TSS.nochr.c.bed


#  Breast

super=/home/user/resources/human/histones/breast/HMEC.bed

time superenhancers breast $super c.bed B11-27235 468LN 468PT

#  Colon

super=/home/user/resources/human/histones/colon/UCSD_Sigmoid_Colon.bed

time superenhancers colon $super c.bed 22A 22B 22C

#  Lung

super=/home/user/resources/human/histones/lung/UCSD_Lung.bed

time superenhancers lung $super c.bed 617N H1437 H1672 H157

#  Brain

super=/home/user/resources/human/histones/brain/BI_Brain_Inferior_Temporal_Lobe.bed

time superenhancers brain $super c.bed G145 W145 SH5YSY U87MG

#  Glioblastoma

super=/home/user/resources/human/histones/brain/glioma/u87.bed

time superenhancers glioma $super c.bed W145 U87MG

#  Colon Cancer

super=/home/user/resources/human/histones/colon/tumor/HCT-116.bed

time superenhancers colontumor $super c.bed 22A 22B 22C





############################
####    Chromatin State
############################

# create directory
mkdir -p chromatineState/

# create file codes
echo NY108 Hepg2 >> chromatineState/codes.txt
echo Laia Gm12878 >> chromatineState/codes.txt
echo B11-27235 Hmec >> chromatineState/codes.txt

# for each sample in chromatineState/codes.txt

while read line
do

# get line elements as array
arr=($line)

# set id and file names
bsname=${arr[0]}
chsname=${arr[1]}

common=$bsname.common.hmr.bed
specific=$bsname.hmr.specific.bed
all=CpG_context_$bsname.hmr.bed

# per region
for i in `seq 1 15`
do
# set chromatine state file
chsfile=/home/user/resources/human/encode/wgEncodeBroadHmm${chsname}HMM.${i}_*.bed

# get chromatine state regions with HMRs
bedtools intersect -u -a $chsfile -b $common > chromatineState/$bsname.common.hmr.in.$i.chromatineState.bed
bedtools intersect -u -a $chsfile -b $specific > chromatineState/$bsname.specific.hmr.in.$i.chromatineState.bed
bedtools intersect -u -a $chsfile -b $all > chromatineState/$bsname.all.hmr.in.$i.chromatineState.bed

# count number of regions
Ncommon=`cat chromatineState/$bsname.common.hmr.in.$i.chromatineState.bed | wc -l`
Nspecific=`cat chromatineState/$bsname.specific.hmr.in.$i.chromatineState.bed | wc -l`
Nall=`cat chromatineState/$bsname.all.hmr.in.$i.chromatineState.bed | wc -l`

NN=`cat $chsfile | wc -l`

# split output
echo $i $Nall $Ncommon $Nspecific $NN
done > chromatineState/$bsname.chrState.with.hmr.txt


# per hmr
# set initial tmp files
echo -e "1\t0\t1" > tmp.common
echo -e "1\t0\t1" > tmp.specific
echo -e "1\t0\t1" > tmp.all

for i in `seq 1 15`
do
# set chromatine state file
chsfile=/home/user/resources/human/encode/wgEncodeBroadHmm${chsname}HMM.${i}_*.bed

# get HMRs with chromatine state regions (respecting hierarchy)
bedtools intersect -v -a $common -b tmp.common | bedtools intersect -u -a - -b $chsfile > chromatineState/$i.chromatineState.in.$bsname.common.hmr.bed
bedtools intersect -v -a $specific -b tmp.specific | bedtools intersect -u -a - -b $chsfile > chromatineState/$i.chromatineState.in.$bsname.specific.hmr.bed
bedtools intersect -v -a $all -b tmp.all | bedtools intersect -u -a - -b $chsfile > chromatineState/$i.chromatineState.in.$bsname.all.hmr.bed

# count number of HMRs
Ncommon=`cat chromatineState/$i.chromatineState.in.$bsname.common.hmr.bed | wc -l`
Nspecific=`cat chromatineState/$i.chromatineState.in.$bsname.specific.hmr.bed | wc -l`
Nall=`cat chromatineState/$i.chromatineState.in.$bsname.all.hmr.bed | wc -l`

# set tmp files
cut -f 1-3 chromatineState/$i.chromatineState.in.$bsname.common.hmr.bed >> tmp.common
cut -f 1-3 chromatineState/$i.chromatineState.in.$bsname.specific.hmr.bed >> tmp.specific
cut -f 1-3 chromatineState/$i.chromatineState.in.$bsname.all.hmr.bed >> tmp.all

# split output
echo $i $Nall $Ncommon $Nspecific
done > chromatineState/$bsname.hmr.with.chrState.txt

# append total number of HMRs
echo 0 `cat $all |wc -l` `cat $common |wc -l` `cat $specific |wc -l` >> chromatineState/$bsname.hmr.with.chrState.txt

# remove tmp files
rm tmp.common tmp.all tmp.specific

done < chromatineState/codes.txt


#### shuffle HMR to make some enrichment tests

# create directory
mkdir -p chromatineState/shuffle

# for each sample with chromatine states info
while read line
do

# get line elements as array
arr=($line)

# set id and file names
bsname=${arr[0]}
chsname=${arr[1]}

echo
echo $bsname
echo

time for i in `seq 1 1000`
do
mktemp
done | parallel -j 40 "scripts/chromatine.state.shuffleHMR.sh $bsname $chsname {} chromatineState/shuffle/"

echo
echo Done
echo

done < chromatineState/codes.txt




##############################
####    Annotate 450k manifest
##############################

# create a bed file with the CG array probe positions
awk -F "," 'BEGIN{for(i=1; i<=22; i++){a[i] = i}; a["X"]="X"}NR>8{OFS="\t"; if($12 in a) print $12, $13, $13 + 1, $2}' /home/user/resources/human/illumina-methyl-450k-manifest.csv | bedtools sort -i - >\
 /home/user/resources/human/450k.manifest.1to22andX.onlypos.bed
bgzip -c /home/user/resources/human/450k.manifest.1to22andX.onlypos.bed > /home/user/resources/human/450k.manifest.1to22andX.onlypos.bed.gz
tabix -p bed /home/user/resources/human/450k.manifest.1to22andX.onlypos.bed.gz

# get normal codes
awk '$1 != "W145" && $1 != "PL5" && $1 != "103" && $1 != "NB" && $1 != "Julia" && substr($3, 1, 3) == "Nor" {OFS="\t"; print $2, $1}' ../samples.info.txt  > ../normals.txt

# get CG probes overlapping specific HMRs per normal sample
while read line
do
arr=($line)
tissue=${arr[0]}
sample=${arr[1]}
bedtools intersect -wo -a /home/user/resources/human/450k.manifest.1to22andX.onlypos.bed.gz -b $sample.hmr.specific.bed | cut -f 1-4,8 > /home/user/resources/human/450k.probes.$tissue.specific.bed
done < ../normals.txt

# get CG probes overlapping common HMRs
awk '{OFS="\t"; print $0, "HYPO"NR}' normal.common.hmr.merged100b.bed | bedtools intersect -wo -a /home/user/resources/human/450k.manifest.1to22andX.onlypos.bed.gz -b - >\
 /home/user/resources/human/450k.probes.common.hmr.bed


##################
####    450k array
##################

##  Get selected array samples
scripts/select.array.samples.r

# get number of HMR with at least one probe
bedtools intersect -wao -a normal.all.hmr.bed\
 -b /home/user/resources/human/450k.manifest.1to22andX.onlypos.bed |\
 awk '{a[$1"\t"$2"\t"$3] += $NF}END{for(i in a) print i "\t" a[i]}' | bedtools sort -i - >\
 normal.all.hmr.450kprobes.bed

awk '{if($NF > 0) a += 1}END{print a/NR}' normal.all.hmr.450kprobes.bed

################
##  genomestudio
################

# get 450k betas
scripts/get.normal.array.genomestudio.r
scripts/get.cancer.array.genomestudio.r

# compress and index bed files in /home/user(bsdata/array/tissue
for i in /home/user/bsdata/array/tissue/*bed
do
echo "bgzip -c $i > $i.gz; tabix -p bed $i.gz"
done | parallel --gnu

# transform coordinates to match bs-data
for i in ../array/450k.betas.*.genomestudio.bed; do new=`echo $i | sed 's/450k.betas/450k.betas.plusone/1'` ; awk '{OFS="\t"; $2 = $2 + 1; $3 = $3 + 1; print }' $i > $new; done


# get average methylation per common HMR
########################################

region=normal.common.hmr.merged100b.bed

for i in ../array/450k.betas.plusone*genomestudio.bed
do
name=`basename $i .bed | sed 's/450k.betas.plusone.//1'`
echo "grep -v NA ../array/450k.betas.plusone.$name.bed | bedtools intersect -wo -a $region -b - | bedtools groupby -g 1,2,3 -c 8,8 -o mean,count > ../array/common.hmr.$name.bed"
done  | parallel --gnu -j 40

# merge all tissue specific HMRs (Adding a tissue id)
# awk '{split(FILENAME, a, "."); $4 = a[1]"_"$4; OFS="\t"; print}' *.hmr.specific.bed | bedtools sort -i - > normals.hmr.specific.bed

# get average methylation per tissue specific HMR
########################################

region=normals.hmr.specific.bed

for i in ../array/450k.betas.plusone*genomestudio.bed
do
name=`basename $i .bed | sed 's/450k.betas.plusone.//1'`
echo "grep -v NA ../array/450k.betas.plusone.$name.bed | bedtools intersect -wo -a $region -b - | bedtools groupby -g 1,2,3,4 -c 16,16 -o mean,count > ../array/specific.hmr.$name.bed"
done  | parallel --gnu -j 40

# tissue specific shrinked HMR

region=normals.shrinkedhmr.specific.bed

for i in ../array/450k.betas.plusone*genomestudio.bed
do
name=`basename $i .bed | sed 's/450k.betas.plusone.//1'`
echo "grep -v NA ../array/450k.betas.plusone.$name.bed | bedtools intersect -wo -a $region -b - | bedtools groupby -g 1,2,3,4 -c 16,16 -o mean,count > ../array/specific.shrinkedhmr.$name.bed"
done  | parallel --gnu -j 40


####    get average methylation per common HMR
##  Loose: more than 0 probes per HMR
for i in /home/user/bsdata/array/tissue/450k.betas*genomestudio.bed.gz
do
name=`basename $i .genomestudio.bed.gz | sed 's/450k.betas.//1'`
echo $name
done |\
parallel --gnu -j 40 \
 "tabix -B /home/user/bsdata/array/tissue/450k.betas.{}.genomestudio.bed.gz  normal.common.hmr.merged100b.ID.bed | bedtools intersect -wao -a normal.common.hmr.merged100b.ID.bed -b - | scripts/aggregate.per.bed.awk -v TH=0 -v K=4 -v FIRST=9 > ../array/common.{}.loose.bed"

##  Mid: more than 1 probes per HMR
for i in /home/user/bsdata/array/tissue/450k.betas*genomestudio.bed.gz
do
name=`basename $i .genomestudio.bed.gz | sed 's/450k.betas.//1'`
echo $name
done |\
parallel --gnu -j 40 \
 "tabix -B /home/user/bsdata/array/tissue/450k.betas.{}.genomestudio.bed.gz  normal.common.hmr.merged100b.ID.bed | bedtools intersect -wao -a normal.common.hmr.merged100b.ID.bed -b - | scripts/aggregate.per.bed.awk -v TH=1 -v K=4 -v FIRST=9 > ../array/common.{}.mid.bed"

##  Strict: more than 4 probes per HMR
for i in /home/user/bsdata/array/tissue/450k.betas*genomestudio.bed.gz
do
name=`basename $i .genomestudio.bed.gz | sed 's/450k.betas.//1'`
echo $name
done |\
parallel --gnu -j 40 \
 "tabix -B /home/user/bsdata/array/tissue/450k.betas.{}.genomestudio.bed.gz  normal.common.hmr.merged100b.ID.bed | bedtools intersect -wao -a normal.common.hmr.merged100b.ID.bed -b - | scripts/aggregate.per.bed.awk -v TH=4 -v K=4 -v FIRST=9 > ../array/common.{}.strict.bed"

####    get average methylation per tissue specific HMR
##  Loose: more than 0 probes per HMR
for i in /home/user/bsdata/array/tissue/450k.betas*genomestudio.bed.gz
do
name=`basename $i .genomestudio.bed.gz | sed 's/450k.betas.//1'`
echo $name
done |\
parallel --gnu -j 40 \
 "tabix -B /home/user/bsdata/array/tissue/450k.betas.{}.genomestudio.bed.gz  normal.specific.hmr.ID.bed | bedtools intersect -wao -a normal.specific.hmr.ID.bed -b - | scripts/aggregate.per.bed.awk -v TH=0 -v K=4 -v FIRST=9 > ../array/specific.{}.loose.bed"

##  Mid: more than 1 probes per HMR
for i in /home/user/bsdata/array/tissue/450k.betas*genomestudio.bed.gz
do
name=`basename $i .genomestudio.bed.gz | sed 's/450k.betas.//1'`
echo $name
done |\
parallel --gnu -j 40 \
 "tabix -B /home/user/bsdata/array/tissue/450k.betas.{}.genomestudio.bed.gz  normal.specific.hmr.ID.bed | bedtools intersect -wao -a normal.specific.hmr.ID.bed -b - | scripts/aggregate.per.bed.awk -v TH=1 -v K=4 -v FIRST=9 > ../array/specific.{}.mid.bed"

##  Strict: more than 4 probes per HMR
for i in /home/user/bsdata/array/tissue/450k.betas*genomestudio.bed.gz
do
name=`basename $i .genomestudio.bed.gz | sed 's/450k.betas.//1'`
echo $name
done |\
parallel --gnu -j 40 \
 "tabix -B /home/user/bsdata/array/tissue/450k.betas.{}.genomestudio.bed.gz  normal.specific.hmr.ID.bed | bedtools intersect -wao -a normal.specific.hmr.ID.bed -b - | scripts/aggregate.per.bed.awk -v TH=4 -v K=4 -v FIRST=9 > ../array/specific.{}.strict.bed"

# arrange results
scripts/normals.hmr.array.genomestudio.r

###   get proportion of HMR with average beta below threshold

# set threshold
thres=40

for i in   ../array/specific.hmr.*.bed; do echo $i; done |\
 parallel --gnu -j 40 "awk -v p=$thres '{if(\$(NF - 1) != \"NA\"){split(\$4, a, \"_\"); if(\$(NF - 1) < p/100) {c[a[1]] += 1}; n[a[1]] += 1}}END{OFS=\"\t\"; split(FILENAME, b, \".\"); for(i in n){ print b[5], i,  0 + c[i],  0 + n[i]}}' {}" >\
 specific.hmr.array.p_below_$thres.tab

for i in   ../array/common.hmr.*.bed; do echo $i; done |\
 parallel --gnu -j 40 "awk -v p=$thres '{if(\$(NF - 1) != \"NA\"){if(\$(NF - 1) < p/100) {c += 1}; n += 1}}END{OFS=\"\t\"; split(FILENAME, b, \".\"); print b[5], 0 + c, 0 + n}' {}" >\
 common.hmr.array.p_below_$thres.tab

for i in   ../array/specific.hmr.*.bed; do echo $i; done |\
 parallel --gnu -j 40 "awk -v p=$thres '{if(\$(NF - 1) != \"NA\" && \$NF > 1){split(\$4, a, \"_\"); if(\$(NF - 1) < p/100) {c[a[1]] += 1}; n[a[1]] += 1}}END{OFS=\"\t\"; split(FILENAME, b, \".\"); for(i in n){ print b[5], i,  0 + c[i],  0 + n[i]}}' {}" >\
 specific.hmr.array.p_below_$thres.2ormoreprobes.tab

for i in   ../array/common.hmr.*.bed; do echo $i; done |\
 parallel --gnu -j 40 "awk -v p=$thres '{if(\$(NF - 1) != \"NA\" && \$NF > 1){if(\$(NF - 1) < p/100) {c += 1}; n += 1}}END{OFS=\"\t\"; split(FILENAME, b, \".\"); print b[5], 0 + c, 0 + n}' {}" >\
 common.hmr.array.p_below_$thres.2ormoreprobes.tab

for i in   ../array/specific.hmr.*.bed; do echo $i; done |\
 parallel --gnu -j 40 "awk -v p=$thres '{if(\$(NF - 1) != \"NA\" && \$NF > 4){split(\$4, a, \"_\"); if(\$(NF - 1) < p/100) {c[a[1]] += 1}; n[a[1]] += 1}}END{OFS=\"\t\"; split(FILENAME, b, \".\"); for(i in n){ print b[5], i,  0 + c[i],  0 + n[i]}}' {}" >\
 specific.hmr.array.p_below_$thres.5ormoreprobes.tab

for i in   ../array/common.hmr.*.bed; do echo $i; done |\
 parallel --gnu -j 40 "awk -v p=$thres '{if(\$(NF - 1) != \"NA\" && \$NF > 4){if(\$(NF - 1) < p/100) {c += 1}; n += 1}}END{OFS=\"\t\"; split(FILENAME, b, \".\"); print b[5], 0 + c, 0 + n}' {}" >\
 common.hmr.array.p_below_$thres.5ormoreprobes.tab

# set threshold
thres=25

for i in   ../array/specific.hmr.*.bed; do echo $i; done |\
 parallel --gnu -j 40 "awk -v p=$thres '{if(\$(NF - 1) != \"NA\"){split(\$4, a, \"_\"); if(\$(NF - 1) < p/100) {c[a[1]] += 1}; n[a[1]] += 1}}END{OFS=\"\t\"; split(FILENAME, b, \".\"); for(i in n){ print b[5], i, 0 + c[i], 0 + n[i]}}' {}" >\
 specific.hmr.array.p_below_$thres.tab

for i in   ../array/common.hmr.*.bed; do echo $i; done |\
 parallel --gnu -j 40 "awk -v p=$thres '{if(\$(NF - 1) != \"NA\"){if(\$(NF - 1) < p/100) {c += 1}; n += 1}}END{OFS=\"\t\"; split(FILENAME, b, \".\"); print b[5], 0 + c, 0 + n}' {}" >\
 common.hmr.array.p_below_$thres.tab

for i in   ../array/specific.hmr.*.bed; do echo $i; done |\
 parallel --gnu -j 40 "awk -v p=$thres '{if(\$(NF - 1) != \"NA\" && \$NF > 1){split(\$4, a, \"_\"); if(\$(NF - 1) < p/100) {c[a[1]] += 1}; n[a[1]] += 1}}END{OFS=\"\t\"; split(FILENAME, b, \".\"); for(i in n){ print b[5], i, 0 + c[i], 0 + n[i]}}' {}" >\
 specific.hmr.array.p_below_$thres.2ormoreprobes.tab

for i in   ../array/common.hmr.*.bed; do echo $i; done |\
 parallel --gnu -j 40 "awk -v p=$thres '{if(\$(NF - 1) != \"NA\" && \$NF > 1){if(\$(NF - 1) < p/100) {c += 1}; n += 1}}END{OFS=\"\t\"; split(FILENAME, b, \".\"); print b[5], 0 + c, 0 + n}' {}" >\
 common.hmr.array.p_below_$thres.2ormoreprobes.tab

for i in   ../array/specific.hmr.*.bed; do echo $i; done |\
 parallel --gnu -j 40 "awk -v p=$thres '{if(\$(NF - 1) != \"NA\" && \$NF > 4){split(\$4, a, \"_\"); if(\$(NF - 1) < p/100) {c[a[1]] += 1}; n[a[1]] += 1}}END{OFS=\"\t\"; split(FILENAME, b, \".\"); for(i in n){ print b[5], i, 0 + c[i], 0 + n[i]}}' {}" >\
 specific.hmr.array.p_below_$thres.5ormoreprobes.tab

for i in   ../array/common.hmr.*.bed; do echo $i; done |\
 parallel --gnu -j 40 "awk -v p=$thres '{if(\$(NF - 1) != \"NA\" && \$NF > 4){if(\$(NF - 1) < p/100) {c += 1}; n += 1}}END{OFS=\"\t\"; split(FILENAME, b, \".\"); print b[5], 0 + c, 0 + n}' {}" >\
 common.hmr.array.p_below_$thres.5ormoreprobes.tab



###   Get proportion of specific shrinked HMR with average beta below threshold

# set threshold
thres=40

for i in   ../array/specific.shrinkedhmr.*.bed; do echo $i; done |\
 parallel --gnu -j 40 "awk -v p=$thres '{if(\$(NF - 1) != \"NA\"){split(\$4, a, \"_\"); if(\$(NF - 1) < p/100) {c[a[1]] += 1}; n[a[1]] += 1}}END{OFS=\"\t\"; split(FILENAME, b, \".\"); for(i in n){ print b[5], i,  0 + c[i],  0 + n[i]}}' {}" >\
 specific.shrinkedhmr.array.p_below_$thres.tab

for i in   ../array/specific.shrinkedhmr.*.bed; do echo $i; done |\
 parallel --gnu -j 40 "awk -v p=$thres '{if(\$(NF - 1) != \"NA\" && \$NF > 1){split(\$4, a, \"_\"); if(\$(NF - 1) < p/100) {c[a[1]] += 1}; n[a[1]] += 1}}END{OFS=\"\t\"; split(FILENAME, b, \".\"); for(i in n){ print b[5], i,  0 + c[i],  0 + n[i]}}' {}" >\
 specific.shrinkedhmr.array.p_below_$thres.2ormoreprobes.tab

for i in   ../array/specific.shrinkedhmr.*.bed; do echo $i; done |\
 parallel --gnu -j 40 "awk -v p=$thres '{if(\$(NF - 1) != \"NA\" && \$NF > 4){split(\$4, a, \"_\"); if(\$(NF - 1) < p/100) {c[a[1]] += 1}; n[a[1]] += 1}}END{OFS=\"\t\"; split(FILENAME, b, \".\"); for(i in n){ print b[5], i,  0 + c[i],  0 + n[i]}}' {}" >\
 specific.shrinkedhmr.array.p_below_$thres.5ormoreprobes.tab

# set threshold
thres=25

for i in   ../array/specific.shrinkedhmr.*.bed; do echo $i; done |\
 parallel --gnu -j 40 "awk -v p=$thres '{if(\$(NF - 1) != \"NA\"){split(\$4, a, \"_\"); if(\$(NF - 1) < p/100) {c[a[1]] += 1}; n[a[1]] += 1}}END{OFS=\"\t\"; split(FILENAME, b, \".\"); for(i in n){ print b[5], i, 0 + c[i], 0 + n[i]}}' {}" >\
 specific.shrinkedhmr.array.p_below_$thres.tab

for i in   ../array/specific.shrinkedhmr.*.bed; do echo $i; done |\
 parallel --gnu -j 40 "awk -v p=$thres '{if(\$(NF - 1) != \"NA\" && \$NF > 1){split(\$4, a, \"_\"); if(\$(NF - 1) < p/100) {c[a[1]] += 1}; n[a[1]] += 1}}END{OFS=\"\t\"; split(FILENAME, b, \".\"); for(i in n){ print b[5], i, 0 + c[i], 0 + n[i]}}' {}" >\
 specific.shrinkedhmr.array.p_below_$thres.2ormoreprobes.tab

for i in   ../array/specific.shrinkedhmr.*.bed; do echo $i; done |\
 parallel --gnu -j 40 "awk -v p=$thres '{if(\$(NF - 1) != \"NA\" && \$NF > 4){split(\$4, a, \"_\"); if(\$(NF - 1) < p/100) {c[a[1]] += 1}; n[a[1]] += 1}}END{OFS=\"\t\"; split(FILENAME, b, \".\"); for(i in n){ print b[5], i, 0 + c[i], 0 + n[i]}}' {}" >\
 specific.shrinkedhmr.array.p_below_$thres.5ormoreprobes.tab


###   ROC curves

# SPECIFIC HMRs: sensitivity and specificity for different thresholds
for th in `seq 0 1 100`
do
for i in   ../array/specific.hmr.*.genomestudio.bed; do echo $i; done | parallel --gnu -j 40 "awk -v p=$th '{if(\$(NF - 1) != \"NA\" && \$NF > 4){split(\$4, a, \"_\"); if(\$(NF - 1) < p/100) {c[a[1]] += 1}; n[a[1]] += 1}}END{OFS=\"\t\"; split(FILENAME, b,      \".\"); for(i in n){ print b[5], i,  0 + c[i],  0 + n[i]}}' {}" | awk -v p=$th 'BEGIN{while((getline < "/home/user/bsdata/normal_tissue_curated.txt") > 0) {code[$1] = $2}}{n[code[$1]"\t"$2] += $3; t[code[$1]"\t"$2] += $4}END{      OFS="\t"; for(i in n) print i, p, n[i], t[i]}'
done > specific.array.ROC.genomestudio.tab

# SPECIFIC HMRs: sensitivity and specificity for different thresholds (include all HMR with at least one probe)
for th in `seq 0 1 100`
do
for i in   ../array/specific.hmr.*.genomestudio.bed; do echo $i; done | parallel --gnu -j 40 "awk -v p=$th '{if(\$(NF - 1) != \"NA\" && \$NF > 0){split(\$4, a, \"_\"); if(\$(NF - 1) < p/100) {c[a[1]] += 1}; n[a[1]] += 1}}END{OFS=\"\t\"; split(FILENAME, b, \".\");
for(i in n){ print b[5], i,  0 + c[i],  0 + n[i]}}' {}" | awk -v p=$th 'BEGIN{while((getline < "/home/user/bsdata/normal_tissue_curated.txt") > 0) {code[$1] = $2}}{n[code[$1]"\t"$2] += $3; t[code[$1]"\t"$2] += $4}END{
OFS="\t"; for(i in n) print i, p, n[i], t[i]}'
done > specific.array.loose.ROC.genomestudio.tab

# SPECIFIC HMRs: sensitivity and specificity for different thresholds (include all HMR with at least 2 probes)
for th in `seq 0 1 100`
do
for i in   ../array/specific.hmr.*.genomestudio.bed; do echo $i; done | parallel --gnu -j 40 "awk -v p=$th '{if(\$(NF - 1) != \"NA\" && \$NF > 1){split(\$4, a, \"_\"); if(\$(NF - 1) < p/100) {c[a[1]] += 1}; n[a[1]] += 1}}END{OFS=\"\t\"; split(FILENAME, b,   \".\");
for(i in n){ print b[5], i,  0 + c[i],  0 + n[i]}}' {}" | awk -v p=$th 'BEGIN{while((getline < "/home/user/bsdata/normal_tissue_curated.txt") > 0) {code[$1] = $2}}{n[code[$1]"\t"$2] += $3; t[code[$1]"\t"$2] += $4}END{
OFS="\t"; for(i in n) print i, p, n[i], t[i]}'
done > specific.array.mid.ROC.genomestudio.tab

# COMMON HMRs: sensitivity for different thresholds
for th in `seq 0 1 100`
do
for i in   ../array/common.hmr.*.genomestudio.bed; do echo $i; done |\
 parallel --gnu -j 40 "awk -v p=$th '{if(\$(NF - 1) != \"NA\" && \$NF > 4){if(\$(NF - 1) < p/100) {c += 1}; n += 1}}END{OFS=\"\t\"; split(FILENAME, b, \".\"); print b[5], p, 0 + c, 0 + n}' {}"
done >  common.array.ROC.genomestudio.tab

# COMMON HMRs: sensitivity for different thresholds
for th in `seq 0 1 100`
do
for i in   ../array/common.hmr.*.genomestudio.bed; do echo $i; done |\
 parallel --gnu -j 40 "awk -v p=$th '{if(\$(NF - 1) != \"NA\" && \$NF > 1){if(\$(NF - 1) < p/100) {c += 1}; n += 1}}END{OFS=\"\t\"; split(FILENAME, b, \".\"); print b[5], p, 0 + c, 0 + n}' {}"
done >  common.array.mid.ROC.genomestudio.tab

# COMMON HMRs: sensitivity for different thresholds
for th in `seq 0 1 100`
do
for i in   ../array/common.hmr.*.genomestudio.bed; do echo $i; done |\
 parallel --gnu -j 40 "awk -v p=$th '{if(\$(NF - 1) != \"NA\" && \$NF > 0){if(\$(NF - 1) < p/100) {c += 1}; n += 1}}END{OFS=\"\t\"; split(FILENAME, b, \".\"); print b[5], p, 0 + c, 0 + n}' {}"
done >  common.array.loose.ROC.genomestudio.tab


####    Clustering of normal array samples based on HMRs
########################################################

# prepare all tissue specific hmr together
awk '{OFS="\t"; $4 = $4"."$1; print $1, $2, $3, $4}' normals.hmr.specific.bed > all.hmr.unsorted.bed
awk '{OFS="\t"; print $1, $2, $3, "common_"NR}' normal.common.hmr.merged100b.bed >> all.hmr.unsorted.bed
bedtools sort -i all.hmr.unsorted.bed > all.hmr.bed

# get array probes of tissue specific HMRs
bedtools intersect -wo -b all.hmr.bed -a /home/user/resources/human/450k.manifest.1to22andX.onlypos.bed | cut -f 4,8 > all.hmr.cgs.tab

# add HMR id
awk '{OFS="\t"; print $0, "HYPO_"NR}' normal.common.hmr.merged100b.bed > normal.common.hmr.merged100b.withid.bed

# get array probes of common HMRs
bedtools intersect -wo -b all.hmr.bed -a /home/user/resources/human/450k.manifest.1to22andX.onlypos.bed | cut -f 4,8 > all.hmr.cgs.tab
grep common all.hmr.cgs.tab > common.hmr.cgs.tab
 
###############################################
####    Validation of t-HMR with cancer samples
###############################################

# Blood
refsample=Laia.hmr.specific.bed

grep Blood ../samples.info.txt | cut -d " " -f 1 |\
 parallel --gnu -j 10 bedtools intersect -c -a $refsample -b CpG_context_{}.hmr.bed |\
 awk -v n=11 '{key=$1; for(i=2; i<=n; i++) key = key"\t"$i; out[key] += 0 + ($NF!=0)}END{for(i in out) print i "\t" out[i]}' | bedtools sort -i - >\
 Laia.hmr.specific.seqvalidation.bed 

# Brain
refsample=G145.hmr.specific.bed

grep Brain ../samples.info.txt | cut -d " " -f 1 |\
 parallel --gnu -j 10 bedtools intersect -c -a $refsample -b CpG_context_{}.hmr.bed |\
 awk -v n=11 '{key=$1; for(i=2; i<=n; i++) key = key"\t"$i; out[key] += 0 + ($NF!=0)}END{for(i in out) print i "\t" out[i]}' | bedtools sort -i - >\
 G145.hmr.specific.seqvalidation.bed 

# Breast
refsample=B11-27235.hmr.specific.bed

grep Breast ../samples.info.txt | cut -d " " -f 1 |\
 parallel --gnu -j 10 bedtools intersect -c -a $refsample -b CpG_context_{}.hmr.bed |\
 awk -v n=11 '{key=$1; for(i=2; i<=n; i++) key = key"\t"$i; out[key] += 0 + ($NF!=0)}END{for(i in out) print i "\t" out[i]}' | bedtools sort -i - >\
 B11-27235.hmr.specific.seqvalidation.bed 

# Colon
refsample=22A.hmr.specific.bed

grep Colon ../samples.info.txt | cut -d " " -f 1 |\
 parallel --gnu -j 10 bedtools intersect -c -a $refsample -b CpG_context_{}.hmr.bed |\
 awk -v n=11 '{key=$1; for(i=2; i<=n; i++) key = key"\t"$i; out[key] += 0 + ($NF!=0)}END{for(i in out) print i "\t" out[i]}' | bedtools sort -i - >\
 22A.hmr.specific.seqvalidation.bed 

# Liver
refsample=NY108.hmr.specific.bed

grep Liver ../samples.info.txt | cut -d " " -f 1 |\
 parallel --gnu -j 10 bedtools intersect -c -a $refsample -b CpG_context_{}.hmr.bed |\
 awk -v n=11 '{key=$1; for(i=2; i<=n; i++) key = key"\t"$i; out[key] += 0 + ($NF!=0)}END{for(i in out) print i "\t" out[i]}' | bedtools sort -i - >\
 NY108.hmr.specific.seqvalidation.bed 

# Lung
refsample=617N.hmr.specific.bed

grep Lung ../samples.info.txt | cut -d " " -f 1 |\
 parallel --gnu -j 10 bedtools intersect -c -a $refsample -b CpG_context_{}.hmr.bed |\
 awk -v n=11 '{key=$1; for(i=2; i<=n; i++) key = key"\t"$i; out[key] += 0 + ($NF!=0)}END{for(i in out) print i "\t" out[i]}' | bedtools sort -i - >\
 617N.hmr.specific.seqvalidation.bed 

# Prostate
refsample=CP10.hmr.specific.bed

grep Prostate ../samples.info.txt | cut -d " " -f 1 |\
 parallel --gnu -j 10 bedtools intersect -c -a $refsample -b CpG_context_{}.hmr.bed |\
 awk -v n=11 '{key=$1; for(i=2; i<=n; i++) key = key"\t"$i; out[key] += 0 + ($NF!=0)}END{for(i in out) print i "\t" out[i]}' | bedtools sort -i - >\
 CP10.hmr.specific.seqvalidation.bed 



#############################################################
####    Annotate HMRs with coding genes at and near promoters
#############################################################

# get common HMRs that overlap promoters
coding=/home/user/resources/human/gencode.v16.alltranscripts.coding.prom.nochr.bed

bedtools intersect -u -a $coding -b normal.common.hmr.merged100b.bed > genomicFeatures/promoters.overlapping.common.hmr.bed

# get c-HMR with promoter annotation
awk '{OFS="\t"; print $0, "common_"NR}' normal.common.hmr.merged100b.bed | bedtools intersect -wo -b $coding -a - |\
 awk '{key=$1"\t"$2"\t"$3"\t"$4; trans[key] = trans[key]";"$9; gene[key] = gene[key]";"$10; cons[key] = cons[key]";"$11}END{OFS="\t"; for(j in gene) print j, trans[j], gene[j], cons[j]}' |\
 bedtools sort -i - | sed 's/\t;/\t/g' > normal.common.hmr.merged100b.atpromoters.bed

# get c-HMR with gene annotation for those not overlapping promoters
awk '{OFS="\t"; print $0, "common_"NR}' normal.common.hmr.merged100b.bed | bedtools closest -D "b" -b $coding -a - |\
 awk '$12 != 0 {key=$1"\t"$2"\t"$3"\t"$4; trans[key] = trans[key]";"$9; gene[key] = gene[key]";"$10; cons[key] = cons[key]";"$11; dis[key] = dis[key]";"$12}END{OFS="\t"; for(j in gene) print j, trans[j], gene[j], cons[j], dis[j]}' |\
 bedtools sort -i - | sed 's/\t;/\t/g' > normal.common.hmr.merged100b.nearpromoters.bed

# get 450k probes for the c-HMR with promoter annotation
tabix -B /home/user/resources/human/450k.manifest.1to22andX.onlypos.bed.gz normal.common.hmr.merged100b.atpromoters.bed |\
 bedtools intersect -wo -a normal.common.hmr.merged100b.atpromoters.bed -b - | cut -f 4,11 >\
 normal.common.hmr.merged100b.atpromoters.450kprobes.txt

# get 450k probes for the c-HMR near promoter annotation
tabix -B /home/user/resources/human/450k.manifest.1to22andX.onlypos.bed.gz normal.common.hmr.merged100b.nearpromoters.bed |\
 bedtools intersect -wo -a normal.common.hmr.merged100b.nearpromoters.bed -b - | cut -f 4,11 >\
 normal.common.hmr.merged100b.nearpromoters.450kprobes.txt

# get t-HMR with promoter annotation
cut -f 1-4 normal.specific.hmr.bed | bedtools intersect -wo -b $coding -a - |\
 awk '{key=$1"\t"$2"\t"$3"\t"$4; trans[key] = trans[key]";"$9; gene[key] = gene[key]";"$10; cons[key] = cons[key]";"$11}END{OFS="\t"; for(j in gene) print j, trans[j], gene[j], cons[j]}' |\
 bedtools sort -i - | sed 's/\t;/\t/g' > normal.specific.hmr.atpromoters.bed

# get t-HMR with gene annotation for those not overlapping promoters
cut -f 1-4 normal.specific.hmr.bed | bedtools closest -D "b" -b $coding -a - |\
 awk '$12 != 0 {key=$1"\t"$2"\t"$3"\t"$4; trans[key] = trans[key]";"$9; gene[key] = gene[key]";"$10; cons[key] = cons[key]";"$11; dis[key] = dis[key]";"$12}END{OFS="\t"; for(j in gene) print j, trans[j], gene[j], cons[j], dis[j]}' |\
 bedtools sort -i - | sed 's/\t;/\t/g' > normal.specific.hmr.nearpromoters.bed


####################
####    Conservation
####################

# get average phastCons per HMR
# placental
time awk '{print$1":"$2"-"$3}' all.hmr.bed | parallel -j40 "tabix /home/user/resources/human/phastCons.placental.bed.gz {} | awk -v region={} 'BEGIN{p=0}{p += \$(NF)}END{print region, p, NR}'" | sed 's/[:\ \-]/\t/g' | bedtools sort -i >  all.hmr.phastCons.placental.txt

# primate
time awk '{print$1":"$2"-"$3}' all.hmr.bed | parallel -j40 "tabix /home/user/resources/human/phastCons.primate.bed.gz {} | awk -v region={} 'BEGIN{p=0}{p += \$(NF)}END{print region, p, NR}'" | sed 's/[:\ \-]/\t/g' | bedtools sort -i >  all.hmr.phastCons.primate.txt

# vertebrate
time awk '{print$1":"$2"-"$3}' all.hmr.bed | parallel -j40 "tabix /home/user/resources/human/phastCons.vertebrate.bed.gz {} | awk -v region={} 'BEGIN{p=0}{p += \$(NF)}END{print region, p, NR}'" | sed 's/[:\ \-]/\t/g' | bedtools sort -i >  all.hmr.phastCons.vertebrate.txt

# get average phyloP per HMR
# placental
time awk '{print$1":"$2"-"$3}' all.hmr.bed | parallel -j40 "tabix /home/user/resources/human/phyloP.placental.bed.gz {} | awk -v region={} 'BEGIN{p=0}{p += \$(NF)}END{gsub(\"-\", \"\t\", region); print region, p, NR}'" | sed 's/[:\ ]/\t/g' | bedtools sort -i >  all.hmr.phyloP.placental.txt

# primate
time awk '{print$1":"$2"-"$3}' all.hmr.bed | parallel -j40 "tabix /home/user/resources/human/phyloP.primate.bed.gz {} | awk -v region={} 'BEGIN{p=0}{p += \$(NF)}END{gsub(\"-\", \"\t\", region); print region, p, NR}'" |\
 sed 's/[:\ ]/\t/g' | bedtools sort -i >  all.hmr.phyloP.primate.txt

# vertebrate
time awk '{print$1":"$2"-"$3}' all.hmr.bed | parallel -j40 "tabix /home/user/resources/human/phyloP.vertebrate.bed.gz {} | awk -v region={} 'BEGIN{p=0}{p += \$(NF)}END{gsub(\"-\", \"\t\", region); print region, p, NR}'" | sed 's/[:\ ]/\t/g' | bedtools sort -i >  all.hmr.phyloP.vertebrate.txt


###############################################################
####    genes with c-HMR (@ prom) gaining methylation in cancer
###############################################################

# get number of probes per HMR
commonHMR=/home/user/bsdata/hmr/normal.common.hmr.merged100b.ID.bed
manifest=/home/user/resources/human/450k.manifest.1to22andX.onlypos.bed.gz

tabix -B $manifest $commonHMR |\
 bedtools intersect -wo -a $commonHMR -b - |\
 bedtools groupby -g 4 -c 8 -o count >\
 /home/user/bsdata/hmr/normal.common.hmr.merged100b.ID.numberof450kprobes.txt

# get sample names
for i in ../array/tissue/450k.betas*bed.gz
do
name=`basename $i .genomestudio.bed.gz | sed 's/450k.betas.//1'`
zcat $i | head -n 1 | cut -f 5- | awk -v name=$name 'BEGIN{RS="\t"; ORS="\n"; OFS="\t"}{print name, NR, $1}'
done > ../array/tissue/sample.names

##  HMRs overlapping promoters
# create outdir
mkdir -p methgain.genomestudio.at

# count number of cancer samples that present an average GpG probe methylation above .33 for each gene/c-HMR
for i in /home/user/bsdata/array/tissue/450k.betas.cancer*genomestudio.bed.gz
do
name=`basename $i .genomestudio.bed.gz | sed 's/450k.betas.cancer.//1'`
tabix -B $i normal.common.hmr.merged100b.atpromoters.bed |\
 bedtools intersect -wao -a normal.common.hmr.merged100b.atpromoters.bed -b - |\
 awk -v thres=.33\
 'BEGIN{OFS="\t"; out=""}{key=$1; for(i=2; i<= 7; i ++) key = key "\t" $i; if(key != oldkey && NR > 1) {if(n[12] > 0 && m[12] != "NA"){out = m[12]/n[12]}else{out = "NA"}for(j = 13; j <= (NF - 1); j++) {if(n[j] > 0 && m[j] != "NA") {out = out "\t" m[j]/n[j]}else{out = out "\t" "NA"}}; print oldkey, out; split("", n);  split("", m)}; for(j = 12; j <= (NF - 1); j++) {if($j != "NA" && $j != "." && $j != -1){m[j] += $j; n[j] += 1}}; oldkey=key}' > methgain.genomestudio.at/common.cancer.$name.gain.bed
done

# count number of normal samples that present an average GpG probe methylation above .33 for each gene/c-HMR
for i in `find /home/user/bsdata/array/tissue/* -name "*genomestudio.bed.gz" -a ! -name "*cancer*"`
do
name=`basename $i genomestudio.bed.gz | sed 's/450k.betas.//1'`

tabix -B $i normal.common.hmr.merged100b.atpromoters.bed |\
 bedtools intersect -wao -a normal.common.hmr.merged100b.atpromoters.bed -b - |\
 awk -v thres=.33\
 'BEGIN{OFS="\t"; out=""}{key=$1; for(i=2; i<= 7; i ++) key = key "\t" $i; if(key != oldkey && NR > 1) {if(n[12] > 0 && m[12] != "NA"){out = m[12]/n[12]}else{out = "NA"}for(j = 13; j <= (NF - 1); j++) {if(n[j] > 0 && m[j] != "NA") {out = out "\t" m[j]/n[j]}else{out = out "\t" "NA"}}; print oldkey, out; split("", n);  split("", m)}; for(j = 12; j <= (NF - 1); j++) {if($j != "NA" && $j != "." && $j != -1){m[j] += $j; n[j] += 1}}; oldkey=key}' > methgain.genomestudio.at/common.normal.$name.gain.bed

done


##  HMRs close to but not overlapping promoters
# create outdir
mkdir -p methgain.genomestudio.near

# count number of cancer samples that present an average GpG probe methylation above .33 for each gene/c-HMR
for i in /home/user/bsdata/array/tissue/450k.betas.cancer*genomestudio.bed.gz
do
name=`basename $i .genomestudio.bed.gz | sed 's/450k.betas.cancer.//1'`
tabix -B $i normal.common.hmr.merged100b.nearpromoters.bed |\
 bedtools intersect -wao -a normal.common.hmr.merged100b.nearpromoters.bed -b - |\
 awk\
 'BEGIN{OFS="\t"; out=""}{key=$1; for(i=2; i<= 8; i ++) key = key "\t" $i; if(key != oldkey && NR > 1) {if(n[13] > 0 && m[13] != "NA"){out = m[13]/n[13]}else{out = "NA"}for(j = 14; j <= (NF - 1); j++) {if(n[j] > 0 && m[j] != "NA") {out = out "\t" m[j]/n[j]}else{out = out "\t" "NA"}}; print oldkey, out; split("", n);  split("", m)}; for(j = 13; j <= (NF - 1); j++) {if($j != "NA" && $j != "." && $j != -1){m[j] += $j; n[j] += 1}}; oldkey=key}' > methgain.genomestudio.near/common.cancer.$name.gain.bed
done

# count number of normal samples that present an average GpG probe methylation above .33 for each gene/c-HMR
for i in `find /home/user/bsdata/array/tissue/* -name "*genomestudio.bed.gz" -a ! -name "*cancer*"`
do
name=`basename $i .genomestudio.bed.gz | sed 's/450k.betas.//1'`

tabix -B $i normal.common.hmr.merged100b.nearpromoters.bed |\
 bedtools intersect -wao -a normal.common.hmr.merged100b.nearpromoters.bed -b - |\
 awk\
 'BEGIN{OFS="\t"; out=""}{key=$1; for(i=2; i<= 8; i ++) key = key "\t" $i; if(key != oldkey && NR > 1) {if(n[13] > 0 && m[13] != "NA"){out = m[13]/n[13]}else{out = "NA"}for(j = 14; j <= (NF - 1); j++) {if(n[j] > 0 && m[j] != "NA") {out = out "\t" m[j]/n[j]}else{out = out "\t" "NA"}}; print oldkey, out; split("", n);  split("", m)}; for(j = 13; j <= (NF - 1); j++) {if($j != "NA" && $j != "." && $j != -1){m[j] += $j; n[j] += 1}}; oldkey=key}' > methgain.genomestudio.near/common.normal.$name.gain.bed

done



###############################################################
####    genes with t-HMR (@ prom) gaining methylation in cancer
###############################################################

# HMRs at promoters
# create outdir
mkdir -p methgain.genomestudio.at

# count number of cancer samples that present an average GpG probe methylation above .33 for each gene/c-HMR
for i in /home/user/bsdata/array/tissue/450k.betas.cancer*genomestudio.bed.gz
do
name=`basename $i .genomestudio.bed.gz | sed 's/450k.betas.cancer.//1'`
tabix -B $i normal.specific.hmr.nearpromoters.bed |\
 bedtools intersect -wao -a normal.specific.hmr.nearpromoters.bed -b - |\
 awk -v thres=.33\
 'BEGIN{OFS="\t"; out=""}{key=$1; for(i=2; i<= 7; i ++) key = key "\t" $i; if(key != oldkey && NR > 1) {if(n[12] > 0 && m[12] != "NA"){out = m[12]/n[12]}else{out = "NA"}for(j = 13; j <= (NF - 1); j++) {if(n[j] > 0 && m[j] != "NA") {out = out "\t" m[j]/n[j]}else{out = out "\t" "NA"}}; print oldkey, out; split("", n);  split("", m)}; for(j = 12; j <= (NF - 1); j++) {if($j != "NA" && $j != "." && $j != -1){m[j] += $j; n[j] += 1}}; oldkey=key}' > methgain.genomestudio.at/specific.cancer.$name.gain.bed
done

# count number of normal samples that present an average GpG probe methylation above .33 for each gene/c-HMR
for i in `find /home/user/bsdata/array/tissue/* -name "*genomestudio.bed.gz" -a ! -name "*cancer*"`
do
name=`basename $i genomestudio.bed.gz | sed 's/450k.betas.//1'`

tabix -B $i normal.specific.hmr.nearpromoters.bed |\
 bedtools intersect -wao -a normal.specific.hmr.nearpromoters.bed -b - |\
 awk -v thres=.33\
 'BEGIN{OFS="\t"; out=""}{key=$1; for(i=2; i<= 7; i ++) key = key "\t" $i; if(key != oldkey && NR > 1) {if(n[12] > 0 && m[12] != "NA"){out = m[12]/n[12]}else{out = "NA"}for(j = 13; j <= (NF - 1); j++) {if(n[j] > 0 && m[j] != "NA") {out = out "\t" m[j]/n[j]}else{out = out "\t" "NA"}}; print oldkey, out; split("", n);  split("", m)}; for(j = 12; j <= (NF - 1); j++) {if($j != "NA" && $j != "." && $j != -1){m[j] += $j; n[j] += 1}}; oldkey=key}' > methgain.genomestudio.at/specific.normal.$name.gain.bed

done


# HMRs near promoters
# create outdir
mkdir -p methgain.genomestudio.near

# count number of cancer samples that present an average GpG probe methylation above .33 for each gene/c-HMR
for i in /home/user/bsdata/array/tissue/450k.betas.cancer*genomestudio.bed.gz
do
name=`basename $i .genomestudio.bed.gz | sed 's/450k.betas.cancer.//1'`
tabix -B $i normal.specific.hmr.nearpromoters.bed |\
 bedtools intersect -wao -a normal.specific.hmr.nearpromoters.bed -b - |\
 awk\
 'BEGIN{OFS="\t"; out=""}{key=$1; for(i=2; i<= 8; i ++) key = key "\t" $i; if(key != oldkey && NR > 1) {if(n[13] > 0 && m[13] != "NA"){out = m[13]/n[13]}else{out = "NA"}for(j = 14; j <= (NF - 1); j++) {if(n[j] > 0 && m[j] != "NA") {out = out "\t" m[j]/n[j]}else{out = out "\t" "NA"}}; print oldkey, out; split("", n);  split("", m)}; for(j = 13; j <= (NF - 1); j++) {if($j != "NA" && $j != "." && $j != -1){m[j] += $j; n[j] += 1}}; oldkey=key}' > methgain.genomestudio.near/specific.cancer.$name.gain.bed
done

# count number of normal samples that present an average GpG probe methylation above .33 for each gene/c-HMR
for i in `find /home/user/bsdata/array/tissue/* -name "*genomestudio.bed.gz" -a ! -name "*cancer*"`
do
name=`basename $i .genomestudio.bed.gz | sed 's/450k.betas.//1'`

tabix -B $i normal.specific.hmr.nearpromoters.bed |\
 bedtools intersect -wao -a normal.specific.hmr.nearpromoters.bed -b - |\
 awk\
 'BEGIN{OFS="\t"; out=""}{key=$1; for(i=2; i<= 8; i ++) key = key "\t" $i; if(key != oldkey && NR > 1) {if(n[13] > 0 && m[13] != "NA"){out = m[13]/n[13]}else{out = "NA"}for(j = 14; j <= (NF - 1); j++) {if(n[j] > 0 && m[j] != "NA") {out = out "\t" m[j]/n[j]}else{out = out "\t" "NA"}}; print oldkey, out; split("", n);  split("", m)}; for(j = 13; j <= (NF - 1); j++) {if($j != "NA" && $j != "." && $j != -1){m[j] += $j; n[j] += 1}}; oldkey=key}' > methgain.genomestudio.near/specific.normal.$name.gain.bed

done


###########################
####    Smiley plots
###########################


####    Check c-HMR shape with array data (smiley plots)

# make output directory
mkdir -p methgain.at/shape

# count number of cancer samples that present an average GpG probe methylation above .33 for each gene/c-HMR
for i in /home/user/bsdata/array/tissue/450k.betas.cancer*bed.gz
do
name=`basename $i .bed.gz | sed 's/450k.betas.cancer.//1'`
tabix -B $i normal.common.hmr.merged100b.atpromoters.bed |\
 bedtools intersect -wao -a normal.common.hmr.merged100b.atpromoters.bed -b - |\
 awk\
 '{size=$3-$2; mid=($2 + $3)/2; dis=int(2*1000*($9-mid)/size + .5)/1000; for(j = 12; j <= (NF - 1); j++) {if($j != "NA" && $j != "." && $j != -1){m[dis"\t"j] += $j; n[dis"\t"j] += 1}}}END{OFS="\t"; for(i in n) print i, m[i], n[i]}' >\
 methgain.at/shape/common.cancer.$name.shape
done

# count number of normal samples that present an average GpG probe methylation above .33 for each gene/c-HMR
for i in `find /home/user/bsdata/array/tissue/* -name "*gz" -a ! -name "*cancer*"`
do
name=`basename $i .bed.gz | sed 's/450k.betas.//1'`

tabix -B $i normal.common.hmr.merged100b.atpromoters.bed |\
 bedtools intersect -wao -a normal.common.hmr.merged100b.atpromoters.bed -b - |\
 awk\
 '{size=$3-$2; mid=($2 + $3)/2; dis=int(2*1000*($9-mid)/size + .5)/1000; for(j = 12; j <= (NF - 1); j++) {if($j != "NA" && $j != "." && $j != -1){m[dis"\t"j] += $j; n[dis"\t"j] += 1}}}END{OFS="\t"; for(i in n) print i, m[i], n[i]}' >\
 methgain.at/shape/common.normal.$name.shape

done


# make output directory
mkdir -p methgain.near/shape

# count number of cancer samples that present an average GpG probe methylation above .33 for each gene/c-HMR
for i in /home/user/bsdata/array/tissue/450k.betas.cancer*bed.gz
do
name=`basename $i .bed.gz | sed 's/450k.betas.cancer.//1'`
tabix -B $i normal.common.hmr.merged100b.nearpromoters.bed |\
 bedtools intersect -wao -a normal.common.hmr.merged100b.nearpromoters.bed -b - |\
 awk\
 '{if(($8 >= -50000) && ($8 <= 50000)){size=$3-$2; mid=($2 + $3)/2; dis=int(2*1000*($10-mid)/size + .5)/1000; for(j = 13; j <= (NF - 1); j++) {if($j != "NA" && $j != "." && $j != -1){m[dis"\t"j] += $j; n[dis"\t"j] += 1}}}}END{OFS="\t";   for(i in n) print i, m[i], n[i]}' >\
 methgain.near/shape/common.cancer.$name.shape
done

# count number of normal samples that present an average GpG probe methylation above .33 for each gene/c-HMR
for i in `find /home/user/bsdata/array/tissue/* -name "*gz" -a ! -name "*cancer*"`
do
name=`basename $i .bed.gz | sed 's/450k.betas.//1'`

tabix -B $i normal.common.hmr.merged100b.nearpromoters.bed |\
 bedtools intersect -wao -a normal.common.hmr.merged100b.nearpromoters.bed -b - |\
 awk\
 '{if(($8 >= -50000) && ($8 <= 50000)){size=$3-$2; mid=($2 + $3)/2; dis=int(2*1000*($10-mid)/size + .5)/1000; for(j = 13; j <= (NF - 1); j++) {if($j != "NA" && $j != "." && $j != -1){m[dis"\t"j] += $j; n[dis"\t"j] += 1}}}}END{OFS="\t";   for(i in n) print i, m[i], n[i]}' >\
 methgain.near/shape/common.normal.$name.shape

done

###########################
####    Repetitive regions
###########################

# set repetitive regions
rmsk=/home/user/resources/human/ucsc.hg19.rmsk.familysure.nochr.bed.gz

# set output directory
outdir=/home/user/bsdata/hmr/rmsk
mkdir -p $outdir

# for each sample
for i in /home/user/bsdata/CpG_context_*smooth.bed.gz
do
# get sample name
name=`basename $i .smooth.bed.gz | sed 's/CpG_context_//1'`

# for each autosome
for chr in `seq 1 22`
do
# compute total U and M reads per rmsk class and chromosome
echo "tabix -B $i <(tabix $rmsk $chr) | bedtools intersect -wo -a <(tabix $rmsk $chr) -b - | awk '{key=\$1\"\t\"\$6; u[key] += \$12; m[key] += \$13; s[key] += \$14; n[key] += 1}END{OFS=\"\t\"; for(i in n) print i,u[i],m[i],s[i], n[i]}'"
done | parallel --gnu > $outdir/$name.tab
done


#######################################
####    TCGA expression and methylation
#######################################

####    Prepare data

##  methylation arrays

# create directory
mkdir -p /home/user/bsdata/tcga_meth

# arrange methylation arrays
scripts/arrange.tcga.metharrays.r

# compress and index
for i in /home/user/bsdata/tcga_meth/*bed
do
	echo "bgzip -c $i > $i.gz; tabix -p bed $i.gz"
done | parallel --gnu

##  RNAseq expression

# create directory
mkdir -p /home/user/bsdata/tcga_expr

# arrange expression values
while read line
do
# get info
lineARR=($line)
tissue=${lineARR[0]}
file=${lineARR[1]}
sample=${lineARR[2]}

# initialize file if needed
if [[ ! -a /home/user/bsdata/tcga_expr/"$tissue".rsem ]]
then
cut -f 1 /home/user/TCGA/$tissue/RNASeqV2/*/*/$file >\
 /home/user/bsdata/tcga_expr/$tissue.rsem
fi

# add a column with the sample name (passing by a tmp file)
tmp=`mktemp`

paste /home/user/bsdata/tcga_expr/$tissue.rsem\
 <(cut -f 2 /home/user/TCGA/$tissue/*/*/*/$file | tail -n +2 | sed '1 i '"$sample"'') >\
 $tmp

mv $tmp /home/user/bsdata/tcga_expr/$tissue.rsem

done < <(cat /home/user/TCGA/all.samples.withIDs.txt /home/user/TCGA/new.samples.withIDs.txt)

##  miRNA expression

# create directory
mkdir -p /home/user/bsdata/tcga_mirnaexpr

# arrange expression values
awk '$1=="miRNASeq"{OFS="\t"; split(FILENAME, a, "/"); print a[6], $7, $6, $3}'\
 /home/user/TCGA/BRCA/miRNA/file_manifest.txt\
 /home/user/TCGA/LUAD/miRNA/file_manifest.txt |\
 grep mirna > /home/user/bsdata/tcga_mirnaexpr/samples.tab

while read line
do
# get info
lineARR=($line)
tissue=${lineARR[0]}
file=${lineARR[1]}
sample=${lineARR[2]}
platform=${lineARR[3]}

# initialize file if needed
if [[ ! -a /home/user/bsdata/tcga_mirnaexpr/"$tissue".rsem ]]
then
cut -f 1 /home/user/TCGA/$tissue/miRNA/*/*/*/$file >\
 /home/user/bsdata/tcga_mirnaexpr/$tissue.rsem
fi

# add a column with the sample name (passing by a tmp file)
tmp=`mktemp`

paste /home/user/bsdata/tcga_mirnaexpr/$tissue.rsem\
 <(cut -f 3 /home/user/TCGA/$tissue/miRNA/*/*/*/$file | tail -n +2 | sed '1 i '"$sample"'') >\
 $tmp

mv $tmp /home/user/bsdata/tcga_mirnaexpr/$tissue.rsem

done < /home/user/bsdata/tcga_mirnaexpr/samples.tab

## Cross expression with methylation

# create output directory
mkdir -p /home/user/bsdata/tcga_res

# gather info
# gene promoters
for i in ARID1A PIK3R1 FANCC MSH2 XRCC3 TGFB2
do
scripts/get.gene.methylationANDexpression.fromTCGA.sh $i\
 /home/user/bsdata/tcga_res/promoter\
 /home/user/resources/human/gencode.v16.alltranscripts.coding.prom.nochr.bed\
 /home/user/bsdata/tcga_expr
done

# common HMRs at promoters
for i in ARID1A PIK3R1 FANCC MSH2 XRCC3 TGFB2
do
scripts/get.gene.methylationANDexpression.fromTCGA.sh $i\
 /home/user/bsdata/tcga_res/commonHMRat\
 /home/user/bsdata/hmr/normal.common.hmr.merged100b.atpromoters.bed\
 /home/user/bsdata/tcga_expr
done

# common HMRs near promoters
for i in ARID1A PIK3R1 FANCC MSH2 XRCC3 TGFB2
do
scripts/get.gene.methylationANDexpression.fromTCGA.sh $i\
 /home/user/bsdata/tcga_res/commonHMRnear\
 /home/user/bsdata/hmr/normal.common.hmr.merged100b.nearpromoters.bed\
 /home/user/bsdata/tcga_expr
done

# tissue specific HMRs at promoters
for i in ARID1A PIK3R1 FANCC MSH2 XRCC3 TGFB2
do
scripts/get.gene.methylationANDexpression.fromTCGA.sh $i\
 /home/user/bsdata/tcga_res/specificHMRat\
 /home/user/bsdata/hmr/normal.specific.hmr.atpromoters.bed\
 /home/user/bsdata/tcga_expr
done

# tissue specific HMRs near promoters
for i in ARID1A PIK3R1 FANCC MSH2 XRCC3 TGFB2
do
scripts/get.gene.methylationANDexpression.fromTCGA.sh $i\
 /home/user/bsdata/tcga_res/specificHMRnear\
 /home/user/bsdata/hmr/normal.specific.hmr.nearpromoters.bed\
 /home/user/bsdata/tcga_expr
done

# now some custom regions & genes
scripts/get.gene.methylationANDexpression.fromTCGA.sh TGFB2\
 /home/user/bsdata/tcga_res/customsuperHMR\
 <(echo -e "1\t218336382\t218349940\tbla\tsuperHMR\tTGFB2")\
 /home/user/bsdata/tcga_expr

# the same but one plot per CG
scripts/get.gene.methylationANDexpression.fromTCGA.sh TGFB2\
 /home/user/bsdata/tcga_res/customsuperHMRpercg\
 <(tabix /home/user/resources/human/450k.manifest.1to22andX.onlypos.bed.gz 1:218336382-218349940 | awk '{OFS="\t"; print $1, $2, $3, "bla", $4, "TGFB2"}')\
 /home/user/bsdata/tcga_expr

# the same but enlarging the region 100 % per side
scripts/get.gene.methylationANDexpression.fromTCGA.sh TGFB2\
 /home/user/bsdata/tcga_res/customsuperHMRenlarged\
 <(echo -e "1\t218336382\t218349940\tbla\tsuperHMR\tTGFB2" | awk '{OFS="\t"; size=$3 - $2; $2=$2 - size; $3=$3 + size; print}')\
 /home/user/bsdata/tcga_expr

# custom super enhancers miRNA
for i in hsa-let-7a-3 hsa-let-7b
do
scripts/get.gene.methylationANDexpression.fromTCGA.sh $i\
 /home/user/bsdata/tcga_res/customsuperENHANCER\
 <(echo -e "22\t46467114\t46488203\tbla\tsuperENHANCER\t$i")\
 /home/user/bsdata/tcga_mirnaexpr
done


# create table with tissue correspondencies
echo BRCA `grep -i breast /home/user/bsdata/samples.info.txt | grep -i normal | cut -d " " -f 1` >>\
 /home/user/bsdata/bsseq.tcga.correspondence.txt
echo COAD `grep -i colon /home/user/bsdata/samples.info.txt | grep -i normal | cut -d " " -f 1` >>\
 /home/user/bsdata/bsseq.tcga.correspondence.txt
echo LUAD `grep -i lung /home/user/bsdata/samples.info.txt | grep -i normal | cut -d " " -f 1` >>\
 /home/user/bsdata/bsseq.tcga.correspondence.txt
echo LIHC `grep -i liver /home/user/bsdata/samples.info.txt | grep -i normal | cut -d " " -f 1` >>\
 /home/user/bsdata/bsseq.tcga.correspondence.txt


# take corresponding common HMRs per tissue
for i in ARID1A PIK3R1 FANCC MSH2 XRCC3 TGFB2
do
scripts/get.gene.methylationANDexpression.fromTCGA.commonHMRpertissue.sh $i\
 /home/user/bsdata/tcga_res/commonHMRatcorrespondingtissue\
 /home/user/bsdata/bsseq.tcga.correspondence.txt\
 /home/user/resources/human/gencode.v16.alltranscripts.coding.prom.nochr.bed
done

#################################################
####    Methylation gain of c-HMR (array samples)
#################################################

# set region file and out dir
region=/home/user/bsdata/hmr/normal.common.hmr.merged100b.ID.bed
outdir=/home/user/bsdata/hmr/methgain

# create outdir
mkdir -p $outdir

# count number of cancer samples that present an average GpG probe methylation above .33 for each gene/c-HMR
for i in /home/user/bsdata/array/tissue/450k.betas.cancer*bed.gz
do

name=`basename $i .bed.gz | sed 's/450k.betas.cancer.//1'`

echo "tabix -B $i $region | bedtools intersect -wao -a $region -b - | awk -v thres=.33 'BEGIN{OFS=\"\t\"; out=0; N=0}{key=\$1; for(i=2; i<= 4; i ++) key = key \"\t\" \$i; if(key != oldkey) {for(j = 9; j <= NF; j++) {if(n[j] > 0) {         if(m[j]/n[j] > thres){ out += 1}; N += 1}}; if(NR > 1) {print oldkey, out, N}; out = 0; N = 0; split(\"\", n); split(\"\", m)}; for(j = 9; j <= NF; j++) {if(\$j != \"NA\"){m[j] += \$j; n[j] += 1}}; oldkey=key}END{print oldkey, out, N}' >  $outdir/common.cancer.$name.gain.bed"

done | parallel --gnu 


# count number of normal samples that present an average GpG probe methylation above .33 for each gene/c-HMR
for i in `find /home/user/bsdata/array/tissue/* -name "*gz" -a ! -name "*cancer*"`
do
name=`basename $i .bed.gz | sed 's/450k.betas.//1'`

echo "tabix -B $i $region | bedtools intersect -wao -a $region -b - | awk -v thres=.33 'BEGIN{OFS=\"\t\"; out=0; N=0}{key=\$1; for(i=2; i<= 4; i ++) key = key \"\t\" \$i; if(key != oldkey) {for(j = 9; j <= NF; j++) {if(n[j] > 0) {         if(m[j]/n[j] > thres){ out += 1}; N += 1}}; if(NR > 1) {print oldkey, out, N}; out = 0; N = 0; split(\"\", n); split(\"\", m)}; for(j = 9; j <= NF; j++) {if(\$j != \"NA\"){m[j] += \$j; n[j] += 1}}; oldkey=key}END{print oldkey, out, N}' > $outdir/common.normal.$name.gain.bed"

done | parallel --gnu 
