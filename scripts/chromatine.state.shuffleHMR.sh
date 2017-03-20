#!/bin/bash

if [ $# -lt 4 ]
then
    echo
    echo "  Takes a sample name, a chromatineHMM name and a tmp filep refix ..."
    echo "  ... and writes a cuple of files at \$outdir with overlap numbers"
	echo
	echo "usage:"
	echo
	echo "     chromatine.shtate.shuffleHMR.sh <sample_name> <chromatineHMM_names> <tmp_file_name> <out_dir>"
	echo
exit 1
fi

# get arguments
bsname=$1
chsname=$2
seed=$3
outdir=$4

# get seed sufix
seedsuf=`basename $seed | sed 's/tmp.//1'`

# set shuffled hmr file names
common=$seed.$bsname.common.bed
specific=$seed.$bsname.specific.bed
all=$seed.$bsname.all.bed

# get suffled HMR files
bedtools shuffle -i $bsname.common.hmr.bed -g /home/user/resources/human/hg19.chrsizes | bedtools sort -i - > $common
bedtools shuffle -i $bsname.hmr.specific.bed -g /home/user/resources/human/hg19.chrsizes | bedtools sort -i -  > $specific
bedtools shuffle -i CpG_context_$bsname.hmr.bed -g /home/user/resources/human/hg19.chrsizes | bedtools sort -i -  > $all

# per region
for i in `seq 1 15`
do
# set chromatine state file
chsfile=/home/user/resources/human/encode/wgEncodeBroadHmm${chsname}HMM.${i}_*.bed

# count number of regions
# get chromatine state regions with HMRs
Ncommon=`bedtools intersect -u -a $chsfile -b $common | wc -l`
Nspecific=`bedtools intersect -u -a $chsfile -b $specific | wc -l`
Nall=`bedtools intersect -u -a $chsfile -b $all | wc -l `

NN=`cat $chsfile | wc -l`

# split output
echo $i $Nall $Ncommon $Nspecific $NN
done > $outdir/$bsname.chrState.with.hmr.$seedsuf.txt

# per hmr
# set initial tmp files
echo -e "1\t0\t1" > $seed.$bsname.common.tmp
echo -e "1\t0\t1" > $seed.$bsname.specific.tmp
echo -e "1\t0\t1" > $seed.$bsname.all.tmp

for i in `seq 1 15`
do
# set chromatine state file
chsfile=/home/user/resources/human/encode/wgEncodeBroadHmm${chsname}HMM.${i}_*.bed

# count number of HMRs (and set up tmp files)
Ncommon=`bedtools intersect -v -a $common -b $seed.$bsname.common.tmp | bedtools intersect -u -a - -b $chsfile | cut -f 1-3 | tee -a $seed.$bsname.common.tmp | wc -l`
Nspecific=`bedtools intersect -v -a $specific -b $seed.$bsname.specific.tmp | bedtools intersect -u -a - -b $chsfile | cut -f 1-3 | tee -a $seed.$bsname.specific.tmp | wc -l`
Nall=`bedtools intersect -v -a $all -b $seed.$bsname.all.tmp | bedtools intersect -u -a - -b $chsfile | cut -f 1-3 | tee -a $seed.$bsname.all.tmp | wc -l`

# split output
echo $i $Nall $Ncommon $Nspecific
done > $outdir/$bsname.hmr.with.chrState.$seedsuf.txt

# append total number of HMRs
echo 0 `cat $all |wc -l` `cat $common |wc -l` `cat $specific |wc -l` >> $outdir/$bsname.hmr.with.chrState.$seedsuf.txt

# remove tmp files
rm $seed.$bsname.common.tmp $seed.$bsname.all.tmp $seed.$bsname.specific.tmp
rm $common $specific $all $seed
