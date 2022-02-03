#!/bin/bash
jacusa_out=$1 
genome=$2
GRCh38_96=$3
path=$4
sites=$5

echo "Preprocessing JACUSA2 output ..."

mkdir -p $path

if [[ $sites = "None" ]]
then
    sort -k5,5gr  $jacusa_out |awk '($13=="A" && $12!~/Y/)' > $path/jacusa_out.txt
    cut -f 1-6 $path/jacusa_out.txt | awk '{print $1"\t"$2"\t"$3"\tCand\t0\t"$6;}' |sort -u > $path/allSites.bed
    grep -v pcDNA5 $path/allSites.bed > $path/allSites2.bed
    bedtools slop -i $path/allSites2.bed -g $genome -b 2 > $path/allSitesExt2.bed 
    bedtools intersect -filenames -wa -wb -a $path/allSitesExt2.bed -b $jacusa_out > $path/call2_SitesExt2.bed
    bedtools getfasta -fi $GRCh38_96 -bed $path/allSitesExt2.bed -s -tab > $path/checkMotif_reformat.txt
else
    bedtools intersect -filenames -wa -wb -a $sites -b $jacusa_out > $path/call2_SitesExt2.bed
    bedtools getfasta -fi $GRCh38_96 -bed $sites -s -tab > $path/checkMotif_reformat.txt
fi
#reformat
cat $path/call2_SitesExt2.bed |perl scripts/get_features.pl > $path/data_reformat.txt
