#!/bin/bash
#Integrate results
#working directory is /prj/MiMB_book_chapter_Amina_Isabel/Nanopore/HEK293/JACUSA2

WT_vs_KO=$1 
# WT_vs_IVT=$2 
# KO_vs_IVT=$3 
genome=$2
GRCh38_96=$3
path=$4
echo "Preprocessing JACUSA2 output ..."

mkdir -p $path
sort -k5,5gr  $WT_vs_KO |awk '($13=="A" && $12!~/Y/)' > $path/WT_vs_KO_call2_result.txt
# sort -k5,5gr  $WT_vs_IVT |awk '($13=="A" && $12!~/Y/)' >  $path/WT_vs_IVT_call2_result.txt
# sort -k5,5gr  $KO_vs_IVT |awk '($13=="A" && $12!~/Y/)' >  $path/KO_vs_IVT_call2_result.txt

# echo $path 
cut -f 1-6 $path/WT_vs_KO_call2_result.txt | awk '{print $1"\t"$2"\t"$3"\tCand\t0\t"$6;}' |sort -u > $path/allSites.bed
# # echo $path
grep -v pcDNA5 $path/allSites.bed > $path/allSites2.bed
# echo $genome
bedtools slop -i $path/allSites2.bed -g $genome -b 2 > $path/allSitesExt2.bed
# echo $path

# #call2
INP="$WT_vs_KO"
bedtools intersect -filenames -loj -a $path/allSitesExt2.bed -b ${INP} > $path/call2_SitesExt2.bed

#reformat
cat $path/call2_SitesExt2.bed |perl HEK293_data_prep.pl > $path/data_reformat.txt

#5mer
# srun bedtools getfasta -fi $GRCh38_96 -bed $path/allSitesExt2.bed -name -s -tab > $path/checkMotif.txt

#sed -e 's/Cand:://' checkMotif.txt |sed -e 's/\-/_/' |awk '{gsub("[\(\)]","");print $1,$2;}' > checkMotif_reformat.txt

# cat $path/checkMotif.txt | perl -e 'while(<>){s/Cand:://; 's/\-/_/'; 's/[\(\)]//g'; print;}' >  $path/checkMotif_reformat.txt

#*************************************BEGIN. Modified part *********************************************
bedtools getfasta -fi $GRCh38_96 -bed $path/allSitesExt2.bed -s -tab > $path/checkMotif_reformat.txt
#*************************************END. Modified part *********************************************

#ok
##
