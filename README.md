# MiMB_JACUSA2_chapter
JACUSA2 - A framework for the mapping of RNA modifications
We provide a pipeline for rRNA modification analysis of Nanopore sequencing reads based on multiple conditions and replicates.

Pairwise condition analysis of rRNA modification is carried out using JACUSA2 tool with call-2 mode. The output is a BED file with many information including Mismatch, Insertion, and Deletion scores reflecting variant discrimination. This output file will be used for the downstream analysis of rRNA modification detection and its evaluation.

# Dependencies and versions
Software | Version 
--- | ---
Minimap2 | 2.22
samtools | 1.12
openjdk | 11.0.13
R | 4.0.5
PERL | 5.28.1
bedtools | 2.29.2
ggplot2 | 3.3.5
NMF | 0.23.0
pROC | 1.18.0
