# MiMB_JACUSA2_chapter
JACUSA2 - A framework for the mapping of RNA modifications

We provide a pipeline for rRNA modification analysis of Nanopore sequencing reads based on multiple conditions and replicates.

# General description
- JACUSA2 can detect RNA modification rapidly.
- JACUSA2 detect modification based on mapping characteristics, it requires BAM files as input.
- JACUSA2 can be run in single mode (call-1) or paired mode (call-2). For sensitivity, we adopt paired mode.
- JACUSA2 output can be used to predict RNA modification.
# Download
Download JACUSA2 JACUSA2 [jar](https://github.com/dieterich-lab/JACUSA2) file.

Download the analysis scripts to your cluster by the following command.

```
git clone https://github.com/dieterich-lab/MiMB_JACUSA2_chapter.git
```
# Usage
The following pipeline is used to predict m6A modification from nanopore RNA direct sequencing data. The benchmark obtained from [PRJEB40872](https://www.ebi.ac.uk/ena/browser/view/PRJEB40872?show=reads) is composed of two samples from two conditions: wild-type cells (modified RNAs) and Mettl3 knockout cells (unmodified RNAs) with two replicates (2 and 3). The analysis is validated against reported m6A sites in the three miCLIP-based studies Bouliaset al. [2019], Koh et al. [2019], Körtel et al. [2021].

## Preprocessing of Direct RNAseq
- Base call the ionic current signal stored in FAST5 file using Guppy basecaller.
```
guppy_basecaller --compress_fastq -i path_to_fast5 -s path_to_output -c config_file.cfg --cpu_threads_per_caller 14 --num_callers 1
```

The inputs are: the path to FAST5 files, the output folder, and the config file or the flowcell/kit combination. The output is FASTQ files that can
be compressed using the option ”–compress fastq”. Set the number of threads ”cpu threads per caller” and the number of parallel basecallers ”num caller” according to your resources. Additional details can be found [here](https://github.com/metagenomics/denbi-nanopore-training/blob/master/docs/basecalling/basecalling.rst).

- Align reads to the transcriptome using the following Minimap2 command. 
```
minimap2 -t 5 -ax map-ont -uf -k14 -p 1.0 -N 100 reference.fasta Reads.fastq |samtools view -bS > mapping.bam	
```
The inputs are the FASTQ files "Reads.fastq" and the reference sequence "reference.fasta". The output is a SAM file that is converted to BAM file "mapping.bam" using samtools command.  Use the default setting for Direct RNA-seq ”-ax map-ont”, ”-uf” to force the alignment to the forward strand of the reference, and a small k-mer size ”-k [=14]” to enhance sensitivity. Also, set  ”-N [=100]” to output primary alignments and up to "100” top secondary alignments if the ratio of their chaining scores compared to the corresponding primary alignments is equal to ”-p [=1]”. Check Minimap2 [manual](https://github.com/lh3/minimap2) for further details.

Sorte BAM file using the following command:
```
samtools sorte mapping.bam mapping.sorted.bam  
```

Create an index of the BAM file
```
samtools index mapping.sorted.bam
```  
## Detect RNA modification
Run JACUSA2 call-2. Make sure that you set the path to the jar file.
```  
java -jar JACUSA2.jar call-2 -q 1 -m 1 -c 4 -p 10 -D -I -a Y -P1 FR-SECONDSTRAND -P2 FR-SECONDSTRAND -r WT_vs_KO_call2_result.out HEK293T-WT-rep2.bam,HEK293T-WT-rep3.bam	HEK293T-KO-rep2.bam, HEK293T-KO-rep3.bam
```
The inputs are BAM files "HEK293T-WT-rep2.bam" and "HEK293T-WT-rep3.bam	" for wild-type replicates and "HEK293T-KO-rep2.bam" and "HEK293T-KO-rep3.bam" for the Knock-out replicates. Make sure to set "-P1" and "-P2" according to the corresponding library. The mismach score is produce by default but you need to add "-D" and "-I" to output the deletion and insertion scores. Filter reads according with base calling quality "-q [$>1$]", mapping quality "-m [$>1$] and read coverage "-c [$>4$]". Plus, consider the filter feature "-a [=Y]" to exclude sites within homopolymer regions to improve sensitivity. The output -r "WT_vs_KO_call2_result.out" consists of a read error profile where the format is a combination of BED6 with JACUSA2 call-2 specific columns and common info columns: info, filter, and ref. The number of threads can be customized via the parameter "-p". Check JACUSA2 [manual](https://github.com/dieterich-lab/JACUSA2) for more details.
```
bash README_processing.sh WT_vs_KO_call2_result.out hg38.genome GRCh38_96.fa path_to_output.
```
```
Rscript HEK293_data_prep_step2.R path_to_output miCLIP_union.bed
```
```
Rscript HEK293_data_prep_step3.R path_to_output miCLIP_union.bed
```  
# Dependencies and versions
Software | Version 
--- | ---
Minimap2 | 2.22
samtools | 1.12
openjdk | 11.0.13
R | 4.0.5
PERL | 5.28.1
bedtools | 2.29.2

R Package | Version
--- | ---
ggplot2 | 3.3.5
NMF | 0.23.0
pROC | 1.18.0
