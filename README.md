# MiMB_JACUSA2_chapter
JACUSA2 - A framework for the mapping of RNA modifications

We provide a pipeline for rRNA modification analysis of Nanopore sequencing reads based on multiple conditions and replicates.

# General description
- JACUSA2 can detect RNA modification rapidly.
- JACUSA2 detect modification based on mapping characteristics, it requires BAM files as input.
- JACUSA2 can be run in single mode (call-1) or paired mode (call-2). For sensitivity, we adopt paired mode.
- JACUSA2 output can be used to predict RNA modifications.
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

- Align reads to the transcriptome using the following Minimap2 command. To reduce the indexing time of the human genome, save the index with the option ”-d” before the mapping and use the index instead of the reference file in the minimap2 command line
```
minimap2 -d reference.mmi reference.fa
minimap2 -t 5 --MD -ax splice --junc-bonus 1 -k14 --secondary=no --junc-bed final_annotation_96.bed -ub reference.mmi Reads.fastq.gz |samtools view -bS > mapping.bam
```
The inputs are the FASTQ files "Reads.fastq.gz" and the reference sequence "reference.fasta". The output is a SAM file that is converted to BAM file "mapping.bam" using samtools command.  The setting is used for spliced alignments of the dierct RNA sequencing :
  - "-ax splice --junc-bed annotation.bed --junc-bonus INT" where annotation.bed is the file contaning  the splice junctions and INT is the bonus score. The BED file can be generated using the following command: 
 ```
 paftools.js gff2bed annotation.gtf > annotation.bed
  ```
  - ”-ub” to align reads to both strands of the reference, 
  - small k-mer size ”-k [=14]” to enhance sensitivity.

We recommend outputting primary alignments ”–secondary=no”. ’–MD’ parameter is used to add the reference sequence information to the alignment; this is recommended for the downstream analysis. Check Minimap2 [manual](https://github.com/lh3/minimap2) for further details.

## Detect RNA modification
We provide a snakemake pipeline for JACUSA2 variant calling using call2 method and downstream analysis for the detection of modification patterns and predict modified sites. The pipeline is composed of many rules and require setting diffrent parameters.

Be aware to set all parameters before running the pipeline. 
describe all parameters ...............................................................................
The inputs are BAM files "HEK293T-WT-rep2.bam" and "HEK293T-WT-rep3.bam	" for wild-type replicates and "HEK293T-KO-rep2.bam" and "HEK293T-KO-rep3.bam" for the Knock-out replicates. Make sure to set "-P1" and "-P2" according to the corresponding library. The mismach score is produce by default but you need to add "-D" and "-I" to output the deletion and insertion scores. Filter reads according with base calling quality "-q [$>1$]", mapping quality "-m [$>1$] and read coverage "-c [$>4$]". Plus, consider the filter feature "-a [=Y]" to exclude sites within homopolymer regions to improve sensitivity. The output -r "WT_vs_KO_call2_result.out" consists of a read error profile where the format is a combination of BED6 with JACUSA2 call-2 specific columns and common info columns: info, filter, and ref. The number of threads can be customized via the parameter "-p". Check JACUSA2 [manual](https://github.com/dieterich-lab/JACUSA2) for more details.
.....................................................

- Run JACUSA2 call-2 rule to call variants using paired conditions. Make sure that you set the path to the jar file.
```  
srun snakemake --cores all jacusa2_call2
```
- Run get_features rule to preprocess JACUSA2 call2 output and extract features.
```
$ srun snakemake --cores all get_features
```
- Run get_pattern rule to lear patterns representing m6A modification
```
$ srun snakemake --cores all get_pattern
```
- Run predict_modification rule to predict modified sites
```
$ srun snakemake --cores all predict_modification
```  
Note that the rules are linked so that the workflow are determined from top (e.g. predict modification) to bottom (e.g. sort bam) and
executed accordingly from bottom to top 4. Therefore, running ”predict modification” rule leads to excuting all rules in its pipeline.

# Dependencies and versions
Software | Version 
--- | ---
Minimap2 | 2.22
samtools | 1.12
openjdk | 11.0.13
R | 4.0.5
PERL | 5.28.1
bedtools | 2.29.2
snakemake | 6.8.1

R Package | Version
--- | ---
ggplot2 | 3.3.5
NMF | 0.23.0
pROC | 1.18.0
