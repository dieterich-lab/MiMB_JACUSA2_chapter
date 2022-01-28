# MiMB_JACUSA2_chapter
JACUSA2 - A framework for the mapping of RNA modifications

We provide a pipeline for RNA modification detection of direct nanopore sequencing based on multiple conditions and replicates.

# General description
- JACUSA2 can detect RNA modification rapidly.
- JACUSA2 detects modification based on mapping characteristics, it requires BAM files as input.
- JACUSA2 can be run in single-mode (call-1) or paired mode (call-2). For sensitivity, we adopt paired mode.
- JACUSA2 output can be used to predict RNA modifications.

# Installation
We recommend installing software dependencies via `Conda` on Linux. You can find Miniconda installation instructions for Linux [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html).
Make sure you install the [Miniconda Python3 distribution](https://docs.conda.io/en/latest/miniconda.html#linux-installers).
```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```
For performance and compatibility reasons you should install `Mamba` via conda to install Snakemake. See [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) for more details.
```
conda install -c conda-forge mamba
```
Once you have installed Conda and Mamba, you can download the Snakemake pipeline and the example datasets.
```
git clone https://github.com/dieterich-lab/MiMB_JACUSA2_chapter.git
cd Nanopore_HEK293
```

Then, you install the required packages after creating an environment with Snakemake installed `environment.yml` file.
```
mamba env create -f environment.yaml
conda activate JACUSA2pipeline_env
```
Before executing the Snakemake workflow, download JACUSA2 [jar](https://github.com/dieterich-lab/JACUSA2) file and make sure that you set the path `jar` in the config file.

# Usage
The following pipeline is used to predict m6A modification from nanopore RNA direct sequencing data. The benchmark obtained from [PRJEB40872](https://www.ebi.ac.uk/ena/browser/view/PRJEB40872?show=reads) is composed of two samples from two conditions: wild-type cells (modified RNAs) and Mettl3 knockout cells (unmodified RNAs) with two replicates (2 and 3). The analysis is validated against reported m6A sites in the three miCLIP-based studies Bouliaset al. [2019], Koh et al. [2019], Körtel et al. [2021]. [We limited the analysis to the set of sites in the 'data/regions.bed' file]

## Preprocessing of Direct RNAseq
1. Base call the ionic current signal stored in FAST5 file using Guppy basecaller.
```
guppy_basecaller --compress_fastq -i path_to_fast5 -s path_to_output -c config_file.cfg --cpu_threads_per_caller 14 --num_callers 1
```

The inputs are the path to FAST5 files, the output folder, and the config file or the flowcell/kit combination. The output is FASTQ files that can
be compressed using the option ”–compress fastq”. Set the number of threads ”cpu threads per caller” and the number of parallel basecallers ”num caller” according to your resources. Additional details can be found [here](https://github.com/metagenomics/denbi-nanopore-training/blob/master/docs/basecalling/basecalling.rst).

2. Align reads to the transcriptome using the following Minimap2 command. To reduce the indexing time of the human genome, save the index with the option ”-d” before the mapping and use the index instead of the reference file in the minimap2 command line
```
minimap2 -d reference.mmi reference.fa
minimap2 -t 5 --MD -ax splice --junc-bonus 1 -k14 --secondary=no --junc-bed final_annotation_96.bed -ub reference.mmi Reads.fastq.gz |samtools view -bS > mapping.bam
```
The inputs are the FASTQ files "Reads.fastq.gz" and the reference sequence "reference.fasta". The output is a SAM file that is converted to BAM file "mapping.bam" using samtools command.  The setting is used for spliced alignments of the direct RNA sequencing :
  * "-ax splice --junc-bed annotation.bed --junc-bonus INT" where annotation.bed is the file containing the splice junctions and INT is the bonus score. The BED file can be generated using the following command: 
 ```
 paftools.js gff2bed annotation.gtf > annotation.bed
  ```
  * ”-ub” to align reads to both strands of the reference, 
  * small k-mer size ”-k [=14]” to enhance sensitivity.

We recommend outputting primary alignments ”–secondary=no”. ’–MD’ parameter is used to add the reference sequence information to the alignment; this is recommended for the downstream analysis. Check Minimap2 [manual](https://github.com/lh3/minimap2) for further details.

## Detect RNA modification
We provide a snakemake pipeline for JACUSA2 variant calling using call2 method and downstream analysis for the detection of modification patterns and predict modified sites. The pipeline is composed of many rules and requires setting different parameters.

- Be aware to set all parameters before running the pipeline. 

      * General
            * label: 'HEK293_WT_KO' label of the analysis
            * jar : 'JACUSA_v2.0.2-RC.jar'  #path to JACUISA2 JAR file 
            * path_out: './output' # path to the output directory, if it doesn't exist it will be created 
            * path_inp: './data' # path to the directory containing inputs - all input files are relative to this directory
            * reference : 'GRCh38_96.fa' # path to reference squence 
            * modified_sites: 'miCLIP_union.bed' #BED6 file containing known modified sites where 'name' refers to the annotation of the position. useful for learning patterns (training and test set).
            * chr_size: "hg38.genome"  #file contaning size of chromosomes (Chromosome     | size )
            * regions: "regions.bed" # BED6 file contaning set of 5-mer (NNANN) to analyze, if ="", all 5-mers (NNANN) will be considered.
            * data: a dictionary of two keys (cond1, cond2) referring to the paired conditions inputs. The value is the list of replicates names without ".bam" extension.
              * cond1: ["HEK293T-WT-rep2","HEK293T-WT-rep3"]
              * cond2: ["HEK293T-KO-rep2","HEK293T-KO-rep3"]
      * jacusa_params: a dictionary where keys refer to parameters (e.g. p: 16 to set the number of threads to 16). Please use "" if no value is affected to the parameter. We use the following parameters:
            * P1: 'FR-SECONDSTRAND'  # Mandatory parameters referring to the library of the first condition sample.
            * P2: 'FR-SECONDSTRAND'  # Mandatory parameters referring to the library of the second condition sample.
            * m: 1  # filter reads by mapping quality
            * q: 1  # filter reads by base calling quality
            * c: 4  # filter reads by coverage
            * a: 'Y'  # Mandatory parameters to filter sites within the holy-polymer regions.
            * p: 16    # parameter to customize the number of threads
            * D: ''  # Mandatory parameter to output deletion score.
            * I: ''  # Mandatory parameter to output insertion score.
      * pattern_params:       # specify patterns and their combinations to be used, please use "" if no value is affected to the field.
            * internal_pattern: "Boulias,Koertel,Koh" # specify the annotation of the set of modified sites to be used as a training set. in case you use an external pattern put "". 
            * external_pattern: ""  # path to an external pattern in case you don't use internal_pattern, else put ""
            * combined_patterns: #patterns to combine, add as many combinations as you want as a [key(any name): value (pattern number)] combination.
                      pt1: [1,2,4,6]  
                      pt2: [1,2,3,4,6]


Please check JACUSA2 [manual](https://github.com/dieterich-lab/JACUSA2) for more details on how to use parameters.

- Run JACUSA2_call2 rule to call variants using paired conditions.
```  
srun snakemake --cores all jacusa2_call2
```
The output is a file called "Cond1vsCond2Call2.out" under "./output/jacusa/label[HEK293_WT_KO]/" and filtered bam file under "./output/bam/label[HEK293_WT_KO]/".
- Run get_features rule to preprocess JACUSA2 call2 output and extract features.
```
$ srun snakemake --cores all get_features
```
The output is an R object "features.rds" under "./output/analysis/label[HEK293_WT_KO]/features/".
- Run get_pattern rule to learn patterns representing m6A modification.
```
$ srun snakemake --cores all get_pattern
```
The output is an R object "NMF.rds" containing the factorization result, including basis and coefficient matrices, plus, plots showing the rank selection result. The output is under "./output/analysis/label[HEK293_WT_KO]/pattern/". Implicitly, training and test set files (resp. train_features.rds, test_features.rds"  under are created and, subsequently, used for the learning model.

For the testing example, the prediction.csv is supposed to contain 1905 sites.

- Run visualize_pattern rule to predict modified sites
```
$ srun snakemake --cores all visualize_pattern
```  
The output is a set of figures representing barplots for the produced patterns, in addition to the pattern scoring barplots and heatmap of NMF resulting matrices, The output can be found under "./output/analysis/label[HEK293_WT_KO]/pattern/viz/".

For the testing example, the scoring pattern will look like the following barplots.

<p align="center">
  <img src="https://github.com/dieterich-lab/MiMB_JACUSA2_chapter/blob/amina/img/pattern_scores.png" width="300">
</p>

The combination of patterns representing more than 80% will look like this:

<p align="center">
  <img src="https://github.com/dieterich-lab/MiMB_JACUSA2_chapter/blob/amina/img/barplot_NMF.png" width="300">
</p>

- Run predict_modification rule to predict modified sites
```
$ srun snakemake --cores all predict_modification
```  
The output is a BED6 file(s) contaning score of the selected pattern(s) for the test set under "./output/analysis/label[HEK293_WT_KO]/prediction/". andthe corresponding eCDF (empirical cumulative distribution) and PPV (positive predictive values) plots.

For the testing example, the eCDF will look like the following figure: 

<p align="center">
  <img src="https://github.com/dieterich-lab/MiMB_JACUSA2_chapter/blob/amina/img/Pattern_ecdf.png" width="500">
</p>

Note that rules are linked so that the workflow is determined from top (e.g. predict modification) to bottom (e.g. sort bam) and
executed accordingly from bottom to top. Therefore, running ”predict_modification” rule leads to executing all rules in its pipeline.


<p align="center">
  <img src="https://github.com/dieterich-lab/MiMB_JACUSA2_chapter/blob/amina/img/snakemake.png" width="500">
</p>

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
stringr | 1.4.0
