# MiMB_JACUSA2_chapter
JACUSA2 - A framework for the mapping of RNA modifications

We provide a Snakemake pipeline for RNA modification detection from nanopore Direct RNA-seq using multiple conditions and replicates.

# General description
- The pipeline is based on JACUSA2 variant calling, where the output of JACUSA2 is analyzed to detect RNA modification patterns and predict modified sites.
- JACUSA2 detects modification based on mapping characteristics, it requires BAM files as input.
- JACUSA2 can be run in single-mode (call-1) or paired mode (call-2). For sensitivity, we adopt paired mode.
- All data used in our manuscript can be found at  [Zenodo](https://doi.org/10.5281/zenodo.5924995).

# Installation
We recommend installing software dependencies via `Conda` on Linux. You can find Miniconda installation instructions for Linux [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html).
Make sure you install the [Miniconda Python3 distribution](https://docs.conda.io/en/latest/miniconda.html#linux-installers).
```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```
In the Conda installer you will have to accept the license and select an installation directory.
When asked also select to run 'conda init'. After the installer completed you can open a shell
to get the basic Conda setup.

For performance and compatibility reasons you should install `Mamba` via conda to install Snakemake. See [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) for more details.
```
conda install -c conda-forge mamba
```
Once you have installed Conda and Mamba, you can download the Snakemake pipeline.
```
git clone https://github.com/dieterich-lab/MiMB_JACUSA2_chapter.git
cd MiMB_JACUSA2_chapter/Nanopore_HEK293
```

Then, install the required packages after creating an environment with the downloaded `environment.yml` file.
```
mamba env create -f environment.yaml
conda activate JACUSA2pipeline_env
```
Before executing the Snakemake workflow, we recommend downloading JACUSA2 [jar](https://github.com/dieterich-lab/JACUSA2) file and make sure that you set the path `jar` in the config file. In case the jar file is not set, "JACUSA_v2.0.1.jar" will be automatically downloaded and used.

# Usage
The following protocol describes how to predict m6A modification from nanopore direct RNA-seq data. The benchmark is obtained from [PRJEB40872](https://www.ebi.ac.uk/ena/browser/view/PRJEB40872?show=reads) and is composed of two samples from two conditions: wild-type cells (modified RNAs) and Mettl3 knockout cells (unmodified RNAs) with two replicates (2 and 3). The analysis is validated against reported m6A sites in the three miCLIP-based studies Bouliaset al. [2019], Koh et al. [2019], Körtel et al. [2021]. [Following the use case 1 of the manuscript, we limited the analysis to the set of sites reported in 'data/selected_regions.bed' file]

We present hereafter the main steps to produce BAM files which would be the input of our pipeline.

## Preprocessing of Direct RNAseq
1. Base call the ionic current signal stored in FAST5 file using Guppy basecaller.
```
guppy_basecaller --compress_fastq -i path_to_fast5 -s path_to_output -c config_file.cfg --cpu_threads_per_caller 14 --num_callers 1
```

The inputs are the path to FAST5 files, the output folder, and the config file or the flowcell/kit combination. The output is FASTQ files that can
be compressed using the option `–compress fastq`. Set the number of threads `cpu threads per caller` and the number of parallel basecallers `num caller` according to your resources. Additional details can be found [here](https://github.com/metagenomics/denbi-nanopore-training/blob/master/docs/basecalling/basecalling.rst).

2. Align reads to the transcriptome using the following Minimap2 command. To reduce the indexing time of the human genome, save the index with the option `-d` before the mapping and use the index instead of the reference file in the minimap2 command line
```
minimap2 -d reference.mmi reference.fa
minimap2 -t 5 --MD -ax splice --junc-bonus 1 -k14 --secondary=no --junc-bed final_annotation_96.bed -ub reference.mmi Reads.fastq.gz |samtools view -bS > mapping.bam
```
The inputs are the FASTQ files `Reads.fastq.gz` and the reference sequence `reference.fasta`. The output is a SAM file that is converted to BAM file `mapping.bam` using samtools command.  The setting is used for spliced alignments of the direct RNA sequencing :
  * `-ax splice --junc-bed annotation.bed --junc-bonus INT` where annotation.bed is the file containing the splice junctions and INT is the bonus score. The BED file can be generated using the following command: 
 ```
 paftools.js gff2bed annotation.gtf > annotation.bed
  ```
  * `-ub` to align reads to both strands of the reference, 
  * small k-mer size `-k [=14]` to enhance sensitivity.

We recommend outputting primary alignments `–secondary=no`. `–MD` parameter is used to add the reference sequence information to the alignment; this is recommended for the downstream analysis. Check Minimap2 [manual](https://github.com/lh3/minimap2) for further details.

## Detect RNA modification
For the testing example, we consider BAM files of wild-type cells ("HEK293T-WT-rep2.bam","HEK293T-WT-rep3.bam") and knockout cells ("HEK293T-KO-rep2.bam","HEK293T-KO-rep3.bam") from [Zenodo](https://doi.org/10.5281/zenodo.5924995) as our inputs. 

The pipeline is composed of many target rules (fig. 4) and requires setting different parameters.

- Be aware to set all parameters before running the pipeline. Please put all required data in `data` folder.

      label: 'HEK293_WT_KO' #label of the analysis
      jar : 'JACUSA_v2.0.2-RC.jar'  #path to JACUISA2 JAR file 
      path_out: 'output' # path to the output directory, if it doesn't exist it will be created 
      path_inp: 'data' # path to the directory containing inputs - all input files are relative to this directory
      reference : 'data/GRCh38_96.fa' # path to reference squence 
      modified_sites: 'data/miCLIP_union.bed' #BED6 file containing known modified sites where 'name' refers to the annotation of the position. Useful for learning patterns (training and test set).
      chr_size: "data/hg38.size"  #file contaning size of chromosomes (Chromosome     | size )
      regions: "data/selected_regions.bed" # BED6 file contaning set of 5-mer (NNANN) to analyze, if ="", all 5-mers (NNANN) will be considered.
      data: # a dictionary of two keys (cond1, cond2) referring to the paired conditions inputs. The value is the list of replicates names without ".bam" extension.
           cond1: ["HEK293T-WT-rep2","HEK293T-WT-rep3"]
           cond2: ["HEK293T-KO-rep2","HEK293T-KO-rep3"]
      jacusa_params: # dictionary where keys refer to parameters (e.g. [p: 16] to set the number of threads to 16). Please use "" if no value is affected to the parameter. We use the following parameters:
           P1: 'FR-SECONDSTRAND'  # mandatory parameter referring to the library of the first condition sample.
           P2: 'FR-SECONDSTRAND'  # mandatory parameter referring to the library of the second condition sample.
           m: 1  # filter reads by mapping quality
           q: 1  # filter reads by base calling quality
           c: 4  # filter reads by coverage
           a: 'Y'  # recommended parameters to filter sites within the holy-polymer regions.
           p: 16    # parameter to customize the number of threads
           D: ''  # mandatory parameter to output deletion score.
           I: ''  # mandatory parameter to output insertion score.
      pattern_params:       # specify patterns and their combinations to be used, please use "" if no value is affected to the field.
           internal_pattern: "Boulias,Koertel,Koh" # specify the annotation of the set of modified sites to be used as a training set. in case you use an external pattern put "". 
           external_pattern: ""  # path to an external pattern in case you don't use internal_pattern, else put ""
           combined_patterns: # patterns to combine, add as many combinations as you want as a [key(any name): value (pattern numbers)] combination.
                      pt1: [1,2,4,6]  
                      pt2: [1,2,3,4,6]


Please check JACUSA2 [manual](https://github.com/dieterich-lab/JACUSA2) for more details on how to use JACUSA2 parameters.

- Run `JACUSA2_call2` target to call variants using paired conditions.
```  
$ snakemake --cores all jacusa2_call2
```
The output is a file called `Cond1vsCond2Call2.out` under `./output/{label}/jacusa` and filtered BAM files under `./output/{label}/bam`.

- Run `get_features` target to preprocess JACUSA2 call2 output and extract features.
```
$ snakemake --cores all get_features
```
The output is an R object `features.rds` under `./output/{label}/features/`.

- Run `get_pattern` target to learn patterns representing m6A modification.
```
$ snakemake --cores all get_pattern
```
The output is an R object `NMF.rds` containing the factorization result, including basis and coefficient matrices, plus, plots showing the rank selection result, pattern scoring barplot and heatmap of NMF resulting matrices. The output is under `./output/{label}/pattern/`. Implicitly, training and test set files (resp. `train_features.rds`, `test_features.rds`) are created under `features` folder and, subsequently, used for the learning model.

For the testing example, `train_features.rds` is supposed to contain 1905 sites. The scoring of patterns will look like the following barplot.

<p align="center">
  <img src="https://github.com/dieterich-lab/MiMB_JACUSA2_chapter/blob/main/img/pattern_scores.png?raw=true" width="300">
</p>
<p align="center"> 
  <em>Figure 1: Membership score of resulted patterns</em>
</p>

- Run `visualize_pattern` target to visualize patterns and their combinations.
```
$ snakemake --cores all visualize_pattern
```  
The output is a set of figures representing barplots for the produced patterns. The outputs can be found under `./output/{label}/pattern/viz/`.

For the testing example, the combination of patterns representing more than 80% of the training set will look like the following barplot:

<p align="center">
  <img src="https://github.com/dieterich-lab/MiMB_JACUSA2_chapter/blob/main/img/barplot_NMF.png?raw=true" width="300">
</p>
<p align="center"> 
  <em>Figure 2: Combination of patterns representing 80% of the training set</em>
</p>

- Run `predict_modification` target to predict modified sites
```
$ snakemake --cores all predict_modification
```  
The output is a BED6 file(s) containing scores of the selected pattern(s) for the test set under `./output/{label}/prediction/` and the corresponding eCDF (empirical cumulative distribution) and PPV (positive predictive values) plots.

For the testing example, the eCDF will look like the following figure: 

<p align="center">
  <img src="https://github.com/dieterich-lab/MiMB_JACUSA2_chapter/blob/main/img/Pattern_ecdf.png?raw=true" width="500">
</p>
<p align="center"> 
  <em>Figure 3: eCDF of estimated scores of sites from diffrent miCLIP categories and non miCLIP sites</em>
</p>

Note that rules are linked so that the workflow is determined from top (e.g. predict modification) to bottom (e.g. sort bam) and
executed accordingly from bottom to top. Therefore, running ”predict_modification” rule leads to executing all rules on its pipeline. Check [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/) for further details.


<p align="center">
  <img src="https://github.com/dieterich-lab/MiMB_JACUSA2_chapter/blob/main/img/snakemake.png?raw=true" width="500">
</p> 
<p align="center"> 
  <em>Figure 4: Snakemake JACUSA2-based pipeline</em>
</p>

# Output files
Once the pipeline has run successfully you should expect the following files in the output directory:
*   **`bam/`:**  This will be skipped if you start from JACUSA2 output
    *   `{mapping BAM file name}.sorted.bam` - sorted BAM file
    *   `{mapping BAM file name}.filtered.bam` - filtered BAM file
    *   `{mapping BAM file name}.filtered.bam.bai` - index of the filtered BAM file
*   **`jacusa/`:**
    *   `Cond1vsCond2Call2.out` - JACUSA2 call-2 output
*   **`features/`:**
    *   `features.rds` - tabular features (target sites (here, 'A') x 15 features)
    *   `train_features.rds` - training set to extract RNA modification patterns
    *   `test_features.rds`- test set
*   **`patterns/`:**
    *   `NMF.rds` - NMF factorization (R object)
    *   `asses_NMF_1.pdf` - NMF rank survey
    *   `asses_NMF_2.pdf` - silhouette + cophenetic correlation result
    *  `pattern_scores.pdf` - barplot of the membership score of patterns based on basis matrix
    *  `NMF_matrices.pdf` - heatmaps for basis and coefficient matrices of the NMF result
 
    *   **`viz/`:**
        *  `pattern_{pattern number}_barplot_NMF.pdf` - barplot of patterns (by position and score type) from coefficient matrix   
*   **`prediction/`:**
    *   `pattern{pattern number}_prediction.bed` - prediction scores from selected patterns
    *   `pattern{pattern number}_ppv.pdf` - PPV plot
    *   `pattern{pattern number}_ecdf.pdf` - eCDF plot
    
# Dependencies and versions

All dependencies and packages are added to the `environment.yaml` file.

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
