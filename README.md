# Towards a consensus microRNA signature of primary and metastatic colorectal cancer
A reproducible pipeline for the integration of public repositories of smallRNA-seq data. Runs the expertly developed miRTrace quality control and contamination detection pipeline, and manually curated MirGeneDB.org repository of bona fide miRNA genes as reference. Outputs both fasta files, sam files as well as count matrix for ease of integration of downstream analysis, including data exploration such as UMAP, and differential expression analysis with DESeq2.

![Pipeline Flowchart](/images/Pipeline_Flowchart2.png)

Quality Control and read processing with miRTrace
![miRTrace logo](/images/miRTrace.png)

> Kang W; Eldfjell Y; Fromm B; Estivill X; Biryukova I; Friedländer MR, 2018. miRTrace reveals the organismal origins of microRNA sequencing data. Genome Biol 19(1):213 [miRTrace](https://github.com/friedlanderlab/mirtrace)



## Getting started

### Dependencies
```
bowtie=1.2.3
biopython
samtools
mirtrace
fastx-toolkit
bioconductor-rsubread
```
### Installing
```
conda env create -f environment.yml
```

### Usage
```
bash run_pipeline.sh <sample_info.tsv>

```

### Tutorial
```
# First create conda environment with (should only be done once):
conda env create -f environment.yml

# then activate environment with (needs to be done for each new terminal session):
conda activate miRNA_pipeline

# To try the pipeline, use the tutorial fastq files and tutorial_sample_info.tsv with:
bash run_pipeline.sh tutorial_sample_info.tsv

```
