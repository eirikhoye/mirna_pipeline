# Towards a consensus microRNA signature of primary and metastatic colorectal cancer
A reproducible pipeline for the integration of public repositories of smallRNA-seq data. Runs the expertly developed miRTrace quality control and contamination detection pipeline, and manually curated MirGeneDB.org repository of bona fide miRNA genes as reference. Outputs both fasta files, sam files as well as count matrix for ease of integration of downstream analysis, including data exploration such as UMAP, and differential expression analysis with DESeq2.

![Pipeline Flowchart](/images/Pipeline_Flowchart.png)


## Getting started

### Dependencies
bowtie=1.2.3
biopython
samtools
mirtrace
fastx-toolkit
bioconductor-rsubread

### Installing
```
conda env create -f environment.yml
```

### Usage
```
bash scripts/create_mirtrace_config.sh <path/to/fastq_files> <name> <adapter_sequence>
bash scripts/run_mirtrace.sh <path/to/config_file>
bash scripts/bowtie_mirgenedb_human_fastafiles.sh
Rscript scripts/featureCounts.R

```
