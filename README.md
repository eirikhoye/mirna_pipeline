# Analysis pipeline for "A microRNA Signature of Metastatic Colorectal Cancer"
A reproducible analysis pipeline for the manuscript "A microRNA Signature of Metastatic Colorectal Cancer" by Høye et al.

Integrates both public repositories and own smallRNA-seq data. Quality control with the expertly developed [miRTrace](https://github.com/friedlanderlab/mirtrace) quality control and contamination detection pipeline. Data processing and read alignment with [miRge3.0](https://github.com/mhalushka/miRge3.0), using our curated database [MirGeneDB](https://mirgenedb.org/) as a high quality miRNA gene reference. miRge3.0 count matrixe used for downstream analysis, including data exploration such with UMAP, differential expression analysis with DESeq2 and inference of bulk tissue cell composition with known cell marker miRNA.

![Pipeline Flowchart](/images/Pipeline_Flowchart2.png)

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
