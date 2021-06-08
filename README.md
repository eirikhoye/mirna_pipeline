## Analysis pipeline for "A microRNA Signature of Metastatic Colorectal Cancer"
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
mirge3.0
```
### Installing
```
# install miRTrace
conda env create -f environment.yml

# Install miRge3.0
# Read instructions here for relevant OS: https://mirge3.readthedocs.io/en/latest/installation.html
# and make sure dependencies are installed (python=3.7 and r)

# Then install mirge3.0 in a new conda environment to avoid conflicts:
conda create --name mirge3
conda activate mirge3
conda install -c bioconda mirge3
```

### Usage
```
# create miRTrace report:
bash run_mirtrace.sh config

# run mirge3.0
bash run_mirge3.sh filepaths.txt

```

## Tutorial
### miRTrace
Assessing read quality with miRTrace
```
# First create conda environment with (should only be done once):
conda env create -f environment.yml

# then activate environment with (needs to be done for each new terminal session):
conda activate miRNA_pipeline

# Set up a config file like this for the tutorial files in this repository:
data/fastq/sub_M12_1.fq.gz,M12,TGGAATTC
data/fastq/sub_M18.fq.gz,M18,TGGAATTC
data/fastq/sub_M19.fq.gz,M19,TGGAATTC
data/fastq/sub_SRR1646483.fastq.gz,SRR1646483,CGCCTTGG
data/fastq/sub_SRR837850.fastq.gz,SRR837826,TCGTATGC
# (The last column is adapter sequence used for sequencing.)

# Create miRTrace QC and contamination report:
bash run_mirtrace.sh config

# This will create a report in data/mirtrace_out/<date_time>/mirtrace-report.html
```

Check to see if the quality of reads are above thresholds:
It is important for miRNA data that the sequencing has been size selected, otherwise reads from other RNA types will be present.
![Read Length](/images/mirtrace-length-plot.png)

Above 4 out of 5 samples had acceptable read length distribution, but SRR837826 has a large number of reads in the piRNA range (~ 30nt) and none in the miRNA range (~18-26nt). 

smallRNA-seq datasets should have a majority of reads in the miRNA length range:

![RNA Type](/images/mirtrace-rnatype-plot.png)

Not surprisingly, in the 5th sample, which did not have any reads in the miRNA length range, no miRNA reads were detected.

Sequencing data should also be checked for contamination. miRTrace can detect which clade a miRNA read maps to:

![Read Contamination](/images/mirtrace-contamination-plot.png)

Here we also see that the 4th sample has a major contamination of Bird/Reptile miRNA reads. Shuch a dataset should certainly not be used in analysing human disease!

### miRge3.0

Using miRge3.0 to align reads to MirGeneDB2.0
```
# start conda environment
conda activate mirge3

# create a text file with paths to fastq files, for example:
data/fastq/sub_M12_1.fq.gz
data/fastq/sub_M18.fq.gz
data/fastq/sub_M19.fq.gz
data/fastq/sub_SRR1646483.fastq.gz
data/fastq/sub_SRR837850.fastq.gz

# run the mirge3.0 with, for instance, the script run_mirge3.sh, 
# or read documentation for more customisation: https://mirge3.readthedocs.io/en/latest/
# miRge3.0 -s $1 -lib miRge3_Lib -on human -db mirgenedb -o output_dir -tcf -cpu 4 -a illumina

bash scripts/run_mirge3.sh filepaths.txt

# output
counts: data/mirge3_output/<date_time>/miR.RPM.csv
RPM: data/mirge3_output/<date_time>/miR.Counts.csv
```




