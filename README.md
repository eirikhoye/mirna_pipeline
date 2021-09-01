## Analysis pipeline for "A microRNA Signature of Metastatic Colorectal Cancer"
A reproducible analysis pipeline for the manuscript "A microRNA Signature of Metastatic Colorectal Cancer" by Høye et al.

Integrates both public repositories and own smallRNA-seq data. Quality control with the expertly developed [miRTrace](https://github.com/friedlanderlab/mirtrace) quality control and contamination detection pipeline. Data processing and read alignment with [miRge3.0](https://github.com/mhalushka/miRge3.0), using our curated database [MirGeneDB](https://mirgenedb.org/) as a high quality miRNA gene reference. miRge3.0 count matrixe used for downstream analysis, including data exploration such with UMAP, differential expression analysis with DESeq2 and inference of bulk tissue cell composition with known cell marker miRNA.![Pipeline Flowchart](/images/Pipeline_Flowchart2.png)

### References:

miRTrace [Genome Biology (2018)](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1588-9). 
> Kang W; Eldfjell Y; Fromm B; Estivill X; Biryukova I; Friedländer MR, 2018. miRTrace reveals the organismal origins of microRNA sequencing data. Genome Biol 19(1):213

miRge3.0 [bioRxiv (2021)](https://www.biorxiv.org/content/10.1101/2021.01.18.427129v1). 
> Patil HP; Halushka MR miRge3.0: a comprehensive microRNA and tRF sequencing analysis pipeline. bioRxiv 2021.01.18.427129

### Set up singularity container
```
# on a host with singularity installed, first clone the repository
git clone https://github.com/eirikhoye/mirna_pipeline
cd mirna_pipeline/

# Configure miRTrace container
sudo singularity build singularity/miRTrace.simg singularity/miRTrace.recipe

# Alternatively, install miRTrace by following instructions here:
https://github.com/friedlanderlab/mirtrace

# Configure miRge3.0 container
sudo singularity build singularity/mirge3.simg singularity/mirge3.recipe

# Alternatively, install miRge3.0 by reading instructions here for relevant OS: 
https://mirge3.readthedocs.io/en/latest/installation.html
# and make sure dependencies are installed (python=3.7 and r) and bowtie, samtools and RNAfold
```
### Usage
```
# running miRTrace singularity container
dt=`date '+%d%m%Y_%H%M%S'`

sudo singularity exec --bind /path/to/project_folder_on_host:/mnt /path/to/miRTrace.simg mirtrace qc --species hsa --custom-db-folder /mnt/custom_databases/ --config /mnt/config -o /mnt/mirtrace_out/"$dt"

# note, to work with your data in singularity you must mount the directory it is stored in (here called "project_folder_on_host") on you host to a path inside the singularity container, here /mnt, with the syntax: /path/to/host_dir:/path/to/singularity_dir/. This way we can run singularity on files stored on our system. See https://sylabs.io/guides/3.0/user-guide/bind_paths_and_mounts.html for details.



# running miRge3.0
sudo singularity exec --bind /path/to/project_folder_on_Host:/mnt /path/to/mirge3.simg miRge3.0 -s /mnt/filepaths.txt -lib /mnt/miRge3_Lib -on human -db mirgenedb -o /mnt/mirge_output -tcf -cpu 4 -a illumina




```
## Tutorial
Lets give an example using toy fastq files fastq_toy directory. First lets do QC on these samples with miRTrace. miRTrace takes raw FASTQ files as input and outputs nicely formatted QC reports, and will also assess potential contamination! 

miRTrace requires a config .csv file with paths, name and library adapter sequence, of the format:
```
/mnt/fastq_toy/mLi_1.fq,mLi_1,TGGAATTC
/mnt/fastq_toy/mLi_2.fq,mLi_2,TGGAATTC
/mnt/fastq_toy/mLi_3.fq,mLi_3,TGGAATTC
/mnt/fastq_toy/mLu_1.fq,mLu_1,TGGAATTC
/mnt/fastq_toy/mLu_2.fq,mLu_2,TGGAATTC
/mnt/fastq_toy/mLu_3.fq,mLu_3,TGGAATTC
/mnt/fastq_toy/nCR_1.fq,nCR_1,TGGAATTC
/mnt/fastq_toy/nCR_2.fq,nCR_2,TGGAATTC
/mnt/fastq_toy/nCR_3.fq,nCR_3,TGGAATTC
/mnt/fastq_toy/nLi_1.fq,nLi_1,TGGAATTC
/mnt/fastq_toy/nLi_2.fq,nLi_2,TGGAATTC
/mnt/fastq_toy/nLi_3.fq,nLi_3,TGGAATTC
/mnt/fastq_toy/nLu_1.fq,nLu_1,TGGAATTC
/mnt/fastq_toy/nLu_2.fq,nLu_2,TGGAATTC
/mnt/fastq_toy/nLu_3.fq,nLu_3,TGGAATTC
/mnt/fastq_toy/pCRC_1.fq,pCRC_1,TGGAATTC
/mnt/fastq_toy/pCRC_2.fq,pCRC_2,TGGAATTC
/mnt/fastq_toy/pCRC_3.fq,pCRC_3,TGGAATTC
/mnt/fastq_toy/contaminated_1.fq,contaminated_1,CGCCTTGGCCGTA
/mnt/fastq_toy/contaminated_2.fq,contaminated_2,CGCCTTGGCCGTA
/mnt/fastq_toy/contaminated_3.fq,contaminated_3,CGCCTTGGCCGTA
/mnt/fastq_toy/read_len_1.fq,read_len_1,TCGTATGC
/mnt/fastq_toy/read_len_2.fq,read_len_2,TCGTATGC
/mnt/fastq_toy/read_len_3.fq,read_len_3,TCGTATGC

# note we set /mnt as prefix to path as we defined our project folder as /mnt
```






















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
# and make sure dependencies are installed (python=3.7 and r) and bowtie, samtools and RNAfold

conda install -c bioconda mirge3

# miRge3 is also available as a docker container:
https://quay.io/repository/biocontainers/mirge3?tab=tags

```

### Usage
```
# create miRTrace report:
bash run_mirtrace.sh config

# run mirge3.0
bash run_mirge3.sh filepaths.txt

```

## Tutorial

Getting started using the five tutorial FASTQ files found in data/fastq/
```
data/fastq/sub_M12_1.fq.gz
data/fastq/sub_M18.fq.gz
data/fastq/sub_M19.fq.gz
data/fastq/sub_SRR1646483.fastq.gz
data/fastq/sub_SRR837850.fastq.gz
```
First, clone the repository:
```
git clone https://github.com/eirikhoye/mirna_pipeline
```


### miRTrace
Assessing read quality with miRTrace
```
# First activate mirtrace environment:
conda activate mirtrace

# Set up a config file, like similar to example below: (first column is path to file, second is sample name, third is adapter sequence)
data/fastq/sub_M12_1.fq.gz,M12,TGGAATTC
data/fastq/sub_M18.fq.gz,M18,TGGAATTC
data/fastq/sub_M19.fq.gz,M19,TGGAATTC
data/fastq/sub_SRR1646483.fastq.gz,SRR1646483,CGCCTTGG
data/fastq/sub_SRR837850.fastq.gz,SRR837826,TCGTATGC

# Run miRTrace to make QC and contamination report:
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

Align reads to MirGeneDB2.0 and create count matrix with miRge3
```
# create a text file with paths to fastq files, for example:
data/fastq/sub_M12_1.fq.gz
data/fastq/sub_M18.fq.gz
data/fastq/sub_M19.fq.gz
data/fastq/sub_SRR1646483.fastq.gz
data/fastq/sub_SRR837850.fastq.gz

# run the mirge3.0 with scripts/run_mirge3.sh and filepaths.txt for tutorial FASTQ
bash scripts/run_mirge3.sh filepaths.txt

# or read documentation for more customisation: https://mirge3.readthedocs.io/en/latest/

# RPM/count matrix is in:
counts: data/mirge3_output/<date_time>/miR.RPM.csv
RPM: data/mirge3_output/<date_time>/miR.Counts.csv
```




