15 February, 2021

``` r
# import libraries
library("tidyverse")
library("DESeq2")
library("knitr")
library('kableExtra')
library("RColorBrewer")
library("pheatmap")
#library("scatterplot3d")
#library("gplots")
library("dplyr")
#library("rBLAST")
#library("Biostrings")
library("BiocParallel")
library("broman")
library("gridExtra")
library("ggpubr")
library("FactoMineR")
library("factoextra")
```

``` r
register(MulticoreParam(12))
```

``` r
"
Define threshold for signature miRNA, including effect size, 
significance therhold, and expression size
"
```

    ## [1] "\nDefine threshold for signature miRNA, including effect size, \nsignificance therhold, and expression size\n"

``` r
lfc.Threshold <- 0.5849625
rpm.Threshold <- 100
p.Threshold   <- 0.05
```

# Import MirGeneDB metadata

``` r
MirGeneDB_info <- read_delim('/Users/eirikhoy/Dropbox/projects/comet_analysis/data/hsa_MirGeneDB_to_miRBase.csv', delim = ';')
```

    ## Parsed with column specification:
    ## cols(
    ##   MirGeneDB_ID = col_character(),
    ##   MiRBase_ID = col_character(),
    ##   Family = col_character(),
    ##   Seed = col_character(),
    ##   `5p accession` = col_character(),
    ##   `3p accession` = col_character(),
    ##   Chromosome = col_character(),
    ##   Start = col_double(),
    ##   End = col_double(),
    ##   Strand = col_character(),
    ##   `Node of origin (locus)` = col_character(),
    ##   `Node of origin (family)` = col_character(),
    ##   `3' NTU` = col_character(),
    ##   ` UG ` = col_double(),
    ##   UGUG = col_double(),
    ##   CNNC = col_double()
    ## )

``` r
MirGeneDB_info <- MirGeneDB_info %>% filter(!grepl("-v[2-9]", MirGeneDB_ID)) # keep only -v1
MirGeneDB_info$MirGeneDB_ID <- str_replace_all(MirGeneDB_info$MirGeneDB_ID, "-v1", "")
```

# Functions

``` r
DeseqObject <- function(DESIGN, countdata, coldata, consensus="None", sample_type="None", Ref) {
  "
  Function to create DESeq2 object
  "
  
  dds <- DESeqDataSetFromMatrix(countData = countdata,
                                colData = coldata,
                                design = as.formula(paste("~", DESIGN)))
    # Kick out non-consensus samples
  if (!(consensus == "None")) {
    dds <- dds[, dds$paper %in% consensus]
  }

  # Kick out samples that are not bulk tissue
  if (!(sample_type == "None")) {
    dds <- dds[, dds$sample_type == sample_type]
  }

  dds$tissue.type <- relevel(dds$tissue.type, ref=ref)
  dds$tissue.type <- droplevels(dds$tissue.type)
  
  dds <- DESeq(dds,
               parallel=TRUE,
               BPPARAM=MulticoreParam(3)
               )

  return(dds)
}
```

``` r
DeseqResult <- function(dds, column, coef, tissue_type_A, tissue_type_B, 
                        lfc.Threshold, rpm.Threshold,
                        norm_adj_up       = "None",
                        norm_adj_down     = "None",
                        pCRC_adj_up   = "None",
                        pCRC_adj_down = "None"){
  "
  Function to return results from DESeq2 for different conditions, 
  including control for normal adjacent tissue, if available
  "
    p.threshold   <- p.Threshold
    lfc.threshold <- lfc.Threshold
    rpm.threshold <- rpm.Threshold
    samples_tissue_type_A <- colData(dds)[, column] == tissue_type_A
    samples_tissue_type_B <- colData(dds)[, column] == tissue_type_B
  res <- results(dds, name = coef, alpha = p.threshold)
  res <- lfcShrink(dds, coef=coef, res=res)
    rpm <- t(t(counts(dds)) / colSums(counts(dds))) * 1000000
    sig <- rownames(res[(abs(res$log2FoldChange) > lfc.threshold) &
                        (res$padj < p.threshold) &
                        !is.na(res$padj), ])
    sig <- sig[ (rowMeans(rpm[sig, samples_tissue_type_A]) > rpm.threshold) |
                (rowMeans(rpm[sig, samples_tissue_type_B]) > rpm.threshold)
                ]
    res_sig    <- res[sig, ]
    up_mirna   <- rownames(res_sig[res_sig$log2FoldChange > lfc.Threshold, ])
    down_mirna <- rownames(res_sig[res_sig$log2FoldChange < -lfc.Threshold, ])
  
    if (!(norm_adj_up == "None")) {
      up_mirna <- setdiff(up_mirna, norm_adj_up)
    }
    if (!(norm_adj_down == "None")) {
      down_mirna <- setdiff(down_mirna, norm_adj_down)
    }
    if (!(pCRC_adj_up == "None")) {
      up_mirna <- setdiff(up_mirna, pCRC_adj_up)
    }
    if (!(pCRC_adj_down == "None")) {
      down_mirna <- setdiff(down_mirna, pCRC_adj_down)
    }
        
    return_list <- list("rpm" = rpm, "res" = res, "sig"=sig, "res_sig"=res_sig, 
                        "down_mirna"=down_mirna, "up_mirna"=up_mirna)
    
    return(return_list)

}
```

``` r
SigList <- function(res, dds, tissue_type_A, tissue_type_B, coef,
                            norm_adj_up, norm_adj_down, 
                            pCRC_adj_up, pCRC_adj_down){
  "
  Function to create annotated lists of signature miRNA
  Return will print upregulated or downregulated miRNA, 
  by printing <signature_list>$up_mirna
  or          <signature_list>$down_mirna
  "
  group_A_rpm <- rowMeans(res$rpm[res$sig, dds$tissue.type == tissue_type_A])
  group_A_rpm_std <- rowSds(res$rpm[res$sig, dds$tissue.type == tissue_type_A])
  group_B_rpm <- rowMeans(res$rpm[res$sig, dds$tissue.type == tissue_type_B])
  group_B_rpm_std <- rowSds(res$rpm[res$sig, dds$tissue.type == tissue_type_B])
  lfc.deseq2  <- res$res[res$sig, ]$log2FoldChange
  lfcSE.deseq2<- res$res[res$sig, ]$lfcSE
  neg.log.10.adj.p <- format(res$res[res$sig, ]$padj, digits=3)
  signature_mirna <- res$sig
  sig_list <- dplyr::tibble(signature_mirna, lfc.deseq2, lfcSE.deseq2,
                            group_A_rpm, #group_A_rpm_std,
                            group_B_rpm, #group_B_rpm_std,
                            neg.log.10.adj.p)
  sig_list$signature_sub <- str_replace_all(signature_mirna, "/.*", "") %>% str_replace_all(., c("_5p" = "", "_3p" = ""))
  sig_list <- left_join(sig_list, MirGeneDB_info, by=c("signature_sub" = "MirGeneDB_ID"))
  
  # create list of upregulated mirna
  up_mirna <- sig_list %>%
    filter(lfc.deseq2 > lfc.Threshold) %>% 
    
    # Annotate which miRNA are cell markers
    mutate(
      cell_marker = ifelse(signature_mirna %in% names(cell_spec_dict_inv), cell_spec_dict_inv[signature_mirna], '')) %>%
    mutate(
      cell_marker = cell_spec(cell_marker, color = ifelse(cell_marker != '', 'white', 'black'),
                              background = ifelse(cell_marker != '', 'blue', 'white'),
                              bold = ifelse(cell_marker != '', F, F)))

    # Annotate which miRNA are in normal_adjacent
    if (norm_adj_up != "None") {
      up_mirna <- up_mirna %>%
        mutate(
          norm_adj = ifelse(signature_mirna %in% norm_adj_up, 'yes', '')) %>%
        mutate(
          norm_adj = cell_spec(norm_adj, color = ifelse(norm_adj == 'yes', 'white', 'black'),
                               background = ifelse(norm_adj == 'yes', 'black', 'white'),
                               bold = ifelse(norm_adj == 'yes', F, F))
        )
    }
    else up_mirna$norm_adj <- "na"
    
  # Annotate which miRNA are in pCRC_adjacent
    if (pCRC_adj_up != "None") {
      up_mirna <- up_mirna %>%
        mutate(
          pCRC_adj = ifelse(signature_mirna %in% pCRC_adj_up, 'yes', '')) %>%
        mutate(
          pCRC_adj = cell_spec(pCRC_adj, color = ifelse(pCRC_adj == 'yes', 'white', 'black'),
                               background = ifelse(pCRC_adj == 'yes', 'black', 'white'),
                               bold = ifelse(pCRC_adj == 'yes', F, F))
        )
    }
    else up_mirna$pCRC_adj <- "na"

    # number of upregulated miRNA
    number_upregulated <- dim(up_mirna)[1]
    
    # select only relevant rows
    up_mirna <- up_mirna %>% select(signature_mirna, lfc.deseq2, lfcSE.deseq2, neg.log.10.adj.p, 
                                    group_A_rpm,  #group_A_rpm_std,
                                    group_B_rpm, #group_B_rpm_std,
                                    MiRBase_ID, Family, Seed, Chromosome,
                                    cell_marker, norm_adj, pCRC_adj)
    
    
  # Create kable list with annotations
    up_mirna <- up_mirna %>%
      arrange(-lfc.deseq2) %>%
      arrange(desc(cell_marker)) %>%
      arrange(pCRC_adj) %>%
      arrange(norm_adj) %>%
      kable(col.names = c("miRNA", "LFC", "lfcSE", "FDR", 
                          paste('RPM', tissue_type_A), #paste('std', tissue_type_A), 
                          paste('RPM', tissue_type_B), #paste('std', tissue_type_B), 
                          "miRBase_ID", "Family", "Seed", "Chr",
                          "Cell-Type Specific", 'Norm Background', 'pCRC Background'),
            escape = F, booktabs = F, caption = paste("Upregulated in ", coef),
            digits = c(0, 2, 2, 3, 0, 0, 2, 3, 0, 0, 0, 0, 0, 0, 0)) %>%
      kable_styling(bootstrap_options = c("striped", "hover", "condensed"), full_width = T, 
                  fixed_thead = list(enabled = T)) %>%
      scroll_box(width = "2000px")
    
 
  # create list of downregulated miRNA
  down_mirna <- sig_list %>%
    filter(lfc.deseq2 < -lfc.Threshold) %>% 
    
    # Annotate which miRNA are cell markers
    mutate(
      cell_marker = ifelse(signature_mirna %in% names(cell_spec_dict_inv), cell_spec_dict_inv[signature_mirna], '')) %>%
    mutate(
      cell_marker = cell_spec(cell_marker, color = ifelse(cell_marker != '', 'white', 'black'),
                              background = ifelse(cell_marker != '', 'blue', 'white'),
                              bold = ifelse(cell_marker != '', F, F)))

    # Annotate which miRNA are in normal_adjacent
    if (norm_adj_down != "None") {
      down_mirna <- down_mirna %>%
        mutate(
          norm_adj = ifelse(signature_mirna %in% norm_adj_down, 'yes', '')) %>%
        mutate(
          norm_adj = cell_spec(norm_adj, color = ifelse(norm_adj == 'yes', 'white', 'black'),
                               background = ifelse(norm_adj == 'yes', 'black', 'white'),
                               bold = ifelse(norm_adj == 'yes', F, F))
        )
    }
    else down_mirna$norm_adj <- "na"

    # Annotate which miRNA are in pCRC_adjacent
    if (pCRC_adj_down != "None") {
      down_mirna <- down_mirna %>%
        mutate(
          pCRC_adj = ifelse(signature_mirna %in% pCRC_adj_down, 'yes', '')) %>%
        mutate(
          pCRC_adj = cell_spec(pCRC_adj, color = ifelse(pCRC_adj == 'yes', 'white', 'black'),
                               background = ifelse(pCRC_adj == 'yes', 'black', 'white'),
                               bold = ifelse(pCRC_adj == 'yes', F, F))
        )
    }
    else down_mirna$pCRC_adj <- "na"
  
    # number of upregulated miRNA
    number_downregulated <- dim(down_mirna)[1]

    down_mirna <- down_mirna %>% select(signature_mirna, lfc.deseq2, lfcSE.deseq2, neg.log.10.adj.p,
                                        group_A_rpm,  #group_A_rpm_std,
                                        group_B_rpm, #group_B_rpm_std,
                                        MiRBase_ID, Family, Seed, Chromosome,
                                        cell_marker, norm_adj, pCRC_adj)    
    
  # Create kable list with annotations    
    down_mirna <- down_mirna %>%
      arrange(lfc.deseq2) %>%
      arrange(desc(cell_marker)) %>%
      arrange(pCRC_adj) %>%
      arrange(norm_adj) %>%
      kable(col.names = c("miRNA", "LFC", "lfcSE", "FDR",
                          paste('RPM', tissue_type_A), #paste('std', tissue_type_A), 
                          paste('RPM', tissue_type_B), #paste('std', tissue_type_B), 
                          "miRBase_ID", "Family", "Seed", "Chr",
                           "Cell-Type Specific", 'Norm Background', 'pCRC Background'),
            escape = F, booktabs = F, caption = paste("Downregulated in ", coef),
            digits = c(0, 2, 2, 3, 0, 0, 2, 3, 0, 0, 0, 0, 0, 0, 0)) %>%
      kable_styling(bootstrap_options = c("striped", "hover", "condensed"), full_width = T, 
                  fixed_thead = list(enabled = T)) %>%
       scroll_box(width = "2000px")

  # Function return is to print kable, either upregulated or downregulated miRNA
  return_list = list("up_mirna" = up_mirna, "down_mirna" = down_mirna, 
                     "number_upregulated" = number_upregulated, 
                     "number_downregulated" = number_downregulated)
  return(return_list)
}    
```

``` r
# Read the sample information into a data frame
sampleinfo <- read_delim("/Users/eirikhoy/Dropbox/projects/comet_analysis/data/sample_info_v9.csv", delim=';')
```

    ## Parsed with column specification:
    ## cols(
    ##   filename = col_character(),
    ##   paper = col_character(),
    ##   sample_name = col_character(),
    ##   type.tissue = col_character(),
    ##   type = col_character(),
    ##   tissue = col_character(),
    ##   paper_sample_name = col_character(),
    ##   `3p-adapter` = col_character(),
    ##   qc_report = col_character(),
    ##   malignant = col_character(),
    ##   new_old = col_character()
    ## )

``` r
sampleinfo <- sampleinfo %>%
  filter(qc_report == 'keep')
sampleinfo$filename <- str_remove(sampleinfo$filename, '.fasta.fas.gz.bam')
sampleinfo$filename <- str_remove(sampleinfo$filename, '.fasta.bam')
sampleinfo$filename <- str_replace_all(sampleinfo$filename, pattern = '\\.', replacement = '_')
sampleinfo$filename <- str_replace_all(sampleinfo$filename, pattern = '-', replacement = '_')
#sampleinfo$filename <- str_replace_all(sampleinfo$filename, pattern = '__', replacement = '_')

# Read the data into R
#seqdata <- read_delim("/Users/eirikhoy/Dropbox/projects/comet_analysis/data/count_matrix_08.12.20.csv", delim = ';')
seqdata_1 <- read_delim("/Users/eirikhoy/Dropbox/projects/mirge3/output_dir/miRge.2021-01-19_10-03-25/miR.Counts.csv", delim = ',')
```

    ## Parsed with column specification:
    ## cols(
    ##   .default = col_double(),
    ##   miRNA = col_character()
    ## )

    ## See spec(...) for full column specifications.

``` r
#seqdata_2 <- read_delim("/Users/eirikhoy/Dropbox/projects/mirge3/output_dir/miRge.2021-01-19_12-04-26/miR.Counts.csv", delim = ',')
seqdata_3 <- read_delim("/Users/eirikhoy/Dropbox/projects/mirge3/output_dir/miRge.2021-01-19_15-00-38/miR.Counts.csv", delim = ',')
```

    ## Parsed with column specification:
    ## cols(
    ##   miRNA = col_character(),
    ##   SRR1273998 = col_double(),
    ##   SRR1273999 = col_double(),
    ##   SRR1274000 = col_double(),
    ##   SRR1274001 = col_double()
    ## )

``` r
seqdata_4 <- read_delim("/Users/eirikhoy/Dropbox/projects/mirge3/output_dir/miRge.2021-01-19_15-33-44/miR.Counts.csv", delim = ',')
```

    ## Parsed with column specification:
    ## cols(
    ##   .default = col_double(),
    ##   miRNA = col_character()
    ## )
    ## See spec(...) for full column specifications.

``` r
seqdata_5 <- read_delim("/Users/eirikhoy/Dropbox/projects/mirge3/output_dir/miRge.2021-01-19_17-33-30/miR.Counts.csv", delim = ',')
```

    ## Parsed with column specification:
    ## cols(
    ##   .default = col_double(),
    ##   miRNA = col_character()
    ## )
    ## See spec(...) for full column specifications.

``` r
seqdata_6 <- read_delim("/Users/eirikhoy/Dropbox/projects/mirge3/output_dir/miRge.2021-02-04_08-16-37/miR.Counts.csv", delim = ',')
```

    ## Parsed with column specification:
    ## cols(
    ##   .default = col_double(),
    ##   miRNA = col_character()
    ## )
    ## See spec(...) for full column specifications.

``` r
seqdata_7 <- read_delim("/Users/eirikhoy/Dropbox/projects/mirge3/output_dir/miRge.2021-02-04_09-09-56/miR.Counts.csv", delim = ',')
```

    ## Parsed with column specification:
    ## cols(
    ##   .default = col_double(),
    ##   miRNA = col_character()
    ## )
    ## See spec(...) for full column specifications.

``` r
seqdata_8 <- read_delim("/Users/eirikhoy/Dropbox/projects/mirge3/output_dir/miRge.2021-02-04_12-11-38/miR.Counts.csv", delim = ',')
```

    ## Parsed with column specification:
    ## cols(
    ##   .default = col_double(),
    ##   miRNA = col_character()
    ## )
    ## See spec(...) for full column specifications.

``` r
seqdata <- inner_join(inner_join(inner_join(inner_join(inner_join(inner_join(seqdata_1, seqdata_3), seqdata_4), seqdata_5), seqdata_6), seqdata_7), seqdata_8)
```

    ## Joining, by = "miRNA"

    ## Joining, by = "miRNA"
    ## Joining, by = "miRNA"
    ## Joining, by = "miRNA"
    ## Joining, by = "miRNA"
    ## Joining, by = "miRNA"

``` r
colnames(seqdata) <- str_replace_all(colnames(seqdata), pattern='-', replacement = '_')
#colnames(seqdata) <- str_replace_all(colnames(seqdata), pattern='__', replacement = '_')

seqdata <- seqdata[- grep("\\*", seqdata$miRNA), ]
seqdata <- seqdata %>% filter(str_detect(miRNA , "chr", negate = TRUE))
#colnames(seqdata)[2:length(colnames(seqdata))] <- sampleinfo$sample

# Format the data
countdata <- seqdata %>%
  column_to_rownames("miRNA") %>%
  #rename_all(str_remove, ".bam") %>%
  select(sampleinfo$filename) %>%
  as.matrix()

# List consensus samples
consensus <- c("neerincx", "fromm", "schee", "selitsky")

# create the design formula
sampleinfo$tissue.type <- as.factor(paste(sampleinfo$type, sampleinfo$tissue, sep="."))
sampleinfo$type <- as.factor(sampleinfo$type)
design <- as.formula(~ tissue.type)
```

# Differential Expression

``` r
# Make a named list of signature miRNA
dict_sig_mirna <- c()
res_dict <- list()
```

``` r
ref <- 'normal.colorect'
dds <- DeseqObject(design, countdata, sampleinfo, "None", "None", ref)
```

``` r
# #datasets in total
dim(dds[, colData(dds)$type.tissue == 'pCRC'])
```

    ## [1] 389 120

``` r
dim(dds[, colData(dds)$type.tissue == 'mLi'])
```

    ## [1] 389  35

``` r
dim(dds[, colData(dds)$type.tissue == 'mLu'])
```

    ## [1] 389  28

``` r
dim(dds[, colData(dds)$type.tissue == 'nCR'])
```

    ## [1] 389  25

``` r
dim(dds[, colData(dds)$type.tissue == 'nLi'])
```

    ## [1] 389  20

``` r
dim(dds[, colData(dds)$type.tissue == 'nLu'])
```

    ## [1] 389  10

``` r
dim(dds[, colData(dds)$type.tissue == 'PM'])
```

    ## [1] 389  30

``` r
# #datasets for Fromm
dim(dds[, colData(dds)$type.tissue == 'pCRC' & colData(dds)$paper == 'fromm'])
```

    ## [1] 389   3

``` r
dim(dds[, colData(dds)$type.tissue == 'mLi' & colData(dds)$paper == 'fromm'])
```

    ## [1] 389  19

``` r
dim(dds[, colData(dds)$type.tissue == 'mLu' & colData(dds)$paper == 'fromm'])
```

    ## [1] 389  24

``` r
dim(dds[, colData(dds)$type.tissue == 'nCR' & colData(dds)$paper == 'fromm'])
```

    ## [1] 389   3

``` r
dim(dds[, colData(dds)$type.tissue == 'nLi' & colData(dds)$paper == 'fromm'])
```

    ## [1] 389   8

``` r
dim(dds[, colData(dds)$type.tissue == 'nLu' & colData(dds)$paper == 'fromm'])
```

    ## [1] 389   7

``` r
dim(dds[, colData(dds)$type.tissue == 'PM' & colData(dds)$paper == 'fromm'])
```

    ## [1] 389  18

``` r
# #datasets for Schee
dim(dds[, colData(dds)$type.tissue == 'pCRC' & colData(dds)$paper == 'schee'])
```

    ## [1] 389  83

``` r
dim(dds[, colData(dds)$type.tissue == 'mLi' & colData(dds)$paper == 'schee'])
```

    ## [1] 389   0

``` r
dim(dds[, colData(dds)$type.tissue == 'mLu' & colData(dds)$paper == 'schee'])
```

    ## [1] 389   0

``` r
dim(dds[, colData(dds)$type.tissue == 'nCR' & colData(dds)$paper == 'schee'])
```

    ## [1] 389   0

``` r
dim(dds[, colData(dds)$type.tissue == 'nLi' & colData(dds)$paper == 'schee'])
```

    ## [1] 389   0

``` r
dim(dds[, colData(dds)$type.tissue == 'nLu' & colData(dds)$paper == 'schee'])
```

    ## [1] 389   0

``` r
dim(dds[, colData(dds)$type.tissue == 'PM' & colData(dds)$paper == 'schee'])
```

    ## [1] 389   0

``` r
# #datasets for Schee
dim(dds[, colData(dds)$type.tissue == 'pCRC' & colData(dds)$paper == 'neerincx'])
```

    ## [1] 389  34

``` r
dim(dds[, colData(dds)$type.tissue == 'mLi' & colData(dds)$paper == 'neerincx'])
```

    ## [1] 389  16

``` r
dim(dds[, colData(dds)$type.tissue == 'mLu' & colData(dds)$paper == 'neerincx'])
```

    ## [1] 389   4

``` r
dim(dds[, colData(dds)$type.tissue == 'nCR' & colData(dds)$paper == 'neerincx'])
```

    ## [1] 389  22

``` r
dim(dds[, colData(dds)$type.tissue == 'nLi' & colData(dds)$paper == 'neerincx'])
```

    ## [1] 389   9

``` r
dim(dds[, colData(dds)$type.tissue == 'nLu' & colData(dds)$paper == 'neerincx'])
```

    ## [1] 389   3

``` r
dim(dds[, colData(dds)$type.tissue == 'PM' & colData(dds)$paper == 'neerincx'])
```

    ## [1] 389  12

``` r
# #datasets for Schee
dim(dds[, colData(dds)$type.tissue == 'pCRC' & colData(dds)$paper == 'selitsky'])
```

    ## [1] 389   0

``` r
dim(dds[, colData(dds)$type.tissue == 'mLi' & colData(dds)$paper == 'selitsky'])
```

    ## [1] 389   0

``` r
dim(dds[, colData(dds)$type.tissue == 'mLu' & colData(dds)$paper == 'selitsky'])
```

    ## [1] 389   0

``` r
dim(dds[, colData(dds)$type.tissue == 'nCR' & colData(dds)$paper == 'selitsky'])
```

    ## [1] 389   0

``` r
dim(dds[, colData(dds)$type.tissue == 'nLi' & colData(dds)$paper == 'selitsky'])
```

    ## [1] 389   3

``` r
dim(dds[, colData(dds)$type.tissue == 'nLu' & colData(dds)$paper == 'selitsky'])
```

    ## [1] 389   0

``` r
dim(dds[, colData(dds)$type.tissue == 'PM' & colData(dds)$paper == 'selitsky'])
```

    ## [1] 389   0

``` r
## Plot dispersion estimates
plotDispEsts(dds)
```

![](comet_analysis_report_automated_no_lfcthreshold_08.12.20_res_alpha_0.05_with_LFC_shrinkage_LFC_0.58_mirge3_7.0_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

McCall, Matthew N; Kim, Min-Sik; Adil, Mohammed; Patil, Arun H; Lu, Yin;
Mitchell, Christopher J; Leal-Rojas, Pamela; Xu, Jinchong; Kumar, Manoj;
Dawson, Valina L; Dawson, Ted M; Baras, Alexander S; Rosenberg, Avi Z;
Arking, Dan E; Burns, Kathleen H; Pandey, Akhilesh; Halushka, Marc K
Toward the human cellular microRNAome Genome Res. October 2017

``` r
cell_spec_dict <- list(
  "CD14+ Monocyte"   = c("Hsa-Mir-15-P1a_5p","Hsa-Mir-15-P1b_5p", "Hsa-Mir-17-P1a_5p/P1b_5p"),
  
  "Dendritic Cell"   = c("Hsa-Mir-146-P2_5p", "Hsa-Mir-342_3p", "Hsa-Mir-142_3p",
                         "Hsa-Mir-223_3p"),
  
  "Endothelial Cell" = c("Hsa-Mir-126_5p"),
  
  "Epithelial Cell"  = c("Hsa-Mir-8-P2a_3p", "Hsa-Mir-8-P2b_3p", "Hsa0Mir-205-P1_5p",
                        "Hsa-Mir-192-P1_5p/P2_5p", "Hsa-Mir-375_3p"),
  
  "Islet Cell"       = c("Hsa-Mir-375_3p", "Hsa-Mir-154-P7_5p", "Hsa-Mir-7-P1_5p/P2_5p/P3_5p"),
  
  "Lymphocyte"       = c("Hsa-Mir-146-P2_5p", "Hsa-Mir-342_3p", "Hsa-Mir-150_5p",
                        "Hsa-Mir-155_5p"),
  
  "Macrophage"       = c("Hsa-Mir-342_3p", "Hsa-Mir-142_3p", "Hsa-Mir-223_3p", 
                         "Hsa-Mir-155_5p", "Hsa-Mir-24-P1_3p/P2_3p",
                         "Hsa-Mir-185_5p"),
  
  "Melanocyte"       = c("Hsa-Mir-185_5p", "Hsa-Mir-204-P2_5p"),
  
  "Mesenchymal"      = c("Hsa-Mir-185_5p", "Hsa-Mir-143_3p", "Hsa-Mir-145_5p"),
  
  "Neural"           = c("Hsa-Mir-375_3p", "Hsa-Mir-154-P7_5p", "Hsa-Mir-7-P1_5p/P2_5p/P3_5p",
                         "Hsa-Mir-128-P1_3p/P2_3p", "Hsa-Mir-129-P1_5p/P2_5p",
                         "Hsa-Mir-9-P1_5p/P2_5p/P3_5p","Hsa-Mir-430-P2_3p",
                         "Hsa-Mir-430-P4_3p"),
  
  "Platelet"         = c("Hsa-Mir-126_5p", "Hsa-Mir-486_5p"),
  
  "Red Blood Cell"   = c("Hsa-Mir-486_5p", "Hsa-Mir-451_5p", "Hsa-Mir-144_5p"),
  
  "Retinal Epithelial Cell" = c("Hsa-Mir-204-P1_5p", "Hsa-Mir-204-P2_5p", "Hsa-Mir-335_5p"),
  
  "Skeletal Myocyte" = c("Hsa-Mir-1-P1_3p/P2_3p", "Hsa-Mir-133-P1_3p/P2_3p/P3_3p"),
  
  "Stem Cell"        = c("Hsa-Mir-430-P2_3p", "Hsa-Mir-430-P4_3p", "Hsa-Mir-133-P1_3p/P2_3p/P3_3p"),
  
  "Hepatocyte"        = c("Hsa-Mir-122_5p")
  )

cell_spec_dict_inv <- topGO::inverseList(cell_spec_dict)
```

    ## 

Pritchard, C C; Kroh, E; Wood, B; Arroyo, J D; Dougherty, K J; Miyaji, M
M; Tait, J F; Tewari, M Blood Cell Origin of Circulating MicroRNAs: A
Cautionary Note for Cancer Biomarker Studies Cancer Prevention Research
2012

``` r
blood.cell.mirna <- c("Hsa-Mir-223_3p",
                      "Hsa-Mir-15-P2a_5p",
                      "Hsa-Mir-15-P2b_5p",
                      "Hsa-Mir-126-v1_3p",
                      "Hsa-Mir-142-v1_3p",
                      "Hsa-Mir-21_5p",
                      "Hsa-Mir-24-P1_3p",
                      "Hsa-Mir-24-P2_3p",
                      "Hsa-Mir-19-P2a_3p",
                      "Hsa-Mir-19-P2b_3p",
                      "Hsa-Mir-103-P1_3p",
                      "Hsa-Mir-103-P2_3p",
                      "Hsa-Let-7-P1a_5p",
                      "Hsa-Let-7-P2a1_5p",
                      "Hsa-Let-7-P2a2_5p",
                      "Hsa-Mir-451_5p",
                      "Hsa-Mir-92-P1a_3p",
                      "Hsa-Mir-92-P1b_3p",
                      "Hsa-Mir-17-P1b_5p",
                      "Hsa-Mir-19-P1_3p",
                      "Hsa-Mir-30-P2c_5p",
                      "Hsa-Mir-17-P1a",
                      "Hsa-Mir-15-P1b_5p",
                      "Hsa-Mir-103-P3_3p",
                      "Hsa-Let-7-P2a3_5p",
                      "Hsa-Let-7-P2b1_5p",
                      "Hsa-Mir-221-P1_3p",
                      "Hsa-Mir-221-P2_3p",
                      "Hsa-Mir-17-P1c_5p",
                      "Hsa-Mir-30-P2a_5p",
                      "Hsa-Mir-30-P2b_5p",
                      "Hsa-Mir-28-P2_5p",
                      "Hsa-Mir-30-P1b_5p",
                      "Hsa-Mir-30-P1c_5p",
                      "Hsa-Mir-486_5p",
                      "Hsa-Mir-92-P2c_3p",
                      "Hsa-Mir-181-P1a_5p",
                      "Hsa-Mir-181-P1b_5p",
                      "Hsa-Mir-146-P1_5p",
                      "Hsa-Let-7-P2c1_5p",
                      "Hsa-Mir-197_3p",
                      "Hsa-Mir-17-P3c_5p",
                      "Hsa-Mir-17-P3c_3p",
                      "Hsa-Mir-148-P3_3p",
                      "Hsa-Mir-766_3p",
                      "Hsa-Mir-17-P3b_5p",
                      "Hsa-Mir-328_3p",
                      "Hsa-Mir-574_3p",
                      "Hsa-Mir-155_5p",
                      "Hsa-Mir-425_5p",
                      "Hsa-Mir-148-P1_3p",
                      "Hsa-Mir-29-P1a_3p",
                      "Hsa-Mir-8-P2b_3p",
                      "Hsa-Mir-92-P1c_3p",
                      "Hsa-Mir-192-P2_5p",
                      "Hsa-Mir-362-P2-v1_3p",
                      "Hsa-Mir-362-P5_5p"
                      )
```

## nCR vs nLi

``` r
column='tissue.type'
tissue_type_A <- 'normal.liver'
tissue_type_B <- 'normal.colorect'
norm_adj_up   = "None"
norm_adj_down = "None"
pCRC_adj_up   = "None"
pCRC_adj_down = "None"

coef <- paste(column, tissue_type_A, 'vs', tissue_type_B, sep='_')
res <- DeseqResult(dds, column, coef, tissue_type_A, tissue_type_B,
                   lfc.Threshold, rpm.Threshold,
                   norm_adj_up,
                   norm_adj_down)
dict_sig_mirna[paste(coef, "up",   sep='_')] <- list(res$up_mirna)
dict_sig_mirna[paste(coef, "down", sep='_')] <- list(res$down_mirna)
res_res <- res$res
res_dict[coef] <- res_res
plotMA(res$res, alpha=0.05)
```

![](comet_analysis_report_automated_no_lfcthreshold_08.12.20_res_alpha_0.05_with_LFC_shrinkage_LFC_0.58_mirge3_7.0_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

``` r
# Plot volcano plot
VolcanoPlot(res$res, coef, res$sig,
            res$up_mirna, res$down_mirna,
            norm_adj_up, norm_adj_down,
            pCRC_adj_up, pCRC_adj_down)
```

![](comet_analysis_report_automated_no_lfcthreshold_08.12.20_res_alpha_0.05_with_LFC_shrinkage_LFC_0.58_mirge3_7.0_files/figure-gfm/unnamed-chunk-15-2.png)<!-- -->

``` r
ExpressionPlot(res$res, res$rpm, coef, res$sig,
               tissue_type_A, tissue_type_B,
               res$up_mirna, res$down_mirna,
               norm_adj_up, norm_adj_down,
               pCRC_adj_up, pCRC_adj_down)
```

![](comet_analysis_report_automated_no_lfcthreshold_08.12.20_res_alpha_0.05_with_LFC_shrinkage_LFC_0.58_mirge3_7.0_files/figure-gfm/unnamed-chunk-15-3.png)<!-- -->

``` r
signature_mirnas <- SigList(res, dds, tissue_type_A, tissue_type_B, coef,
                            norm_adj_up, norm_adj_down, 
                            pCRC_adj_up, pCRC_adj_down)
# Print list upregulated miRNA
signature_mirnas$up_mirna
```

<div style="border: 1px solid #ddd; padding: 5px; overflow-x: scroll; width:2000px; ">

<table class="table table-striped table-hover table-condensed" style="margin-left: auto; margin-right: auto;">

<caption>

Upregulated in tissue.type\_normal.liver\_vs\_normal.colorect

</caption>

<thead>

<tr>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

miRNA

</th>

<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">

LFC

</th>

<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">

lfcSE

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

FDR

</th>

<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">

RPM normal.liver

</th>

<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">

RPM normal.colorect

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

miRBase\_ID

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

Family

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

Seed

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

Chr

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

Cell-Type Specific

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

Norm Background

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

pCRC Background

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

Hsa-Mir-204-P1\_5p

</td>

<td style="text-align:right;">

1.56

</td>

<td style="text-align:right;">

0.51

</td>

<td style="text-align:left;">

4.46e-03

</td>

<td style="text-align:right;">

190

</td>

<td style="text-align:right;">

43

</td>

<td style="text-align:left;">

hsa-mir-204

</td>

<td style="text-align:left;">

MIR-204

</td>

<td style="text-align:left;">

UCCCUUU

</td>

<td style="text-align:left;">

chr9

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">Retinal
Epithelial Cell</span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-335\_5p

</td>

<td style="text-align:right;">

1.18

</td>

<td style="text-align:right;">

0.17

</td>

<td style="text-align:left;">

2.56e-11

</td>

<td style="text-align:right;">

205

</td>

<td style="text-align:right;">

79

</td>

<td style="text-align:left;">

hsa-mir-335

</td>

<td style="text-align:left;">

MIR-335

</td>

<td style="text-align:left;">

CAAGAGC

</td>

<td style="text-align:left;">

chr7

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">Retinal
Epithelial Cell</span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-144\_5p

</td>

<td style="text-align:right;">

1.53

</td>

<td style="text-align:right;">

0.33

</td>

<td style="text-align:left;">

5.21e-05

</td>

<td style="text-align:right;">

206

</td>

<td style="text-align:right;">

65

</td>

<td style="text-align:left;">

hsa-mir-144

</td>

<td style="text-align:left;">

MIR-144

</td>

<td style="text-align:left;">

GAUAUCA

</td>

<td style="text-align:left;">

chr17

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">Red
Blood Cell</span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-128-P1\_3p/P2\_3p

</td>

<td style="text-align:right;">

0.65

</td>

<td style="text-align:right;">

0.16

</td>

<td style="text-align:left;">

2.62e-04

</td>

<td style="text-align:right;">

194

</td>

<td style="text-align:right;">

110

</td>

<td style="text-align:left;">

hsa-mir-128-1

</td>

<td style="text-align:left;">

MIR-128

</td>

<td style="text-align:left;">

CACAGUG

</td>

<td style="text-align:left;">

chr2

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">Neural</span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-122\_5p

</td>

<td style="text-align:right;">

10.58

</td>

<td style="text-align:right;">

0.78

</td>

<td style="text-align:left;">

7.97e-28

</td>

<td style="text-align:right;">

148842

</td>

<td style="text-align:right;">

26

</td>

<td style="text-align:left;">

hsa-mir-122

</td>

<td style="text-align:left;">

MIR-122

</td>

<td style="text-align:left;">

GGAGUGU

</td>

<td style="text-align:left;">

chr18

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">Hepatocyte</span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-486\_5p

</td>

<td style="text-align:right;">

1.34

</td>

<td style="text-align:right;">

0.36

</td>

<td style="text-align:left;">

2.99e-03

</td>

<td style="text-align:right;">

10617

</td>

<td style="text-align:right;">

3961

</td>

<td style="text-align:left;">

hsa-mir-486-1

</td>

<td style="text-align:left;">

MIR-486

</td>

<td style="text-align:left;">

CCUGUAC

</td>

<td style="text-align:left;">

chr8

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">c(“Platelet”,
“Red Blood Cell”)</span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-126\_5p

</td>

<td style="text-align:right;">

1.25

</td>

<td style="text-align:right;">

0.21

</td>

<td style="text-align:left;">

1.87e-08

</td>

<td style="text-align:right;">

11148

</td>

<td style="text-align:right;">

4230

</td>

<td style="text-align:left;">

hsa-mir-126

</td>

<td style="text-align:left;">

MIR-126

</td>

<td style="text-align:left;">

AUUAUUA

</td>

<td style="text-align:left;">

chr9

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">c(“Endothelial
Cell”, “Platelet”)</span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-885\_5p

</td>

<td style="text-align:right;">

8.77

</td>

<td style="text-align:right;">

0.70

</td>

<td style="text-align:left;">

8.60e-24

</td>

<td style="text-align:right;">

510

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:left;">

hsa-mir-885

</td>

<td style="text-align:left;">

MIR-885

</td>

<td style="text-align:left;">

CCAUUAC

</td>

<td style="text-align:left;">

chr3

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-483\_5p

</td>

<td style="text-align:right;">

4.03

</td>

<td style="text-align:right;">

0.71

</td>

<td style="text-align:left;">

1.04e-10

</td>

<td style="text-align:right;">

115

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

hsa-mir-483

</td>

<td style="text-align:left;">

MIR-483

</td>

<td style="text-align:left;">

AGACGGG

</td>

<td style="text-align:left;">

chr11

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Let-7-P1c\_5p

</td>

<td style="text-align:right;">

3.18

</td>

<td style="text-align:right;">

0.28

</td>

<td style="text-align:left;">

1.60e-27

</td>

<td style="text-align:right;">

4945

</td>

<td style="text-align:right;">

466

</td>

<td style="text-align:left;">

hsa-let-7c

</td>

<td style="text-align:left;">

LET-7

</td>

<td style="text-align:left;">

GAGGUAG

</td>

<td style="text-align:left;">

chr21

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-10-P2c\_5p

</td>

<td style="text-align:right;">

3.08

</td>

<td style="text-align:right;">

0.33

</td>

<td style="text-align:left;">

6.34e-19

</td>

<td style="text-align:right;">

1202

</td>

<td style="text-align:right;">

114

</td>

<td style="text-align:left;">

hsa-mir-99a

</td>

<td style="text-align:left;">

MIR-10

</td>

<td style="text-align:left;">

ACCCGUA

</td>

<td style="text-align:left;">

chr21

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-455\_5p

</td>

<td style="text-align:right;">

2.48

</td>

<td style="text-align:right;">

0.18

</td>

<td style="text-align:left;">

9.25e-41

</td>

<td style="text-align:right;">

279

</td>

<td style="text-align:right;">

44

</td>

<td style="text-align:left;">

hsa-mir-455

</td>

<td style="text-align:left;">

MIR-455

</td>

<td style="text-align:left;">

AUGUGCC

</td>

<td style="text-align:left;">

chr9

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-193-P2a\_3p/P2b\_3p

</td>

<td style="text-align:right;">

2.37

</td>

<td style="text-align:right;">

0.23

</td>

<td style="text-align:left;">

5.29e-22

</td>

<td style="text-align:right;">

219

</td>

<td style="text-align:right;">

37

</td>

<td style="text-align:left;">

hsa-mir-365b

</td>

<td style="text-align:left;">

MIR-193

</td>

<td style="text-align:left;">

AAUGCCC

</td>

<td style="text-align:left;">

chr17

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-193-P1b\_3p

</td>

<td style="text-align:right;">

2.34

</td>

<td style="text-align:right;">

0.23

</td>

<td style="text-align:left;">

2.49e-22

</td>

<td style="text-align:right;">

514

</td>

<td style="text-align:right;">

90

</td>

<td style="text-align:left;">

hsa-mir-193b

</td>

<td style="text-align:left;">

MIR-193

</td>

<td style="text-align:left;">

ACUGGCC

</td>

<td style="text-align:left;">

chr16

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-15-P1d\_5p

</td>

<td style="text-align:right;">

2.26

</td>

<td style="text-align:right;">

0.30

</td>

<td style="text-align:left;">

5.34e-14

</td>

<td style="text-align:right;">

247

</td>

<td style="text-align:right;">

40

</td>

<td style="text-align:left;">

hsa-mir-424

</td>

<td style="text-align:left;">

MIR-15

</td>

<td style="text-align:left;">

AGCAGCA

</td>

<td style="text-align:left;">

chrX

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-193-P1a\_5p

</td>

<td style="text-align:right;">

2.14

</td>

<td style="text-align:right;">

0.28

</td>

<td style="text-align:left;">

3.66e-12

</td>

<td style="text-align:right;">

298

</td>

<td style="text-align:right;">

57

</td>

<td style="text-align:left;">

hsa-mir-193a

</td>

<td style="text-align:left;">

MIR-193

</td>

<td style="text-align:left;">

GGGUCUU

</td>

<td style="text-align:left;">

chr17

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-10-P3c\_5p

</td>

<td style="text-align:right;">

2.12

</td>

<td style="text-align:right;">

0.30

</td>

<td style="text-align:left;">

6.48e-11

</td>

<td style="text-align:right;">

1620

</td>

<td style="text-align:right;">

326

</td>

<td style="text-align:left;">

hsa-mir-125b-2

</td>

<td style="text-align:left;">

MIR-10

</td>

<td style="text-align:left;">

CCCUGAG

</td>

<td style="text-align:left;">

chr21

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-139\_5p

</td>

<td style="text-align:right;">

2.00

</td>

<td style="text-align:right;">

0.26

</td>

<td style="text-align:left;">

1.06e-11

</td>

<td style="text-align:right;">

174

</td>

<td style="text-align:right;">

41

</td>

<td style="text-align:left;">

hsa-mir-139

</td>

<td style="text-align:left;">

MIR-139

</td>

<td style="text-align:left;">

CUACAGU

</td>

<td style="text-align:left;">

chr11

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-148-P1\_3p

</td>

<td style="text-align:right;">

2.00

</td>

<td style="text-align:right;">

0.24

</td>

<td style="text-align:left;">

8.34e-15

</td>

<td style="text-align:right;">

105086

</td>

<td style="text-align:right;">

21730

</td>

<td style="text-align:left;">

hsa-mir-148a

</td>

<td style="text-align:left;">

MIR-148

</td>

<td style="text-align:left;">

CAGUGCA

</td>

<td style="text-align:left;">

chr7

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-10-P3a\_5p

</td>

<td style="text-align:right;">

1.93

</td>

<td style="text-align:right;">

0.31

</td>

<td style="text-align:left;">

1.24e-08

</td>

<td style="text-align:right;">

611

</td>

<td style="text-align:right;">

137

</td>

<td style="text-align:left;">

hsa-mir-125b-1

</td>

<td style="text-align:left;">

MIR-10

</td>

<td style="text-align:left;">

CCCUGAG

</td>

<td style="text-align:left;">

chr11

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-423\_5p

</td>

<td style="text-align:right;">

1.62

</td>

<td style="text-align:right;">

0.24

</td>

<td style="text-align:left;">

1.23e-09

</td>

<td style="text-align:right;">

1198

</td>

<td style="text-align:right;">

345

</td>

<td style="text-align:left;">

hsa-mir-423

</td>

<td style="text-align:left;">

MIR-423

</td>

<td style="text-align:left;">

GAGGGGC

</td>

<td style="text-align:left;">

chr17

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-92-P1a\_3p/P1b\_3p

</td>

<td style="text-align:right;">

1.60

</td>

<td style="text-align:right;">

0.22

</td>

<td style="text-align:left;">

1.98e-12

</td>

<td style="text-align:right;">

51173

</td>

<td style="text-align:right;">

13849

</td>

<td style="text-align:left;">

hsa-mir-92a-1

</td>

<td style="text-align:left;">

MIR-92

</td>

<td style="text-align:left;">

AUUGCAC

</td>

<td style="text-align:left;">

chr13

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-101-P1\_3p/P2\_3p

</td>

<td style="text-align:right;">

1.58

</td>

<td style="text-align:right;">

0.18

</td>

<td style="text-align:left;">

3.44e-17

</td>

<td style="text-align:right;">

14649

</td>

<td style="text-align:right;">

4197

</td>

<td style="text-align:left;">

hsa-mir-101-1

</td>

<td style="text-align:left;">

MIR-101

</td>

<td style="text-align:left;">

UACAGUA

</td>

<td style="text-align:left;">

chr1

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-22-P1a\_3p

</td>

<td style="text-align:right;">

1.58

</td>

<td style="text-align:right;">

0.17

</td>

<td style="text-align:left;">

1.76e-18

</td>

<td style="text-align:right;">

78932

</td>

<td style="text-align:right;">

23453

</td>

<td style="text-align:left;">

hsa-mir-22

</td>

<td style="text-align:left;">

MIR-22

</td>

<td style="text-align:left;">

AGCUGCC

</td>

<td style="text-align:left;">

chr17

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-340\_5p

</td>

<td style="text-align:right;">

1.42

</td>

<td style="text-align:right;">

0.16

</td>

<td style="text-align:left;">

3.45e-18

</td>

<td style="text-align:right;">

1030

</td>

<td style="text-align:right;">

337

</td>

<td style="text-align:left;">

hsa-mir-340

</td>

<td style="text-align:left;">

MIR-340

</td>

<td style="text-align:left;">

UAUAAAG

</td>

<td style="text-align:left;">

chr5

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-30-P1a\_5p

</td>

<td style="text-align:right;">

1.32

</td>

<td style="text-align:right;">

0.23

</td>

<td style="text-align:left;">

1.07e-07

</td>

<td style="text-align:right;">

15175

</td>

<td style="text-align:right;">

5150

</td>

<td style="text-align:left;">

hsa-mir-30a

</td>

<td style="text-align:left;">

MIR-30

</td>

<td style="text-align:left;">

GUAAACA

</td>

<td style="text-align:left;">

chr6

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-744\_5p

</td>

<td style="text-align:right;">

1.31

</td>

<td style="text-align:right;">

0.21

</td>

<td style="text-align:left;">

1.28e-08

</td>

<td style="text-align:right;">

165

</td>

<td style="text-align:right;">

59

</td>

<td style="text-align:left;">

hsa-mir-744

</td>

<td style="text-align:left;">

MIR-744

</td>

<td style="text-align:left;">

GCGGGGC

</td>

<td style="text-align:left;">

chr17

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-574\_3p

</td>

<td style="text-align:right;">

1.21

</td>

<td style="text-align:right;">

0.19

</td>

<td style="text-align:left;">

1.27e-08

</td>

<td style="text-align:right;">

745

</td>

<td style="text-align:right;">

305

</td>

<td style="text-align:left;">

hsa-mir-574

</td>

<td style="text-align:left;">

MIR-574

</td>

<td style="text-align:left;">

ACGCUCA

</td>

<td style="text-align:left;">

chr4

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-130-P1a\_3p

</td>

<td style="text-align:right;">

1.09

</td>

<td style="text-align:right;">

0.20

</td>

<td style="text-align:left;">

4.15e-07

</td>

<td style="text-align:right;">

907

</td>

<td style="text-align:right;">

379

</td>

<td style="text-align:left;">

hsa-mir-130a

</td>

<td style="text-align:left;">

MIR-130

</td>

<td style="text-align:left;">

AGUGCAA

</td>

<td style="text-align:left;">

chr11

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-10-P2a\_5p

</td>

<td style="text-align:right;">

1.07

</td>

<td style="text-align:right;">

0.32

</td>

<td style="text-align:left;">

3.39e-03

</td>

<td style="text-align:right;">

2634

</td>

<td style="text-align:right;">

986

</td>

<td style="text-align:left;">

hsa-mir-100

</td>

<td style="text-align:left;">

MIR-10

</td>

<td style="text-align:left;">

ACCCGUA

</td>

<td style="text-align:left;">

chr11

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-30-P2a\_5p/P2b\_5p/P2c\_5p

</td>

<td style="text-align:right;">

0.99

</td>

<td style="text-align:right;">

0.15

</td>

<td style="text-align:left;">

4.69e-10

</td>

<td style="text-align:right;">

5891

</td>

<td style="text-align:right;">

2677

</td>

<td style="text-align:left;">

hsa-mir-30c-2

</td>

<td style="text-align:left;">

MIR-30

</td>

<td style="text-align:left;">

GUAAACA

</td>

<td style="text-align:left;">

chr6

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-154-P23\_3p

</td>

<td style="text-align:right;">

0.89

</td>

<td style="text-align:right;">

0.20

</td>

<td style="text-align:left;">

5.21e-05

</td>

<td style="text-align:right;">

222

</td>

<td style="text-align:right;">

108

</td>

<td style="text-align:left;">

hsa-mir-654

</td>

<td style="text-align:left;">

MIR-154

</td>

<td style="text-align:left;">

AUGUCUG

</td>

<td style="text-align:left;">

chr14

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-197\_3p

</td>

<td style="text-align:right;">

0.79

</td>

<td style="text-align:right;">

0.19

</td>

<td style="text-align:left;">

2.06e-04

</td>

<td style="text-align:right;">

344

</td>

<td style="text-align:right;">

180

</td>

<td style="text-align:left;">

hsa-mir-197

</td>

<td style="text-align:left;">

MIR-197

</td>

<td style="text-align:left;">

UCACCAC

</td>

<td style="text-align:left;">

chr1

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-214\_3p

</td>

<td style="text-align:right;">

0.73

</td>

<td style="text-align:right;">

0.23

</td>

<td style="text-align:left;">

5.91e-03

</td>

<td style="text-align:right;">

137

</td>

<td style="text-align:right;">

74

</td>

<td style="text-align:left;">

hsa-mir-214

</td>

<td style="text-align:left;">

MIR-214

</td>

<td style="text-align:left;">

CAGCAGG

</td>

<td style="text-align:left;">

chr1

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-30-P1b\_5p

</td>

<td style="text-align:right;">

0.71

</td>

<td style="text-align:right;">

0.15

</td>

<td style="text-align:left;">

7.77e-06

</td>

<td style="text-align:right;">

11240

</td>

<td style="text-align:right;">

6046

</td>

<td style="text-align:left;">

hsa-mir-30e

</td>

<td style="text-align:left;">

MIR-30

</td>

<td style="text-align:left;">

GUAAACA

</td>

<td style="text-align:left;">

chr1

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-30-P1c\_5p

</td>

<td style="text-align:right;">

0.70

</td>

<td style="text-align:right;">

0.16

</td>

<td style="text-align:left;">

7.79e-05

</td>

<td style="text-align:right;">

13680

</td>

<td style="text-align:right;">

7243

</td>

<td style="text-align:left;">

hsa-mir-30d

</td>

<td style="text-align:left;">

MIR-30

</td>

<td style="text-align:left;">

GUAAACA

</td>

<td style="text-align:left;">

chr8

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

</tbody>

</table>

</div>

``` r
# Number of upregulated miRNA
signature_mirnas$number_upregulated
```

    ## [1] 36

``` r
# Print list downregulated miRNA
signature_mirnas$down_mirna
```

<div style="border: 1px solid #ddd; padding: 5px; overflow-x: scroll; width:2000px; ">

<table class="table table-striped table-hover table-condensed" style="margin-left: auto; margin-right: auto;">

<caption>

Downregulated in tissue.type\_normal.liver\_vs\_normal.colorect

</caption>

<thead>

<tr>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

miRNA

</th>

<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">

LFC

</th>

<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">

lfcSE

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

FDR

</th>

<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">

RPM normal.liver

</th>

<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">

RPM normal.colorect

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

miRBase\_ID

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

Family

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

Seed

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

Chr

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

Cell-Type Specific

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

Norm Background

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

pCRC Background

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

Hsa-Mir-145\_5p

</td>

<td style="text-align:right;">

\-2.51

</td>

<td style="text-align:right;">

0.32

</td>

<td style="text-align:left;">

5.55e-15

</td>

<td style="text-align:right;">

484

</td>

<td style="text-align:right;">

2832

</td>

<td style="text-align:left;">

hsa-mir-145

</td>

<td style="text-align:left;">

MIR-145

</td>

<td style="text-align:left;">

UCCAGUU

</td>

<td style="text-align:left;">

chr5

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">Mesenchymal</span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-143\_3p

</td>

<td style="text-align:right;">

\-2.35

</td>

<td style="text-align:right;">

0.27

</td>

<td style="text-align:left;">

1.80e-17

</td>

<td style="text-align:right;">

37572

</td>

<td style="text-align:right;">

164812

</td>

<td style="text-align:left;">

hsa-mir-143

</td>

<td style="text-align:left;">

MIR-143

</td>

<td style="text-align:left;">

GAGAUGA

</td>

<td style="text-align:left;">

chr5

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">Mesenchymal</span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-24-P1\_3p/P2\_3p

</td>

<td style="text-align:right;">

\-0.63

</td>

<td style="text-align:right;">

0.20

</td>

<td style="text-align:left;">

5.98e-03

</td>

<td style="text-align:right;">

319

</td>

<td style="text-align:right;">

413

</td>

<td style="text-align:left;">

hsa-mir-24-2

</td>

<td style="text-align:left;">

MIR-24

</td>

<td style="text-align:left;">

GGCUCAG

</td>

<td style="text-align:left;">

chr19

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">Macrophage</span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-8-P2b\_3p

</td>

<td style="text-align:right;">

\-6.38

</td>

<td style="text-align:right;">

0.23

</td>

<td style="text-align:left;">

7.72e-155

</td>

<td style="text-align:right;">

31

</td>

<td style="text-align:right;">

2444

</td>

<td style="text-align:left;">

hsa-mir-200c

</td>

<td style="text-align:left;">

MIR-8

</td>

<td style="text-align:left;">

AAUACUG

</td>

<td style="text-align:left;">

chr12

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">Epithelial
Cell</span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-8-P2a\_3p

</td>

<td style="text-align:right;">

\-4.47

</td>

<td style="text-align:right;">

0.24

</td>

<td style="text-align:left;">

4.08e-76

</td>

<td style="text-align:right;">

369

</td>

<td style="text-align:right;">

7741

</td>

<td style="text-align:left;">

hsa-mir-200b

</td>

<td style="text-align:left;">

MIR-8

</td>

<td style="text-align:left;">

AAUACUG

</td>

<td style="text-align:left;">

chr1

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">Epithelial
Cell</span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-192-P1\_5p/P2\_5p

</td>

<td style="text-align:right;">

\-0.72

</td>

<td style="text-align:right;">

0.32

</td>

<td style="text-align:left;">

1.03e-02

</td>

<td style="text-align:right;">

29546

</td>

<td style="text-align:right;">

46010

</td>

<td style="text-align:left;">

hsa-mir-192

</td>

<td style="text-align:left;">

MIR-192

</td>

<td style="text-align:left;">

UGACCUA

</td>

<td style="text-align:left;">

chr11

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">Epithelial
Cell</span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-15-P1b\_5p

</td>

<td style="text-align:right;">

\-0.68

</td>

<td style="text-align:right;">

0.16

</td>

<td style="text-align:left;">

1.18e-04

</td>

<td style="text-align:right;">

114

</td>

<td style="text-align:right;">

168

</td>

<td style="text-align:left;">

hsa-mir-15b

</td>

<td style="text-align:left;">

MIR-15

</td>

<td style="text-align:left;">

AGCAGCA

</td>

<td style="text-align:left;">

chr3

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">CD14+
Monocyte</span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-133-P1\_3p/P2\_3p/P3\_3p

</td>

<td style="text-align:right;">

\-4.10

</td>

<td style="text-align:right;">

0.38

</td>

<td style="text-align:left;">

9.37e-29

</td>

<td style="text-align:right;">

22

</td>

<td style="text-align:right;">

460

</td>

<td style="text-align:left;">

hsa-mir-133a-2

</td>

<td style="text-align:left;">

MIR-133

</td>

<td style="text-align:left;">

UUGGUCC

</td>

<td style="text-align:left;">

chr20

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">c(“Skeletal
Myocyte”, “Stem Cell”)</span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-155\_5p

</td>

<td style="text-align:right;">

\-1.51

</td>

<td style="text-align:right;">

0.26

</td>

<td style="text-align:left;">

8.41e-08

</td>

<td style="text-align:right;">

113

</td>

<td style="text-align:right;">

295

</td>

<td style="text-align:left;">

hsa-mir-155

</td>

<td style="text-align:left;">

MIR-155

</td>

<td style="text-align:left;">

UAAUGCU

</td>

<td style="text-align:left;">

chr21

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">c(“Lymphocyte”,
“Macrophage”)</span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-375\_3p

</td>

<td style="text-align:right;">

\-2.65

</td>

<td style="text-align:right;">

0.37

</td>

<td style="text-align:left;">

1.69e-12

</td>

<td style="text-align:right;">

2553

</td>

<td style="text-align:right;">

17215

</td>

<td style="text-align:left;">

hsa-mir-375

</td>

<td style="text-align:left;">

MIR-375

</td>

<td style="text-align:left;">

UUGUUCG

</td>

<td style="text-align:left;">

chr2

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">c(“Epithelial
Cell”, “Islet Cell”, “Neural”)</span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-196-P1\_5p/P2\_5p

</td>

<td style="text-align:right;">

\-6.76

</td>

<td style="text-align:right;">

0.34

</td>

<td style="text-align:left;">

3.97e-79

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

250

</td>

<td style="text-align:left;">

hsa-mir-196a-1

</td>

<td style="text-align:left;">

MIR-196

</td>

<td style="text-align:left;">

AGGUAGU

</td>

<td style="text-align:left;">

chr17

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-147\_3p

</td>

<td style="text-align:right;">

\-6.57

</td>

<td style="text-align:right;">

0.38

</td>

<td style="text-align:left;">

2.57e-64

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

110

</td>

<td style="text-align:left;">

hsa-mir-147b

</td>

<td style="text-align:left;">

MIR-147

</td>

<td style="text-align:left;">

UGUGCGG

</td>

<td style="text-align:left;">

chr15

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-577\_5p

</td>

<td style="text-align:right;">

\-6.02

</td>

<td style="text-align:right;">

0.34

</td>

<td style="text-align:left;">

1.65e-65

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

138

</td>

<td style="text-align:left;">

hsa-mir-577

</td>

<td style="text-align:left;">

MIR-577

</td>

<td style="text-align:left;">

UAGAUAA

</td>

<td style="text-align:left;">

chr4

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-8-P1b\_3p

</td>

<td style="text-align:right;">

\-6.00

</td>

<td style="text-align:right;">

0.30

</td>

<td style="text-align:left;">

3.97e-79

</td>

<td style="text-align:right;">

154

</td>

<td style="text-align:right;">

9180

</td>

<td style="text-align:left;">

hsa-mir-141

</td>

<td style="text-align:left;">

MIR-8

</td>

<td style="text-align:left;">

AACACUG

</td>

<td style="text-align:left;">

chr12

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-196-P3\_5p

</td>

<td style="text-align:right;">

\-5.53

</td>

<td style="text-align:right;">

0.37

</td>

<td style="text-align:left;">

2.24e-42

</td>

<td style="text-align:right;">

9

</td>

<td style="text-align:right;">

400

</td>

<td style="text-align:left;">

hsa-mir-196b

</td>

<td style="text-align:left;">

MIR-196

</td>

<td style="text-align:left;">

AGGUAGU

</td>

<td style="text-align:left;">

chr7

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-10-P1b\_5p

</td>

<td style="text-align:right;">

\-4.80

</td>

<td style="text-align:right;">

0.36

</td>

<td style="text-align:left;">

6.31e-38

</td>

<td style="text-align:right;">

2704

</td>

<td style="text-align:right;">

75067

</td>

<td style="text-align:left;">

hsa-mir-10b

</td>

<td style="text-align:left;">

MIR-10

</td>

<td style="text-align:left;">

ACCCUGU

</td>

<td style="text-align:left;">

chr2

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-8-P3a\_3p

</td>

<td style="text-align:right;">

\-4.32

</td>

<td style="text-align:right;">

0.24

</td>

<td style="text-align:left;">

2.00e-68

</td>

<td style="text-align:right;">

68

</td>

<td style="text-align:right;">

1282

</td>

<td style="text-align:left;">

hsa-mir-429

</td>

<td style="text-align:left;">

MIR-8

</td>

<td style="text-align:left;">

AAUACUG

</td>

<td style="text-align:left;">

chr1

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-8-P1a\_3p

</td>

<td style="text-align:right;">

\-4.10

</td>

<td style="text-align:right;">

0.24

</td>

<td style="text-align:left;">

4.38e-64

</td>

<td style="text-align:right;">

108

</td>

<td style="text-align:right;">

1754

</td>

<td style="text-align:left;">

hsa-mir-200a

</td>

<td style="text-align:left;">

MIR-8

</td>

<td style="text-align:left;">

AACACUG

</td>

<td style="text-align:left;">

chr1

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-190-P1\_5p

</td>

<td style="text-align:right;">

\-3.71

</td>

<td style="text-align:right;">

0.26

</td>

<td style="text-align:left;">

4.48e-45

</td>

<td style="text-align:right;">

23

</td>

<td style="text-align:right;">

294

</td>

<td style="text-align:left;">

hsa-mir-190a

</td>

<td style="text-align:left;">

MIR-190

</td>

<td style="text-align:left;">

GAUAUGU

</td>

<td style="text-align:left;">

chr15

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-96-P3\_5p

</td>

<td style="text-align:right;">

\-3.57

</td>

<td style="text-align:right;">

0.32

</td>

<td style="text-align:left;">

1.72e-22

</td>

<td style="text-align:right;">

27

</td>

<td style="text-align:right;">

224

</td>

<td style="text-align:left;">

hsa-mir-183

</td>

<td style="text-align:left;">

MIR-96

</td>

<td style="text-align:left;">

AUGGCAC

</td>

<td style="text-align:left;">

chr7

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-203\_3p

</td>

<td style="text-align:right;">

\-2.68

</td>

<td style="text-align:right;">

0.28

</td>

<td style="text-align:left;">

6.56e-18

</td>

<td style="text-align:right;">

171

</td>

<td style="text-align:right;">

922

</td>

<td style="text-align:left;">

hsa-mir-203a

</td>

<td style="text-align:left;">

MIR-203

</td>

<td style="text-align:left;">

UGAAAUG

</td>

<td style="text-align:left;">

chr14

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-96-P2\_5p

</td>

<td style="text-align:right;">

\-2.65

</td>

<td style="text-align:right;">

0.24

</td>

<td style="text-align:left;">

2.99e-23

</td>

<td style="text-align:right;">

370

</td>

<td style="text-align:right;">

1867

</td>

<td style="text-align:left;">

hsa-mir-182

</td>

<td style="text-align:left;">

MIR-96

</td>

<td style="text-align:left;">

UUGGCAA

</td>

<td style="text-align:left;">

chr7

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-221-P1\_3p

</td>

<td style="text-align:right;">

\-2.30

</td>

<td style="text-align:right;">

0.19

</td>

<td style="text-align:left;">

5.81e-32

</td>

<td style="text-align:right;">

217

</td>

<td style="text-align:right;">

964

</td>

<td style="text-align:left;">

hsa-mir-221

</td>

<td style="text-align:left;">

MIR-221

</td>

<td style="text-align:left;">

GCUACAU

</td>

<td style="text-align:left;">

chrX

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-338-P1\_3p

</td>

<td style="text-align:right;">

\-1.89

</td>

<td style="text-align:right;">

0.31

</td>

<td style="text-align:left;">

1.21e-08

</td>

<td style="text-align:right;">

48

</td>

<td style="text-align:right;">

166

</td>

<td style="text-align:left;">

hsa-mir-338

</td>

<td style="text-align:left;">

MIR-338

</td>

<td style="text-align:left;">

CCAGCAU

</td>

<td style="text-align:left;">

chr17

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-221-P2\_3p

</td>

<td style="text-align:right;">

\-1.87

</td>

<td style="text-align:right;">

0.20

</td>

<td style="text-align:left;">

9.25e-18

</td>

<td style="text-align:right;">

173

</td>

<td style="text-align:right;">

574

</td>

<td style="text-align:left;">

hsa-mir-222

</td>

<td style="text-align:left;">

MIR-221

</td>

<td style="text-align:left;">

GCUACAU

</td>

<td style="text-align:left;">

chrX

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-10-P1a\_5p

</td>

<td style="text-align:right;">

\-1.77

</td>

<td style="text-align:right;">

0.30

</td>

<td style="text-align:left;">

1.66e-07

</td>

<td style="text-align:right;">

25828

</td>

<td style="text-align:right;">

70250

</td>

<td style="text-align:left;">

hsa-mir-10a

</td>

<td style="text-align:left;">

MIR-10

</td>

<td style="text-align:left;">

ACCCUGU

</td>

<td style="text-align:left;">

chr17

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-146-P1\_5p

</td>

<td style="text-align:right;">

\-1.69

</td>

<td style="text-align:right;">

0.34

</td>

<td style="text-align:left;">

9.05e-06

</td>

<td style="text-align:right;">

351

</td>

<td style="text-align:right;">

891

</td>

<td style="text-align:left;">

hsa-mir-146a

</td>

<td style="text-align:left;">

MIR-146

</td>

<td style="text-align:left;">

GAGAACU

</td>

<td style="text-align:left;">

chr5

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-92-P1c\_3p

</td>

<td style="text-align:right;">

\-1.32

</td>

<td style="text-align:right;">

0.26

</td>

<td style="text-align:left;">

4.00e-06

</td>

<td style="text-align:right;">

303

</td>

<td style="text-align:right;">

634

</td>

<td style="text-align:left;">

hsa-mir-92b

</td>

<td style="text-align:left;">

MIR-92

</td>

<td style="text-align:left;">

AUUGCAC

</td>

<td style="text-align:left;">

chr1

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-181-P1c\_5p

</td>

<td style="text-align:right;">

\-1.25

</td>

<td style="text-align:right;">

0.21

</td>

<td style="text-align:left;">

6.34e-08

</td>

<td style="text-align:right;">

261

</td>

<td style="text-align:right;">

527

</td>

<td style="text-align:left;">

hsa-mir-181c

</td>

<td style="text-align:left;">

MIR-181

</td>

<td style="text-align:left;">

ACAUUCA

</td>

<td style="text-align:left;">

chr19

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-425\_5p

</td>

<td style="text-align:right;">

\-1.16

</td>

<td style="text-align:right;">

0.19

</td>

<td style="text-align:left;">

6.36e-09

</td>

<td style="text-align:right;">

247

</td>

<td style="text-align:right;">

499

</td>

<td style="text-align:left;">

hsa-mir-425

</td>

<td style="text-align:left;">

MIR-425

</td>

<td style="text-align:left;">

AUGACAC

</td>

<td style="text-align:left;">

chr3

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-362-P3\_3p

</td>

<td style="text-align:right;">

\-1.15

</td>

<td style="text-align:right;">

0.36

</td>

<td style="text-align:left;">

2.88e-03

</td>

<td style="text-align:right;">

60

</td>

<td style="text-align:right;">

103

</td>

<td style="text-align:left;">

hsa-mir-501

</td>

<td style="text-align:left;">

MIR-362

</td>

<td style="text-align:left;">

AUGCACC

</td>

<td style="text-align:left;">

chrX

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-192-P1\_5p

</td>

<td style="text-align:right;">

\-1.10

</td>

<td style="text-align:right;">

0.27

</td>

<td style="text-align:left;">

7.22e-05

</td>

<td style="text-align:right;">

109276

</td>

<td style="text-align:right;">

212293

</td>

<td style="text-align:left;">

hsa-mir-192

</td>

<td style="text-align:left;">

MIR-192

</td>

<td style="text-align:left;">

UGACCUA

</td>

<td style="text-align:left;">

chr11

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-10-P2b\_5p

</td>

<td style="text-align:right;">

\-1.08

</td>

<td style="text-align:right;">

0.36

</td>

<td style="text-align:left;">

6.81e-03

</td>

<td style="text-align:right;">

1478

</td>

<td style="text-align:right;">

2499

</td>

<td style="text-align:left;">

hsa-mir-99b

</td>

<td style="text-align:left;">

MIR-10

</td>

<td style="text-align:left;">

ACCCGUA

</td>

<td style="text-align:left;">

chr19

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-210\_3p

</td>

<td style="text-align:right;">

\-1.07

</td>

<td style="text-align:right;">

0.28

</td>

<td style="text-align:left;">

1.53e-03

</td>

<td style="text-align:right;">

106

</td>

<td style="text-align:right;">

185

</td>

<td style="text-align:left;">

hsa-mir-210

</td>

<td style="text-align:left;">

MIR-210

</td>

<td style="text-align:left;">

UGUGCGU

</td>

<td style="text-align:left;">

chr11

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-132-P1\_3p

</td>

<td style="text-align:right;">

\-1.03

</td>

<td style="text-align:right;">

0.20

</td>

<td style="text-align:left;">

2.09e-06

</td>

<td style="text-align:right;">

61

</td>

<td style="text-align:right;">

117

</td>

<td style="text-align:left;">

hsa-mir-132

</td>

<td style="text-align:left;">

MIR-132

</td>

<td style="text-align:left;">

AACAGUC

</td>

<td style="text-align:left;">

chr17

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-194-P1\_5p/P2\_5p

</td>

<td style="text-align:right;">

\-0.94

</td>

<td style="text-align:right;">

0.27

</td>

<td style="text-align:left;">

5.19e-04

</td>

<td style="text-align:right;">

6323

</td>

<td style="text-align:right;">

11619

</td>

<td style="text-align:left;">

hsa-mir-194-2

</td>

<td style="text-align:left;">

MIR-194

</td>

<td style="text-align:left;">

GUAACAG

</td>

<td style="text-align:left;">

chr11

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-130-P4a\_3p

</td>

<td style="text-align:right;">

\-0.84

</td>

<td style="text-align:right;">

0.13

</td>

<td style="text-align:left;">

2.77e-10

</td>

<td style="text-align:right;">

70

</td>

<td style="text-align:right;">

111

</td>

<td style="text-align:left;">

hsa-mir-454

</td>

<td style="text-align:left;">

MIR-130

</td>

<td style="text-align:left;">

AGUGCAA

</td>

<td style="text-align:left;">

chr17

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-188-P2\_5p

</td>

<td style="text-align:right;">

\-0.82

</td>

<td style="text-align:right;">

0.19

</td>

<td style="text-align:left;">

9.56e-05

</td>

<td style="text-align:right;">

279

</td>

<td style="text-align:right;">

427

</td>

<td style="text-align:left;">

hsa-mir-532

</td>

<td style="text-align:left;">

MIR-188

</td>

<td style="text-align:left;">

AUGCCUU

</td>

<td style="text-align:left;">

chrX

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-191\_5p

</td>

<td style="text-align:right;">

\-0.78

</td>

<td style="text-align:right;">

0.23

</td>

<td style="text-align:left;">

1.91e-03

</td>

<td style="text-align:right;">

10721

</td>

<td style="text-align:right;">

14608

</td>

<td style="text-align:left;">

hsa-mir-191

</td>

<td style="text-align:left;">

MIR-191

</td>

<td style="text-align:left;">

AACGGAA

</td>

<td style="text-align:left;">

chr3

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-21\_5p

</td>

<td style="text-align:right;">

\-0.76

</td>

<td style="text-align:right;">

0.17

</td>

<td style="text-align:left;">

1.54e-04

</td>

<td style="text-align:right;">

19240

</td>

<td style="text-align:right;">

28173

</td>

<td style="text-align:left;">

hsa-mir-21

</td>

<td style="text-align:left;">

MIR-21

</td>

<td style="text-align:left;">

AGCUUAU

</td>

<td style="text-align:left;">

chr17

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-130-P2a\_3p

</td>

<td style="text-align:right;">

\-0.71

</td>

<td style="text-align:right;">

0.20

</td>

<td style="text-align:left;">

3.57e-03

</td>

<td style="text-align:right;">

92

</td>

<td style="text-align:right;">

131

</td>

<td style="text-align:left;">

hsa-mir-301a

</td>

<td style="text-align:left;">

MIR-130

</td>

<td style="text-align:left;">

AGUGCAA

</td>

<td style="text-align:left;">

chr17

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-1307\_3p

</td>

<td style="text-align:right;">

\-0.66

</td>

<td style="text-align:right;">

0.19

</td>

<td style="text-align:left;">

1.89e-03

</td>

<td style="text-align:right;">

89

</td>

<td style="text-align:right;">

127

</td>

<td style="text-align:left;">

hsa-mir-1307

</td>

<td style="text-align:left;">

MIR-1307

</td>

<td style="text-align:left;">

CGACCGG

</td>

<td style="text-align:left;">

chr10

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-362-P2\_3p/P4\_3p

</td>

<td style="text-align:right;">

\-0.61

</td>

<td style="text-align:right;">

0.20

</td>

<td style="text-align:left;">

7.29e-03

</td>

<td style="text-align:right;">

211

</td>

<td style="text-align:right;">

268

</td>

<td style="text-align:left;">

hsa-mir-500a

</td>

<td style="text-align:left;">

MIR-362

</td>

<td style="text-align:left;">

UGCACCU

</td>

<td style="text-align:left;">

chrX

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-1307\_5p

</td>

<td style="text-align:right;">

\-0.61

</td>

<td style="text-align:right;">

0.26

</td>

<td style="text-align:left;">

3.36e-02

</td>

<td style="text-align:right;">

502

</td>

<td style="text-align:right;">

713

</td>

<td style="text-align:left;">

hsa-mir-1307

</td>

<td style="text-align:left;">

MIR-1307

</td>

<td style="text-align:left;">

CGACCGG

</td>

<td style="text-align:left;">

chr10

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

</tbody>

</table>

</div>

``` r
# Number of downregulated miRNA
signature_mirnas$number_downregulated
```

    ## [1] 44

## nCR vs nLu

``` r
column='tissue.type'
tissue_type_A <- 'normal.lung'
tissue_type_B <- 'normal.colorect'
norm_adj_up       = "None"
norm_adj_down     = "None"
pCRC_adj_up   = "None"
pCRC_adj_down = "None"

coef <- paste(column, tissue_type_A, 'vs', tissue_type_B, sep='_')
res <- DeseqResult(dds, column, coef, tissue_type_A, tissue_type_B,
                   lfc.Threshold, rpm.Threshold,
                   norm_adj_up,
                   norm_adj_down)
dict_sig_mirna[paste(coef, "up",   sep='_')] <- list(res$up_mirna)
dict_sig_mirna[paste(coef, "down", sep='_')] <- list(res$down_mirna)
res_res <- res$res
res_dict[coef] <- res_res
plotMA(res$res, alpha=0.05)
```

![](comet_analysis_report_automated_no_lfcthreshold_08.12.20_res_alpha_0.05_with_LFC_shrinkage_LFC_0.58_mirge3_7.0_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

``` r
# Plot volcano plot
VolcanoPlot(res$res, coef, res$sig,
            res$up_mirna, res$down_mirna,
            norm_adj_up, norm_adj_down,
            pCRC_adj_up, pCRC_adj_down)
```

![](comet_analysis_report_automated_no_lfcthreshold_08.12.20_res_alpha_0.05_with_LFC_shrinkage_LFC_0.58_mirge3_7.0_files/figure-gfm/unnamed-chunk-16-2.png)<!-- -->

``` r
ExpressionPlot(res$res, res$rpm, coef, res$sig,
               tissue_type_A, tissue_type_B,
               res$up_mirna, res$down_mirna,
               norm_adj_up, norm_adj_down,
               pCRC_adj_up, pCRC_adj_down)
```

![](comet_analysis_report_automated_no_lfcthreshold_08.12.20_res_alpha_0.05_with_LFC_shrinkage_LFC_0.58_mirge3_7.0_files/figure-gfm/unnamed-chunk-16-3.png)<!-- -->

``` r
signature_mirnas <- SigList(res, dds, tissue_type_A, tissue_type_B, coef,
                            norm_adj_up, norm_adj_down, 
                            pCRC_adj_up, pCRC_adj_down)
# Print list upregulated miRNA
signature_mirnas$up_mirna
```

<div style="border: 1px solid #ddd; padding: 5px; overflow-x: scroll; width:2000px; ">

<table class="table table-striped table-hover table-condensed" style="margin-left: auto; margin-right: auto;">

<caption>

Upregulated in tissue.type\_normal.lung\_vs\_normal.colorect

</caption>

<thead>

<tr>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

miRNA

</th>

<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">

LFC

</th>

<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">

lfcSE

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

FDR

</th>

<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">

RPM normal.lung

</th>

<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">

RPM normal.colorect

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

miRBase\_ID

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

Family

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

Seed

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

Chr

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

Cell-Type Specific

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

Norm Background

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

pCRC Background

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

Hsa-Mir-335\_5p

</td>

<td style="text-align:right;">

1.48

</td>

<td style="text-align:right;">

0.21

</td>

<td style="text-align:left;">

3.73e-11

</td>

<td style="text-align:right;">

300

</td>

<td style="text-align:right;">

79

</td>

<td style="text-align:left;">

hsa-mir-335

</td>

<td style="text-align:left;">

MIR-335

</td>

<td style="text-align:left;">

CAAGAGC

</td>

<td style="text-align:left;">

chr7

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">Retinal
Epithelial Cell</span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-144\_5p

</td>

<td style="text-align:right;">

2.32

</td>

<td style="text-align:right;">

0.41

</td>

<td style="text-align:left;">

5.46e-07

</td>

<td style="text-align:right;">

451

</td>

<td style="text-align:right;">

65

</td>

<td style="text-align:left;">

hsa-mir-144

</td>

<td style="text-align:left;">

MIR-144

</td>

<td style="text-align:left;">

GAUAUCA

</td>

<td style="text-align:left;">

chr17

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">Red
Blood Cell</span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-451\_5p

</td>

<td style="text-align:right;">

1.90

</td>

<td style="text-align:right;">

0.47

</td>

<td style="text-align:left;">

1.06e-03

</td>

<td style="text-align:right;">

12496

</td>

<td style="text-align:right;">

2815

</td>

<td style="text-align:left;">

hsa-mir-451a

</td>

<td style="text-align:left;">

MIR-451

</td>

<td style="text-align:left;">

AACCGUU

</td>

<td style="text-align:left;">

chr17

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">Red
Blood Cell</span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-24-P1\_3p/P2\_3p

</td>

<td style="text-align:right;">

0.86

</td>

<td style="text-align:right;">

0.25

</td>

<td style="text-align:left;">

1.88e-03

</td>

<td style="text-align:right;">

1084

</td>

<td style="text-align:right;">

413

</td>

<td style="text-align:left;">

hsa-mir-24-2

</td>

<td style="text-align:left;">

MIR-24

</td>

<td style="text-align:left;">

GGCUCAG

</td>

<td style="text-align:left;">

chr19

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">Macrophage</span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-486\_5p

</td>

<td style="text-align:right;">

2.04

</td>

<td style="text-align:right;">

0.45

</td>

<td style="text-align:left;">

1.33e-04

</td>

<td style="text-align:right;">

22081

</td>

<td style="text-align:right;">

3961

</td>

<td style="text-align:left;">

hsa-mir-486-1

</td>

<td style="text-align:left;">

MIR-486

</td>

<td style="text-align:left;">

CCUGUAC

</td>

<td style="text-align:left;">

chr8

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">c(“Platelet”,
“Red Blood Cell”)</span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-126\_5p

</td>

<td style="text-align:right;">

2.57

</td>

<td style="text-align:right;">

0.26

</td>

<td style="text-align:left;">

4.06e-21

</td>

<td style="text-align:right;">

33620

</td>

<td style="text-align:right;">

4230

</td>

<td style="text-align:left;">

hsa-mir-126

</td>

<td style="text-align:left;">

MIR-126

</td>

<td style="text-align:left;">

AUUAUUA

</td>

<td style="text-align:left;">

chr9

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">c(“Endothelial
Cell”, “Platelet”)</span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-146-P2\_5p

</td>

<td style="text-align:right;">

1.65

</td>

<td style="text-align:right;">

0.41

</td>

<td style="text-align:left;">

1.94e-04

</td>

<td style="text-align:right;">

17435

</td>

<td style="text-align:right;">

3259

</td>

<td style="text-align:left;">

hsa-mir-146b

</td>

<td style="text-align:left;">

MIR-146

</td>

<td style="text-align:left;">

GAGAACU

</td>

<td style="text-align:left;">

chr10

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">c(“Dendritic
Cell”, “Lymphocyte”)</span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-342\_3p

</td>

<td style="text-align:right;">

1.00

</td>

<td style="text-align:right;">

0.31

</td>

<td style="text-align:left;">

6.94e-03

</td>

<td style="text-align:right;">

681

</td>

<td style="text-align:right;">

260

</td>

<td style="text-align:left;">

hsa-mir-342

</td>

<td style="text-align:left;">

MIR-342

</td>

<td style="text-align:left;">

CUCACAC

</td>

<td style="text-align:left;">

chr14

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">c(“Dendritic
Cell”, “Lymphocyte”, “Macrophage”)</span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-34-P2b\_5p

</td>

<td style="text-align:right;">

5.29

</td>

<td style="text-align:right;">

0.43

</td>

<td style="text-align:left;">

1.18e-33

</td>

<td style="text-align:right;">

748

</td>

<td style="text-align:right;">

11

</td>

<td style="text-align:left;">

hsa-mir-34c

</td>

<td style="text-align:left;">

MIR-34

</td>

<td style="text-align:left;">

GGCAGUG

</td>

<td style="text-align:left;">

chr11

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-34-P2a\_5p

</td>

<td style="text-align:right;">

4.98

</td>

<td style="text-align:right;">

0.49

</td>

<td style="text-align:left;">

3.17e-22

</td>

<td style="text-align:right;">

116

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

hsa-mir-34b

</td>

<td style="text-align:left;">

MIR-34

</td>

<td style="text-align:left;">

GGCAGUG

</td>

<td style="text-align:left;">

chr11

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-184\_3p

</td>

<td style="text-align:right;">

4.75

</td>

<td style="text-align:right;">

0.92

</td>

<td style="text-align:left;">

5.98e-06

</td>

<td style="text-align:right;">

111

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

hsa-mir-184

</td>

<td style="text-align:left;">

MIR-184

</td>

<td style="text-align:left;">

GGACGGA

</td>

<td style="text-align:left;">

chr15

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-30-P1a\_5p

</td>

<td style="text-align:right;">

3.04

</td>

<td style="text-align:right;">

0.28

</td>

<td style="text-align:left;">

6.17e-24

</td>

<td style="text-align:right;">

60860

</td>

<td style="text-align:right;">

5150

</td>

<td style="text-align:left;">

hsa-mir-30a

</td>

<td style="text-align:left;">

MIR-30

</td>

<td style="text-align:left;">

GUAAACA

</td>

<td style="text-align:left;">

chr6

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-218-P1\_5p/P2\_5p

</td>

<td style="text-align:right;">

2.39

</td>

<td style="text-align:right;">

0.33

</td>

<td style="text-align:left;">

7.14e-11

</td>

<td style="text-align:right;">

374

</td>

<td style="text-align:right;">

52

</td>

<td style="text-align:left;">

hsa-mir-218-1

</td>

<td style="text-align:left;">

MIR-218

</td>

<td style="text-align:left;">

UGUGCUU

</td>

<td style="text-align:left;">

chr4

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-10-P2c\_5p

</td>

<td style="text-align:right;">

2.23

</td>

<td style="text-align:right;">

0.41

</td>

<td style="text-align:left;">

4.40e-07

</td>

<td style="text-align:right;">

804

</td>

<td style="text-align:right;">

114

</td>

<td style="text-align:left;">

hsa-mir-99a

</td>

<td style="text-align:left;">

MIR-10

</td>

<td style="text-align:left;">

ACCCGUA

</td>

<td style="text-align:left;">

chr21

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-181-P1a\_5p/P1b\_5p

</td>

<td style="text-align:right;">

2.19

</td>

<td style="text-align:right;">

0.24

</td>

<td style="text-align:left;">

1.16e-18

</td>

<td style="text-align:right;">

80063

</td>

<td style="text-align:right;">

12321

</td>

<td style="text-align:left;">

hsa-mir-181a-1

</td>

<td style="text-align:left;">

MIR-181

</td>

<td style="text-align:left;">

ACAUUCA

</td>

<td style="text-align:left;">

chr1

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-181-P2a\_5p/P2b\_5p

</td>

<td style="text-align:right;">

2.06

</td>

<td style="text-align:right;">

0.23

</td>

<td style="text-align:left;">

2.25e-18

</td>

<td style="text-align:right;">

3799

</td>

<td style="text-align:right;">

650

</td>

<td style="text-align:left;">

hsa-mir-181b-1

</td>

<td style="text-align:left;">

MIR-181

</td>

<td style="text-align:left;">

ACAUUCA

</td>

<td style="text-align:left;">

chr1

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Let-7-P1c\_5p

</td>

<td style="text-align:right;">

1.97

</td>

<td style="text-align:right;">

0.35

</td>

<td style="text-align:left;">

1.57e-07

</td>

<td style="text-align:right;">

2606

</td>

<td style="text-align:right;">

466

</td>

<td style="text-align:left;">

hsa-let-7c

</td>

<td style="text-align:left;">

LET-7

</td>

<td style="text-align:left;">

GAGGUAG

</td>

<td style="text-align:left;">

chr21

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-130-P1a\_3p

</td>

<td style="text-align:right;">

1.92

</td>

<td style="text-align:right;">

0.25

</td>

<td style="text-align:left;">

5.89e-13

</td>

<td style="text-align:right;">

1941

</td>

<td style="text-align:right;">

379

</td>

<td style="text-align:left;">

hsa-mir-130a

</td>

<td style="text-align:left;">

MIR-130

</td>

<td style="text-align:left;">

AGUGCAA

</td>

<td style="text-align:left;">

chr11

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-30-P1c\_5p

</td>

<td style="text-align:right;">

1.74

</td>

<td style="text-align:right;">

0.20

</td>

<td style="text-align:left;">

4.37e-16

</td>

<td style="text-align:right;">

34088

</td>

<td style="text-align:right;">

7243

</td>

<td style="text-align:left;">

hsa-mir-30d

</td>

<td style="text-align:left;">

MIR-30

</td>

<td style="text-align:left;">

GUAAACA

</td>

<td style="text-align:left;">

chr8

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-101-P1\_3p/P2\_3p

</td>

<td style="text-align:right;">

1.45

</td>

<td style="text-align:right;">

0.22

</td>

<td style="text-align:left;">

1.42e-09

</td>

<td style="text-align:right;">

16084

</td>

<td style="text-align:right;">

4197

</td>

<td style="text-align:left;">

hsa-mir-101-1

</td>

<td style="text-align:left;">

MIR-101

</td>

<td style="text-align:left;">

UACAGUA

</td>

<td style="text-align:left;">

chr1

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-338-P1\_3p

</td>

<td style="text-align:right;">

1.45

</td>

<td style="text-align:right;">

0.39

</td>

<td style="text-align:left;">

2.15e-03

</td>

<td style="text-align:right;">

605

</td>

<td style="text-align:right;">

166

</td>

<td style="text-align:left;">

hsa-mir-338

</td>

<td style="text-align:left;">

MIR-338

</td>

<td style="text-align:left;">

CCAGCAU

</td>

<td style="text-align:left;">

chr17

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-181-P2c\_5p

</td>

<td style="text-align:right;">

1.42

</td>

<td style="text-align:right;">

0.28

</td>

<td style="text-align:left;">

1.14e-06

</td>

<td style="text-align:right;">

200

</td>

<td style="text-align:right;">

51

</td>

<td style="text-align:left;">

hsa-mir-181d

</td>

<td style="text-align:left;">

MIR-181

</td>

<td style="text-align:left;">

ACAUUCA

</td>

<td style="text-align:left;">

chr19

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-10-P2a\_5p

</td>

<td style="text-align:right;">

1.40

</td>

<td style="text-align:right;">

0.40

</td>

<td style="text-align:left;">

2.08e-03

</td>

<td style="text-align:right;">

4056

</td>

<td style="text-align:right;">

986

</td>

<td style="text-align:left;">

hsa-mir-100

</td>

<td style="text-align:left;">

MIR-10

</td>

<td style="text-align:left;">

ACCCGUA

</td>

<td style="text-align:left;">

chr11

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-10-P3a\_5p

</td>

<td style="text-align:right;">

1.37

</td>

<td style="text-align:right;">

0.38

</td>

<td style="text-align:left;">

1.96e-03

</td>

<td style="text-align:right;">

507

</td>

<td style="text-align:right;">

137

</td>

<td style="text-align:left;">

hsa-mir-125b-1

</td>

<td style="text-align:left;">

MIR-10

</td>

<td style="text-align:left;">

CCCUGAG

</td>

<td style="text-align:left;">

chr11

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-10-P3b\_5p

</td>

<td style="text-align:right;">

1.34

</td>

<td style="text-align:right;">

0.34

</td>

<td style="text-align:left;">

7.55e-04

</td>

<td style="text-align:right;">

8561

</td>

<td style="text-align:right;">

2392

</td>

<td style="text-align:left;">

hsa-mir-125a

</td>

<td style="text-align:left;">

MIR-10

</td>

<td style="text-align:left;">

CCCUGAG

</td>

<td style="text-align:left;">

chr19

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-181-P1c\_5p

</td>

<td style="text-align:right;">

1.25

</td>

<td style="text-align:right;">

0.26

</td>

<td style="text-align:left;">

5.14e-06

</td>

<td style="text-align:right;">

1790

</td>

<td style="text-align:right;">

527

</td>

<td style="text-align:left;">

hsa-mir-181c

</td>

<td style="text-align:left;">

MIR-181

</td>

<td style="text-align:left;">

ACAUUCA

</td>

<td style="text-align:left;">

chr19

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-140\_3p

</td>

<td style="text-align:right;">

1.25

</td>

<td style="text-align:right;">

0.19

</td>

<td style="text-align:left;">

1.88e-09

</td>

<td style="text-align:right;">

2195

</td>

<td style="text-align:right;">

686

</td>

<td style="text-align:left;">

hsa-mir-140

</td>

<td style="text-align:left;">

MIR-140

</td>

<td style="text-align:left;">

CCACAGG

</td>

<td style="text-align:left;">

chr16

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Let-7-P2b1\_5p

</td>

<td style="text-align:right;">

1.10

</td>

<td style="text-align:right;">

0.15

</td>

<td style="text-align:left;">

1.06e-11

</td>

<td style="text-align:right;">

6253

</td>

<td style="text-align:right;">

2110

</td>

<td style="text-align:left;">

hsa-let-7f-1

</td>

<td style="text-align:left;">

LET-7

</td>

<td style="text-align:left;">

GAGGUAG

</td>

<td style="text-align:left;">

chr9

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-30-P2a\_5p/P2b\_5p/P2c\_5p

</td>

<td style="text-align:right;">

1.08

</td>

<td style="text-align:right;">

0.19

</td>

<td style="text-align:left;">

9.04e-08

</td>

<td style="text-align:right;">

7603

</td>

<td style="text-align:right;">

2677

</td>

<td style="text-align:left;">

hsa-mir-30c-2

</td>

<td style="text-align:left;">

MIR-30

</td>

<td style="text-align:left;">

GUAAACA

</td>

<td style="text-align:left;">

chr6

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-221-P1\_3p

</td>

<td style="text-align:right;">

0.95

</td>

<td style="text-align:right;">

0.23

</td>

<td style="text-align:left;">

1.40e-04

</td>

<td style="text-align:right;">

2502

</td>

<td style="text-align:right;">

964

</td>

<td style="text-align:left;">

hsa-mir-221

</td>

<td style="text-align:left;">

MIR-221

</td>

<td style="text-align:left;">

GCUACAU

</td>

<td style="text-align:left;">

chrX

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-92-P1c\_3p

</td>

<td style="text-align:right;">

0.90

</td>

<td style="text-align:right;">

0.32

</td>

<td style="text-align:left;">

1.56e-02

</td>

<td style="text-align:right;">

1737

</td>

<td style="text-align:right;">

634

</td>

<td style="text-align:left;">

hsa-mir-92b

</td>

<td style="text-align:left;">

MIR-92

</td>

<td style="text-align:left;">

AUUGCAC

</td>

<td style="text-align:left;">

chr1

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-15-P2c\_5p

</td>

<td style="text-align:right;">

0.86

</td>

<td style="text-align:right;">

0.30

</td>

<td style="text-align:left;">

2.04e-02

</td>

<td style="text-align:right;">

1552

</td>

<td style="text-align:right;">

672

</td>

<td style="text-align:left;">

hsa-mir-195

</td>

<td style="text-align:left;">

MIR-15

</td>

<td style="text-align:left;">

AGCAGCA

</td>

<td style="text-align:left;">

chr17

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-374-P2\_5p

</td>

<td style="text-align:right;">

0.80

</td>

<td style="text-align:right;">

0.25

</td>

<td style="text-align:left;">

3.59e-03

</td>

<td style="text-align:right;">

148

</td>

<td style="text-align:right;">

65

</td>

<td style="text-align:left;">

hsa-mir-374b

</td>

<td style="text-align:left;">

MIR-374

</td>

<td style="text-align:left;">

UAUAAUA

</td>

<td style="text-align:left;">

chrX

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-221-P2\_3p

</td>

<td style="text-align:right;">

0.77

</td>

<td style="text-align:right;">

0.25

</td>

<td style="text-align:left;">

5.86e-03

</td>

<td style="text-align:right;">

1346

</td>

<td style="text-align:right;">

574

</td>

<td style="text-align:left;">

hsa-mir-222

</td>

<td style="text-align:left;">

MIR-221

</td>

<td style="text-align:left;">

GCUACAU

</td>

<td style="text-align:left;">

chrX

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Let-7-P2c2\_5p

</td>

<td style="text-align:right;">

0.76

</td>

<td style="text-align:right;">

0.21

</td>

<td style="text-align:left;">

1.16e-03

</td>

<td style="text-align:right;">

5088

</td>

<td style="text-align:right;">

2218

</td>

<td style="text-align:left;">

hsa-let-7i

</td>

<td style="text-align:left;">

LET-7

</td>

<td style="text-align:left;">

GAGGUAG

</td>

<td style="text-align:left;">

chr12

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-130-P2a\_3p

</td>

<td style="text-align:right;">

0.75

</td>

<td style="text-align:right;">

0.25

</td>

<td style="text-align:left;">

6.48e-03

</td>

<td style="text-align:right;">

305

</td>

<td style="text-align:right;">

131

</td>

<td style="text-align:left;">

hsa-mir-301a

</td>

<td style="text-align:left;">

MIR-130

</td>

<td style="text-align:left;">

AGUGCAA

</td>

<td style="text-align:left;">

chr17

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-652\_3p

</td>

<td style="text-align:right;">

0.71

</td>

<td style="text-align:right;">

0.20

</td>

<td style="text-align:left;">

1.38e-03

</td>

<td style="text-align:right;">

187

</td>

<td style="text-align:right;">

86

</td>

<td style="text-align:left;">

hsa-mir-652

</td>

<td style="text-align:left;">

MIR-652

</td>

<td style="text-align:left;">

AUGGCGC

</td>

<td style="text-align:left;">

chrX

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-26-P1\_5p/P2\_5p

</td>

<td style="text-align:right;">

0.68

</td>

<td style="text-align:right;">

0.17

</td>

<td style="text-align:left;">

3.54e-04

</td>

<td style="text-align:right;">

123947

</td>

<td style="text-align:right;">

56268

</td>

<td style="text-align:left;">

hsa-mir-26b

</td>

<td style="text-align:left;">

MIR-26

</td>

<td style="text-align:left;">

UCAAGUA

</td>

<td style="text-align:left;">

chr2

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-28-P2\_3p

</td>

<td style="text-align:right;">

0.63

</td>

<td style="text-align:right;">

0.21

</td>

<td style="text-align:left;">

9.02e-03

</td>

<td style="text-align:right;">

6007

</td>

<td style="text-align:right;">

2638

</td>

<td style="text-align:left;">

hsa-mir-151a

</td>

<td style="text-align:left;">

MIR-28

</td>

<td style="text-align:left;">

CGAGGAG

</td>

<td style="text-align:left;">

chr8

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

</tbody>

</table>

</div>

``` r
# Number of upregulated miRNA
signature_mirnas$number_upregulated
```

    ## [1] 39

``` r
# Print list downregulated miRNA
signature_mirnas$down_mirna
```

<div style="border: 1px solid #ddd; padding: 5px; overflow-x: scroll; width:2000px; ">

<table class="table table-striped table-hover table-condensed" style="margin-left: auto; margin-right: auto;">

<caption>

Downregulated in tissue.type\_normal.lung\_vs\_normal.colorect

</caption>

<thead>

<tr>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

miRNA

</th>

<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">

LFC

</th>

<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">

lfcSE

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

FDR

</th>

<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">

RPM normal.lung

</th>

<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">

RPM normal.colorect

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

miRBase\_ID

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

Family

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

Seed

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

Chr

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

Cell-Type Specific

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

Norm Background

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

pCRC Background

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

Hsa-Mir-143\_3p

</td>

<td style="text-align:right;">

\-0.81

</td>

<td style="text-align:right;">

0.34

</td>

<td style="text-align:left;">

2.00e-02

</td>

<td style="text-align:right;">

133600

</td>

<td style="text-align:right;">

164812

</td>

<td style="text-align:left;">

hsa-mir-143

</td>

<td style="text-align:left;">

MIR-143

</td>

<td style="text-align:left;">

GAGAUGA

</td>

<td style="text-align:left;">

chr5

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">Mesenchymal</span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-192-P1\_5p/P2\_5p

</td>

<td style="text-align:right;">

\-8.51

</td>

<td style="text-align:right;">

0.40

</td>

<td style="text-align:left;">

3.64e-97

</td>

<td style="text-align:right;">

142

</td>

<td style="text-align:right;">

46010

</td>

<td style="text-align:left;">

hsa-mir-192

</td>

<td style="text-align:left;">

MIR-192

</td>

<td style="text-align:left;">

UGACCUA

</td>

<td style="text-align:left;">

chr11

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">Epithelial
Cell</span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-8-P2a\_3p

</td>

<td style="text-align:right;">

\-4.00

</td>

<td style="text-align:right;">

0.29

</td>

<td style="text-align:left;">

6.74e-40

</td>

<td style="text-align:right;">

618

</td>

<td style="text-align:right;">

7741

</td>

<td style="text-align:left;">

hsa-mir-200b

</td>

<td style="text-align:left;">

MIR-8

</td>

<td style="text-align:left;">

AAUACUG

</td>

<td style="text-align:left;">

chr1

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">Epithelial
Cell</span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-8-P2b\_3p

</td>

<td style="text-align:right;">

\-2.06

</td>

<td style="text-align:right;">

0.29

</td>

<td style="text-align:left;">

1.04e-11

</td>

<td style="text-align:right;">

766

</td>

<td style="text-align:right;">

2444

</td>

<td style="text-align:left;">

hsa-mir-200c

</td>

<td style="text-align:left;">

MIR-8

</td>

<td style="text-align:left;">

AAUACUG

</td>

<td style="text-align:left;">

chr12

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">Epithelial
Cell</span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-17-P1a\_5p/P1b\_5p

</td>

<td style="text-align:right;">

\-0.77

</td>

<td style="text-align:right;">

0.24

</td>

<td style="text-align:left;">

9.07e-03

</td>

<td style="text-align:right;">

184

</td>

<td style="text-align:right;">

225

</td>

<td style="text-align:left;">

hsa-mir-17

</td>

<td style="text-align:left;">

MIR-17

</td>

<td style="text-align:left;">

AAAGUGC

</td>

<td style="text-align:left;">

chr13

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">CD14+
Monocyte</span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-133-P1\_3p/P2\_3p/P3\_3p

</td>

<td style="text-align:right;">

\-2.33

</td>

<td style="text-align:right;">

0.47

</td>

<td style="text-align:left;">

9.21e-08

</td>

<td style="text-align:right;">

92

</td>

<td style="text-align:right;">

460

</td>

<td style="text-align:left;">

hsa-mir-133a-2

</td>

<td style="text-align:left;">

MIR-133

</td>

<td style="text-align:left;">

UUGGUCC

</td>

<td style="text-align:left;">

chr20

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">c(“Skeletal
Myocyte”, “Stem Cell”)</span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-375\_3p

</td>

<td style="text-align:right;">

\-3.01

</td>

<td style="text-align:right;">

0.47

</td>

<td style="text-align:left;">

2.15e-10

</td>

<td style="text-align:right;">

2346

</td>

<td style="text-align:right;">

17215

</td>

<td style="text-align:left;">

hsa-mir-375

</td>

<td style="text-align:left;">

MIR-375

</td>

<td style="text-align:left;">

UUGUUCG

</td>

<td style="text-align:left;">

chr2

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">c(“Epithelial
Cell”, “Islet Cell”, “Neural”)</span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-194-P1\_5p/P2\_5p

</td>

<td style="text-align:right;">

\-7.89

</td>

<td style="text-align:right;">

0.34

</td>

<td style="text-align:left;">

2.93e-114

</td>

<td style="text-align:right;">

58

</td>

<td style="text-align:right;">

11619

</td>

<td style="text-align:left;">

hsa-mir-194-2

</td>

<td style="text-align:left;">

MIR-194

</td>

<td style="text-align:left;">

GUAACAG

</td>

<td style="text-align:left;">

chr11

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-192-P1\_5p

</td>

<td style="text-align:right;">

\-7.28

</td>

<td style="text-align:right;">

0.34

</td>

<td style="text-align:left;">

3.81e-95

</td>

<td style="text-align:right;">

1684

</td>

<td style="text-align:right;">

212293

</td>

<td style="text-align:left;">

hsa-mir-192

</td>

<td style="text-align:left;">

MIR-192

</td>

<td style="text-align:left;">

UGACCUA

</td>

<td style="text-align:left;">

chr11

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-196-P1\_5p/P2\_5p

</td>

<td style="text-align:right;">

\-6.68

</td>

<td style="text-align:right;">

0.41

</td>

<td style="text-align:left;">

1.42e-54

</td>

<td style="text-align:right;">

3

</td>

<td style="text-align:right;">

250

</td>

<td style="text-align:left;">

hsa-mir-196a-1

</td>

<td style="text-align:left;">

MIR-196

</td>

<td style="text-align:left;">

AGGUAGU

</td>

<td style="text-align:left;">

chr17

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-196-P3\_5p

</td>

<td style="text-align:right;">

\-5.95

</td>

<td style="text-align:right;">

0.47

</td>

<td style="text-align:left;">

9.59e-33

</td>

<td style="text-align:right;">

7

</td>

<td style="text-align:right;">

400

</td>

<td style="text-align:left;">

hsa-mir-196b

</td>

<td style="text-align:left;">

MIR-196

</td>

<td style="text-align:left;">

AGGUAGU

</td>

<td style="text-align:left;">

chr7

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-577\_5p

</td>

<td style="text-align:right;">

\-5.95

</td>

<td style="text-align:right;">

0.41

</td>

<td style="text-align:left;">

4.39e-45

</td>

<td style="text-align:right;">

3

</td>

<td style="text-align:right;">

138

</td>

<td style="text-align:left;">

hsa-mir-577

</td>

<td style="text-align:left;">

MIR-577

</td>

<td style="text-align:left;">

UAGAUAA

</td>

<td style="text-align:left;">

chr4

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-147\_3p

</td>

<td style="text-align:right;">

\-5.82

</td>

<td style="text-align:right;">

0.43

</td>

<td style="text-align:left;">

5.54e-40

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

110

</td>

<td style="text-align:left;">

hsa-mir-147b

</td>

<td style="text-align:left;">

MIR-147

</td>

<td style="text-align:left;">

UGUGCGG

</td>

<td style="text-align:left;">

chr15

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-190-P1\_5p

</td>

<td style="text-align:right;">

\-4.03

</td>

<td style="text-align:right;">

0.32

</td>

<td style="text-align:left;">

1.35e-34

</td>

<td style="text-align:right;">

22

</td>

<td style="text-align:right;">

294

</td>

<td style="text-align:left;">

hsa-mir-190a

</td>

<td style="text-align:left;">

MIR-190

</td>

<td style="text-align:left;">

GAUAUGU

</td>

<td style="text-align:left;">

chr15

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-378\_3p

</td>

<td style="text-align:right;">

\-3.17

</td>

<td style="text-align:right;">

0.24

</td>

<td style="text-align:left;">

2.10e-38

</td>

<td style="text-align:right;">

2189

</td>

<td style="text-align:right;">

14834

</td>

<td style="text-align:left;">

hsa-mir-378a

</td>

<td style="text-align:left;">

MIR-378

</td>

<td style="text-align:left;">

CUGGACU

</td>

<td style="text-align:left;">

chr5

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-127\_3p

</td>

<td style="text-align:right;">

\-2.96

</td>

<td style="text-align:right;">

0.40

</td>

<td style="text-align:left;">

4.55e-12

</td>

<td style="text-align:right;">

327

</td>

<td style="text-align:right;">

1967

</td>

<td style="text-align:left;">

hsa-mir-127

</td>

<td style="text-align:left;">

MIR-127

</td>

<td style="text-align:left;">

CGGAUCC

</td>

<td style="text-align:left;">

chr14

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-10-P1b\_5p

</td>

<td style="text-align:right;">

\-2.67

</td>

<td style="text-align:right;">

0.45

</td>

<td style="text-align:left;">

5.72e-09

</td>

<td style="text-align:right;">

14125

</td>

<td style="text-align:right;">

75067

</td>

<td style="text-align:left;">

hsa-mir-10b

</td>

<td style="text-align:left;">

MIR-10

</td>

<td style="text-align:left;">

ACCCUGU

</td>

<td style="text-align:left;">

chr2

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-8-P1a\_3p

</td>

<td style="text-align:right;">

\-2.66

</td>

<td style="text-align:right;">

0.29

</td>

<td style="text-align:left;">

3.55e-18

</td>

<td style="text-align:right;">

353

</td>

<td style="text-align:right;">

1754

</td>

<td style="text-align:left;">

hsa-mir-200a

</td>

<td style="text-align:left;">

MIR-8

</td>

<td style="text-align:left;">

AACACUG

</td>

<td style="text-align:left;">

chr1

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-8-P3a\_3p

</td>

<td style="text-align:right;">

\-2.65

</td>

<td style="text-align:right;">

0.30

</td>

<td style="text-align:left;">

1.86e-17

</td>

<td style="text-align:right;">

262

</td>

<td style="text-align:right;">

1282

</td>

<td style="text-align:left;">

hsa-mir-429

</td>

<td style="text-align:left;">

MIR-8

</td>

<td style="text-align:left;">

AAUACUG

</td>

<td style="text-align:left;">

chr1

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-154-P23\_3p

</td>

<td style="text-align:right;">

\-1.79

</td>

<td style="text-align:right;">

0.25

</td>

<td style="text-align:left;">

3.73e-11

</td>

<td style="text-align:right;">

41

</td>

<td style="text-align:right;">

108

</td>

<td style="text-align:left;">

hsa-mir-654

</td>

<td style="text-align:left;">

MIR-154

</td>

<td style="text-align:left;">

AUGUCUG

</td>

<td style="text-align:left;">

chr14

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-425\_5p

</td>

<td style="text-align:right;">

\-1.66

</td>

<td style="text-align:right;">

0.23

</td>

<td style="text-align:left;">

3.12e-11

</td>

<td style="text-align:right;">

213

</td>

<td style="text-align:right;">

499

</td>

<td style="text-align:left;">

hsa-mir-425

</td>

<td style="text-align:left;">

MIR-425

</td>

<td style="text-align:left;">

AUGACAC

</td>

<td style="text-align:left;">

chr3

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-154-P9\_3p

</td>

<td style="text-align:right;">

\-1.37

</td>

<td style="text-align:right;">

0.26

</td>

<td style="text-align:left;">

6.98e-07

</td>

<td style="text-align:right;">

135

</td>

<td style="text-align:right;">

265

</td>

<td style="text-align:left;">

hsa-mir-381

</td>

<td style="text-align:left;">

MIR-154

</td>

<td style="text-align:left;">

AUACAAG

</td>

<td style="text-align:left;">

chr14

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-154-P13\_5p

</td>

<td style="text-align:right;">

\-1.35

</td>

<td style="text-align:right;">

0.27

</td>

<td style="text-align:left;">

3.42e-06

</td>

<td style="text-align:right;">

161

</td>

<td style="text-align:right;">

311

</td>

<td style="text-align:left;">

hsa-mir-411

</td>

<td style="text-align:left;">

MIR-154

</td>

<td style="text-align:left;">

AGUAGAC

</td>

<td style="text-align:left;">

chr14

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-28-P1\_3p

</td>

<td style="text-align:right;">

\-1.24

</td>

<td style="text-align:right;">

0.23

</td>

<td style="text-align:left;">

3.62e-07

</td>

<td style="text-align:right;">

2385

</td>

<td style="text-align:right;">

3956

</td>

<td style="text-align:left;">

hsa-mir-28

</td>

<td style="text-align:left;">

MIR-28

</td>

<td style="text-align:left;">

ACUAGAU

</td>

<td style="text-align:left;">

chr3

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-17-P3a\_5p

</td>

<td style="text-align:right;">

\-1.23

</td>

<td style="text-align:right;">

0.29

</td>

<td style="text-align:left;">

2.53e-04

</td>

<td style="text-align:right;">

194

</td>

<td style="text-align:right;">

322

</td>

<td style="text-align:left;">

hsa-mir-20a

</td>

<td style="text-align:left;">

MIR-17

</td>

<td style="text-align:left;">

AAAGUGC

</td>

<td style="text-align:left;">

chr13

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-1307\_3p

</td>

<td style="text-align:right;">

\-1.09

</td>

<td style="text-align:right;">

0.24

</td>

<td style="text-align:left;">

2.47e-05

</td>

<td style="text-align:right;">

80

</td>

<td style="text-align:right;">

127

</td>

<td style="text-align:left;">

hsa-mir-1307

</td>

<td style="text-align:left;">

MIR-1307

</td>

<td style="text-align:left;">

CGACCGG

</td>

<td style="text-align:left;">

chr10

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-8-P1b\_3p

</td>

<td style="text-align:right;">

\-1.06

</td>

<td style="text-align:right;">

0.38

</td>

<td style="text-align:left;">

1.07e-02

</td>

<td style="text-align:right;">

6145

</td>

<td style="text-align:right;">

9180

</td>

<td style="text-align:left;">

hsa-mir-141

</td>

<td style="text-align:left;">

MIR-8

</td>

<td style="text-align:left;">

AACACUG

</td>

<td style="text-align:left;">

chr12

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-1307\_5p

</td>

<td style="text-align:right;">

\-1.03

</td>

<td style="text-align:right;">

0.32

</td>

<td style="text-align:left;">

4.10e-03

</td>

<td style="text-align:right;">

444

</td>

<td style="text-align:right;">

713

</td>

<td style="text-align:left;">

hsa-mir-1307

</td>

<td style="text-align:left;">

MIR-1307

</td>

<td style="text-align:left;">

CGACCGG

</td>

<td style="text-align:left;">

chr10

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-362-P3\_3p

</td>

<td style="text-align:right;">

\-1.00

</td>

<td style="text-align:right;">

0.45

</td>

<td style="text-align:left;">

4.09e-02

</td>

<td style="text-align:right;">

82

</td>

<td style="text-align:right;">

103

</td>

<td style="text-align:left;">

hsa-mir-501

</td>

<td style="text-align:left;">

MIR-362

</td>

<td style="text-align:left;">

AUGCACC

</td>

<td style="text-align:left;">

chrX

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-210\_3p

</td>

<td style="text-align:right;">

\-0.93

</td>

<td style="text-align:right;">

0.35

</td>

<td style="text-align:left;">

3.36e-02

</td>

<td style="text-align:right;">

138

</td>

<td style="text-align:right;">

185

</td>

<td style="text-align:left;">

hsa-mir-210

</td>

<td style="text-align:left;">

MIR-210

</td>

<td style="text-align:left;">

UGUGCGU

</td>

<td style="text-align:left;">

chr11

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-136\_3p

</td>

<td style="text-align:right;">

\-0.89

</td>

<td style="text-align:right;">

0.27

</td>

<td style="text-align:left;">

2.95e-03

</td>

<td style="text-align:right;">

88

</td>

<td style="text-align:right;">

127

</td>

<td style="text-align:left;">

hsa-mir-136

</td>

<td style="text-align:left;">

MIR-136

</td>

<td style="text-align:left;">

AUCAUCG

</td>

<td style="text-align:left;">

chr14

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-21\_5p

</td>

<td style="text-align:right;">

\-0.86

</td>

<td style="text-align:right;">

0.22

</td>

<td style="text-align:left;">

6.21e-04

</td>

<td style="text-align:right;">

21349

</td>

<td style="text-align:right;">

28173

</td>

<td style="text-align:left;">

hsa-mir-21

</td>

<td style="text-align:left;">

MIR-21

</td>

<td style="text-align:left;">

AGCUUAU

</td>

<td style="text-align:left;">

chr17

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-30-P1b\_5p

</td>

<td style="text-align:right;">

\-0.86

</td>

<td style="text-align:right;">

0.18

</td>

<td style="text-align:left;">

8.68e-06

</td>

<td style="text-align:right;">

4513

</td>

<td style="text-align:right;">

6046

</td>

<td style="text-align:left;">

hsa-mir-30e

</td>

<td style="text-align:left;">

MIR-30

</td>

<td style="text-align:left;">

GUAAACA

</td>

<td style="text-align:left;">

chr1

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-191\_5p

</td>

<td style="text-align:right;">

\-0.82

</td>

<td style="text-align:right;">

0.28

</td>

<td style="text-align:left;">

1.02e-02

</td>

<td style="text-align:right;">

12597

</td>

<td style="text-align:right;">

14608

</td>

<td style="text-align:left;">

hsa-mir-191

</td>

<td style="text-align:left;">

MIR-191

</td>

<td style="text-align:left;">

AACGGAA

</td>

<td style="text-align:left;">

chr3

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-148-P1\_3p

</td>

<td style="text-align:right;">

\-0.72

</td>

<td style="text-align:right;">

0.30

</td>

<td style="text-align:left;">

4.78e-02

</td>

<td style="text-align:right;">

18874

</td>

<td style="text-align:right;">

21730

</td>

<td style="text-align:left;">

hsa-mir-148a

</td>

<td style="text-align:left;">

MIR-148

</td>

<td style="text-align:left;">

CAGUGCA

</td>

<td style="text-align:left;">

chr7

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

</tbody>

</table>

</div>

``` r
# Number of downregulated miRNA
signature_mirnas$number_downregulated
```

    ## [1] 35

## nCR vs pCRC

``` r
column='tissue.type'
tissue_type_A <- 'tumor.colorect'
tissue_type_B <- 'normal.colorect'
norm_adj_up       = "None"
norm_adj_down     = "None"
pCRC_adj_up   = "None"
pCRC_adj_down = "None"

coef <- paste(column, tissue_type_A, 'vs', tissue_type_B, sep='_')
res <- DeseqResult(dds, column, coef, tissue_type_A, tissue_type_B,
                   lfc.Threshold, rpm.Threshold,
                   norm_adj_up,
                   norm_adj_down)
dict_sig_mirna[paste(coef, "up",   sep='_')] <- list(res$up_mirna)
dict_sig_mirna[paste(coef, "down", sep='_')] <- list(res$down_mirna)
res_res <- res$res
res_dict[coef] <- res_res
plotMA(res$res, alpha=0.05)
```

![](comet_analysis_report_automated_no_lfcthreshold_08.12.20_res_alpha_0.05_with_LFC_shrinkage_LFC_0.58_mirge3_7.0_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

``` r
# Plot volcano plot
VolcanoPlot(res$res, coef, res$sig,
            res$up_mirna, res$down_mirna,
            norm_adj_up, norm_adj_down,
            pCRC_adj_up, pCRC_adj_down)
```

![](comet_analysis_report_automated_no_lfcthreshold_08.12.20_res_alpha_0.05_with_LFC_shrinkage_LFC_0.58_mirge3_7.0_files/figure-gfm/unnamed-chunk-17-2.png)<!-- -->

``` r
ExpressionPlot(res$res, res$rpm, coef, res$sig,
               tissue_type_A, tissue_type_B,
               res$up_mirna, res$down_mirna,
               norm_adj_up, norm_adj_down,
               pCRC_adj_up, pCRC_adj_down)
```

![](comet_analysis_report_automated_no_lfcthreshold_08.12.20_res_alpha_0.05_with_LFC_shrinkage_LFC_0.58_mirge3_7.0_files/figure-gfm/unnamed-chunk-17-3.png)<!-- -->

``` r
signature_mirnas <- SigList(res, dds, tissue_type_A, tissue_type_B, coef,
                            norm_adj_up, norm_adj_down, 
                            pCRC_adj_up, pCRC_adj_down)
# Print list upregulated miRNA
signature_mirnas$up_mirna
```

<div style="border: 1px solid #ddd; padding: 5px; overflow-x: scroll; width:2000px; ">

<table class="table table-striped table-hover table-condensed" style="margin-left: auto; margin-right: auto;">

<caption>

Upregulated in tissue.type\_tumor.colorect\_vs\_normal.colorect

</caption>

<thead>

<tr>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

miRNA

</th>

<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">

LFC

</th>

<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">

lfcSE

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

FDR

</th>

<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">

RPM tumor.colorect

</th>

<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">

RPM normal.colorect

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

miRBase\_ID

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

Family

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

Seed

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

Chr

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

Cell-Type Specific

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

Norm Background

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

pCRC Background

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

Hsa-Mir-17-P1a\_5p/P1b\_5p

</td>

<td style="text-align:right;">

1.42

</td>

<td style="text-align:right;">

0.14

</td>

<td style="text-align:left;">

1.30e-22

</td>

<td style="text-align:right;">

900

</td>

<td style="text-align:right;">

225

</td>

<td style="text-align:left;">

hsa-mir-17

</td>

<td style="text-align:left;">

MIR-17

</td>

<td style="text-align:left;">

AAAGUGC

</td>

<td style="text-align:left;">

chr13

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">CD14+
Monocyte</span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-7-P1\_5p/P2\_5p/P3\_5p

</td>

<td style="text-align:right;">

2.32

</td>

<td style="text-align:right;">

0.25

</td>

<td style="text-align:left;">

8.96e-18

</td>

<td style="text-align:right;">

199

</td>

<td style="text-align:right;">

22

</td>

<td style="text-align:left;">

hsa-mir-7-1

</td>

<td style="text-align:left;">

MIR-7

</td>

<td style="text-align:left;">

GGAAGAC

</td>

<td style="text-align:left;">

chr9

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">c(“Islet
Cell”, “Neural”)</span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-223\_3p

</td>

<td style="text-align:right;">

1.22

</td>

<td style="text-align:right;">

0.24

</td>

<td style="text-align:left;">

2.50e-06

</td>

<td style="text-align:right;">

619

</td>

<td style="text-align:right;">

177

</td>

<td style="text-align:left;">

hsa-mir-223

</td>

<td style="text-align:left;">

MIR-223

</td>

<td style="text-align:left;">

GUCAGUU

</td>

<td style="text-align:left;">

chrX

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">c(“Dendritic
Cell”, “Macrophage”)</span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-31\_5p

</td>

<td style="text-align:right;">

4.32

</td>

<td style="text-align:right;">

0.33

</td>

<td style="text-align:left;">

1.40e-44

</td>

<td style="text-align:right;">

627

</td>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

hsa-mir-31

</td>

<td style="text-align:left;">

MIR-31

</td>

<td style="text-align:left;">

GGCAAGA

</td>

<td style="text-align:left;">

chr9

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-135-P3\_5p

</td>

<td style="text-align:right;">

4.10

</td>

<td style="text-align:right;">

0.23

</td>

<td style="text-align:left;">

1.68e-69

</td>

<td style="text-align:right;">

155

</td>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

hsa-mir-135b

</td>

<td style="text-align:left;">

MIR-135

</td>

<td style="text-align:left;">

AUGGCUU

</td>

<td style="text-align:left;">

chr1

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-224\_5p

</td>

<td style="text-align:right;">

2.48

</td>

<td style="text-align:right;">

0.19

</td>

<td style="text-align:left;">

2.72e-36

</td>

<td style="text-align:right;">

395

</td>

<td style="text-align:right;">

43

</td>

<td style="text-align:left;">

hsa-mir-224

</td>

<td style="text-align:left;">

MIR-224

</td>

<td style="text-align:left;">

AAGUCAC

</td>

<td style="text-align:left;">

chrX

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-584\_5p

</td>

<td style="text-align:right;">

1.94

</td>

<td style="text-align:right;">

0.21

</td>

<td style="text-align:left;">

2.16e-19

</td>

<td style="text-align:right;">

110

</td>

<td style="text-align:right;">

17

</td>

<td style="text-align:left;">

hsa-mir-584

</td>

<td style="text-align:left;">

MIR-584

</td>

<td style="text-align:left;">

UAUGGUU

</td>

<td style="text-align:left;">

chr5

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-15-P1d\_5p

</td>

<td style="text-align:right;">

1.87

</td>

<td style="text-align:right;">

0.21

</td>

<td style="text-align:left;">

1.15e-17

</td>

<td style="text-align:right;">

235

</td>

<td style="text-align:right;">

40

</td>

<td style="text-align:left;">

hsa-mir-424

</td>

<td style="text-align:left;">

MIR-15

</td>

<td style="text-align:left;">

AGCAGCA

</td>

<td style="text-align:left;">

chrX

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-96-P2\_5p

</td>

<td style="text-align:right;">

1.84

</td>

<td style="text-align:right;">

0.17

</td>

<td style="text-align:left;">

6.44e-24

</td>

<td style="text-align:right;">

10417

</td>

<td style="text-align:right;">

1867

</td>

<td style="text-align:left;">

hsa-mir-182

</td>

<td style="text-align:left;">

MIR-96

</td>

<td style="text-align:left;">

UUGGCAA

</td>

<td style="text-align:left;">

chr7

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-17-P3a\_5p

</td>

<td style="text-align:right;">

1.62

</td>

<td style="text-align:right;">

0.16

</td>

<td style="text-align:left;">

2.85e-21

</td>

<td style="text-align:right;">

1513

</td>

<td style="text-align:right;">

322

</td>

<td style="text-align:left;">

hsa-mir-20a

</td>

<td style="text-align:left;">

MIR-17

</td>

<td style="text-align:left;">

AAAGUGC

</td>

<td style="text-align:left;">

chr13

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-96-P3\_5p

</td>

<td style="text-align:right;">

1.53

</td>

<td style="text-align:right;">

0.23

</td>

<td style="text-align:left;">

1.96e-10

</td>

<td style="text-align:right;">

1063

</td>

<td style="text-align:right;">

224

</td>

<td style="text-align:left;">

hsa-mir-183

</td>

<td style="text-align:left;">

MIR-96

</td>

<td style="text-align:left;">

AUGGCAC

</td>

<td style="text-align:left;">

chr7

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-19-P1\_3p

</td>

<td style="text-align:right;">

1.47

</td>

<td style="text-align:right;">

0.17

</td>

<td style="text-align:left;">

1.43e-17

</td>

<td style="text-align:right;">

379

</td>

<td style="text-align:right;">

87

</td>

<td style="text-align:left;">

hsa-mir-19a

</td>

<td style="text-align:left;">

MIR-19

</td>

<td style="text-align:left;">

GUGCAAA

</td>

<td style="text-align:left;">

chr13

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-21\_5p

</td>

<td style="text-align:right;">

1.35

</td>

<td style="text-align:right;">

0.13

</td>

<td style="text-align:left;">

1.01e-24

</td>

<td style="text-align:right;">

105735

</td>

<td style="text-align:right;">

28173

</td>

<td style="text-align:left;">

hsa-mir-21

</td>

<td style="text-align:left;">

MIR-21

</td>

<td style="text-align:left;">

AGCUUAU

</td>

<td style="text-align:left;">

chr17

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-17-P2a\_5p

</td>

<td style="text-align:right;">

1.33

</td>

<td style="text-align:right;">

0.19

</td>

<td style="text-align:left;">

6.62e-11

</td>

<td style="text-align:right;">

168

</td>

<td style="text-align:right;">

44

</td>

<td style="text-align:left;">

hsa-mir-18a

</td>

<td style="text-align:left;">

MIR-17

</td>

<td style="text-align:left;">

AAGGUGC

</td>

<td style="text-align:left;">

chr13

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-181-P2c\_5p

</td>

<td style="text-align:right;">

1.33

</td>

<td style="text-align:right;">

0.16

</td>

<td style="text-align:left;">

1.40e-15

</td>

<td style="text-align:right;">

196

</td>

<td style="text-align:right;">

51

</td>

<td style="text-align:left;">

hsa-mir-181d

</td>

<td style="text-align:left;">

MIR-181

</td>

<td style="text-align:left;">

ACAUUCA

</td>

<td style="text-align:left;">

chr19

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-95-P2\_3p

</td>

<td style="text-align:right;">

1.17

</td>

<td style="text-align:right;">

0.14

</td>

<td style="text-align:left;">

1.40e-15

</td>

<td style="text-align:right;">

142

</td>

<td style="text-align:right;">

41

</td>

<td style="text-align:left;">

hsa-mir-421

</td>

<td style="text-align:left;">

MIR-95

</td>

<td style="text-align:left;">

UCAACAG

</td>

<td style="text-align:left;">

chrX

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-19-P2a\_3p/P2b\_3p

</td>

<td style="text-align:right;">

1.14

</td>

<td style="text-align:right;">

0.16

</td>

<td style="text-align:left;">

1.59e-12

</td>

<td style="text-align:right;">

1297

</td>

<td style="text-align:right;">

389

</td>

<td style="text-align:left;">

hsa-mir-19b-1

</td>

<td style="text-align:left;">

MIR-19

</td>

<td style="text-align:left;">

GUGCAAA

</td>

<td style="text-align:left;">

chr13

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-181-P1c\_5p

</td>

<td style="text-align:right;">

1.11

</td>

<td style="text-align:right;">

0.15

</td>

<td style="text-align:left;">

1.06e-12

</td>

<td style="text-align:right;">

1724

</td>

<td style="text-align:right;">

527

</td>

<td style="text-align:left;">

hsa-mir-181c

</td>

<td style="text-align:left;">

MIR-181

</td>

<td style="text-align:left;">

ACAUUCA

</td>

<td style="text-align:left;">

chr19

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-29-P2a\_3p/P2b\_3p

</td>

<td style="text-align:right;">

1.08

</td>

<td style="text-align:right;">

0.17

</td>

<td style="text-align:left;">

5.43e-10

</td>

<td style="text-align:right;">

380

</td>

<td style="text-align:right;">

124

</td>

<td style="text-align:left;">

hsa-mir-29b-1

</td>

<td style="text-align:left;">

MIR-29

</td>

<td style="text-align:left;">

AGCACCA

</td>

<td style="text-align:left;">

chr7

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-130-P2a\_3p

</td>

<td style="text-align:right;">

1.05

</td>

<td style="text-align:right;">

0.15

</td>

<td style="text-align:left;">

8.09e-12

</td>

<td style="text-align:right;">

409

</td>

<td style="text-align:right;">

131

</td>

<td style="text-align:left;">

hsa-mir-301a

</td>

<td style="text-align:left;">

MIR-130

</td>

<td style="text-align:left;">

AGUGCAA

</td>

<td style="text-align:left;">

chr17

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-17-P3c\_5p

</td>

<td style="text-align:right;">

0.97

</td>

<td style="text-align:right;">

0.12

</td>

<td style="text-align:left;">

1.00e-15

</td>

<td style="text-align:right;">

370

</td>

<td style="text-align:right;">

130

</td>

<td style="text-align:left;">

hsa-mir-106b

</td>

<td style="text-align:left;">

MIR-17

</td>

<td style="text-align:left;">

AAAGUGC

</td>

<td style="text-align:left;">

chr7

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-221-P2\_3p

</td>

<td style="text-align:right;">

0.94

</td>

<td style="text-align:right;">

0.15

</td>

<td style="text-align:left;">

1.32e-09

</td>

<td style="text-align:right;">

1633

</td>

<td style="text-align:right;">

574

</td>

<td style="text-align:left;">

hsa-mir-222

</td>

<td style="text-align:left;">

MIR-221

</td>

<td style="text-align:left;">

GCUACAU

</td>

<td style="text-align:left;">

chrX

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-92-P1a\_3p/P1b\_3p

</td>

<td style="text-align:right;">

0.90

</td>

<td style="text-align:right;">

0.16

</td>

<td style="text-align:left;">

5.66e-08

</td>

<td style="text-align:right;">

39320

</td>

<td style="text-align:right;">

13849

</td>

<td style="text-align:left;">

hsa-mir-92a-1

</td>

<td style="text-align:left;">

MIR-92

</td>

<td style="text-align:left;">

AUUGCAC

</td>

<td style="text-align:left;">

chr13

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-17-P1c\_5p

</td>

<td style="text-align:right;">

0.88

</td>

<td style="text-align:right;">

0.10

</td>

<td style="text-align:left;">

2.69e-17

</td>

<td style="text-align:right;">

2287

</td>

<td style="text-align:right;">

843

</td>

<td style="text-align:left;">

hsa-mir-93

</td>

<td style="text-align:left;">

MIR-17

</td>

<td style="text-align:left;">

AAAGUGC

</td>

<td style="text-align:left;">

chr7

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-221-P1\_3p

</td>

<td style="text-align:right;">

0.83

</td>

<td style="text-align:right;">

0.13

</td>

<td style="text-align:left;">

4.23e-09

</td>

<td style="text-align:right;">

2529

</td>

<td style="text-align:right;">

964

</td>

<td style="text-align:left;">

hsa-mir-221

</td>

<td style="text-align:left;">

MIR-221

</td>

<td style="text-align:left;">

GCUACAU

</td>

<td style="text-align:left;">

chrX

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-203\_3p

</td>

<td style="text-align:right;">

0.74

</td>

<td style="text-align:right;">

0.20

</td>

<td style="text-align:left;">

9.13e-04

</td>

<td style="text-align:right;">

2307

</td>

<td style="text-align:right;">

922

</td>

<td style="text-align:left;">

hsa-mir-203a

</td>

<td style="text-align:left;">

MIR-203

</td>

<td style="text-align:left;">

UGAAAUG

</td>

<td style="text-align:left;">

chr14

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Let-7-P2c2\_5p

</td>

<td style="text-align:right;">

0.74

</td>

<td style="text-align:right;">

0.12

</td>

<td style="text-align:left;">

1.50e-08

</td>

<td style="text-align:right;">

5432

</td>

<td style="text-align:right;">

2218

</td>

<td style="text-align:left;">

hsa-let-7i

</td>

<td style="text-align:left;">

LET-7

</td>

<td style="text-align:left;">

GAGGUAG

</td>

<td style="text-align:left;">

chr12

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Let-7-P2b3\_5p

</td>

<td style="text-align:right;">

0.68

</td>

<td style="text-align:right;">

0.14

</td>

<td style="text-align:left;">

4.33e-06

</td>

<td style="text-align:right;">

2149

</td>

<td style="text-align:right;">

937

</td>

<td style="text-align:left;">

hsa-mir-98

</td>

<td style="text-align:left;">

LET-7

</td>

<td style="text-align:left;">

GAGGUAG

</td>

<td style="text-align:left;">

chrX

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-29-P1a\_3p

</td>

<td style="text-align:right;">

0.62

</td>

<td style="text-align:right;">

0.12

</td>

<td style="text-align:left;">

1.40e-06

</td>

<td style="text-align:right;">

3667

</td>

<td style="text-align:right;">

1692

</td>

<td style="text-align:left;">

hsa-mir-29a

</td>

<td style="text-align:left;">

MIR-29

</td>

<td style="text-align:left;">

AGCACCA

</td>

<td style="text-align:left;">

chr7

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-92-P2c\_3p

</td>

<td style="text-align:right;">

0.62

</td>

<td style="text-align:right;">

0.10

</td>

<td style="text-align:left;">

4.70e-09

</td>

<td style="text-align:right;">

3274

</td>

<td style="text-align:right;">

1443

</td>

<td style="text-align:left;">

hsa-mir-25

</td>

<td style="text-align:left;">

MIR-92

</td>

<td style="text-align:left;">

AUUGCAC

</td>

<td style="text-align:left;">

chr7

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-92-P1c\_3p

</td>

<td style="text-align:right;">

0.60

</td>

<td style="text-align:right;">

0.18

</td>

<td style="text-align:left;">

4.90e-03

</td>

<td style="text-align:right;">

1398

</td>

<td style="text-align:right;">

634

</td>

<td style="text-align:left;">

hsa-mir-92b

</td>

<td style="text-align:left;">

MIR-92

</td>

<td style="text-align:left;">

AUUGCAC

</td>

<td style="text-align:left;">

chr1

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-769\_5p

</td>

<td style="text-align:right;">

0.59

</td>

<td style="text-align:right;">

0.11

</td>

<td style="text-align:left;">

1.18e-06

</td>

<td style="text-align:right;">

439

</td>

<td style="text-align:right;">

195

</td>

<td style="text-align:left;">

hsa-mir-769

</td>

<td style="text-align:left;">

MIR-769

</td>

<td style="text-align:left;">

GAGACCU

</td>

<td style="text-align:left;">

chr19

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

</tbody>

</table>

</div>

``` r
# Number of upregulated miRNA
signature_mirnas$number_upregulated
```

    ## [1] 32

``` r
# Print list downregulated miRNA
signature_mirnas$down_mirna
```

<div style="border: 1px solid #ddd; padding: 5px; overflow-x: scroll; width:2000px; ">

<table class="table table-striped table-hover table-condensed" style="margin-left: auto; margin-right: auto;">

<caption>

Downregulated in tissue.type\_tumor.colorect\_vs\_normal.colorect

</caption>

<thead>

<tr>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

miRNA

</th>

<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">

LFC

</th>

<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">

lfcSE

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

FDR

</th>

<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">

RPM tumor.colorect

</th>

<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">

RPM normal.colorect

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

miRBase\_ID

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

Family

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

Seed

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

Chr

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

Cell-Type Specific

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

Norm Background

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

pCRC Background

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

Hsa-Mir-451\_5p

</td>

<td style="text-align:right;">

\-1.19

</td>

<td style="text-align:right;">

0.26

</td>

<td style="text-align:left;">

3.69e-05

</td>

<td style="text-align:right;">

1565

</td>

<td style="text-align:right;">

2815

</td>

<td style="text-align:left;">

hsa-mir-451a

</td>

<td style="text-align:left;">

MIR-451

</td>

<td style="text-align:left;">

AACCGUU

</td>

<td style="text-align:left;">

chr17

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">Red
Blood Cell</span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-145\_5p

</td>

<td style="text-align:right;">

\-1.97

</td>

<td style="text-align:right;">

0.22

</td>

<td style="text-align:left;">

1.63e-17

</td>

<td style="text-align:right;">

840

</td>

<td style="text-align:right;">

2832

</td>

<td style="text-align:left;">

hsa-mir-145

</td>

<td style="text-align:left;">

MIR-145

</td>

<td style="text-align:left;">

UCCAGUU

</td>

<td style="text-align:left;">

chr5

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">Mesenchymal</span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-143\_3p

</td>

<td style="text-align:right;">

\-0.60

</td>

<td style="text-align:right;">

0.19

</td>

<td style="text-align:left;">

1.39e-03

</td>

<td style="text-align:right;">

148431

</td>

<td style="text-align:right;">

164812

</td>

<td style="text-align:left;">

hsa-mir-143

</td>

<td style="text-align:left;">

MIR-143

</td>

<td style="text-align:left;">

GAGAUGA

</td>

<td style="text-align:left;">

chr5

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">Mesenchymal</span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-150\_5p

</td>

<td style="text-align:right;">

\-1.36

</td>

<td style="text-align:right;">

0.24

</td>

<td style="text-align:left;">

3.58e-07

</td>

<td style="text-align:right;">

280

</td>

<td style="text-align:right;">

586

</td>

<td style="text-align:left;">

hsa-mir-150

</td>

<td style="text-align:left;">

MIR-150

</td>

<td style="text-align:left;">

CUCCCAA

</td>

<td style="text-align:left;">

chr19

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">Lymphocyte</span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-192-P1\_5p/P2\_5p

</td>

<td style="text-align:right;">

\-2.05

</td>

<td style="text-align:right;">

0.22

</td>

<td style="text-align:left;">

3.04e-19

</td>

<td style="text-align:right;">

13717

</td>

<td style="text-align:right;">

46010

</td>

<td style="text-align:left;">

hsa-mir-192

</td>

<td style="text-align:left;">

MIR-192

</td>

<td style="text-align:left;">

UGACCUA

</td>

<td style="text-align:left;">

chr11

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">Epithelial
Cell</span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-133-P1\_3p/P2\_3p/P3\_3p

</td>

<td style="text-align:right;">

\-1.98

</td>

<td style="text-align:right;">

0.26

</td>

<td style="text-align:left;">

6.98e-16

</td>

<td style="text-align:right;">

116

</td>

<td style="text-align:right;">

460

</td>

<td style="text-align:left;">

hsa-mir-133a-2

</td>

<td style="text-align:left;">

MIR-133

</td>

<td style="text-align:left;">

UUGGUCC

</td>

<td style="text-align:left;">

chr20

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">c(“Skeletal
Myocyte”, “Stem Cell”)</span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-486\_5p

</td>

<td style="text-align:right;">

\-1.18

</td>

<td style="text-align:right;">

0.25

</td>

<td style="text-align:left;">

2.02e-05

</td>

<td style="text-align:right;">

2013

</td>

<td style="text-align:right;">

3961

</td>

<td style="text-align:left;">

hsa-mir-486-1

</td>

<td style="text-align:left;">

MIR-486

</td>

<td style="text-align:left;">

CCUGUAC

</td>

<td style="text-align:left;">

chr8

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">c(“Platelet”,
“Red Blood Cell”)</span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-375\_3p

</td>

<td style="text-align:right;">

\-1.86

</td>

<td style="text-align:right;">

0.26

</td>

<td style="text-align:left;">

2.92e-12

</td>

<td style="text-align:right;">

5263

</td>

<td style="text-align:right;">

17215

</td>

<td style="text-align:left;">

hsa-mir-375

</td>

<td style="text-align:left;">

MIR-375

</td>

<td style="text-align:left;">

UUGUUCG

</td>

<td style="text-align:left;">

chr2

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">c(“Epithelial
Cell”, “Islet Cell”, “Neural”)</span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-126\_5p

</td>

<td style="text-align:right;">

\-0.72

</td>

<td style="text-align:right;">

0.15

</td>

<td style="text-align:left;">

1.61e-05

</td>

<td style="text-align:right;">

3455

</td>

<td style="text-align:right;">

4230

</td>

<td style="text-align:left;">

hsa-mir-126

</td>

<td style="text-align:left;">

MIR-126

</td>

<td style="text-align:left;">

AUUAUUA

</td>

<td style="text-align:left;">

chr9

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">c(“Endothelial
Cell”, “Platelet”)</span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-342\_3p

</td>

<td style="text-align:right;">

\-1.26

</td>

<td style="text-align:right;">

0.18

</td>

<td style="text-align:left;">

1.59e-10

</td>

<td style="text-align:right;">

145

</td>

<td style="text-align:right;">

260

</td>

<td style="text-align:left;">

hsa-mir-342

</td>

<td style="text-align:left;">

MIR-342

</td>

<td style="text-align:left;">

CUCACAC

</td>

<td style="text-align:left;">

chr14

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">c(“Dendritic
Cell”, “Lymphocyte”, “Macrophage”)</span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-147\_3p

</td>

<td style="text-align:right;">

\-1.92

</td>

<td style="text-align:right;">

0.23

</td>

<td style="text-align:left;">

8.04e-17

</td>

<td style="text-align:right;">

33

</td>

<td style="text-align:right;">

110

</td>

<td style="text-align:left;">

hsa-mir-147b

</td>

<td style="text-align:left;">

MIR-147

</td>

<td style="text-align:left;">

UGUGCGG

</td>

<td style="text-align:left;">

chr15

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-15-P2c\_5p

</td>

<td style="text-align:right;">

\-1.76

</td>

<td style="text-align:right;">

0.17

</td>

<td style="text-align:left;">

6.42e-23

</td>

<td style="text-align:right;">

257

</td>

<td style="text-align:right;">

672

</td>

<td style="text-align:left;">

hsa-mir-195

</td>

<td style="text-align:left;">

MIR-15

</td>

<td style="text-align:left;">

AGCAGCA

</td>

<td style="text-align:left;">

chr17

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-378\_3p

</td>

<td style="text-align:right;">

\-1.72

</td>

<td style="text-align:right;">

0.14

</td>

<td style="text-align:left;">

8.99e-34

</td>

<td style="text-align:right;">

6426

</td>

<td style="text-align:right;">

14834

</td>

<td style="text-align:left;">

hsa-mir-378a

</td>

<td style="text-align:left;">

MIR-378

</td>

<td style="text-align:left;">

CUGGACU

</td>

<td style="text-align:left;">

chr5

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-15-P1c\_5p

</td>

<td style="text-align:right;">

\-1.62

</td>

<td style="text-align:right;">

0.16

</td>

<td style="text-align:left;">

5.86e-22

</td>

<td style="text-align:right;">

143

</td>

<td style="text-align:right;">

337

</td>

<td style="text-align:left;">

hsa-mir-497

</td>

<td style="text-align:left;">

MIR-15

</td>

<td style="text-align:left;">

AGCAGCA

</td>

<td style="text-align:left;">

chr17

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-190-P1\_5p

</td>

<td style="text-align:right;">

\-1.47

</td>

<td style="text-align:right;">

0.18

</td>

<td style="text-align:left;">

1.43e-15

</td>

<td style="text-align:right;">

140

</td>

<td style="text-align:right;">

294

</td>

<td style="text-align:left;">

hsa-mir-190a

</td>

<td style="text-align:left;">

MIR-190

</td>

<td style="text-align:left;">

GAUAUGU

</td>

<td style="text-align:left;">

chr15

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-194-P1\_5p/P2\_5p

</td>

<td style="text-align:right;">

\-1.31

</td>

<td style="text-align:right;">

0.19

</td>

<td style="text-align:left;">

3.24e-11

</td>

<td style="text-align:right;">

6167

</td>

<td style="text-align:right;">

11619

</td>

<td style="text-align:left;">

hsa-mir-194-2

</td>

<td style="text-align:left;">

MIR-194

</td>

<td style="text-align:left;">

GUAACAG

</td>

<td style="text-align:left;">

chr11

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-30-P1a\_5p

</td>

<td style="text-align:right;">

\-1.21

</td>

<td style="text-align:right;">

0.16

</td>

<td style="text-align:left;">

7.30e-12

</td>

<td style="text-align:right;">

3165

</td>

<td style="text-align:right;">

5150

</td>

<td style="text-align:left;">

hsa-mir-30a

</td>

<td style="text-align:left;">

MIR-30

</td>

<td style="text-align:left;">

GUAAACA

</td>

<td style="text-align:left;">

chr6

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-26-P3\_5p

</td>

<td style="text-align:right;">

\-1.21

</td>

<td style="text-align:right;">

0.14

</td>

<td style="text-align:left;">

3.12e-16

</td>

<td style="text-align:right;">

3467

</td>

<td style="text-align:right;">

5898

</td>

<td style="text-align:left;">

hsa-mir-26a-1

</td>

<td style="text-align:left;">

MIR-26

</td>

<td style="text-align:left;">

UCAAGUA

</td>

<td style="text-align:left;">

chr3

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-338-P1\_3p

</td>

<td style="text-align:right;">

\-1.18

</td>

<td style="text-align:right;">

0.22

</td>

<td style="text-align:left;">

8.59e-07

</td>

<td style="text-align:right;">

98

</td>

<td style="text-align:right;">

166

</td>

<td style="text-align:left;">

hsa-mir-338

</td>

<td style="text-align:left;">

MIR-338

</td>

<td style="text-align:left;">

CCAGCAU

</td>

<td style="text-align:left;">

chr17

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-26-P1\_5p/P2\_5p

</td>

<td style="text-align:right;">

\-1.15

</td>

<td style="text-align:right;">

0.10

</td>

<td style="text-align:left;">

4.50e-29

</td>

<td style="text-align:right;">

35584

</td>

<td style="text-align:right;">

56268

</td>

<td style="text-align:left;">

hsa-mir-26b

</td>

<td style="text-align:left;">

MIR-26

</td>

<td style="text-align:left;">

UCAAGUA

</td>

<td style="text-align:left;">

chr2

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-10-P1b\_5p

</td>

<td style="text-align:right;">

\-1.11

</td>

<td style="text-align:right;">

0.25

</td>

<td style="text-align:left;">

4.63e-06

</td>

<td style="text-align:right;">

39286

</td>

<td style="text-align:right;">

75067

</td>

<td style="text-align:left;">

hsa-mir-10b

</td>

<td style="text-align:left;">

MIR-10

</td>

<td style="text-align:left;">

ACCCUGU

</td>

<td style="text-align:left;">

chr2

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-192-P1\_5p

</td>

<td style="text-align:right;">

\-1.11

</td>

<td style="text-align:right;">

0.20

</td>

<td style="text-align:left;">

2.37e-08

</td>

<td style="text-align:right;">

132000

</td>

<td style="text-align:right;">

212293

</td>

<td style="text-align:left;">

hsa-mir-192

</td>

<td style="text-align:left;">

MIR-192

</td>

<td style="text-align:left;">

UGACCUA

</td>

<td style="text-align:left;">

chr11

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-148-P2\_3p

</td>

<td style="text-align:right;">

\-1.04

</td>

<td style="text-align:right;">

0.21

</td>

<td style="text-align:left;">

5.54e-06

</td>

<td style="text-align:right;">

137

</td>

<td style="text-align:right;">

210

</td>

<td style="text-align:left;">

hsa-mir-152

</td>

<td style="text-align:left;">

MIR-148

</td>

<td style="text-align:left;">

CAGUGCA

</td>

<td style="text-align:left;">

chr17

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-15-P2a\_5p/P2b\_5p

</td>

<td style="text-align:right;">

\-1.02

</td>

<td style="text-align:right;">

0.12

</td>

<td style="text-align:left;">

8.21e-15

</td>

<td style="text-align:right;">

4522

</td>

<td style="text-align:right;">

6926

</td>

<td style="text-align:left;">

hsa-mir-16-1

</td>

<td style="text-align:left;">

MIR-15

</td>

<td style="text-align:left;">

AGCAGCA

</td>

<td style="text-align:left;">

chr13

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-10-P3b\_5p

</td>

<td style="text-align:right;">

\-0.97

</td>

<td style="text-align:right;">

0.19

</td>

<td style="text-align:left;">

4.09e-06

</td>

<td style="text-align:right;">

1649

</td>

<td style="text-align:right;">

2392

</td>

<td style="text-align:left;">

hsa-mir-125a

</td>

<td style="text-align:left;">

MIR-10

</td>

<td style="text-align:left;">

CCCUGAG

</td>

<td style="text-align:left;">

chr19

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-28-P1\_3p

</td>

<td style="text-align:right;">

\-0.95

</td>

<td style="text-align:right;">

0.13

</td>

<td style="text-align:left;">

6.44e-12

</td>

<td style="text-align:right;">

2915

</td>

<td style="text-align:right;">

3956

</td>

<td style="text-align:left;">

hsa-mir-28

</td>

<td style="text-align:left;">

MIR-28

</td>

<td style="text-align:left;">

ACUAGAU

</td>

<td style="text-align:left;">

chr3

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-1307\_5p

</td>

<td style="text-align:right;">

\-0.95

</td>

<td style="text-align:right;">

0.18

</td>

<td style="text-align:left;">

2.27e-06

</td>

<td style="text-align:right;">

492

</td>

<td style="text-align:right;">

713

</td>

<td style="text-align:left;">

hsa-mir-1307

</td>

<td style="text-align:left;">

MIR-1307

</td>

<td style="text-align:left;">

CGACCGG

</td>

<td style="text-align:left;">

chr10

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-29-P1b\_3p

</td>

<td style="text-align:right;">

\-0.93

</td>

<td style="text-align:right;">

0.15

</td>

<td style="text-align:left;">

1.11e-09

</td>

<td style="text-align:right;">

331

</td>

<td style="text-align:right;">

464

</td>

<td style="text-align:left;">

hsa-mir-29c

</td>

<td style="text-align:left;">

MIR-29

</td>

<td style="text-align:left;">

AGCACCA

</td>

<td style="text-align:left;">

chr1

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-10-P3c\_5p

</td>

<td style="text-align:right;">

\-0.80

</td>

<td style="text-align:right;">

0.21

</td>

<td style="text-align:left;">

1.44e-03

</td>

<td style="text-align:right;">

239

</td>

<td style="text-align:right;">

326

</td>

<td style="text-align:left;">

hsa-mir-125b-2

</td>

<td style="text-align:left;">

MIR-10

</td>

<td style="text-align:left;">

CCCUGAG

</td>

<td style="text-align:left;">

chr21

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-574\_3p

</td>

<td style="text-align:right;">

\-0.78

</td>

<td style="text-align:right;">

0.14

</td>

<td style="text-align:left;">

1.83e-07

</td>

<td style="text-align:right;">

240

</td>

<td style="text-align:right;">

305

</td>

<td style="text-align:left;">

hsa-mir-574

</td>

<td style="text-align:left;">

MIR-574

</td>

<td style="text-align:left;">

ACGCUCA

</td>

<td style="text-align:left;">

chr4

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-8-P1b\_3p

</td>

<td style="text-align:right;">

\-0.69

</td>

<td style="text-align:right;">

0.21

</td>

<td style="text-align:left;">

2.87e-03

</td>

<td style="text-align:right;">

7518

</td>

<td style="text-align:right;">

9180

</td>

<td style="text-align:left;">

hsa-mir-141

</td>

<td style="text-align:left;">

MIR-8

</td>

<td style="text-align:left;">

AACACUG

</td>

<td style="text-align:left;">

chr12

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-191\_5p

</td>

<td style="text-align:right;">

\-0.68

</td>

<td style="text-align:right;">

0.16

</td>

<td style="text-align:left;">

1.46e-04

</td>

<td style="text-align:right;">

13630

</td>

<td style="text-align:right;">

14608

</td>

<td style="text-align:left;">

hsa-mir-191

</td>

<td style="text-align:left;">

MIR-191

</td>

<td style="text-align:left;">

AACGGAA

</td>

<td style="text-align:left;">

chr3

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-362-P3\_3p

</td>

<td style="text-align:right;">

\-0.62

</td>

<td style="text-align:right;">

0.25

</td>

<td style="text-align:left;">

1.92e-02

</td>

<td style="text-align:right;">

108

</td>

<td style="text-align:right;">

103

</td>

<td style="text-align:left;">

hsa-mir-501

</td>

<td style="text-align:left;">

MIR-362

</td>

<td style="text-align:left;">

AUGCACC

</td>

<td style="text-align:left;">

chrX

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-142\_5p

</td>

<td style="text-align:right;">

\-0.60

</td>

<td style="text-align:right;">

0.19

</td>

<td style="text-align:left;">

4.90e-03

</td>

<td style="text-align:right;">

3453

</td>

<td style="text-align:right;">

3817

</td>

<td style="text-align:left;">

hsa-mir-142

</td>

<td style="text-align:left;">

MIR-142

</td>

<td style="text-align:left;">

AUAAAGU

</td>

<td style="text-align:left;">

chr17

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-154-P9\_3p

</td>

<td style="text-align:right;">

\-0.59

</td>

<td style="text-align:right;">

0.15

</td>

<td style="text-align:left;">

1.81e-04

</td>

<td style="text-align:right;">

253

</td>

<td style="text-align:right;">

265

</td>

<td style="text-align:left;">

hsa-mir-381

</td>

<td style="text-align:left;">

MIR-154

</td>

<td style="text-align:left;">

AUACAAG

</td>

<td style="text-align:left;">

chr14

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

</tbody>

</table>

</div>

``` r
# Number of downregulated miRNA
signature_mirnas$number_downregulated
```

    ## [1] 35

``` r
ref <- 'tumor.colorect'
dds <- DeseqObject(design, countdata, sampleinfo, "None", "None", ref)
#
```

``` r
# #datasets in total
dim(dds[, colData(dds)$type.tissue == 'pCRC'])
```

    ## [1] 389 120

``` r
dim(dds[, colData(dds)$type.tissue == 'mLi'])
```

    ## [1] 389  35

``` r
dim(dds[, colData(dds)$type.tissue == 'mLu'])
```

    ## [1] 389  28

``` r
dim(dds[, colData(dds)$type.tissue == 'nCR'])
```

    ## [1] 389  25

``` r
dim(dds[, colData(dds)$type.tissue == 'nLi'])
```

    ## [1] 389  20

``` r
dim(dds[, colData(dds)$type.tissue == 'nLu'])
```

    ## [1] 389  10

``` r
dim(dds[, colData(dds)$type.tissue == 'PM'])
```

    ## [1] 389  30

``` r
# #datasets for Fromm
dim(dds[, colData(dds)$type.tissue == 'pCRC' & colData(dds)$paper == 'fromm'])
```

    ## [1] 389   3

``` r
dim(dds[, colData(dds)$type.tissue == 'mLi' & colData(dds)$paper == 'fromm'])
```

    ## [1] 389  19

``` r
dim(dds[, colData(dds)$type.tissue == 'mLu' & colData(dds)$paper == 'fromm'])
```

    ## [1] 389  24

``` r
dim(dds[, colData(dds)$type.tissue == 'nCR' & colData(dds)$paper == 'fromm'])
```

    ## [1] 389   3

``` r
dim(dds[, colData(dds)$type.tissue == 'nLi' & colData(dds)$paper == 'fromm'])
```

    ## [1] 389   8

``` r
dim(dds[, colData(dds)$type.tissue == 'nLu' & colData(dds)$paper == 'fromm'])
```

    ## [1] 389   7

``` r
dim(dds[, colData(dds)$type.tissue == 'PM' & colData(dds)$paper == 'fromm'])
```

    ## [1] 389  18

``` r
# #datasets for Schee
dim(dds[, colData(dds)$type.tissue == 'pCRC' & colData(dds)$paper == 'schee'])
```

    ## [1] 389  83

``` r
dim(dds[, colData(dds)$type.tissue == 'mLi' & colData(dds)$paper == 'schee'])
```

    ## [1] 389   0

``` r
dim(dds[, colData(dds)$type.tissue == 'mLu' & colData(dds)$paper == 'schee'])
```

    ## [1] 389   0

``` r
dim(dds[, colData(dds)$type.tissue == 'nCR' & colData(dds)$paper == 'schee'])
```

    ## [1] 389   0

``` r
dim(dds[, colData(dds)$type.tissue == 'nLi' & colData(dds)$paper == 'schee'])
```

    ## [1] 389   0

``` r
dim(dds[, colData(dds)$type.tissue == 'nLu' & colData(dds)$paper == 'schee'])
```

    ## [1] 389   0

``` r
dim(dds[, colData(dds)$type.tissue == 'PM' & colData(dds)$paper == 'schee'])
```

    ## [1] 389   0

``` r
# #datasets for Schee
dim(dds[, colData(dds)$type.tissue == 'pCRC' & colData(dds)$paper == 'neerincx'])
```

    ## [1] 389  34

``` r
dim(dds[, colData(dds)$type.tissue == 'mLi' & colData(dds)$paper == 'neerincx'])
```

    ## [1] 389  16

``` r
dim(dds[, colData(dds)$type.tissue == 'mLu' & colData(dds)$paper == 'neerincx'])
```

    ## [1] 389   4

``` r
dim(dds[, colData(dds)$type.tissue == 'nCR' & colData(dds)$paper == 'neerincx'])
```

    ## [1] 389  22

``` r
dim(dds[, colData(dds)$type.tissue == 'nLi' & colData(dds)$paper == 'neerincx'])
```

    ## [1] 389   9

``` r
dim(dds[, colData(dds)$type.tissue == 'nLu' & colData(dds)$paper == 'neerincx'])
```

    ## [1] 389   3

``` r
dim(dds[, colData(dds)$type.tissue == 'PM' & colData(dds)$paper == 'neerincx'])
```

    ## [1] 389  12

``` r
# #datasets for Schee
dim(dds[, colData(dds)$type.tissue == 'pCRC' & colData(dds)$paper == 'selitsky'])
```

    ## [1] 389   0

``` r
dim(dds[, colData(dds)$type.tissue == 'mLi' & colData(dds)$paper == 'selitsky'])
```

    ## [1] 389   0

``` r
dim(dds[, colData(dds)$type.tissue == 'mLu' & colData(dds)$paper == 'selitsky'])
```

    ## [1] 389   0

``` r
dim(dds[, colData(dds)$type.tissue == 'nCR' & colData(dds)$paper == 'selitsky'])
```

    ## [1] 389   0

``` r
dim(dds[, colData(dds)$type.tissue == 'nLi' & colData(dds)$paper == 'selitsky'])
```

    ## [1] 389   3

``` r
dim(dds[, colData(dds)$type.tissue == 'nLu' & colData(dds)$paper == 'selitsky'])
```

    ## [1] 389   0

``` r
dim(dds[, colData(dds)$type.tissue == 'PM' & colData(dds)$paper == 'selitsky'])
```

    ## [1] 389   0

``` r
plotDispEsts(dds)
```

![](comet_analysis_report_automated_no_lfcthreshold_08.12.20_res_alpha_0.05_with_LFC_shrinkage_LFC_0.58_mirge3_7.0_files/figure-gfm/unnamed-chunk-24-1.png)<!-- -->
\#\# pCRC vs nLi

``` r
column='tissue.type'
tissue_type_A <- 'normal.liver'
tissue_type_B <- 'tumor.colorect'
norm_adj_up       = "None"
norm_adj_down     = "None"
pCRC_adj_up   = "None"
pCRC_adj_down = "None"

coef <- paste(column, tissue_type_A, 'vs', tissue_type_B, sep='_')
res <- DeseqResult(dds, column, coef, tissue_type_A, tissue_type_B,
                   lfc.Threshold, rpm.Threshold,
                   norm_adj_up,
                   norm_adj_down)
dict_sig_mirna[paste(coef, "up",   sep='_')] <- list(res$up_mirna)
dict_sig_mirna[paste(coef, "down", sep='_')] <- list(res$down_mirna)
res_res <- res$res
res_dict[coef] <- res_res
plotMA(res$res, alpha=0.05)
```

![](comet_analysis_report_automated_no_lfcthreshold_08.12.20_res_alpha_0.05_with_LFC_shrinkage_LFC_0.58_mirge3_7.0_files/figure-gfm/unnamed-chunk-25-1.png)<!-- -->

``` r
# Plot volcano plot
VolcanoPlot(res$res, coef, res$sig,
            res$up_mirna, res$down_mirna,
            norm_adj_up, norm_adj_down,
            pCRC_adj_up, pCRC_adj_down)
```

![](comet_analysis_report_automated_no_lfcthreshold_08.12.20_res_alpha_0.05_with_LFC_shrinkage_LFC_0.58_mirge3_7.0_files/figure-gfm/unnamed-chunk-25-2.png)<!-- -->

``` r
ExpressionPlot(res$res, res$rpm, coef, res$sig,
               tissue_type_A, tissue_type_B,
               res$up_mirna, res$down_mirna,
               norm_adj_up, norm_adj_down,
               pCRC_adj_up, pCRC_adj_down)
```

![](comet_analysis_report_automated_no_lfcthreshold_08.12.20_res_alpha_0.05_with_LFC_shrinkage_LFC_0.58_mirge3_7.0_files/figure-gfm/unnamed-chunk-25-3.png)<!-- -->

``` r
signature_mirnas <- SigList(res, dds, tissue_type_A, tissue_type_B, coef,
                            norm_adj_up, norm_adj_down, 
                            pCRC_adj_up, pCRC_adj_down)
# Print list upregulated miRNA
signature_mirnas$up_mirna
```

<div style="border: 1px solid #ddd; padding: 5px; overflow-x: scroll; width:2000px; ">

<table class="table table-striped table-hover table-condensed" style="margin-left: auto; margin-right: auto;">

<caption>

Upregulated in tissue.type\_normal.liver\_vs\_tumor.colorect

</caption>

<thead>

<tr>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

miRNA

</th>

<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">

LFC

</th>

<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">

lfcSE

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

FDR

</th>

<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">

RPM normal.liver

</th>

<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">

RPM tumor.colorect

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

miRBase\_ID

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

Family

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

Seed

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

Chr

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

Cell-Type Specific

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

Norm Background

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

pCRC Background

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

Hsa-Mir-204-P1\_5p

</td>

<td style="text-align:right;">

1.03

</td>

<td style="text-align:right;">

0.46

</td>

<td style="text-align:left;">

1.25e-04

</td>

<td style="text-align:right;">

190

</td>

<td style="text-align:right;">

62

</td>

<td style="text-align:left;">

hsa-mir-204

</td>

<td style="text-align:left;">

MIR-204

</td>

<td style="text-align:left;">

UCCCUUU

</td>

<td style="text-align:left;">

chr9

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">Retinal
Epithelial Cell</span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-335\_5p

</td>

<td style="text-align:right;">

0.81

</td>

<td style="text-align:right;">

0.14

</td>

<td style="text-align:left;">

1.89e-08

</td>

<td style="text-align:right;">

205

</td>

<td style="text-align:right;">

150

</td>

<td style="text-align:left;">

hsa-mir-335

</td>

<td style="text-align:left;">

MIR-335

</td>

<td style="text-align:left;">

CAAGAGC

</td>

<td style="text-align:left;">

chr7

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">Retinal
Epithelial Cell</span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-144\_5p

</td>

<td style="text-align:right;">

2.22

</td>

<td style="text-align:right;">

0.28

</td>

<td style="text-align:left;">

3.68e-15

</td>

<td style="text-align:right;">

206

</td>

<td style="text-align:right;">

57

</td>

<td style="text-align:left;">

hsa-mir-144

</td>

<td style="text-align:left;">

MIR-144

</td>

<td style="text-align:left;">

GAUAUCA

</td>

<td style="text-align:left;">

chr17

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">Red
Blood Cell</span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-451\_5p

</td>

<td style="text-align:right;">

1.48

</td>

<td style="text-align:right;">

0.32

</td>

<td style="text-align:left;">

7.08e-06

</td>

<td style="text-align:right;">

3358

</td>

<td style="text-align:right;">

1565

</td>

<td style="text-align:left;">

hsa-mir-451a

</td>

<td style="text-align:left;">

MIR-451

</td>

<td style="text-align:left;">

AACCGUU

</td>

<td style="text-align:left;">

chr17

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">Red
Blood Cell</span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-150\_5p

</td>

<td style="text-align:right;">

2.07

</td>

<td style="text-align:right;">

0.29

</td>

<td style="text-align:left;">

5.24e-13

</td>

<td style="text-align:right;">

1053

</td>

<td style="text-align:right;">

280

</td>

<td style="text-align:left;">

hsa-mir-150

</td>

<td style="text-align:left;">

MIR-150

</td>

<td style="text-align:left;">

CUCCCAA

</td>

<td style="text-align:left;">

chr19

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">Lymphocyte</span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-122\_5p

</td>

<td style="text-align:right;">

11.53

</td>

<td style="text-align:right;">

0.78

</td>

<td style="text-align:left;">

4.04e-55

</td>

<td style="text-align:right;">

148842

</td>

<td style="text-align:right;">

9

</td>

<td style="text-align:left;">

hsa-mir-122

</td>

<td style="text-align:left;">

MIR-122

</td>

<td style="text-align:left;">

GGAGUGU

</td>

<td style="text-align:left;">

chr18

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">Hepatocyte</span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-192-P1\_5p/P2\_5p

</td>

<td style="text-align:right;">

1.31

</td>

<td style="text-align:right;">

0.27

</td>

<td style="text-align:left;">

2.62e-06

</td>

<td style="text-align:right;">

29546

</td>

<td style="text-align:right;">

13717

</td>

<td style="text-align:left;">

hsa-mir-192

</td>

<td style="text-align:left;">

MIR-192

</td>

<td style="text-align:left;">

UGACCUA

</td>

<td style="text-align:left;">

chr11

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">Epithelial
Cell</span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-15-P1a\_5p

</td>

<td style="text-align:right;">

0.62

</td>

<td style="text-align:right;">

0.12

</td>

<td style="text-align:left;">

1.17e-06

</td>

<td style="text-align:right;">

555

</td>

<td style="text-align:right;">

453

</td>

<td style="text-align:left;">

hsa-mir-15a

</td>

<td style="text-align:left;">

MIR-15

</td>

<td style="text-align:left;">

AGCAGCA

</td>

<td style="text-align:left;">

chr13

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">CD14+
Monocyte</span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-486\_5p

</td>

<td style="text-align:right;">

2.48

</td>

<td style="text-align:right;">

0.30

</td>

<td style="text-align:left;">

8.95e-16

</td>

<td style="text-align:right;">

10617

</td>

<td style="text-align:right;">

2013

</td>

<td style="text-align:left;">

hsa-mir-486-1

</td>

<td style="text-align:left;">

MIR-486

</td>

<td style="text-align:left;">

CCUGUAC

</td>

<td style="text-align:left;">

chr8

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">c(“Platelet”,
“Red Blood Cell”)</span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-126\_5p

</td>

<td style="text-align:right;">

1.95

</td>

<td style="text-align:right;">

0.17

</td>

<td style="text-align:left;">

7.32e-30

</td>

<td style="text-align:right;">

11148

</td>

<td style="text-align:right;">

3455

</td>

<td style="text-align:left;">

hsa-mir-126

</td>

<td style="text-align:left;">

MIR-126

</td>

<td style="text-align:left;">

AUUAUUA

</td>

<td style="text-align:left;">

chr9

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">c(“Endothelial
Cell”, “Platelet”)</span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-342\_3p

</td>

<td style="text-align:right;">

1.05

</td>

<td style="text-align:right;">

0.21

</td>

<td style="text-align:left;">

5.60e-07

</td>

<td style="text-align:right;">

254

</td>

<td style="text-align:right;">

145

</td>

<td style="text-align:left;">

hsa-mir-342

</td>

<td style="text-align:left;">

MIR-342

</td>

<td style="text-align:left;">

CUCACAC

</td>

<td style="text-align:left;">

chr14

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">c(“Dendritic
Cell”, “Lymphocyte”, “Macrophage”)</span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-885\_5p

</td>

<td style="text-align:right;">

9.48

</td>

<td style="text-align:right;">

0.68

</td>

<td style="text-align:left;">

8.47e-41

</td>

<td style="text-align:right;">

510

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:left;">

hsa-mir-885

</td>

<td style="text-align:left;">

MIR-885

</td>

<td style="text-align:left;">

CCAUUAC

</td>

<td style="text-align:left;">

chr3

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-139\_5p

</td>

<td style="text-align:right;">

4.20

</td>

<td style="text-align:right;">

0.21

</td>

<td style="text-align:left;">

6.49e-84

</td>

<td style="text-align:right;">

174

</td>

<td style="text-align:right;">

11

</td>

<td style="text-align:left;">

hsa-mir-139

</td>

<td style="text-align:left;">

MIR-139

</td>

<td style="text-align:left;">

CUACAGU

</td>

<td style="text-align:left;">

chr11

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Let-7-P1c\_5p

</td>

<td style="text-align:right;">

3.25

</td>

<td style="text-align:right;">

0.23

</td>

<td style="text-align:left;">

3.14e-43

</td>

<td style="text-align:right;">

4945

</td>

<td style="text-align:right;">

606

</td>

<td style="text-align:left;">

hsa-let-7c

</td>

<td style="text-align:left;">

LET-7

</td>

<td style="text-align:left;">

GAGGUAG

</td>

<td style="text-align:left;">

chr21

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-10-P2c\_5p

</td>

<td style="text-align:right;">

3.24

</td>

<td style="text-align:right;">

0.28

</td>

<td style="text-align:left;">

2.36e-31

</td>

<td style="text-align:right;">

1202

</td>

<td style="text-align:right;">

142

</td>

<td style="text-align:left;">

hsa-mir-99a

</td>

<td style="text-align:left;">

MIR-10

</td>

<td style="text-align:left;">

ACCCGUA

</td>

<td style="text-align:left;">

chr21

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-10-P3c\_5p

</td>

<td style="text-align:right;">

2.87

</td>

<td style="text-align:right;">

0.25

</td>

<td style="text-align:left;">

2.44e-30

</td>

<td style="text-align:right;">

1620

</td>

<td style="text-align:right;">

239

</td>

<td style="text-align:left;">

hsa-mir-125b-2

</td>

<td style="text-align:left;">

MIR-10

</td>

<td style="text-align:left;">

CCCUGAG

</td>

<td style="text-align:left;">

chr21

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-30-P1a\_5p

</td>

<td style="text-align:right;">

2.51

</td>

<td style="text-align:right;">

0.19

</td>

<td style="text-align:left;">

3.24e-40

</td>

<td style="text-align:right;">

15175

</td>

<td style="text-align:right;">

3165

</td>

<td style="text-align:left;">

hsa-mir-30a

</td>

<td style="text-align:left;">

MIR-30

</td>

<td style="text-align:left;">

GUAAACA

</td>

<td style="text-align:left;">

chr6

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-193-P1a\_5p

</td>

<td style="text-align:right;">

2.48

</td>

<td style="text-align:right;">

0.23

</td>

<td style="text-align:left;">

3.82e-25

</td>

<td style="text-align:right;">

298

</td>

<td style="text-align:right;">

65

</td>

<td style="text-align:left;">

hsa-mir-193a

</td>

<td style="text-align:left;">

MIR-193

</td>

<td style="text-align:left;">

GGGUCUU

</td>

<td style="text-align:left;">

chr17

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-193-P2a\_3p/P2b\_3p

</td>

<td style="text-align:right;">

2.33

</td>

<td style="text-align:right;">

0.19

</td>

<td style="text-align:left;">

1.61e-32

</td>

<td style="text-align:right;">

219

</td>

<td style="text-align:right;">

56

</td>

<td style="text-align:left;">

hsa-mir-365b

</td>

<td style="text-align:left;">

MIR-193

</td>

<td style="text-align:left;">

AAUGCCC

</td>

<td style="text-align:left;">

chr17

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-10-P3a\_5p

</td>

<td style="text-align:right;">

2.29

</td>

<td style="text-align:right;">

0.26

</td>

<td style="text-align:left;">

3.08e-18

</td>

<td style="text-align:right;">

611

</td>

<td style="text-align:right;">

145

</td>

<td style="text-align:left;">

hsa-mir-125b-1

</td>

<td style="text-align:left;">

MIR-10

</td>

<td style="text-align:left;">

CCCUGAG

</td>

<td style="text-align:left;">

chr11

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-193-P1b\_3p

</td>

<td style="text-align:right;">

2.17

</td>

<td style="text-align:right;">

0.19

</td>

<td style="text-align:left;">

2.30e-29

</td>

<td style="text-align:right;">

514

</td>

<td style="text-align:right;">

145

</td>

<td style="text-align:left;">

hsa-mir-193b

</td>

<td style="text-align:left;">

MIR-193

</td>

<td style="text-align:left;">

ACUGGCC

</td>

<td style="text-align:left;">

chr16

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-455\_5p

</td>

<td style="text-align:right;">

2.10

</td>

<td style="text-align:right;">

0.15

</td>

<td style="text-align:left;">

6.57e-45

</td>

<td style="text-align:right;">

279

</td>

<td style="text-align:right;">

83

</td>

<td style="text-align:left;">

hsa-mir-455

</td>

<td style="text-align:left;">

MIR-455

</td>

<td style="text-align:left;">

AUGUGCC

</td>

<td style="text-align:left;">

chr9

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-574\_3p

</td>

<td style="text-align:right;">

1.99

</td>

<td style="text-align:right;">

0.16

</td>

<td style="text-align:left;">

1.44e-34

</td>

<td style="text-align:right;">

745

</td>

<td style="text-align:right;">

240

</td>

<td style="text-align:left;">

hsa-mir-574

</td>

<td style="text-align:left;">

MIR-574

</td>

<td style="text-align:left;">

ACGCUCA

</td>

<td style="text-align:left;">

chr4

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-101-P1\_3p/P2\_3p

</td>

<td style="text-align:right;">

1.98

</td>

<td style="text-align:right;">

0.15

</td>

<td style="text-align:left;">

6.40e-41

</td>

<td style="text-align:right;">

14649

</td>

<td style="text-align:right;">

4638

</td>

<td style="text-align:left;">

hsa-mir-101-1

</td>

<td style="text-align:left;">

MIR-101

</td>

<td style="text-align:left;">

UACAGUA

</td>

<td style="text-align:left;">

chr1

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-148-P1\_3p

</td>

<td style="text-align:right;">

1.86

</td>

<td style="text-align:right;">

0.20

</td>

<td style="text-align:left;">

1.75e-19

</td>

<td style="text-align:right;">

105086

</td>

<td style="text-align:right;">

34704

</td>

<td style="text-align:left;">

hsa-mir-148a

</td>

<td style="text-align:left;">

MIR-148

</td>

<td style="text-align:left;">

CAGUGCA

</td>

<td style="text-align:left;">

chr7

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-423\_5p

</td>

<td style="text-align:right;">

1.83

</td>

<td style="text-align:right;">

0.20

</td>

<td style="text-align:left;">

1.11e-18

</td>

<td style="text-align:right;">

1198

</td>

<td style="text-align:right;">

433

</td>

<td style="text-align:left;">

hsa-mir-423

</td>

<td style="text-align:left;">

MIR-423

</td>

<td style="text-align:left;">

GAGGGGC

</td>

<td style="text-align:left;">

chr17

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-22-P1a\_3p

</td>

<td style="text-align:right;">

1.79

</td>

<td style="text-align:right;">

0.14

</td>

<td style="text-align:left;">

1.08e-36

</td>

<td style="text-align:right;">

78932

</td>

<td style="text-align:right;">

29018

</td>

<td style="text-align:left;">

hsa-mir-22

</td>

<td style="text-align:left;">

MIR-22

</td>

<td style="text-align:left;">

AGCUGCC

</td>

<td style="text-align:left;">

chr17

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-744\_5p

</td>

<td style="text-align:right;">

1.79

</td>

<td style="text-align:right;">

0.17

</td>

<td style="text-align:left;">

1.28e-23

</td>

<td style="text-align:right;">

165

</td>

<td style="text-align:right;">

61

</td>

<td style="text-align:left;">

hsa-mir-744

</td>

<td style="text-align:left;">

MIR-744

</td>

<td style="text-align:left;">

GCGGGGC

</td>

<td style="text-align:left;">

chr17

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-30-P2a\_5p/P2b\_5p/P2c\_5p

</td>

<td style="text-align:right;">

1.45

</td>

<td style="text-align:right;">

0.12

</td>

<td style="text-align:left;">

7.16e-32

</td>

<td style="text-align:right;">

5891

</td>

<td style="text-align:right;">

2793

</td>

<td style="text-align:left;">

hsa-mir-30c-2

</td>

<td style="text-align:left;">

MIR-30

</td>

<td style="text-align:left;">

GUAAACA

</td>

<td style="text-align:left;">

chr6

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-130-P1a\_3p

</td>

<td style="text-align:right;">

1.43

</td>

<td style="text-align:right;">

0.16

</td>

<td style="text-align:left;">

1.39e-17

</td>

<td style="text-align:right;">

907

</td>

<td style="text-align:right;">

421

</td>

<td style="text-align:left;">

hsa-mir-130a

</td>

<td style="text-align:left;">

MIR-130

</td>

<td style="text-align:left;">

AGUGCAA

</td>

<td style="text-align:left;">

chr11

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-378\_3p

</td>

<td style="text-align:right;">

1.35

</td>

<td style="text-align:right;">

0.16

</td>

<td style="text-align:left;">

1.07e-16

</td>

<td style="text-align:right;">

12736

</td>

<td style="text-align:right;">

6426

</td>

<td style="text-align:left;">

hsa-mir-378a

</td>

<td style="text-align:left;">

MIR-378

</td>

<td style="text-align:left;">

CUGGACU

</td>

<td style="text-align:left;">

chr5

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-148-P2\_3p

</td>

<td style="text-align:right;">

1.34

</td>

<td style="text-align:right;">

0.25

</td>

<td style="text-align:left;">

1.57e-07

</td>

<td style="text-align:right;">

289

</td>

<td style="text-align:right;">

137

</td>

<td style="text-align:left;">

hsa-mir-152

</td>

<td style="text-align:left;">

MIR-148

</td>

<td style="text-align:left;">

CAGUGCA

</td>

<td style="text-align:left;">

chr17

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-197\_3p

</td>

<td style="text-align:right;">

1.28

</td>

<td style="text-align:right;">

0.15

</td>

<td style="text-align:left;">

5.68e-16

</td>

<td style="text-align:right;">

344

</td>

<td style="text-align:right;">

181

</td>

<td style="text-align:left;">

hsa-mir-197

</td>

<td style="text-align:left;">

MIR-197

</td>

<td style="text-align:left;">

UCACCAC

</td>

<td style="text-align:left;">

chr1

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-15-P1c\_5p

</td>

<td style="text-align:right;">

1.25

</td>

<td style="text-align:right;">

0.19

</td>

<td style="text-align:left;">

3.90e-11

</td>

<td style="text-align:right;">

276

</td>

<td style="text-align:right;">

143

</td>

<td style="text-align:left;">

hsa-mir-497

</td>

<td style="text-align:left;">

MIR-15

</td>

<td style="text-align:left;">

AGCAGCA

</td>

<td style="text-align:left;">

chr17

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-26-P3\_5p

</td>

<td style="text-align:right;">

1.24

</td>

<td style="text-align:right;">

0.16

</td>

<td style="text-align:left;">

2.22e-14

</td>

<td style="text-align:right;">

6787

</td>

<td style="text-align:right;">

3467

</td>

<td style="text-align:left;">

hsa-mir-26a-1

</td>

<td style="text-align:left;">

MIR-26

</td>

<td style="text-align:left;">

UCAAGUA

</td>

<td style="text-align:left;">

chr3

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-15-P2c\_5p

</td>

<td style="text-align:right;">

1.19

</td>

<td style="text-align:right;">

0.20

</td>

<td style="text-align:left;">

1.78e-09

</td>

<td style="text-align:right;">

496

</td>

<td style="text-align:right;">

257

</td>

<td style="text-align:left;">

hsa-mir-195

</td>

<td style="text-align:left;">

MIR-15

</td>

<td style="text-align:left;">

AGCAGCA

</td>

<td style="text-align:left;">

chr17

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-331\_3p

</td>

<td style="text-align:right;">

1.16

</td>

<td style="text-align:right;">

0.21

</td>

<td style="text-align:left;">

6.75e-08

</td>

<td style="text-align:right;">

110

</td>

<td style="text-align:right;">

60

</td>

<td style="text-align:left;">

hsa-mir-331

</td>

<td style="text-align:left;">

MIR-331

</td>

<td style="text-align:left;">

CCCCUGG

</td>

<td style="text-align:left;">

chr12

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-26-P1\_5p/P2\_5p

</td>

<td style="text-align:right;">

1.10

</td>

<td style="text-align:right;">

0.11

</td>

<td style="text-align:left;">

1.13e-22

</td>

<td style="text-align:right;">

62274

</td>

<td style="text-align:right;">

35584

</td>

<td style="text-align:left;">

hsa-mir-26b

</td>

<td style="text-align:left;">

MIR-26

</td>

<td style="text-align:left;">

UCAAGUA

</td>

<td style="text-align:left;">

chr2

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-10-P3b\_5p

</td>

<td style="text-align:right;">

1.09

</td>

<td style="text-align:right;">

0.23

</td>

<td style="text-align:left;">

2.99e-06

</td>

<td style="text-align:right;">

3040

</td>

<td style="text-align:right;">

1649

</td>

<td style="text-align:left;">

hsa-mir-125a

</td>

<td style="text-align:left;">

MIR-10

</td>

<td style="text-align:left;">

CCCUGAG

</td>

<td style="text-align:left;">

chr19

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-30-P1b\_5p

</td>

<td style="text-align:right;">

1.08

</td>

<td style="text-align:right;">

0.12

</td>

<td style="text-align:left;">

5.44e-19

</td>

<td style="text-align:right;">

11240

</td>

<td style="text-align:right;">

6833

</td>

<td style="text-align:left;">

hsa-mir-30e

</td>

<td style="text-align:left;">

MIR-30

</td>

<td style="text-align:left;">

GUAAACA

</td>

<td style="text-align:left;">

chr1

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-340\_5p

</td>

<td style="text-align:right;">

1.05

</td>

<td style="text-align:right;">

0.13

</td>

<td style="text-align:left;">

1.20e-15

</td>

<td style="text-align:right;">

1030

</td>

<td style="text-align:right;">

657

</td>

<td style="text-align:left;">

hsa-mir-340

</td>

<td style="text-align:left;">

MIR-340

</td>

<td style="text-align:left;">

UAUAAAG

</td>

<td style="text-align:left;">

chr5

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-27-P1\_3p/P2\_3p

</td>

<td style="text-align:right;">

1.04

</td>

<td style="text-align:right;">

0.12

</td>

<td style="text-align:left;">

4.13e-17

</td>

<td style="text-align:right;">

50009

</td>

<td style="text-align:right;">

29444

</td>

<td style="text-align:left;">

hsa-mir-27a

</td>

<td style="text-align:left;">

MIR-27

</td>

<td style="text-align:left;">

UCACAGU

</td>

<td style="text-align:left;">

chr19

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-10-P2a\_5p

</td>

<td style="text-align:right;">

1.00

</td>

<td style="text-align:right;">

0.27

</td>

<td style="text-align:left;">

4.80e-04

</td>

<td style="text-align:right;">

2634

</td>

<td style="text-align:right;">

1724

</td>

<td style="text-align:left;">

hsa-mir-100

</td>

<td style="text-align:left;">

MIR-10

</td>

<td style="text-align:left;">

ACCCGUA

</td>

<td style="text-align:left;">

chr11

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-29-P1b\_3p

</td>

<td style="text-align:right;">

0.98

</td>

<td style="text-align:right;">

0.17

</td>

<td style="text-align:left;">

1.55e-08

</td>

<td style="text-align:right;">

516

</td>

<td style="text-align:right;">

331

</td>

<td style="text-align:left;">

hsa-mir-29c

</td>

<td style="text-align:left;">

MIR-29

</td>

<td style="text-align:left;">

AGCACCA

</td>

<td style="text-align:left;">

chr1

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-154-P23\_3p

</td>

<td style="text-align:right;">

0.94

</td>

<td style="text-align:right;">

0.16

</td>

<td style="text-align:left;">

3.85e-08

</td>

<td style="text-align:right;">

222

</td>

<td style="text-align:right;">

148

</td>

<td style="text-align:left;">

hsa-mir-654

</td>

<td style="text-align:left;">

MIR-154

</td>

<td style="text-align:left;">

AUGUCUG

</td>

<td style="text-align:left;">

chr14

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-15-P2a\_5p/P2b\_5p

</td>

<td style="text-align:right;">

0.87

</td>

<td style="text-align:right;">

0.14

</td>

<td style="text-align:left;">

9.10e-10

</td>

<td style="text-align:right;">

6884

</td>

<td style="text-align:right;">

4522

</td>

<td style="text-align:left;">

hsa-mir-16-1

</td>

<td style="text-align:left;">

MIR-15

</td>

<td style="text-align:left;">

AGCAGCA

</td>

<td style="text-align:left;">

chr13

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-30-P1c\_5p

</td>

<td style="text-align:right;">

0.84

</td>

<td style="text-align:right;">

0.13

</td>

<td style="text-align:left;">

9.49e-10

</td>

<td style="text-align:right;">

13680

</td>

<td style="text-align:right;">

9703

</td>

<td style="text-align:left;">

hsa-mir-30d

</td>

<td style="text-align:left;">

MIR-30

</td>

<td style="text-align:left;">

GUAAACA

</td>

<td style="text-align:left;">

chr8

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-154-P13\_5p

</td>

<td style="text-align:right;">

0.78

</td>

<td style="text-align:right;">

0.18

</td>

<td style="text-align:left;">

2.94e-05

</td>

<td style="text-align:right;">

442

</td>

<td style="text-align:right;">

328

</td>

<td style="text-align:left;">

hsa-mir-411

</td>

<td style="text-align:left;">

MIR-154

</td>

<td style="text-align:left;">

AGUAGAC

</td>

<td style="text-align:left;">

chr14

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-28-P2\_5p

</td>

<td style="text-align:right;">

0.76

</td>

<td style="text-align:right;">

0.15

</td>

<td style="text-align:left;">

1.01e-06

</td>

<td style="text-align:right;">

3098

</td>

<td style="text-align:right;">

2244

</td>

<td style="text-align:left;">

hsa-mir-151a

</td>

<td style="text-align:left;">

MIR-28

</td>

<td style="text-align:left;">

CGAGGAG

</td>

<td style="text-align:left;">

chr8

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-92-P1a\_3p/P1b\_3p

</td>

<td style="text-align:right;">

0.71

</td>

<td style="text-align:right;">

0.18

</td>

<td style="text-align:left;">

2.33e-04

</td>

<td style="text-align:right;">

51173

</td>

<td style="text-align:right;">

39320

</td>

<td style="text-align:left;">

hsa-mir-92a-1

</td>

<td style="text-align:left;">

MIR-92

</td>

<td style="text-align:left;">

AUUGCAC

</td>

<td style="text-align:left;">

chr13

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-136\_3p

</td>

<td style="text-align:right;">

0.65

</td>

<td style="text-align:right;">

0.18

</td>

<td style="text-align:left;">

4.46e-04

</td>

<td style="text-align:right;">

208

</td>

<td style="text-align:right;">

173

</td>

<td style="text-align:left;">

hsa-mir-136

</td>

<td style="text-align:left;">

MIR-136

</td>

<td style="text-align:left;">

AUCAUCG

</td>

<td style="text-align:left;">

chr14

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

</tbody>

</table>

</div>

``` r
# Number of upregulated miRNA
signature_mirnas$number_upregulated
```

    ## [1] 51

``` r
# Print list downregulated miRNA
signature_mirnas$down_mirna
```

<div style="border: 1px solid #ddd; padding: 5px; overflow-x: scroll; width:2000px; ">

<table class="table table-striped table-hover table-condensed" style="margin-left: auto; margin-right: auto;">

<caption>

Downregulated in tissue.type\_normal.liver\_vs\_tumor.colorect

</caption>

<thead>

<tr>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

miRNA

</th>

<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">

LFC

</th>

<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">

lfcSE

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

FDR

</th>

<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">

RPM normal.liver

</th>

<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">

RPM tumor.colorect

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

miRBase\_ID

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

Family

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

Seed

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

Chr

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

Cell-Type Specific

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

Norm Background

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

pCRC Background

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

Hsa-Mir-143\_3p

</td>

<td style="text-align:right;">

\-1.72

</td>

<td style="text-align:right;">

0.23

</td>

<td style="text-align:left;">

4.38e-14

</td>

<td style="text-align:right;">

37572

</td>

<td style="text-align:right;">

148431

</td>

<td style="text-align:left;">

hsa-mir-143

</td>

<td style="text-align:left;">

MIR-143

</td>

<td style="text-align:left;">

GAGAUGA

</td>

<td style="text-align:left;">

chr5

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">Mesenchymal</span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-24-P1\_3p/P2\_3p

</td>

<td style="text-align:right;">

\-0.77

</td>

<td style="text-align:right;">

0.16

</td>

<td style="text-align:left;">

8.36e-06

</td>

<td style="text-align:right;">

319

</td>

<td style="text-align:right;">

626

</td>

<td style="text-align:left;">

hsa-mir-24-2

</td>

<td style="text-align:left;">

MIR-24

</td>

<td style="text-align:left;">

GGCUCAG

</td>

<td style="text-align:left;">

chr19

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">Macrophage</span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-8-P2b\_3p

</td>

<td style="text-align:right;">

\-6.07

</td>

<td style="text-align:right;">

0.19

</td>

<td style="text-align:left;">

3.90e-210

</td>

<td style="text-align:right;">

31

</td>

<td style="text-align:right;">

2727

</td>

<td style="text-align:left;">

hsa-mir-200c

</td>

<td style="text-align:left;">

MIR-8

</td>

<td style="text-align:left;">

AAUACUG

</td>

<td style="text-align:left;">

chr12

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">Epithelial
Cell</span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-8-P2a\_3p

</td>

<td style="text-align:right;">

\-4.36

</td>

<td style="text-align:right;">

0.19

</td>

<td style="text-align:left;">

2.71e-109

</td>

<td style="text-align:right;">

369

</td>

<td style="text-align:right;">

10330

</td>

<td style="text-align:left;">

hsa-mir-200b

</td>

<td style="text-align:left;">

MIR-8

</td>

<td style="text-align:left;">

AAUACUG

</td>

<td style="text-align:left;">

chr1

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">Epithelial
Cell</span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-17-P1a\_5p/P1b\_5p

</td>

<td style="text-align:right;">

\-1.32

</td>

<td style="text-align:right;">

0.16

</td>

<td style="text-align:left;">

4.02e-16

</td>

<td style="text-align:right;">

278

</td>

<td style="text-align:right;">

900

</td>

<td style="text-align:left;">

hsa-mir-17

</td>

<td style="text-align:left;">

MIR-17

</td>

<td style="text-align:left;">

AAAGUGC

</td>

<td style="text-align:left;">

chr13

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">CD14+
Monocyte</span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-15-P1b\_5p

</td>

<td style="text-align:right;">

\-0.92

</td>

<td style="text-align:right;">

0.13

</td>

<td style="text-align:left;">

1.30e-11

</td>

<td style="text-align:right;">

114

</td>

<td style="text-align:right;">

273

</td>

<td style="text-align:left;">

hsa-mir-15b

</td>

<td style="text-align:left;">

MIR-15

</td>

<td style="text-align:left;">

AGCAGCA

</td>

<td style="text-align:left;">

chr3

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">CD14+
Monocyte</span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-133-P1\_3p/P2\_3p/P3\_3p

</td>

<td style="text-align:right;">

\-2.03

</td>

<td style="text-align:right;">

0.32

</td>

<td style="text-align:left;">

1.91e-10

</td>

<td style="text-align:right;">

22

</td>

<td style="text-align:right;">

116

</td>

<td style="text-align:left;">

hsa-mir-133a-2

</td>

<td style="text-align:left;">

MIR-133

</td>

<td style="text-align:left;">

UUGGUCC

</td>

<td style="text-align:left;">

chr20

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">c(“Skeletal
Myocyte”, “Stem Cell”)</span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-155\_5p

</td>

<td style="text-align:right;">

\-1.44

</td>

<td style="text-align:right;">

0.21

</td>

<td style="text-align:left;">

1.42e-10

</td>

<td style="text-align:right;">

113

</td>

<td style="text-align:right;">

387

</td>

<td style="text-align:left;">

hsa-mir-155

</td>

<td style="text-align:left;">

MIR-155

</td>

<td style="text-align:left;">

UAAUGCU

</td>

<td style="text-align:left;">

chr21

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">c(“Lymphocyte”,
“Macrophage”)</span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-7-P1\_5p/P2\_5p/P3\_5p

</td>

<td style="text-align:right;">

\-5.51

</td>

<td style="text-align:right;">

0.32

</td>

<td style="text-align:left;">

2.03e-63

</td>

<td style="text-align:right;">

3

</td>

<td style="text-align:right;">

199

</td>

<td style="text-align:left;">

hsa-mir-7-1

</td>

<td style="text-align:left;">

MIR-7

</td>

<td style="text-align:left;">

GGAAGAC

</td>

<td style="text-align:left;">

chr9

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">c(“Islet
Cell”, “Neural”)</span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-375\_3p

</td>

<td style="text-align:right;">

\-0.84

</td>

<td style="text-align:right;">

0.32

</td>

<td style="text-align:left;">

2.21e-02

</td>

<td style="text-align:right;">

2553

</td>

<td style="text-align:right;">

5263

</td>

<td style="text-align:left;">

hsa-mir-375

</td>

<td style="text-align:left;">

MIR-375

</td>

<td style="text-align:left;">

UUGUUCG

</td>

<td style="text-align:left;">

chr2

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">c(“Epithelial
Cell”, “Islet Cell”, “Neural”)</span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-196-P1\_5p/P2\_5p

</td>

<td style="text-align:right;">

\-6.71

</td>

<td style="text-align:right;">

0.30

</td>

<td style="text-align:left;">

2.31e-108

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

319

</td>

<td style="text-align:left;">

hsa-mir-196a-1

</td>

<td style="text-align:left;">

MIR-196

</td>

<td style="text-align:left;">

AGGUAGU

</td>

<td style="text-align:left;">

chr17

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-135-P3\_5p

</td>

<td style="text-align:right;">

\-6.12

</td>

<td style="text-align:right;">

0.31

</td>

<td style="text-align:left;">

3.71e-82

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

155

</td>

<td style="text-align:left;">

hsa-mir-135b

</td>

<td style="text-align:left;">

MIR-135

</td>

<td style="text-align:left;">

AUGGCUU

</td>

<td style="text-align:left;">

chr1

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-196-P3\_5p

</td>

<td style="text-align:right;">

\-6.11

</td>

<td style="text-align:right;">

0.32

</td>

<td style="text-align:left;">

1.17e-77

</td>

<td style="text-align:right;">

9

</td>

<td style="text-align:right;">

828

</td>

<td style="text-align:left;">

hsa-mir-196b

</td>

<td style="text-align:left;">

MIR-196

</td>

<td style="text-align:left;">

AGGUAGU

</td>

<td style="text-align:left;">

chr7

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-31\_5p

</td>

<td style="text-align:right;">

\-5.76

</td>

<td style="text-align:right;">

0.45

</td>

<td style="text-align:left;">

3.20e-37

</td>

<td style="text-align:right;">

7

</td>

<td style="text-align:right;">

627

</td>

<td style="text-align:left;">

hsa-mir-31

</td>

<td style="text-align:left;">

MIR-31

</td>

<td style="text-align:left;">

GGCAAGA

</td>

<td style="text-align:left;">

chr9

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-577\_5p

</td>

<td style="text-align:right;">

\-5.70

</td>

<td style="text-align:right;">

0.29

</td>

<td style="text-align:left;">

1.37e-80

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

142

</td>

<td style="text-align:left;">

hsa-mir-577

</td>

<td style="text-align:left;">

MIR-577

</td>

<td style="text-align:left;">

UAGAUAA

</td>

<td style="text-align:left;">

chr4

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-8-P1b\_3p

</td>

<td style="text-align:right;">

\-5.35

</td>

<td style="text-align:right;">

0.25

</td>

<td style="text-align:left;">

1.29e-93

</td>

<td style="text-align:right;">

154

</td>

<td style="text-align:right;">

7518

</td>

<td style="text-align:left;">

hsa-mir-141

</td>

<td style="text-align:left;">

MIR-8

</td>

<td style="text-align:left;">

AACACUG

</td>

<td style="text-align:left;">

chr12

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-96-P3\_5p

</td>

<td style="text-align:right;">

\-5.08

</td>

<td style="text-align:right;">

0.27

</td>

<td style="text-align:left;">

2.84e-74

</td>

<td style="text-align:right;">

27

</td>

<td style="text-align:right;">

1063

</td>

<td style="text-align:left;">

hsa-mir-183

</td>

<td style="text-align:left;">

MIR-96

</td>

<td style="text-align:left;">

AUGGCAC

</td>

<td style="text-align:left;">

chr7

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-8-P3a\_3p

</td>

<td style="text-align:right;">

\-4.60

</td>

<td style="text-align:right;">

0.20

</td>

<td style="text-align:left;">

5.55e-118

</td>

<td style="text-align:right;">

68

</td>

<td style="text-align:right;">

2303

</td>

<td style="text-align:left;">

hsa-mir-429

</td>

<td style="text-align:left;">

MIR-8

</td>

<td style="text-align:left;">

AAUACUG

</td>

<td style="text-align:left;">

chr1

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-96-P2\_5p

</td>

<td style="text-align:right;">

\-4.47

</td>

<td style="text-align:right;">

0.20

</td>

<td style="text-align:left;">

4.46e-107

</td>

<td style="text-align:right;">

370

</td>

<td style="text-align:right;">

10417

</td>

<td style="text-align:left;">

hsa-mir-182

</td>

<td style="text-align:left;">

MIR-96

</td>

<td style="text-align:left;">

UUGGCAA

</td>

<td style="text-align:left;">

chr7

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-10-P1b\_5p

</td>

<td style="text-align:right;">

\-3.68

</td>

<td style="text-align:right;">

0.31

</td>

<td style="text-align:left;">

1.37e-31

</td>

<td style="text-align:right;">

2704

</td>

<td style="text-align:right;">

39286

</td>

<td style="text-align:left;">

hsa-mir-10b

</td>

<td style="text-align:left;">

MIR-10

</td>

<td style="text-align:left;">

ACCCUGU

</td>

<td style="text-align:left;">

chr2

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-8-P1a\_3p

</td>

<td style="text-align:right;">

\-3.64

</td>

<td style="text-align:right;">

0.19

</td>

<td style="text-align:left;">

9.42e-76

</td>

<td style="text-align:right;">

108

</td>

<td style="text-align:right;">

1772

</td>

<td style="text-align:left;">

hsa-mir-200a

</td>

<td style="text-align:left;">

MIR-8

</td>

<td style="text-align:left;">

AACACUG

</td>

<td style="text-align:left;">

chr1

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-203\_3p

</td>

<td style="text-align:right;">

\-3.42

</td>

<td style="text-align:right;">

0.24

</td>

<td style="text-align:left;">

2.77e-45

</td>

<td style="text-align:right;">

171

</td>

<td style="text-align:right;">

2307

</td>

<td style="text-align:left;">

hsa-mir-203a

</td>

<td style="text-align:left;">

MIR-203

</td>

<td style="text-align:left;">

UGAAAUG

</td>

<td style="text-align:left;">

chr14

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-221-P1\_3p

</td>

<td style="text-align:right;">

\-3.13

</td>

<td style="text-align:right;">

0.15

</td>

<td style="text-align:left;">

5.83e-92

</td>

<td style="text-align:right;">

217

</td>

<td style="text-align:right;">

2529

</td>

<td style="text-align:left;">

hsa-mir-221

</td>

<td style="text-align:left;">

MIR-221

</td>

<td style="text-align:left;">

GCUACAU

</td>

<td style="text-align:left;">

chrX

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-224\_5p

</td>

<td style="text-align:right;">

\-2.99

</td>

<td style="text-align:right;">

0.22

</td>

<td style="text-align:left;">

2.60e-39

</td>

<td style="text-align:right;">

36

</td>

<td style="text-align:right;">

395

</td>

<td style="text-align:left;">

hsa-mir-224

</td>

<td style="text-align:left;">

MIR-224

</td>

<td style="text-align:left;">

AAGUCAC

</td>

<td style="text-align:left;">

chrX

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-221-P2\_3p

</td>

<td style="text-align:right;">

\-2.81

</td>

<td style="text-align:right;">

0.17

</td>

<td style="text-align:left;">

1.32e-61

</td>

<td style="text-align:right;">

173

</td>

<td style="text-align:right;">

1633

</td>

<td style="text-align:left;">

hsa-mir-222

</td>

<td style="text-align:left;">

MIR-221

</td>

<td style="text-align:left;">

GCUACAU

</td>

<td style="text-align:left;">

chrX

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-181-P2c\_5p

</td>

<td style="text-align:right;">

\-2.59

</td>

<td style="text-align:right;">

0.19

</td>

<td style="text-align:left;">

7.58e-42

</td>

<td style="text-align:right;">

26

</td>

<td style="text-align:right;">

196

</td>

<td style="text-align:left;">

hsa-mir-181d

</td>

<td style="text-align:left;">

MIR-181

</td>

<td style="text-align:left;">

ACAUUCA

</td>

<td style="text-align:left;">

chr19

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-181-P1c\_5p

</td>

<td style="text-align:right;">

\-2.35

</td>

<td style="text-align:right;">

0.17

</td>

<td style="text-align:left;">

3.48e-42

</td>

<td style="text-align:right;">

261

</td>

<td style="text-align:right;">

1724

</td>

<td style="text-align:left;">

hsa-mir-181c

</td>

<td style="text-align:left;">

MIR-181

</td>

<td style="text-align:left;">

ACAUUCA

</td>

<td style="text-align:left;">

chr19

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-584\_5p

</td>

<td style="text-align:right;">

\-2.25

</td>

<td style="text-align:right;">

0.25

</td>

<td style="text-align:left;">

4.29e-18

</td>

<td style="text-align:right;">

18

</td>

<td style="text-align:right;">

110

</td>

<td style="text-align:left;">

hsa-mir-584

</td>

<td style="text-align:left;">

MIR-584

</td>

<td style="text-align:left;">

UAUGGUU

</td>

<td style="text-align:left;">

chr5

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-190-P1\_5p

</td>

<td style="text-align:right;">

\-2.24

</td>

<td style="text-align:right;">

0.22

</td>

<td style="text-align:left;">

4.84e-24

</td>

<td style="text-align:right;">

23

</td>

<td style="text-align:right;">

140

</td>

<td style="text-align:left;">

hsa-mir-190a

</td>

<td style="text-align:left;">

MIR-190

</td>

<td style="text-align:left;">

GAUAUGU

</td>

<td style="text-align:left;">

chr15

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-17-P2a\_5p

</td>

<td style="text-align:right;">

\-2.22

</td>

<td style="text-align:right;">

0.23

</td>

<td style="text-align:left;">

1.80e-21

</td>

<td style="text-align:right;">

26

</td>

<td style="text-align:right;">

168

</td>

<td style="text-align:left;">

hsa-mir-18a

</td>

<td style="text-align:left;">

MIR-17

</td>

<td style="text-align:left;">

AAGGUGC

</td>

<td style="text-align:left;">

chr13

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-21\_5p

</td>

<td style="text-align:right;">

\-2.10

</td>

<td style="text-align:right;">

0.14

</td>

<td style="text-align:left;">

5.68e-48

</td>

<td style="text-align:right;">

19240

</td>

<td style="text-align:right;">

105735

</td>

<td style="text-align:left;">

hsa-mir-21

</td>

<td style="text-align:left;">

MIR-21

</td>

<td style="text-align:left;">

AGCUUAU

</td>

<td style="text-align:left;">

chr17

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-92-P1c\_3p

</td>

<td style="text-align:right;">

\-1.90

</td>

<td style="text-align:right;">

0.21

</td>

<td style="text-align:left;">

3.08e-18

</td>

<td style="text-align:right;">

303

</td>

<td style="text-align:right;">

1398

</td>

<td style="text-align:left;">

hsa-mir-92b

</td>

<td style="text-align:left;">

MIR-92

</td>

<td style="text-align:left;">

AUUGCAC

</td>

<td style="text-align:left;">

chr1

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-10-P1a\_5p

</td>

<td style="text-align:right;">

\-1.79

</td>

<td style="text-align:right;">

0.25

</td>

<td style="text-align:left;">

2.08e-11

</td>

<td style="text-align:right;">

25828

</td>

<td style="text-align:right;">

97123

</td>

<td style="text-align:left;">

hsa-mir-10a

</td>

<td style="text-align:left;">

MIR-10

</td>

<td style="text-align:left;">

ACCCUGU

</td>

<td style="text-align:left;">

chr17

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-146-P1\_5p

</td>

<td style="text-align:right;">

\-1.77

</td>

<td style="text-align:right;">

0.29

</td>

<td style="text-align:left;">

4.42e-09

</td>

<td style="text-align:right;">

351

</td>

<td style="text-align:right;">

1453

</td>

<td style="text-align:left;">

hsa-mir-146a

</td>

<td style="text-align:left;">

MIR-146

</td>

<td style="text-align:left;">

GAGAACU

</td>

<td style="text-align:left;">

chr5

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-130-P2a\_3p

</td>

<td style="text-align:right;">

\-1.75

</td>

<td style="text-align:right;">

0.17

</td>

<td style="text-align:left;">

1.02e-24

</td>

<td style="text-align:right;">

92

</td>

<td style="text-align:right;">

409

</td>

<td style="text-align:left;">

hsa-mir-301a

</td>

<td style="text-align:left;">

MIR-130

</td>

<td style="text-align:left;">

AGUGCAA

</td>

<td style="text-align:left;">

chr17

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-29-P2a\_3p/P2b\_3p

</td>

<td style="text-align:right;">

\-1.48

</td>

<td style="text-align:right;">

0.19

</td>

<td style="text-align:left;">

9.97e-14

</td>

<td style="text-align:right;">

107

</td>

<td style="text-align:right;">

380

</td>

<td style="text-align:left;">

hsa-mir-29b-1

</td>

<td style="text-align:left;">

MIR-29

</td>

<td style="text-align:left;">

AGCACCA

</td>

<td style="text-align:left;">

chr7

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-17-P3c\_5p

</td>

<td style="text-align:right;">

\-1.47

</td>

<td style="text-align:right;">

0.13

</td>

<td style="text-align:left;">

8.74e-29

</td>

<td style="text-align:right;">

100

</td>

<td style="text-align:right;">

370

</td>

<td style="text-align:left;">

hsa-mir-106b

</td>

<td style="text-align:left;">

MIR-17

</td>

<td style="text-align:left;">

AAAGUGC

</td>

<td style="text-align:left;">

chr7

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-17-P1c\_5p

</td>

<td style="text-align:right;">

\-1.45

</td>

<td style="text-align:right;">

0.11

</td>

<td style="text-align:left;">

7.11e-37

</td>

<td style="text-align:right;">

647

</td>

<td style="text-align:right;">

2287

</td>

<td style="text-align:left;">

hsa-mir-93

</td>

<td style="text-align:left;">

MIR-17

</td>

<td style="text-align:left;">

AAAGUGC

</td>

<td style="text-align:left;">

chr7

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-95-P2\_3p

</td>

<td style="text-align:right;">

\-1.42

</td>

<td style="text-align:right;">

0.16

</td>

<td style="text-align:left;">

7.76e-18

</td>

<td style="text-align:right;">

40

</td>

<td style="text-align:right;">

142

</td>

<td style="text-align:left;">

hsa-mir-421

</td>

<td style="text-align:left;">

MIR-95

</td>

<td style="text-align:left;">

UCAACAG

</td>

<td style="text-align:left;">

chrX

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-17-P3a\_5p

</td>

<td style="text-align:right;">

\-1.41

</td>

<td style="text-align:right;">

0.19

</td>

<td style="text-align:left;">

1.99e-13

</td>

<td style="text-align:right;">

437

</td>

<td style="text-align:right;">

1513

</td>

<td style="text-align:left;">

hsa-mir-20a

</td>

<td style="text-align:left;">

MIR-17

</td>

<td style="text-align:left;">

AAAGUGC

</td>

<td style="text-align:left;">

chr13

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-188-P2\_5p

</td>

<td style="text-align:right;">

\-1.33

</td>

<td style="text-align:right;">

0.15

</td>

<td style="text-align:left;">

3.22e-17

</td>

<td style="text-align:right;">

279

</td>

<td style="text-align:right;">

914

</td>

<td style="text-align:left;">

hsa-mir-532

</td>

<td style="text-align:left;">

MIR-188

</td>

<td style="text-align:left;">

AUGCCUU

</td>

<td style="text-align:left;">

chrX

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Let-7-P2c2\_5p

</td>

<td style="text-align:right;">

\-1.26

</td>

<td style="text-align:right;">

0.14

</td>

<td style="text-align:left;">

4.74e-19

</td>

<td style="text-align:right;">

1738

</td>

<td style="text-align:right;">

5432

</td>

<td style="text-align:left;">

hsa-let-7i

</td>

<td style="text-align:left;">

LET-7

</td>

<td style="text-align:left;">

GAGGUAG

</td>

<td style="text-align:left;">

chr12

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Let-7-P2b3\_5p

</td>

<td style="text-align:right;">

\-1.22

</td>

<td style="text-align:right;">

0.15

</td>

<td style="text-align:left;">

7.97e-15

</td>

<td style="text-align:right;">

760

</td>

<td style="text-align:right;">

2149

</td>

<td style="text-align:left;">

hsa-mir-98

</td>

<td style="text-align:left;">

LET-7

</td>

<td style="text-align:left;">

GAGGUAG

</td>

<td style="text-align:left;">

chrX

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-210\_3p

</td>

<td style="text-align:right;">

\-1.16

</td>

<td style="text-align:right;">

0.23

</td>

<td style="text-align:left;">

5.10e-06

</td>

<td style="text-align:right;">

106

</td>

<td style="text-align:right;">

291

</td>

<td style="text-align:left;">

hsa-mir-210

</td>

<td style="text-align:left;">

MIR-210

</td>

<td style="text-align:left;">

UGUGCGU

</td>

<td style="text-align:left;">

chr11

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-425\_5p

</td>

<td style="text-align:right;">

\-1.15

</td>

<td style="text-align:right;">

0.15

</td>

<td style="text-align:left;">

2.77e-13

</td>

<td style="text-align:right;">

247

</td>

<td style="text-align:right;">

700

</td>

<td style="text-align:left;">

hsa-mir-425

</td>

<td style="text-align:left;">

MIR-425

</td>

<td style="text-align:left;">

AUGACAC

</td>

<td style="text-align:left;">

chr3

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-19-P1\_3p

</td>

<td style="text-align:right;">

\-1.12

</td>

<td style="text-align:right;">

0.19

</td>

<td style="text-align:left;">

1.87e-08

</td>

<td style="text-align:right;">

134

</td>

<td style="text-align:right;">

379

</td>

<td style="text-align:left;">

hsa-mir-19a

</td>

<td style="text-align:left;">

MIR-19

</td>

<td style="text-align:left;">

GUGCAAA

</td>

<td style="text-align:left;">

chr13

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-132-P1\_3p

</td>

<td style="text-align:right;">

\-1.06

</td>

<td style="text-align:right;">

0.17

</td>

<td style="text-align:left;">

6.29e-10

</td>

<td style="text-align:right;">

61

</td>

<td style="text-align:right;">

165

</td>

<td style="text-align:left;">

hsa-mir-132

</td>

<td style="text-align:left;">

MIR-132

</td>

<td style="text-align:left;">

AACAGUC

</td>

<td style="text-align:left;">

chr17

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-130-P4a\_3p

</td>

<td style="text-align:right;">

\-0.99

</td>

<td style="text-align:right;">

0.10

</td>

<td style="text-align:left;">

5.70e-21

</td>

<td style="text-align:right;">

70

</td>

<td style="text-align:right;">

177

</td>

<td style="text-align:left;">

hsa-mir-454

</td>

<td style="text-align:left;">

MIR-130

</td>

<td style="text-align:left;">

AGUGCAA

</td>

<td style="text-align:left;">

chr17

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-652\_3p

</td>

<td style="text-align:right;">

\-0.91

</td>

<td style="text-align:right;">

0.13

</td>

<td style="text-align:left;">

2.45e-11

</td>

<td style="text-align:right;">

44

</td>

<td style="text-align:right;">

106

</td>

<td style="text-align:left;">

hsa-mir-652

</td>

<td style="text-align:left;">

MIR-652

</td>

<td style="text-align:left;">

AUGGCGC

</td>

<td style="text-align:left;">

chrX

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-769\_5p

</td>

<td style="text-align:right;">

\-0.81

</td>

<td style="text-align:right;">

0.13

</td>

<td style="text-align:left;">

9.49e-10

</td>

<td style="text-align:right;">

190

</td>

<td style="text-align:right;">

439

</td>

<td style="text-align:left;">

hsa-mir-769

</td>

<td style="text-align:left;">

MIR-769

</td>

<td style="text-align:left;">

GAGACCU

</td>

<td style="text-align:left;">

chr19

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-362-P2\_3p/P4\_3p

</td>

<td style="text-align:right;">

\-0.69

</td>

<td style="text-align:right;">

0.17

</td>

<td style="text-align:left;">

7.46e-05

</td>

<td style="text-align:right;">

211

</td>

<td style="text-align:right;">

441

</td>

<td style="text-align:left;">

hsa-mir-500a

</td>

<td style="text-align:left;">

MIR-362

</td>

<td style="text-align:left;">

UGCACCU

</td>

<td style="text-align:left;">

chrX

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-181-P2a\_5p/P2b\_5p

</td>

<td style="text-align:right;">

\-0.67

</td>

<td style="text-align:right;">

0.15

</td>

<td style="text-align:left;">

1.77e-05

</td>

<td style="text-align:right;">

611

</td>

<td style="text-align:right;">

1251

</td>

<td style="text-align:left;">

hsa-mir-181b-1

</td>

<td style="text-align:left;">

MIR-181

</td>

<td style="text-align:left;">

ACAUUCA

</td>

<td style="text-align:left;">

chr1

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-17-P3c\_3p

</td>

<td style="text-align:right;">

\-0.64

</td>

<td style="text-align:right;">

0.14

</td>

<td style="text-align:left;">

5.82e-06

</td>

<td style="text-align:right;">

105

</td>

<td style="text-align:right;">

212

</td>

<td style="text-align:left;">

hsa-mir-106b

</td>

<td style="text-align:left;">

MIR-17

</td>

<td style="text-align:left;">

AAAGUGC

</td>

<td style="text-align:left;">

chr7

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-19-P2a\_3p/P2b\_3p

</td>

<td style="text-align:right;">

\-0.60

</td>

<td style="text-align:right;">

0.18

</td>

<td style="text-align:left;">

1.90e-03

</td>

<td style="text-align:right;">

677

</td>

<td style="text-align:right;">

1297

</td>

<td style="text-align:left;">

hsa-mir-19b-1

</td>

<td style="text-align:left;">

MIR-19

</td>

<td style="text-align:left;">

GUGCAAA

</td>

<td style="text-align:left;">

chr13

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

</tbody>

</table>

</div>

``` r
# Number of downregulated miRNA
signature_mirnas$number_downregulated
```

    ## [1] 54

## pCRC vs nLu

``` r
column='tissue.type'
tissue_type_A <- 'normal.lung'
tissue_type_B <- 'tumor.colorect'
norm_adj_up       = "None"
norm_adj_down     = "None"
pCRC_adj_up   = "None"
pCRC_adj_down = "None"

coef <- paste(column, tissue_type_A, 'vs', tissue_type_B, sep='_')
res <- DeseqResult(dds, column, coef, tissue_type_A, tissue_type_B,
                   lfc.Threshold, rpm.Threshold,
                   norm_adj_up,
                   norm_adj_down)

dict_sig_mirna[paste(coef, "up",   sep='_')] <- list(res$up_mirna)
dict_sig_mirna[paste(coef, "down", sep='_')] <- list(res$down_mirna)
res_res <- res$res
res_dict[coef] <- res_res
plotMA(res$res, alpha=0.05)
```

![](comet_analysis_report_automated_no_lfcthreshold_08.12.20_res_alpha_0.05_with_LFC_shrinkage_LFC_0.58_mirge3_7.0_files/figure-gfm/unnamed-chunk-26-1.png)<!-- -->

``` r
# Plot volcano plot
VolcanoPlot(res$res, coef, res$sig,
            res$up_mirna, res$down_mirna,
            norm_adj_up, norm_adj_down,
            pCRC_adj_up, pCRC_adj_down)
```

![](comet_analysis_report_automated_no_lfcthreshold_08.12.20_res_alpha_0.05_with_LFC_shrinkage_LFC_0.58_mirge3_7.0_files/figure-gfm/unnamed-chunk-26-2.png)<!-- -->

``` r
ExpressionPlot(res$res, res$rpm, coef, res$sig,
               tissue_type_A, tissue_type_B,
               res$up_mirna, res$down_mirna,
               norm_adj_up, norm_adj_down,
               pCRC_adj_up, pCRC_adj_down)
```

![](comet_analysis_report_automated_no_lfcthreshold_08.12.20_res_alpha_0.05_with_LFC_shrinkage_LFC_0.58_mirge3_7.0_files/figure-gfm/unnamed-chunk-26-3.png)<!-- -->

``` r
signature_mirnas <- SigList(res, dds, tissue_type_A, tissue_type_B, coef,
                            norm_adj_up, norm_adj_down, 
                            pCRC_adj_up, pCRC_adj_down)
# Print list upregulated miRNA
signature_mirnas$up_mirna
```

<div style="border: 1px solid #ddd; padding: 5px; overflow-x: scroll; width:2000px; ">

<table class="table table-striped table-hover table-condensed" style="margin-left: auto; margin-right: auto;">

<caption>

Upregulated in tissue.type\_normal.lung\_vs\_tumor.colorect

</caption>

<thead>

<tr>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

miRNA

</th>

<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">

LFC

</th>

<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">

lfcSE

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

FDR

</th>

<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">

RPM normal.lung

</th>

<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">

RPM tumor.colorect

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

miRBase\_ID

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

Family

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

Seed

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

Chr

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

Cell-Type Specific

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

Norm Background

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

pCRC Background

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

Hsa-Mir-335\_5p

</td>

<td style="text-align:right;">

1.10

</td>

<td style="text-align:right;">

0.19

</td>

<td style="text-align:left;">

1.98e-08

</td>

<td style="text-align:right;">

300

</td>

<td style="text-align:right;">

150

</td>

<td style="text-align:left;">

hsa-mir-335

</td>

<td style="text-align:left;">

MIR-335

</td>

<td style="text-align:left;">

CAAGAGC

</td>

<td style="text-align:left;">

chr7

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">Retinal
Epithelial Cell</span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-451\_5p

</td>

<td style="text-align:right;">

3.02

</td>

<td style="text-align:right;">

0.43

</td>

<td style="text-align:left;">

1.09e-11

</td>

<td style="text-align:right;">

12496

</td>

<td style="text-align:right;">

1565

</td>

<td style="text-align:left;">

hsa-mir-451a

</td>

<td style="text-align:left;">

MIR-451

</td>

<td style="text-align:left;">

AACCGUU

</td>

<td style="text-align:left;">

chr17

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">Red
Blood Cell</span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-144\_5p

</td>

<td style="text-align:right;">

2.99

</td>

<td style="text-align:right;">

0.37

</td>

<td style="text-align:left;">

6.58e-15

</td>

<td style="text-align:right;">

451

</td>

<td style="text-align:right;">

57

</td>

<td style="text-align:left;">

hsa-mir-144

</td>

<td style="text-align:left;">

MIR-144

</td>

<td style="text-align:left;">

GAUAUCA

</td>

<td style="text-align:left;">

chr17

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">Red
Blood Cell</span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-145\_5p

</td>

<td style="text-align:right;">

1.49

</td>

<td style="text-align:right;">

0.36

</td>

<td style="text-align:left;">

9.83e-05

</td>

<td style="text-align:right;">

2445

</td>

<td style="text-align:right;">

840

</td>

<td style="text-align:left;">

hsa-mir-145

</td>

<td style="text-align:left;">

MIR-145

</td>

<td style="text-align:left;">

UCCAGUU

</td>

<td style="text-align:left;">

chr5

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">Mesenchymal</span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-24-P1\_3p/P2\_3p

</td>

<td style="text-align:right;">

0.73

</td>

<td style="text-align:right;">

0.22

</td>

<td style="text-align:left;">

2.23e-03

</td>

<td style="text-align:right;">

1084

</td>

<td style="text-align:right;">

626

</td>

<td style="text-align:left;">

hsa-mir-24-2

</td>

<td style="text-align:left;">

MIR-24

</td>

<td style="text-align:left;">

GGCUCAG

</td>

<td style="text-align:left;">

chr19

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">Macrophage</span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-150\_5p

</td>

<td style="text-align:right;">

1.73

</td>

<td style="text-align:right;">

0.39

</td>

<td style="text-align:left;">

1.07e-05

</td>

<td style="text-align:right;">

973

</td>

<td style="text-align:right;">

280

</td>

<td style="text-align:left;">

hsa-mir-150

</td>

<td style="text-align:left;">

MIR-150

</td>

<td style="text-align:left;">

CUCCCAA

</td>

<td style="text-align:left;">

chr19

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">Lymphocyte</span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-15-P1a\_5p

</td>

<td style="text-align:right;">

0.96

</td>

<td style="text-align:right;">

0.17

</td>

<td style="text-align:left;">

2.82e-08

</td>

<td style="text-align:right;">

854

</td>

<td style="text-align:right;">

453

</td>

<td style="text-align:left;">

hsa-mir-15a

</td>

<td style="text-align:left;">

MIR-15

</td>

<td style="text-align:left;">

AGCAGCA

</td>

<td style="text-align:left;">

chr13

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">CD14+
Monocyte</span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-486\_5p

</td>

<td style="text-align:right;">

3.15

</td>

<td style="text-align:right;">

0.41

</td>

<td style="text-align:left;">

6.09e-14

</td>

<td style="text-align:right;">

22081

</td>

<td style="text-align:right;">

2013

</td>

<td style="text-align:left;">

hsa-mir-486-1

</td>

<td style="text-align:left;">

MIR-486

</td>

<td style="text-align:left;">

CCUGUAC

</td>

<td style="text-align:left;">

chr8

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">c(“Platelet”,
“Red Blood Cell”)</span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-126\_5p

</td>

<td style="text-align:right;">

3.27

</td>

<td style="text-align:right;">

0.23

</td>

<td style="text-align:left;">

4.58e-44

</td>

<td style="text-align:right;">

33620

</td>

<td style="text-align:right;">

3455

</td>

<td style="text-align:left;">

hsa-mir-126

</td>

<td style="text-align:left;">

MIR-126

</td>

<td style="text-align:left;">

AUUAUUA

</td>

<td style="text-align:left;">

chr9

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">c(“Endothelial
Cell”, “Platelet”)</span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-146-P2\_5p

</td>

<td style="text-align:right;">

1.39

</td>

<td style="text-align:right;">

0.37

</td>

<td style="text-align:left;">

3.72e-04

</td>

<td style="text-align:right;">

17435

</td>

<td style="text-align:right;">

6279

</td>

<td style="text-align:left;">

hsa-mir-146b

</td>

<td style="text-align:left;">

MIR-146

</td>

<td style="text-align:left;">

GAGAACU

</td>

<td style="text-align:left;">

chr10

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">c(“Dendritic
Cell”, “Lymphocyte”)</span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-342\_3p

</td>

<td style="text-align:right;">

2.20

</td>

<td style="text-align:right;">

0.28

</td>

<td style="text-align:left;">

2.02e-14

</td>

<td style="text-align:right;">

681

</td>

<td style="text-align:right;">

145

</td>

<td style="text-align:left;">

hsa-mir-342

</td>

<td style="text-align:left;">

MIR-342

</td>

<td style="text-align:left;">

CUCACAC

</td>

<td style="text-align:left;">

chr14

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">c(“Dendritic
Cell”, “Lymphocyte”, “Macrophage”)</span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-184\_3p

</td>

<td style="text-align:right;">

4.90

</td>

<td style="text-align:right;">

0.90

</td>

<td style="text-align:left;">

1.36e-07

</td>

<td style="text-align:right;">

111

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

hsa-mir-184

</td>

<td style="text-align:left;">

MIR-184

</td>

<td style="text-align:left;">

GGACGGA

</td>

<td style="text-align:left;">

chr15

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-34-P2a\_5p

</td>

<td style="text-align:right;">

4.41

</td>

<td style="text-align:right;">

0.45

</td>

<td style="text-align:left;">

6.03e-22

</td>

<td style="text-align:right;">

116

</td>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

hsa-mir-34b

</td>

<td style="text-align:left;">

MIR-34

</td>

<td style="text-align:left;">

GGCAGUG

</td>

<td style="text-align:left;">

chr11

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-34-P2b\_5p

</td>

<td style="text-align:right;">

4.40

</td>

<td style="text-align:right;">

0.39

</td>

<td style="text-align:left;">

1.18e-28

</td>

<td style="text-align:right;">

748

</td>

<td style="text-align:right;">

33

</td>

<td style="text-align:left;">

hsa-mir-34c

</td>

<td style="text-align:left;">

MIR-34

</td>

<td style="text-align:left;">

GGCAGUG

</td>

<td style="text-align:left;">

chr11

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-30-P1a\_5p

</td>

<td style="text-align:right;">

4.21

</td>

<td style="text-align:right;">

0.25

</td>

<td style="text-align:left;">

3.44e-60

</td>

<td style="text-align:right;">

60860

</td>

<td style="text-align:right;">

3165

</td>

<td style="text-align:left;">

hsa-mir-30a

</td>

<td style="text-align:left;">

MIR-30

</td>

<td style="text-align:left;">

GUAAACA

</td>

<td style="text-align:left;">

chr6

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-218-P1\_5p/P2\_5p

</td>

<td style="text-align:right;">

3.31

</td>

<td style="text-align:right;">

0.30

</td>

<td style="text-align:left;">

3.30e-27

</td>

<td style="text-align:right;">

374

</td>

<td style="text-align:right;">

37

</td>

<td style="text-align:left;">

hsa-mir-218-1

</td>

<td style="text-align:left;">

MIR-218

</td>

<td style="text-align:left;">

UGUGCUU

</td>

<td style="text-align:left;">

chr4

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-15-P2c\_5p

</td>

<td style="text-align:right;">

2.58

</td>

<td style="text-align:right;">

0.26

</td>

<td style="text-align:left;">

1.28e-21

</td>

<td style="text-align:right;">

1552

</td>

<td style="text-align:right;">

257

</td>

<td style="text-align:left;">

hsa-mir-195

</td>

<td style="text-align:left;">

MIR-15

</td>

<td style="text-align:left;">

AGCAGCA

</td>

<td style="text-align:left;">

chr17

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-338-P1\_3p

</td>

<td style="text-align:right;">

2.57

</td>

<td style="text-align:right;">

0.35

</td>

<td style="text-align:left;">

1.28e-12

</td>

<td style="text-align:right;">

605

</td>

<td style="text-align:right;">

98

</td>

<td style="text-align:left;">

hsa-mir-338

</td>

<td style="text-align:left;">

MIR-338

</td>

<td style="text-align:left;">

CCAGCAU

</td>

<td style="text-align:left;">

chr17

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-10-P2c\_5p

</td>

<td style="text-align:right;">

2.37

</td>

<td style="text-align:right;">

0.37

</td>

<td style="text-align:left;">

4.44e-10

</td>

<td style="text-align:right;">

804

</td>

<td style="text-align:right;">

142

</td>

<td style="text-align:left;">

hsa-mir-99a

</td>

<td style="text-align:left;">

MIR-10

</td>

<td style="text-align:left;">

ACCCGUA

</td>

<td style="text-align:left;">

chr21

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-10-P3b\_5p

</td>

<td style="text-align:right;">

2.28

</td>

<td style="text-align:right;">

0.30

</td>

<td style="text-align:left;">

5.11e-13

</td>

<td style="text-align:right;">

8561

</td>

<td style="text-align:right;">

1649

</td>

<td style="text-align:left;">

hsa-mir-125a

</td>

<td style="text-align:left;">

MIR-10

</td>

<td style="text-align:left;">

CCCUGAG

</td>

<td style="text-align:left;">

chr19

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-130-P1a\_3p

</td>

<td style="text-align:right;">

2.26

</td>

<td style="text-align:right;">

0.22

</td>

<td style="text-align:left;">

6.52e-23

</td>

<td style="text-align:right;">

1941

</td>

<td style="text-align:right;">

421

</td>

<td style="text-align:left;">

hsa-mir-130a

</td>

<td style="text-align:left;">

MIR-130

</td>

<td style="text-align:left;">

AGUGCAA

</td>

<td style="text-align:left;">

chr11

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Let-7-P1c\_5p

</td>

<td style="text-align:right;">

2.02

</td>

<td style="text-align:right;">

0.31

</td>

<td style="text-align:left;">

4.62e-10

</td>

<td style="text-align:right;">

2606

</td>

<td style="text-align:right;">

606

</td>

<td style="text-align:left;">

hsa-let-7c

</td>

<td style="text-align:left;">

LET-7

</td>

<td style="text-align:left;">

GAGGUAG

</td>

<td style="text-align:left;">

chr21

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-181-P1a\_5p/P1b\_5p

</td>

<td style="text-align:right;">

1.97

</td>

<td style="text-align:right;">

0.21

</td>

<td style="text-align:left;">

1.65e-19

</td>

<td style="text-align:right;">

80063

</td>

<td style="text-align:right;">

21233

</td>

<td style="text-align:left;">

hsa-mir-181a-1

</td>

<td style="text-align:left;">

MIR-181

</td>

<td style="text-align:left;">

ACAUUCA

</td>

<td style="text-align:left;">

chr1

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-30-P1c\_5p

</td>

<td style="text-align:right;">

1.87

</td>

<td style="text-align:right;">

0.18

</td>

<td style="text-align:left;">

4.93e-24

</td>

<td style="text-align:right;">

34088

</td>

<td style="text-align:right;">

9703

</td>

<td style="text-align:left;">

hsa-mir-30d

</td>

<td style="text-align:left;">

MIR-30

</td>

<td style="text-align:left;">

GUAAACA

</td>

<td style="text-align:left;">

chr8

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-101-P1\_3p/P2\_3p

</td>

<td style="text-align:right;">

1.85

</td>

<td style="text-align:right;">

0.20

</td>

<td style="text-align:left;">

1.47e-19

</td>

<td style="text-align:right;">

16084

</td>

<td style="text-align:right;">

4638

</td>

<td style="text-align:left;">

hsa-mir-101-1

</td>

<td style="text-align:left;">

MIR-101

</td>

<td style="text-align:left;">

UACAGUA

</td>

<td style="text-align:left;">

chr1

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-26-P1\_5p/P2\_5p

</td>

<td style="text-align:right;">

1.82

</td>

<td style="text-align:right;">

0.15

</td>

<td style="text-align:left;">

1.88e-32

</td>

<td style="text-align:right;">

123947

</td>

<td style="text-align:right;">

35584

</td>

<td style="text-align:left;">

hsa-mir-26b

</td>

<td style="text-align:left;">

MIR-26

</td>

<td style="text-align:left;">

UCAAGUA

</td>

<td style="text-align:left;">

chr2

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-10-P3a\_5p

</td>

<td style="text-align:right;">

1.72

</td>

<td style="text-align:right;">

0.35

</td>

<td style="text-align:left;">

2.27e-06

</td>

<td style="text-align:right;">

507

</td>

<td style="text-align:right;">

145

</td>

<td style="text-align:left;">

hsa-mir-125b-1

</td>

<td style="text-align:left;">

MIR-10

</td>

<td style="text-align:left;">

CCCUGAG

</td>

<td style="text-align:left;">

chr11

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-181-P2a\_5p/P2b\_5p

</td>

<td style="text-align:right;">

1.68

</td>

<td style="text-align:right;">

0.20

</td>

<td style="text-align:left;">

4.86e-16

</td>

<td style="text-align:right;">

3799

</td>

<td style="text-align:right;">

1251

</td>

<td style="text-align:left;">

hsa-mir-181b-1

</td>

<td style="text-align:left;">

MIR-181

</td>

<td style="text-align:left;">

ACAUUCA

</td>

<td style="text-align:left;">

chr1

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-10-P3c\_5p

</td>

<td style="text-align:right;">

1.60

</td>

<td style="text-align:right;">

0.34

</td>

<td style="text-align:left;">

3.99e-06

</td>

<td style="text-align:right;">

822

</td>

<td style="text-align:right;">

239

</td>

<td style="text-align:left;">

hsa-mir-125b-2

</td>

<td style="text-align:left;">

MIR-10

</td>

<td style="text-align:left;">

CCCUGAG

</td>

<td style="text-align:left;">

chr21

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-140\_3p

</td>

<td style="text-align:right;">

1.56

</td>

<td style="text-align:right;">

0.17

</td>

<td style="text-align:left;">

1.01e-18

</td>

<td style="text-align:right;">

2195

</td>

<td style="text-align:right;">

784

</td>

<td style="text-align:left;">

hsa-mir-140

</td>

<td style="text-align:left;">

MIR-140

</td>

<td style="text-align:left;">

CCACAGG

</td>

<td style="text-align:left;">

chr16

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-15-P1c\_5p

</td>

<td style="text-align:right;">

1.56

</td>

<td style="text-align:right;">

0.25

</td>

<td style="text-align:left;">

2.16e-09

</td>

<td style="text-align:right;">

405

</td>

<td style="text-align:right;">

143

</td>

<td style="text-align:left;">

hsa-mir-497

</td>

<td style="text-align:left;">

MIR-15

</td>

<td style="text-align:left;">

AGCAGCA

</td>

<td style="text-align:left;">

chr17

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-30-P2a\_5p/P2b\_5p/P2c\_5p

</td>

<td style="text-align:right;">

1.54

</td>

<td style="text-align:right;">

0.17

</td>

<td style="text-align:left;">

1.46e-19

</td>

<td style="text-align:right;">

7603

</td>

<td style="text-align:right;">

2793

</td>

<td style="text-align:left;">

hsa-mir-30c-2

</td>

<td style="text-align:left;">

MIR-30

</td>

<td style="text-align:left;">

GUAAACA

</td>

<td style="text-align:left;">

chr6

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-92-P2b\_3p

</td>

<td style="text-align:right;">

1.53

</td>

<td style="text-align:right;">

0.36

</td>

<td style="text-align:left;">

8.68e-05

</td>

<td style="text-align:right;">

137

</td>

<td style="text-align:right;">

47

</td>

<td style="text-align:left;">

hsa-mir-363

</td>

<td style="text-align:left;">

MIR-92

</td>

<td style="text-align:left;">

AUUGCAC

</td>

<td style="text-align:left;">

chrX

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-15-P2a\_5p/P2b\_5p

</td>

<td style="text-align:right;">

1.48

</td>

<td style="text-align:right;">

0.19

</td>

<td style="text-align:left;">

3.24e-14

</td>

<td style="text-align:right;">

12836

</td>

<td style="text-align:right;">

4522

</td>

<td style="text-align:left;">

hsa-mir-16-1

</td>

<td style="text-align:left;">

MIR-15

</td>

<td style="text-align:left;">

AGCAGCA

</td>

<td style="text-align:left;">

chr13

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-26-P3\_5p

</td>

<td style="text-align:right;">

1.37

</td>

<td style="text-align:right;">

0.22

</td>

<td style="text-align:left;">

1.19e-09

</td>

<td style="text-align:right;">

9029

</td>

<td style="text-align:right;">

3467

</td>

<td style="text-align:left;">

hsa-mir-26a-1

</td>

<td style="text-align:left;">

MIR-26

</td>

<td style="text-align:left;">

UCAAGUA

</td>

<td style="text-align:left;">

chr3

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-10-P2a\_5p

</td>

<td style="text-align:right;">

1.32

</td>

<td style="text-align:right;">

0.36

</td>

<td style="text-align:left;">

6.59e-04

</td>

<td style="text-align:right;">

4056

</td>

<td style="text-align:right;">

1724

</td>

<td style="text-align:left;">

hsa-mir-100

</td>

<td style="text-align:left;">

MIR-10

</td>

<td style="text-align:left;">

ACCCGUA

</td>

<td style="text-align:left;">

chr11

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Let-7-P2b1\_5p

</td>

<td style="text-align:right;">

1.23

</td>

<td style="text-align:right;">

0.14

</td>

<td style="text-align:left;">

9.88e-19

</td>

<td style="text-align:right;">

6253

</td>

<td style="text-align:right;">

2784

</td>

<td style="text-align:left;">

hsa-let-7f-1

</td>

<td style="text-align:left;">

LET-7

</td>

<td style="text-align:left;">

GAGGUAG

</td>

<td style="text-align:left;">

chr9

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-331\_3p

</td>

<td style="text-align:right;">

1.11

</td>

<td style="text-align:right;">

0.29

</td>

<td style="text-align:left;">

2.07e-04

</td>

<td style="text-align:right;">

128

</td>

<td style="text-align:right;">

60

</td>

<td style="text-align:left;">

hsa-mir-331

</td>

<td style="text-align:left;">

MIR-331

</td>

<td style="text-align:left;">

CCCCUGG

</td>

<td style="text-align:left;">

chr12

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-10-P2b\_5p

</td>

<td style="text-align:right;">

1.06

</td>

<td style="text-align:right;">

0.42

</td>

<td style="text-align:left;">

2.28e-02

</td>

<td style="text-align:right;">

5968

</td>

<td style="text-align:right;">

2594

</td>

<td style="text-align:left;">

hsa-mir-99b

</td>

<td style="text-align:left;">

MIR-10

</td>

<td style="text-align:left;">

ACCCGUA

</td>

<td style="text-align:left;">

chr19

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-148-P2\_3p

</td>

<td style="text-align:right;">

0.98

</td>

<td style="text-align:right;">

0.34

</td>

<td style="text-align:left;">

6.10e-03

</td>

<td style="text-align:right;">

270

</td>

<td style="text-align:right;">

137

</td>

<td style="text-align:left;">

hsa-mir-152

</td>

<td style="text-align:left;">

MIR-148

</td>

<td style="text-align:left;">

CAGUGCA

</td>

<td style="text-align:left;">

chr17

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-29-P1b\_3p

</td>

<td style="text-align:right;">

0.96

</td>

<td style="text-align:right;">

0.23

</td>

<td style="text-align:left;">

6.29e-05

</td>

<td style="text-align:right;">

607

</td>

<td style="text-align:right;">

331

</td>

<td style="text-align:left;">

hsa-mir-29c

</td>

<td style="text-align:left;">

MIR-29

</td>

<td style="text-align:left;">

AGCACCA

</td>

<td style="text-align:left;">

chr1

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-652\_3p

</td>

<td style="text-align:right;">

0.88

</td>

<td style="text-align:right;">

0.18

</td>

<td style="text-align:left;">

1.77e-06

</td>

<td style="text-align:right;">

187

</td>

<td style="text-align:right;">

106

</td>

<td style="text-align:left;">

hsa-mir-652

</td>

<td style="text-align:left;">

MIR-652

</td>

<td style="text-align:left;">

AUGGCGC

</td>

<td style="text-align:left;">

chrX

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-744\_5p

</td>

<td style="text-align:right;">

0.88

</td>

<td style="text-align:right;">

0.24

</td>

<td style="text-align:left;">

5.22e-04

</td>

<td style="text-align:right;">

109

</td>

<td style="text-align:right;">

61

</td>

<td style="text-align:left;">

hsa-mir-744

</td>

<td style="text-align:left;">

MIR-744

</td>

<td style="text-align:left;">

GCGGGGC

</td>

<td style="text-align:left;">

chr17

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Let-7-P2b2\_5p

</td>

<td style="text-align:right;">

0.84

</td>

<td style="text-align:right;">

0.17

</td>

<td style="text-align:left;">

3.85e-06

</td>

<td style="text-align:right;">

6575

</td>

<td style="text-align:right;">

3787

</td>

<td style="text-align:left;">

hsa-let-7b

</td>

<td style="text-align:left;">

LET-7

</td>

<td style="text-align:left;">

GAGGUAG

</td>

<td style="text-align:left;">

chr22

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-27-P1\_3p/P2\_3p

</td>

<td style="text-align:right;">

0.80

</td>

<td style="text-align:right;">

0.16

</td>

<td style="text-align:left;">

4.13e-06

</td>

<td style="text-align:right;">

50854

</td>

<td style="text-align:right;">

29444

</td>

<td style="text-align:left;">

hsa-mir-27a

</td>

<td style="text-align:left;">

MIR-27

</td>

<td style="text-align:left;">

UCACAGU

</td>

<td style="text-align:left;">

chr19

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-374-P2\_5p

</td>

<td style="text-align:right;">

0.74

</td>

<td style="text-align:right;">

0.22

</td>

<td style="text-align:left;">

1.39e-03

</td>

<td style="text-align:right;">

148

</td>

<td style="text-align:right;">

95

</td>

<td style="text-align:left;">

hsa-mir-374b

</td>

<td style="text-align:left;">

MIR-374

</td>

<td style="text-align:left;">

UAUAAUA

</td>

<td style="text-align:left;">

chrX

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-23-P1\_3p/P2\_3p

</td>

<td style="text-align:right;">

0.68

</td>

<td style="text-align:right;">

0.17

</td>

<td style="text-align:left;">

1.35e-04

</td>

<td style="text-align:right;">

5786

</td>

<td style="text-align:right;">

3766

</td>

<td style="text-align:left;">

hsa-mir-23a

</td>

<td style="text-align:left;">

MIR-23

</td>

<td style="text-align:left;">

UCACAUU

</td>

<td style="text-align:left;">

chr19

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-28-P2\_3p

</td>

<td style="text-align:right;">

0.65

</td>

<td style="text-align:right;">

0.19

</td>

<td style="text-align:left;">

1.33e-03

</td>

<td style="text-align:right;">

6007

</td>

<td style="text-align:right;">

3830

</td>

<td style="text-align:left;">

hsa-mir-151a

</td>

<td style="text-align:left;">

MIR-28

</td>

<td style="text-align:left;">

CGAGGAG

</td>

<td style="text-align:left;">

chr8

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

</tbody>

</table>

</div>

``` r
# Number of upregulated miRNA
signature_mirnas$number_upregulated
```

    ## [1] 48

``` r
# Print list downregulated miRNA
signature_mirnas$down_mirna
```

<div style="border: 1px solid #ddd; padding: 5px; overflow-x: scroll; width:2000px; ">

<table class="table table-striped table-hover table-condensed" style="margin-left: auto; margin-right: auto;">

<caption>

Downregulated in tissue.type\_normal.lung\_vs\_tumor.colorect

</caption>

<thead>

<tr>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

miRNA

</th>

<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">

LFC

</th>

<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">

lfcSE

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

FDR

</th>

<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">

RPM normal.lung

</th>

<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">

RPM tumor.colorect

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

miRBase\_ID

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

Family

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

Seed

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

Chr

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

Cell-Type Specific

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

Norm Background

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

pCRC Background

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

Hsa-Mir-128-P1\_3p/P2\_3p

</td>

<td style="text-align:right;">

\-0.86

</td>

<td style="text-align:right;">

0.18

</td>

<td style="text-align:left;">

6.72e-06

</td>

<td style="text-align:right;">

116

</td>

<td style="text-align:right;">

233

</td>

<td style="text-align:left;">

hsa-mir-128-1

</td>

<td style="text-align:left;">

MIR-128

</td>

<td style="text-align:left;">

CACAGUG

</td>

<td style="text-align:left;">

chr2

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">Neural</span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-192-P1\_5p/P2\_5p

</td>

<td style="text-align:right;">

\-6.50

</td>

<td style="text-align:right;">

0.36

</td>

<td style="text-align:left;">

4.31e-69

</td>

<td style="text-align:right;">

142

</td>

<td style="text-align:right;">

13717

</td>

<td style="text-align:left;">

hsa-mir-192

</td>

<td style="text-align:left;">

MIR-192

</td>

<td style="text-align:left;">

UGACCUA

</td>

<td style="text-align:left;">

chr11

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">Epithelial
Cell</span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-8-P2a\_3p

</td>

<td style="text-align:right;">

\-3.88

</td>

<td style="text-align:right;">

0.26

</td>

<td style="text-align:left;">

9.46e-48

</td>

<td style="text-align:right;">

618

</td>

<td style="text-align:right;">

10330

</td>

<td style="text-align:left;">

hsa-mir-200b

</td>

<td style="text-align:left;">

MIR-8

</td>

<td style="text-align:left;">

AAUACUG

</td>

<td style="text-align:left;">

chr1

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">Epithelial
Cell</span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-8-P2b\_3p

</td>

<td style="text-align:right;">

\-1.75

</td>

<td style="text-align:right;">

0.26

</td>

<td style="text-align:left;">

7.14e-11

</td>

<td style="text-align:right;">

766

</td>

<td style="text-align:right;">

2727

</td>

<td style="text-align:left;">

hsa-mir-200c

</td>

<td style="text-align:left;">

MIR-8

</td>

<td style="text-align:left;">

AAUACUG

</td>

<td style="text-align:left;">

chr12

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">Epithelial
Cell</span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-17-P1a\_5p/P1b\_5p

</td>

<td style="text-align:right;">

\-2.17

</td>

<td style="text-align:right;">

0.21

</td>

<td style="text-align:left;">

9.51e-23

</td>

<td style="text-align:right;">

184

</td>

<td style="text-align:right;">

900

</td>

<td style="text-align:left;">

hsa-mir-17

</td>

<td style="text-align:left;">

MIR-17

</td>

<td style="text-align:left;">

AAAGUGC

</td>

<td style="text-align:left;">

chr13

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">CD14+
Monocyte</span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-7-P1\_5p/P2\_5p/P3\_5p

</td>

<td style="text-align:right;">

\-5.51

</td>

<td style="text-align:right;">

0.42

</td>

<td style="text-align:left;">

5.13e-38

</td>

<td style="text-align:right;">

3

</td>

<td style="text-align:right;">

199

</td>

<td style="text-align:left;">

hsa-mir-7-1

</td>

<td style="text-align:left;">

MIR-7

</td>

<td style="text-align:left;">

GGAAGAC

</td>

<td style="text-align:left;">

chr9

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">c(“Islet
Cell”, “Neural”)</span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-375\_3p

</td>

<td style="text-align:right;">

\-1.22

</td>

<td style="text-align:right;">

0.43

</td>

<td style="text-align:left;">

1.11e-02

</td>

<td style="text-align:right;">

2346

</td>

<td style="text-align:right;">

5263

</td>

<td style="text-align:left;">

hsa-mir-375

</td>

<td style="text-align:left;">

MIR-375

</td>

<td style="text-align:left;">

UUGUUCG

</td>

<td style="text-align:left;">

chr2

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">c(“Epithelial
Cell”, “Islet Cell”, “Neural”)</span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-196-P1\_5p/P2\_5p

</td>

<td style="text-align:right;">

\-6.60

</td>

<td style="text-align:right;">

0.38

</td>

<td style="text-align:left;">

7.65e-66

</td>

<td style="text-align:right;">

3

</td>

<td style="text-align:right;">

319

</td>

<td style="text-align:left;">

hsa-mir-196a-1

</td>

<td style="text-align:left;">

MIR-196

</td>

<td style="text-align:left;">

AGGUAGU

</td>

<td style="text-align:left;">

chr17

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-194-P1\_5p/P2\_5p

</td>

<td style="text-align:right;">

\-6.59

</td>

<td style="text-align:right;">

0.31

</td>

<td style="text-align:left;">

7.31e-100

</td>

<td style="text-align:right;">

58

</td>

<td style="text-align:right;">

6167

</td>

<td style="text-align:left;">

hsa-mir-194-2

</td>

<td style="text-align:left;">

MIR-194

</td>

<td style="text-align:left;">

GUAACAG

</td>

<td style="text-align:left;">

chr11

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-196-P3\_5p

</td>

<td style="text-align:right;">

\-6.49

</td>

<td style="text-align:right;">

0.43

</td>

<td style="text-align:left;">

9.67e-50

</td>

<td style="text-align:right;">

7

</td>

<td style="text-align:right;">

828

</td>

<td style="text-align:left;">

hsa-mir-196b

</td>

<td style="text-align:left;">

MIR-196

</td>

<td style="text-align:left;">

AGGUAGU

</td>

<td style="text-align:left;">

chr7

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-192-P1\_5p

</td>

<td style="text-align:right;">

\-6.17

</td>

<td style="text-align:right;">

0.31

</td>

<td style="text-align:left;">

1.45e-85

</td>

<td style="text-align:right;">

1684

</td>

<td style="text-align:right;">

132000

</td>

<td style="text-align:left;">

hsa-mir-192

</td>

<td style="text-align:left;">

MIR-192

</td>

<td style="text-align:left;">

UGACCUA

</td>

<td style="text-align:left;">

chr11

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-577\_5p

</td>

<td style="text-align:right;">

\-5.61

</td>

<td style="text-align:right;">

0.37

</td>

<td style="text-align:left;">

3.74e-49

</td>

<td style="text-align:right;">

3

</td>

<td style="text-align:right;">

142

</td>

<td style="text-align:left;">

hsa-mir-577

</td>

<td style="text-align:left;">

MIR-577

</td>

<td style="text-align:left;">

UAGAUAA

</td>

<td style="text-align:left;">

chr4

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-31\_5p

</td>

<td style="text-align:right;">

\-4.10

</td>

<td style="text-align:right;">

0.59

</td>

<td style="text-align:left;">

3.10e-12

</td>

<td style="text-align:right;">

24

</td>

<td style="text-align:right;">

627

</td>

<td style="text-align:left;">

hsa-mir-31

</td>

<td style="text-align:left;">

MIR-31

</td>

<td style="text-align:left;">

GGCAAGA

</td>

<td style="text-align:left;">

chr9

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-8-P3a\_3p

</td>

<td style="text-align:right;">

\-2.93

</td>

<td style="text-align:right;">

0.27

</td>

<td style="text-align:left;">

4.00e-27

</td>

<td style="text-align:right;">

262

</td>

<td style="text-align:right;">

2303

</td>

<td style="text-align:left;">

hsa-mir-429

</td>

<td style="text-align:left;">

MIR-8

</td>

<td style="text-align:left;">

AAUACUG

</td>

<td style="text-align:left;">

chr1

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-17-P3a\_5p

</td>

<td style="text-align:right;">

\-2.82

</td>

<td style="text-align:right;">

0.26

</td>

<td style="text-align:left;">

3.30e-27

</td>

<td style="text-align:right;">

194

</td>

<td style="text-align:right;">

1513

</td>

<td style="text-align:left;">

hsa-mir-20a

</td>

<td style="text-align:left;">

MIR-17

</td>

<td style="text-align:left;">

AAAGUGC

</td>

<td style="text-align:left;">

chr13

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-135-P3\_5p

</td>

<td style="text-align:right;">

\-2.73

</td>

<td style="text-align:right;">

0.37

</td>

<td style="text-align:left;">

6.36e-13

</td>

<td style="text-align:right;">

20

</td>

<td style="text-align:right;">

155

</td>

<td style="text-align:left;">

hsa-mir-135b

</td>

<td style="text-align:left;">

MIR-135

</td>

<td style="text-align:left;">

AUGGCUU

</td>

<td style="text-align:left;">

chr1

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-127\_3p

</td>

<td style="text-align:right;">

\-2.70

</td>

<td style="text-align:right;">

0.36

</td>

<td style="text-align:left;">

7.71e-13

</td>

<td style="text-align:right;">

327

</td>

<td style="text-align:right;">

2263

</td>

<td style="text-align:left;">

hsa-mir-127

</td>

<td style="text-align:left;">

MIR-127

</td>

<td style="text-align:left;">

CGGAUCC

</td>

<td style="text-align:left;">

chr14

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-190-P1\_5p

</td>

<td style="text-align:right;">

\-2.56

</td>

<td style="text-align:right;">

0.29

</td>

<td style="text-align:left;">

1.90e-17

</td>

<td style="text-align:right;">

22

</td>

<td style="text-align:right;">

140

</td>

<td style="text-align:left;">

hsa-mir-190a

</td>

<td style="text-align:left;">

MIR-190

</td>

<td style="text-align:left;">

GAUAUGU

</td>

<td style="text-align:left;">

chr15

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-224\_5p

</td>

<td style="text-align:right;">

\-2.55

</td>

<td style="text-align:right;">

0.30

</td>

<td style="text-align:left;">

2.30e-16

</td>

<td style="text-align:right;">

59

</td>

<td style="text-align:right;">

395

</td>

<td style="text-align:left;">

hsa-mir-224

</td>

<td style="text-align:left;">

MIR-224

</td>

<td style="text-align:left;">

AAGUCAC

</td>

<td style="text-align:left;">

chrX

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-154-P36\_3p

</td>

<td style="text-align:right;">

\-2.22

</td>

<td style="text-align:right;">

0.23

</td>

<td style="text-align:left;">

2.04e-21

</td>

<td style="text-align:right;">

39

</td>

<td style="text-align:right;">

195

</td>

<td style="text-align:left;">

hsa-mir-409

</td>

<td style="text-align:left;">

MIR-154

</td>

<td style="text-align:left;">

AAUGUUG

</td>

<td style="text-align:left;">

chr14

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-8-P1a\_3p

</td>

<td style="text-align:right;">

\-2.20

</td>

<td style="text-align:right;">

0.26

</td>

<td style="text-align:left;">

7.55e-16

</td>

<td style="text-align:right;">

353

</td>

<td style="text-align:right;">

1772

</td>

<td style="text-align:left;">

hsa-mir-200a

</td>

<td style="text-align:left;">

MIR-8

</td>

<td style="text-align:left;">

AACACUG

</td>

<td style="text-align:left;">

chr1

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-21\_5p

</td>

<td style="text-align:right;">

\-2.20

</td>

<td style="text-align:right;">

0.19

</td>

<td style="text-align:left;">

1.45e-28

</td>

<td style="text-align:right;">

21349

</td>

<td style="text-align:right;">

105735

</td>

<td style="text-align:left;">

hsa-mir-21

</td>

<td style="text-align:left;">

MIR-21

</td>

<td style="text-align:left;">

AGCUUAU

</td>

<td style="text-align:left;">

chr17

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-15-P1d\_5p

</td>

<td style="text-align:right;">

\-2.15

</td>

<td style="text-align:right;">

0.34

</td>

<td style="text-align:left;">

9.48e-10

</td>

<td style="text-align:right;">

46

</td>

<td style="text-align:right;">

235

</td>

<td style="text-align:left;">

hsa-mir-424

</td>

<td style="text-align:left;">

MIR-15

</td>

<td style="text-align:left;">

AGCAGCA

</td>

<td style="text-align:left;">

chrX

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-17-P2a\_5p

</td>

<td style="text-align:right;">

\-2.09

</td>

<td style="text-align:right;">

0.31

</td>

<td style="text-align:left;">

5.05e-11

</td>

<td style="text-align:right;">

33

</td>

<td style="text-align:right;">

168

</td>

<td style="text-align:left;">

hsa-mir-18a

</td>

<td style="text-align:left;">

MIR-17

</td>

<td style="text-align:left;">

AAGGUGC

</td>

<td style="text-align:left;">

chr13

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-130-P1b\_3p

</td>

<td style="text-align:right;">

\-1.94

</td>

<td style="text-align:right;">

0.22

</td>

<td style="text-align:left;">

2.67e-18

</td>

<td style="text-align:right;">

28

</td>

<td style="text-align:right;">

116

</td>

<td style="text-align:left;">

hsa-mir-130b

</td>

<td style="text-align:left;">

MIR-130

</td>

<td style="text-align:left;">

AGUGCAA

</td>

<td style="text-align:left;">

chr22

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-96-P3\_5p

</td>

<td style="text-align:right;">

\-1.94

</td>

<td style="text-align:right;">

0.37

</td>

<td style="text-align:left;">

3.95e-07

</td>

<td style="text-align:right;">

291

</td>

<td style="text-align:right;">

1063

</td>

<td style="text-align:left;">

hsa-mir-183

</td>

<td style="text-align:left;">

MIR-96

</td>

<td style="text-align:left;">

AUGGCAC

</td>

<td style="text-align:left;">

chr7

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-96-P2\_5p

</td>

<td style="text-align:right;">

\-1.91

</td>

<td style="text-align:right;">

0.27

</td>

<td style="text-align:left;">

1.11e-11

</td>

<td style="text-align:right;">

2685

</td>

<td style="text-align:right;">

10417

</td>

<td style="text-align:left;">

hsa-mir-182

</td>

<td style="text-align:left;">

MIR-96

</td>

<td style="text-align:left;">

UUGGCAA

</td>

<td style="text-align:left;">

chr7

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-154-P23\_3p

</td>

<td style="text-align:right;">

\-1.74

</td>

<td style="text-align:right;">

0.22

</td>

<td style="text-align:left;">

1.04e-13

</td>

<td style="text-align:right;">

41

</td>

<td style="text-align:right;">

148

</td>

<td style="text-align:left;">

hsa-mir-654

</td>

<td style="text-align:left;">

MIR-154

</td>

<td style="text-align:left;">

AUGUCUG

</td>

<td style="text-align:left;">

chr14

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-19-P1\_3p

</td>

<td style="text-align:right;">

\-1.65

</td>

<td style="text-align:right;">

0.26

</td>

<td style="text-align:left;">

1.09e-09

</td>

<td style="text-align:right;">

108

</td>

<td style="text-align:right;">

379

</td>

<td style="text-align:left;">

hsa-mir-19a

</td>

<td style="text-align:left;">

MIR-19

</td>

<td style="text-align:left;">

GUGCAAA

</td>

<td style="text-align:left;">

chr13

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-425\_5p

</td>

<td style="text-align:right;">

\-1.65

</td>

<td style="text-align:right;">

0.21

</td>

<td style="text-align:left;">

2.02e-14

</td>

<td style="text-align:right;">

213

</td>

<td style="text-align:right;">

700

</td>

<td style="text-align:left;">

hsa-mir-425

</td>

<td style="text-align:left;">

MIR-425

</td>

<td style="text-align:left;">

AUGACAC

</td>

<td style="text-align:left;">

chr3

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-154-P12\_3p

</td>

<td style="text-align:right;">

\-1.63

</td>

<td style="text-align:right;">

0.24

</td>

<td style="text-align:left;">

1.28e-10

</td>

<td style="text-align:right;">

33

</td>

<td style="text-align:right;">

115

</td>

<td style="text-align:left;">

hsa-mir-410

</td>

<td style="text-align:left;">

MIR-154

</td>

<td style="text-align:left;">

AUAUAAC

</td>

<td style="text-align:left;">

chr14

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-10-P1b\_5p

</td>

<td style="text-align:right;">

\-1.55

</td>

<td style="text-align:right;">

0.41

</td>

<td style="text-align:left;">

4.05e-04

</td>

<td style="text-align:right;">

14125

</td>

<td style="text-align:right;">

39286

</td>

<td style="text-align:left;">

hsa-mir-10b

</td>

<td style="text-align:left;">

MIR-10

</td>

<td style="text-align:left;">

ACCCUGU

</td>

<td style="text-align:left;">

chr2

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-378\_3p

</td>

<td style="text-align:right;">

\-1.45

</td>

<td style="text-align:right;">

0.22

</td>

<td style="text-align:left;">

9.02e-11

</td>

<td style="text-align:right;">

2189

</td>

<td style="text-align:right;">

6426

</td>

<td style="text-align:left;">

hsa-mir-378a

</td>

<td style="text-align:left;">

MIR-378

</td>

<td style="text-align:left;">

CUGGACU

</td>

<td style="text-align:left;">

chr5

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-584\_5p

</td>

<td style="text-align:right;">

\-1.38

</td>

<td style="text-align:right;">

0.34

</td>

<td style="text-align:left;">

1.18e-04

</td>

<td style="text-align:right;">

39

</td>

<td style="text-align:right;">

110

</td>

<td style="text-align:left;">

hsa-mir-584

</td>

<td style="text-align:left;">

MIR-584

</td>

<td style="text-align:left;">

UAUGGUU

</td>

<td style="text-align:left;">

chr5

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-92-P1a\_3p/P1b\_3p

</td>

<td style="text-align:right;">

\-1.27

</td>

<td style="text-align:right;">

0.25

</td>

<td style="text-align:left;">

8.94e-07

</td>

<td style="text-align:right;">

15850

</td>

<td style="text-align:right;">

39320

</td>

<td style="text-align:left;">

hsa-mir-92a-1

</td>

<td style="text-align:left;">

MIR-92

</td>

<td style="text-align:left;">

AUUGCAC

</td>

<td style="text-align:left;">

chr13

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-17-P3c\_5p

</td>

<td style="text-align:right;">

\-1.21

</td>

<td style="text-align:right;">

0.18

</td>

<td style="text-align:left;">

2.57e-11

</td>

<td style="text-align:right;">

143

</td>

<td style="text-align:right;">

370

</td>

<td style="text-align:left;">

hsa-mir-106b

</td>

<td style="text-align:left;">

MIR-17

</td>

<td style="text-align:left;">

AAAGUGC

</td>

<td style="text-align:left;">

chr7

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-17-P1c\_5p

</td>

<td style="text-align:right;">

\-1.06

</td>

<td style="text-align:right;">

0.15

</td>

<td style="text-align:left;">

2.20e-11

</td>

<td style="text-align:right;">

1024

</td>

<td style="text-align:right;">

2287

</td>

<td style="text-align:left;">

hsa-mir-93

</td>

<td style="text-align:left;">

MIR-17

</td>

<td style="text-align:left;">

AAAGUGC

</td>

<td style="text-align:left;">

chr7

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-203\_3p

</td>

<td style="text-align:right;">

\-1.03

</td>

<td style="text-align:right;">

0.32

</td>

<td style="text-align:left;">

2.53e-03

</td>

<td style="text-align:right;">

1068

</td>

<td style="text-align:right;">

2307

</td>

<td style="text-align:left;">

hsa-mir-203a

</td>

<td style="text-align:left;">

MIR-203

</td>

<td style="text-align:left;">

UGAAAUG

</td>

<td style="text-align:left;">

chr14

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-210\_3p

</td>

<td style="text-align:right;">

\-1.01

</td>

<td style="text-align:right;">

0.31

</td>

<td style="text-align:left;">

4.02e-03

</td>

<td style="text-align:right;">

138

</td>

<td style="text-align:right;">

291

</td>

<td style="text-align:left;">

hsa-mir-210

</td>

<td style="text-align:left;">

MIR-210

</td>

<td style="text-align:left;">

UGUGCGU

</td>

<td style="text-align:left;">

chr11

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-95-P2\_3p

</td>

<td style="text-align:right;">

\-0.99

</td>

<td style="text-align:right;">

0.22

</td>

<td style="text-align:left;">

1.35e-05

</td>

<td style="text-align:right;">

65

</td>

<td style="text-align:right;">

142

</td>

<td style="text-align:left;">

hsa-mir-421

</td>

<td style="text-align:left;">

MIR-95

</td>

<td style="text-align:left;">

UCAACAG

</td>

<td style="text-align:left;">

chrX

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-154-P13\_5p

</td>

<td style="text-align:right;">

\-0.93

</td>

<td style="text-align:right;">

0.24

</td>

<td style="text-align:left;">

3.43e-04

</td>

<td style="text-align:right;">

161

</td>

<td style="text-align:right;">

328

</td>

<td style="text-align:left;">

hsa-mir-411

</td>

<td style="text-align:left;">

MIR-154

</td>

<td style="text-align:left;">

AGUAGAC

</td>

<td style="text-align:left;">

chr14

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-193-P1b\_3p

</td>

<td style="text-align:right;">

\-0.92

</td>

<td style="text-align:right;">

0.26

</td>

<td style="text-align:left;">

8.53e-04

</td>

<td style="text-align:right;">

72

</td>

<td style="text-align:right;">

145

</td>

<td style="text-align:left;">

hsa-mir-193b

</td>

<td style="text-align:left;">

MIR-193

</td>

<td style="text-align:left;">

ACUGGCC

</td>

<td style="text-align:left;">

chr16

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-19-P2a\_3p/P2b\_3p

</td>

<td style="text-align:right;">

\-0.91

</td>

<td style="text-align:right;">

0.24

</td>

<td style="text-align:left;">

5.22e-04

</td>

<td style="text-align:right;">

645

</td>

<td style="text-align:right;">

1297

</td>

<td style="text-align:left;">

hsa-mir-19b-1

</td>

<td style="text-align:left;">

MIR-19

</td>

<td style="text-align:left;">

GUGCAAA

</td>

<td style="text-align:left;">

chr13

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-148-P1\_3p

</td>

<td style="text-align:right;">

\-0.86

</td>

<td style="text-align:right;">

0.27

</td>

<td style="text-align:left;">

3.44e-03

</td>

<td style="text-align:right;">

18874

</td>

<td style="text-align:right;">

34704

</td>

<td style="text-align:left;">

hsa-mir-148a

</td>

<td style="text-align:left;">

MIR-148

</td>

<td style="text-align:left;">

CAGUGCA

</td>

<td style="text-align:left;">

chr7

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-146-P1\_5p

</td>

<td style="text-align:right;">

\-0.85

</td>

<td style="text-align:right;">

0.39

</td>

<td style="text-align:left;">

4.76e-02

</td>

<td style="text-align:right;">

786

</td>

<td style="text-align:right;">

1453

</td>

<td style="text-align:left;">

hsa-mir-146a

</td>

<td style="text-align:left;">

MIR-146

</td>

<td style="text-align:left;">

GAGAACU

</td>

<td style="text-align:left;">

chr5

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-136\_3p

</td>

<td style="text-align:right;">

\-0.81

</td>

<td style="text-align:right;">

0.24

</td>

<td style="text-align:left;">

1.54e-03

</td>

<td style="text-align:right;">

88

</td>

<td style="text-align:right;">

173

</td>

<td style="text-align:left;">

hsa-mir-136

</td>

<td style="text-align:left;">

MIR-136

</td>

<td style="text-align:left;">

AUCAUCG

</td>

<td style="text-align:left;">

chr14

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-10-P1a\_5p

</td>

<td style="text-align:right;">

\-0.79

</td>

<td style="text-align:right;">

0.34

</td>

<td style="text-align:left;">

4.36e-02

</td>

<td style="text-align:right;">

61989

</td>

<td style="text-align:right;">

97123

</td>

<td style="text-align:left;">

hsa-mir-10a

</td>

<td style="text-align:left;">

MIR-10

</td>

<td style="text-align:left;">

ACCCUGU

</td>

<td style="text-align:left;">

chr17

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-154-P9\_3p

</td>

<td style="text-align:right;">

\-0.78

</td>

<td style="text-align:right;">

0.23

</td>

<td style="text-align:left;">

1.67e-03

</td>

<td style="text-align:right;">

135

</td>

<td style="text-align:right;">

253

</td>

<td style="text-align:left;">

hsa-mir-381

</td>

<td style="text-align:left;">

MIR-154

</td>

<td style="text-align:left;">

AUACAAG

</td>

<td style="text-align:left;">

chr14

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-214\_3p

</td>

<td style="text-align:right;">

\-0.74

</td>

<td style="text-align:right;">

0.26

</td>

<td style="text-align:left;">

8.05e-03

</td>

<td style="text-align:right;">

77

</td>

<td style="text-align:right;">

137

</td>

<td style="text-align:left;">

hsa-mir-214

</td>

<td style="text-align:left;">

MIR-214

</td>

<td style="text-align:left;">

CAGCAGG

</td>

<td style="text-align:left;">

chr1

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-1307\_3p

</td>

<td style="text-align:right;">

\-0.74

</td>

<td style="text-align:right;">

0.21

</td>

<td style="text-align:left;">

1.32e-03

</td>

<td style="text-align:right;">

80

</td>

<td style="text-align:right;">

139

</td>

<td style="text-align:left;">

hsa-mir-1307

</td>

<td style="text-align:left;">

MIR-1307

</td>

<td style="text-align:left;">

CGACCGG

</td>

<td style="text-align:left;">

chr10

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Let-7-P2b3\_5p

</td>

<td style="text-align:right;">

\-0.72

</td>

<td style="text-align:right;">

0.21

</td>

<td style="text-align:left;">

1.18e-03

</td>

<td style="text-align:right;">

1289

</td>

<td style="text-align:right;">

2149

</td>

<td style="text-align:left;">

hsa-mir-98

</td>

<td style="text-align:left;">

LET-7

</td>

<td style="text-align:left;">

GAGGUAG

</td>

<td style="text-align:left;">

chrX

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-199-P1\_3p/P2\_3p/P3\_3p

</td>

<td style="text-align:right;">

\-0.64

</td>

<td style="text-align:right;">

0.22

</td>

<td style="text-align:left;">

5.97e-03

</td>

<td style="text-align:right;">

3382

</td>

<td style="text-align:right;">

5772

</td>

<td style="text-align:left;">

hsa-mir-199b

</td>

<td style="text-align:left;">

MIR-199

</td>

<td style="text-align:left;">

CAGUAGU

</td>

<td style="text-align:left;">

chr9

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

na

</td>

<td style="text-align:left;">

na

</td>

</tr>

</tbody>

</table>

</div>

``` r
# Number of downregulated miRNA
signature_mirnas$number_downregulated
```

    ## [1] 52

``` r
SubtractLFC <- function(x, y){
  z = x
  if ( is.na(x) | is.na(y) ){ return( z ) }
  else if (sign(x) == sign(y)){ z = x - y }
  if (sign(z) != sign(x)) { z = 0 }
  return( z )
}

SubtractAdjP <- function(x , y, xP, yP){
  z = xP
  if ( is.na(xP) | is.na(yP) ){ return( z ) }
  if ( sign(x) == sign(y) ){ 
    z = (xP + ( 1 - yP )) }
  if (z > 1) {z = 1}
  return(z)
}
```

## pCRC vs mLi

``` r
#pCRC versus liver metastasis, control also with pCRC versus normal liver

column='tissue.type'
tissue_type_A <- 'metastasis.liver'
tissue_type_B <- 'tumor.colorect'
norm_adj_up       = dict_sig_mirna$tissue.type_normal.liver_vs_normal.colorect_up
norm_adj_down     = dict_sig_mirna$tissue.type_normal.liver_vs_normal.colorect_down
pCRC_adj_up   = dict_sig_mirna$tissue.type_normal.liver_vs_tumor.colorect_up
pCRC_adj_down = dict_sig_mirna$tissue.type_normal.liver_vs_tumor.colorect_down
palette <- 'jco'

coef <- paste(column, tissue_type_A, 'vs', tissue_type_B, sep='_')
res <- DeseqResult(dds, column, coef, tissue_type_A, tissue_type_B,
                   lfc.Threshold, rpm.Threshold,
                   norm_adj_up,
                   norm_adj_down,
                   pCRC_adj_up,
                   pCRC_adj_down)

dict_sig_mirna[paste(coef, "up",   sep='_')] <- list(res$up_mirna)
dict_sig_mirna[paste(coef, "down", sep='_')] <- list(res$down_mirna)
res_res <- res$res
res_dict[coef] <- res_res
plotMA(res$res, alpha=0.05)
```

![](comet_analysis_report_automated_no_lfcthreshold_08.12.20_res_alpha_0.05_with_LFC_shrinkage_LFC_0.58_mirge3_7.0_files/figure-gfm/unnamed-chunk-28-1.png)<!-- -->

``` r
# Plot volcano plot
VolcanoPlot(res$res, coef, res$sig,
            res$up_mirna, res$down_mirna,
            norm_adj_up, norm_adj_down,
            pCRC_adj_up, pCRC_adj_down)
```

![](comet_analysis_report_automated_no_lfcthreshold_08.12.20_res_alpha_0.05_with_LFC_shrinkage_LFC_0.58_mirge3_7.0_files/figure-gfm/unnamed-chunk-28-2.png)<!-- -->

``` r
ExpressionPlot(res$res, res$rpm, coef, res$sig,
               tissue_type_A, tissue_type_B,
               res$up_mirna, res$down_mirna,
               norm_adj_up, norm_adj_down,
               pCRC_adj_up, pCRC_adj_down)
```

![](comet_analysis_report_automated_no_lfcthreshold_08.12.20_res_alpha_0.05_with_LFC_shrinkage_LFC_0.58_mirge3_7.0_files/figure-gfm/unnamed-chunk-28-3.png)<!-- -->

``` r
signature_mirnas <- SigList(res, dds, tissue_type_A, tissue_type_B, coef,
                            norm_adj_up, norm_adj_down, 
                            pCRC_adj_up, pCRC_adj_down)
# Print list upregulated miRNA
signature_mirnas$up_mirna
```

<div style="border: 1px solid #ddd; padding: 5px; overflow-x: scroll; width:2000px; ">

<table class="table table-striped table-hover table-condensed" style="margin-left: auto; margin-right: auto;">

<caption>

Upregulated in tissue.type\_metastasis.liver\_vs\_tumor.colorect

</caption>

<thead>

<tr>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

miRNA

</th>

<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">

LFC

</th>

<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">

lfcSE

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

FDR

</th>

<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">

RPM metastasis.liver

</th>

<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">

RPM tumor.colorect

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

miRBase\_ID

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

Family

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

Seed

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

Chr

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

Cell-Type Specific

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

Norm Background

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

pCRC Background

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

Hsa-Mir-210\_3p

</td>

<td style="text-align:right;">

1.26

</td>

<td style="text-align:right;">

0.18

</td>

<td style="text-align:left;">

4.73e-11

</td>

<td style="text-align:right;">

685

</td>

<td style="text-align:right;">

291

</td>

<td style="text-align:left;">

hsa-mir-210

</td>

<td style="text-align:left;">

MIR-210

</td>

<td style="text-align:left;">

UGUGCGU

</td>

<td style="text-align:left;">

chr11

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-592\_5p

</td>

<td style="text-align:right;">

0.98

</td>

<td style="text-align:right;">

0.27

</td>

<td style="text-align:left;">

4.88e-03

</td>

<td style="text-align:right;">

137

</td>

<td style="text-align:right;">

75

</td>

<td style="text-align:left;">

hsa-mir-592

</td>

<td style="text-align:left;">

MIR-592

</td>

<td style="text-align:left;">

UGUGUCA

</td>

<td style="text-align:left;">

chr7

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-10-P1a\_5p

</td>

<td style="text-align:right;">

0.86

</td>

<td style="text-align:right;">

0.19

</td>

<td style="text-align:left;">

8.04e-05

</td>

<td style="text-align:right;">

164925

</td>

<td style="text-align:right;">

97123

</td>

<td style="text-align:left;">

hsa-mir-10a

</td>

<td style="text-align:left;">

MIR-10

</td>

<td style="text-align:left;">

ACCCUGU

</td>

<td style="text-align:left;">

chr17

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-1307\_5p

</td>

<td style="text-align:right;">

0.85

</td>

<td style="text-align:right;">

0.17

</td>

<td style="text-align:left;">

2.99e-06

</td>

<td style="text-align:right;">

887

</td>

<td style="text-align:right;">

492

</td>

<td style="text-align:left;">

hsa-mir-1307

</td>

<td style="text-align:left;">

MIR-1307

</td>

<td style="text-align:left;">

CGACCGG

</td>

<td style="text-align:left;">

chr10

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-1247\_5p

</td>

<td style="text-align:right;">

0.85

</td>

<td style="text-align:right;">

0.28

</td>

<td style="text-align:left;">

2.58e-02

</td>

<td style="text-align:right;">

127

</td>

<td style="text-align:right;">

69

</td>

<td style="text-align:left;">

hsa-mir-1247

</td>

<td style="text-align:left;">

MIR-1247

</td>

<td style="text-align:left;">

CCCGUCC

</td>

<td style="text-align:left;">

chr14

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-191\_5p

</td>

<td style="text-align:right;">

0.74

</td>

<td style="text-align:right;">

0.15

</td>

<td style="text-align:left;">

3.67e-06

</td>

<td style="text-align:right;">

23253

</td>

<td style="text-align:right;">

13630

</td>

<td style="text-align:left;">

hsa-mir-191

</td>

<td style="text-align:left;">

MIR-191

</td>

<td style="text-align:left;">

AACGGAA

</td>

<td style="text-align:left;">

chr3

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-425\_5p

</td>

<td style="text-align:right;">

0.73

</td>

<td style="text-align:right;">

0.12

</td>

<td style="text-align:left;">

2.90e-08

</td>

<td style="text-align:right;">

1117

</td>

<td style="text-align:right;">

700

</td>

<td style="text-align:left;">

hsa-mir-425

</td>

<td style="text-align:left;">

MIR-425

</td>

<td style="text-align:left;">

AUGACAC

</td>

<td style="text-align:left;">

chr3

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-8-P1b\_3p

</td>

<td style="text-align:right;">

0.64

</td>

<td style="text-align:right;">

0.19

</td>

<td style="text-align:left;">

4.07e-03

</td>

<td style="text-align:right;">

11957

</td>

<td style="text-align:right;">

7518

</td>

<td style="text-align:left;">

hsa-mir-141

</td>

<td style="text-align:left;">

MIR-8

</td>

<td style="text-align:left;">

AACACUG

</td>

<td style="text-align:left;">

chr12

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-150\_5p

</td>

<td style="text-align:right;">

1.08

</td>

<td style="text-align:right;">

0.22

</td>

<td style="text-align:left;">

2.99e-06

</td>

<td style="text-align:right;">

721

</td>

<td style="text-align:right;">

280

</td>

<td style="text-align:left;">

hsa-mir-150

</td>

<td style="text-align:left;">

MIR-150

</td>

<td style="text-align:left;">

CUCCCAA

</td>

<td style="text-align:left;">

chr19

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">Lymphocyte</span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-342\_3p

</td>

<td style="text-align:right;">

0.93

</td>

<td style="text-align:right;">

0.16

</td>

<td style="text-align:left;">

5.32e-08

</td>

<td style="text-align:right;">

320

</td>

<td style="text-align:right;">

145

</td>

<td style="text-align:left;">

hsa-mir-342

</td>

<td style="text-align:left;">

MIR-342

</td>

<td style="text-align:left;">

CUCACAC

</td>

<td style="text-align:left;">

chr14

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">c(“Dendritic
Cell”, “Lymphocyte”, “Macrophage”)</span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-331\_3p

</td>

<td style="text-align:right;">

0.75

</td>

<td style="text-align:right;">

0.16

</td>

<td style="text-align:left;">

3.16e-05

</td>

<td style="text-align:right;">

107

</td>

<td style="text-align:right;">

60

</td>

<td style="text-align:left;">

hsa-mir-331

</td>

<td style="text-align:left;">

MIR-331

</td>

<td style="text-align:left;">

CCCCUGG

</td>

<td style="text-align:left;">

chr12

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-15-P2a\_5p/P2b\_5p

</td>

<td style="text-align:right;">

0.68

</td>

<td style="text-align:right;">

0.11

</td>

<td style="text-align:left;">

1.34e-08

</td>

<td style="text-align:right;">

7588

</td>

<td style="text-align:right;">

4522

</td>

<td style="text-align:left;">

hsa-mir-16-1

</td>

<td style="text-align:left;">

MIR-15

</td>

<td style="text-align:left;">

AGCAGCA

</td>

<td style="text-align:left;">

chr13

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-204-P1\_5p

</td>

<td style="text-align:right;">

1.04

</td>

<td style="text-align:right;">

0.32

</td>

<td style="text-align:left;">

1.21e-06

</td>

<td style="text-align:right;">

260

</td>

<td style="text-align:right;">

62

</td>

<td style="text-align:left;">

hsa-mir-204

</td>

<td style="text-align:left;">

MIR-204

</td>

<td style="text-align:left;">

UCCCUUU

</td>

<td style="text-align:left;">

chr9

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">Retinal
Epithelial Cell</span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-335\_5p

</td>

<td style="text-align:right;">

0.78

</td>

<td style="text-align:right;">

0.11

</td>

<td style="text-align:left;">

5.33e-11

</td>

<td style="text-align:right;">

261

</td>

<td style="text-align:right;">

150

</td>

<td style="text-align:left;">

hsa-mir-335

</td>

<td style="text-align:left;">

MIR-335

</td>

<td style="text-align:left;">

CAAGAGC

</td>

<td style="text-align:left;">

chr7

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">Retinal
Epithelial Cell</span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-122\_5p

</td>

<td style="text-align:right;">

5.14

</td>

<td style="text-align:right;">

0.43

</td>

<td style="text-align:left;">

2.78e-29

</td>

<td style="text-align:right;">

3668

</td>

<td style="text-align:right;">

9

</td>

<td style="text-align:left;">

hsa-mir-122

</td>

<td style="text-align:left;">

MIR-122

</td>

<td style="text-align:left;">

GGAGUGU

</td>

<td style="text-align:left;">

chr18

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">Hepatocyte</span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-10-P3c\_5p

</td>

<td style="text-align:right;">

0.71

</td>

<td style="text-align:right;">

0.19

</td>

<td style="text-align:left;">

7.30e-04

</td>

<td style="text-align:right;">

480

</td>

<td style="text-align:right;">

239

</td>

<td style="text-align:left;">

hsa-mir-125b-2

</td>

<td style="text-align:left;">

MIR-10

</td>

<td style="text-align:left;">

CCCUGAG

</td>

<td style="text-align:left;">

chr21

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

</tr>

</tbody>

</table>

</div>

``` r
# Number of upregulated miRNA
signature_mirnas$number_upregulated
```

    ## [1] 16

``` r
# Print list downregulated miRNA
signature_mirnas$down_mirna
```

<div style="border: 1px solid #ddd; padding: 5px; overflow-x: scroll; width:2000px; ">

<table class="table table-striped table-hover table-condensed" style="margin-left: auto; margin-right: auto;">

<caption>

Downregulated in tissue.type\_metastasis.liver\_vs\_tumor.colorect

</caption>

<thead>

<tr>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

miRNA

</th>

<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">

LFC

</th>

<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">

lfcSE

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

FDR

</th>

<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">

RPM metastasis.liver

</th>

<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">

RPM tumor.colorect

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

miRBase\_ID

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

Family

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

Seed

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

Chr

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

Cell-Type Specific

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

Norm Background

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

pCRC Background

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

Hsa-Mir-486\_5p

</td>

<td style="text-align:right;">

\-0.66

</td>

<td style="text-align:right;">

0.23

</td>

<td style="text-align:left;">

3.47e-02

</td>

<td style="text-align:right;">

1489

</td>

<td style="text-align:right;">

2013

</td>

<td style="text-align:left;">

hsa-mir-486-1

</td>

<td style="text-align:left;">

MIR-486

</td>

<td style="text-align:left;">

CCUGUAC

</td>

<td style="text-align:left;">

chr8

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">c(“Platelet”,
“Red Blood Cell”)</span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-31\_5p

</td>

<td style="text-align:right;">

\-2.34

</td>

<td style="text-align:right;">

0.31

</td>

<td style="text-align:left;">

4.05e-13

</td>

<td style="text-align:right;">

88

</td>

<td style="text-align:right;">

627

</td>

<td style="text-align:left;">

hsa-mir-31

</td>

<td style="text-align:left;">

MIR-31

</td>

<td style="text-align:left;">

GGCAAGA

</td>

<td style="text-align:left;">

chr9

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Let-7-P2c2\_5p

</td>

<td style="text-align:right;">

\-0.90

</td>

<td style="text-align:right;">

0.11

</td>

<td style="text-align:left;">

1.69e-14

</td>

<td style="text-align:right;">

2840

</td>

<td style="text-align:right;">

5432

</td>

<td style="text-align:left;">

hsa-let-7i

</td>

<td style="text-align:left;">

LET-7

</td>

<td style="text-align:left;">

GAGGUAG

</td>

<td style="text-align:left;">

chr12

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-143\_3p

</td>

<td style="text-align:right;">

\-1.02

</td>

<td style="text-align:right;">

0.17

</td>

<td style="text-align:left;">

4.70e-08

</td>

<td style="text-align:right;">

72945

</td>

<td style="text-align:right;">

148431

</td>

<td style="text-align:left;">

hsa-mir-143

</td>

<td style="text-align:left;">

MIR-143

</td>

<td style="text-align:left;">

GAGAUGA

</td>

<td style="text-align:left;">

chr5

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">Mesenchymal</span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-133-P1\_3p/P2\_3p/P3\_3p

</td>

<td style="text-align:right;">

\-1.56

</td>

<td style="text-align:right;">

0.24

</td>

<td style="text-align:left;">

3.43e-10

</td>

<td style="text-align:right;">

38

</td>

<td style="text-align:right;">

116

</td>

<td style="text-align:left;">

hsa-mir-133a-2

</td>

<td style="text-align:left;">

MIR-133

</td>

<td style="text-align:left;">

UUGGUCC

</td>

<td style="text-align:left;">

chr20

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">c(“Skeletal
Myocyte”, “Stem Cell”)</span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-10-P1b\_5p

</td>

<td style="text-align:right;">

\-1.62

</td>

<td style="text-align:right;">

0.23

</td>

<td style="text-align:left;">

1.76e-10

</td>

<td style="text-align:right;">

11785

</td>

<td style="text-align:right;">

39286

</td>

<td style="text-align:left;">

hsa-mir-10b

</td>

<td style="text-align:left;">

MIR-10

</td>

<td style="text-align:left;">

ACCCUGU

</td>

<td style="text-align:left;">

chr2

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-92-P1c\_3p

</td>

<td style="text-align:right;">

\-0.68

</td>

<td style="text-align:right;">

0.17

</td>

<td style="text-align:left;">

2.83e-04

</td>

<td style="text-align:right;">

879

</td>

<td style="text-align:right;">

1398

</td>

<td style="text-align:left;">

hsa-mir-92b

</td>

<td style="text-align:left;">

MIR-92

</td>

<td style="text-align:left;">

AUUGCAC

</td>

<td style="text-align:left;">

chr1

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-146-P1\_5p

</td>

<td style="text-align:right;">

\-0.63

</td>

<td style="text-align:right;">

0.22

</td>

<td style="text-align:left;">

1.91e-02

</td>

<td style="text-align:right;">

909

</td>

<td style="text-align:right;">

1453

</td>

<td style="text-align:left;">

hsa-mir-146a

</td>

<td style="text-align:left;">

MIR-146

</td>

<td style="text-align:left;">

GAGAACU

</td>

<td style="text-align:left;">

chr5

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

</tr>

</tbody>

</table>

</div>

``` r
# Number of downregulated miRNA
signature_mirnas$number_downregulated
```

    ## [1] 8

``` r
res_tibble <- res$res
res_tibble$miRNA <- rownames(res_tibble)
res_tibble <- as_tibble(res_tibble)

metslfc <- res_dict$tissue.type_metastasis.liver_vs_tumor.colorect$log2FoldChange
normlfc <- res_dict$tissue.type_normal.liver_vs_normal.colorect$log2FoldChange

res_tibble$LFC_adj_background <- mapply(SubtractLFC, metslfc, normlfc)

metsP <- res_dict$tissue.type_metastasis.liver_vs_tumor.colorect$padj
normP <- res_dict$tissue.type_normal.liver_vs_normal.colorect$padj

res_tibble$padj_subt_normal <- mapply( SubtractAdjP, metslfc, normlfc, metsP, normP )

res_tibble %>% select(miRNA, log2FoldChange, lfcSE, LFC_adj_background, padj_subt_normal, baseMean, stat, pvalue, padj) %>% write_csv(path = '/Users/eirikhoy/Dropbox/projects/comet_analysis/data/Deseq_result_clm_vs_pcrc.csv')
```

## pCRC vs mLu

``` r
#pCRC versus lung metastasis, control also with pCRC versus normal liver

column='tissue.type'
tissue_type_A <- 'metastasis.lung'
tissue_type_B <- 'tumor.colorect'
norm_adj_up       = dict_sig_mirna$tissue.type_normal.lung_vs_normal.colorect_up
norm_adj_down     = dict_sig_mirna$tissue.type_normal.lung_vs_normal.colorect_down
pCRC_adj_up   = dict_sig_mirna$tissue.type_normal.lung_vs_tumor.colorect_up
pCRC_adj_down = dict_sig_mirna$tissue.type_normal.lung_vs_tumor.colorect_down

coef <- paste(column, tissue_type_A, 'vs', tissue_type_B, sep='_')
res <- DeseqResult(dds, column, coef, tissue_type_A, tissue_type_B,
                   lfc.Threshold, rpm.Threshold,
                   norm_adj_up,
                   norm_adj_down,
                   pCRC_adj_up,
                   pCRC_adj_down)
dict_sig_mirna[paste(coef, "up",   sep='_')] <- list(res$up_mirna)
dict_sig_mirna[paste(coef, "down", sep='_')] <- list(res$down_mirna)
res_res <- res$res
res_dict[coef] <- res_res
plotMA(res$res, alpha=0.05)
```

![](comet_analysis_report_automated_no_lfcthreshold_08.12.20_res_alpha_0.05_with_LFC_shrinkage_LFC_0.58_mirge3_7.0_files/figure-gfm/unnamed-chunk-29-1.png)<!-- -->

``` r
# Plot volcano plot
VolcanoPlot(res$res, coef, res$sig,
            res$up_mirna, res$down_mirna,
            norm_adj_up, norm_adj_down,
            pCRC_adj_up, pCRC_adj_down)
```

![](comet_analysis_report_automated_no_lfcthreshold_08.12.20_res_alpha_0.05_with_LFC_shrinkage_LFC_0.58_mirge3_7.0_files/figure-gfm/unnamed-chunk-29-2.png)<!-- -->

``` r
ExpressionPlot(res$res, res$rpm, coef, res$sig,
               tissue_type_A, tissue_type_B,
               res$up_mirna, res$down_mirna,
               norm_adj_up, norm_adj_down,
               pCRC_adj_up, pCRC_adj_down)
```

![](comet_analysis_report_automated_no_lfcthreshold_08.12.20_res_alpha_0.05_with_LFC_shrinkage_LFC_0.58_mirge3_7.0_files/figure-gfm/unnamed-chunk-29-3.png)<!-- -->

``` r
signature_mirnas <- SigList(res, dds, tissue_type_A, tissue_type_B, coef,
                            norm_adj_up, norm_adj_down, 
                            pCRC_adj_up, pCRC_adj_down)
# Print list upregulated miRNA
signature_mirnas$up_mirna
```

<div style="border: 1px solid #ddd; padding: 5px; overflow-x: scroll; width:2000px; ">

<table class="table table-striped table-hover table-condensed" style="margin-left: auto; margin-right: auto;">

<caption>

Upregulated in tissue.type\_metastasis.lung\_vs\_tumor.colorect

</caption>

<thead>

<tr>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

miRNA

</th>

<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">

LFC

</th>

<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">

lfcSE

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

FDR

</th>

<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">

RPM metastasis.lung

</th>

<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">

RPM tumor.colorect

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

miRBase\_ID

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

Family

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

Seed

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

Chr

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

Cell-Type Specific

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

Norm Background

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

pCRC Background

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

Hsa-Mir-155\_5p

</td>

<td style="text-align:right;">

0.76

</td>

<td style="text-align:right;">

0.18

</td>

<td style="text-align:left;">

7.24e-05

</td>

<td style="text-align:right;">

733

</td>

<td style="text-align:right;">

387

</td>

<td style="text-align:left;">

hsa-mir-155

</td>

<td style="text-align:left;">

MIR-155

</td>

<td style="text-align:left;">

UAAUGCU

</td>

<td style="text-align:left;">

chr21

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">c(“Lymphocyte”,
“Macrophage”)</span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-210\_3p

</td>

<td style="text-align:right;">

1.15

</td>

<td style="text-align:right;">

0.19

</td>

<td style="text-align:left;">

1.55e-08

</td>

<td style="text-align:right;">

693

</td>

<td style="text-align:right;">

291

</td>

<td style="text-align:left;">

hsa-mir-210

</td>

<td style="text-align:left;">

MIR-210

</td>

<td style="text-align:left;">

UGUGCGU

</td>

<td style="text-align:left;">

chr11

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-142\_5p

</td>

<td style="text-align:right;">

0.84

</td>

<td style="text-align:right;">

0.18

</td>

<td style="text-align:left;">

2.01e-05

</td>

<td style="text-align:right;">

7117

</td>

<td style="text-align:right;">

3453

</td>

<td style="text-align:left;">

hsa-mir-142

</td>

<td style="text-align:left;">

MIR-142

</td>

<td style="text-align:left;">

AUAAAGU

</td>

<td style="text-align:left;">

chr17

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-19-P2a\_3p/P2b\_3p

</td>

<td style="text-align:right;">

0.80

</td>

<td style="text-align:right;">

0.15

</td>

<td style="text-align:left;">

9.38e-07

</td>

<td style="text-align:right;">

2345

</td>

<td style="text-align:right;">

1297

</td>

<td style="text-align:left;">

hsa-mir-19b-1

</td>

<td style="text-align:left;">

MIR-19

</td>

<td style="text-align:left;">

GUGCAAA

</td>

<td style="text-align:left;">

chr13

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-8-P1b\_3p

</td>

<td style="text-align:right;">

0.75

</td>

<td style="text-align:right;">

0.21

</td>

<td style="text-align:left;">

1.12e-03

</td>

<td style="text-align:right;">

15111

</td>

<td style="text-align:right;">

7518

</td>

<td style="text-align:left;">

hsa-mir-141

</td>

<td style="text-align:left;">

MIR-8

</td>

<td style="text-align:left;">

AACACUG

</td>

<td style="text-align:left;">

chr12

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-191\_5p

</td>

<td style="text-align:right;">

0.74

</td>

<td style="text-align:right;">

0.16

</td>

<td style="text-align:left;">

9.81e-06

</td>

<td style="text-align:right;">

26989

</td>

<td style="text-align:right;">

13630

</td>

<td style="text-align:left;">

hsa-mir-191

</td>

<td style="text-align:left;">

MIR-191

</td>

<td style="text-align:left;">

AACGGAA

</td>

<td style="text-align:left;">

chr3

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-374-P1\_5p

</td>

<td style="text-align:right;">

0.72

</td>

<td style="text-align:right;">

0.12

</td>

<td style="text-align:left;">

9.07e-09

</td>

<td style="text-align:right;">

217

</td>

<td style="text-align:right;">

115

</td>

<td style="text-align:left;">

hsa-mir-374a

</td>

<td style="text-align:left;">

MIR-374

</td>

<td style="text-align:left;">

UAUAAUA

</td>

<td style="text-align:left;">

chrX

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-19-P1\_3p

</td>

<td style="text-align:right;">

0.65

</td>

<td style="text-align:right;">

0.16

</td>

<td style="text-align:left;">

3.42e-04

</td>

<td style="text-align:right;">

597

</td>

<td style="text-align:right;">

379

</td>

<td style="text-align:left;">

hsa-mir-19a

</td>

<td style="text-align:left;">

MIR-19

</td>

<td style="text-align:left;">

GUGCAAA

</td>

<td style="text-align:left;">

chr13

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-145\_5p

</td>

<td style="text-align:right;">

0.67

</td>

<td style="text-align:right;">

0.22

</td>

<td style="text-align:left;">

6.74e-03

</td>

<td style="text-align:right;">

1629

</td>

<td style="text-align:right;">

840

</td>

<td style="text-align:left;">

hsa-mir-145

</td>

<td style="text-align:left;">

MIR-145

</td>

<td style="text-align:left;">

UCCAGUU

</td>

<td style="text-align:left;">

chr5

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">Mesenchymal</span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-150\_5p

</td>

<td style="text-align:right;">

1.14

</td>

<td style="text-align:right;">

0.24

</td>

<td style="text-align:left;">

2.51e-06

</td>

<td style="text-align:right;">

738

</td>

<td style="text-align:right;">

280

</td>

<td style="text-align:left;">

hsa-mir-150

</td>

<td style="text-align:left;">

MIR-150

</td>

<td style="text-align:left;">

CUCCCAA

</td>

<td style="text-align:left;">

chr19

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">Lymphocyte</span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-10-P3c\_5p

</td>

<td style="text-align:right;">

0.85

</td>

<td style="text-align:right;">

0.21

</td>

<td style="text-align:left;">

9.78e-05

</td>

<td style="text-align:right;">

574

</td>

<td style="text-align:right;">

239

</td>

<td style="text-align:left;">

hsa-mir-125b-2

</td>

<td style="text-align:left;">

MIR-10

</td>

<td style="text-align:left;">

CCCUGAG

</td>

<td style="text-align:left;">

chr21

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-29-P1b\_3p

</td>

<td style="text-align:right;">

0.71

</td>

<td style="text-align:right;">

0.14

</td>

<td style="text-align:left;">

3.39e-06

</td>

<td style="text-align:right;">

597

</td>

<td style="text-align:right;">

331

</td>

<td style="text-align:left;">

hsa-mir-29c

</td>

<td style="text-align:left;">

MIR-29

</td>

<td style="text-align:left;">

AGCACCA

</td>

<td style="text-align:left;">

chr1

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-15-P2a\_5p/P2b\_5p

</td>

<td style="text-align:right;">

0.71

</td>

<td style="text-align:right;">

0.12

</td>

<td style="text-align:left;">

1.51e-08

</td>

<td style="text-align:right;">

8258

</td>

<td style="text-align:right;">

4522

</td>

<td style="text-align:left;">

hsa-mir-16-1

</td>

<td style="text-align:left;">

MIR-15

</td>

<td style="text-align:left;">

AGCAGCA

</td>

<td style="text-align:left;">

chr13

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-26-P3\_5p

</td>

<td style="text-align:right;">

0.65

</td>

<td style="text-align:right;">

0.14

</td>

<td style="text-align:left;">

6.79e-06

</td>

<td style="text-align:right;">

6205

</td>

<td style="text-align:right;">

3467

</td>

<td style="text-align:left;">

hsa-mir-26a-1

</td>

<td style="text-align:left;">

MIR-26

</td>

<td style="text-align:left;">

UCAAGUA

</td>

<td style="text-align:left;">

chr3

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-331\_3p

</td>

<td style="text-align:right;">

0.62

</td>

<td style="text-align:right;">

0.18

</td>

<td style="text-align:left;">

1.05e-03

</td>

<td style="text-align:right;">

105

</td>

<td style="text-align:right;">

60

</td>

<td style="text-align:left;">

hsa-mir-331

</td>

<td style="text-align:left;">

MIR-331

</td>

<td style="text-align:left;">

CCCCUGG

</td>

<td style="text-align:left;">

chr12

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-335\_5p

</td>

<td style="text-align:right;">

0.88

</td>

<td style="text-align:right;">

0.12

</td>

<td style="text-align:left;">

1.45e-12

</td>

<td style="text-align:right;">

295

</td>

<td style="text-align:right;">

150

</td>

<td style="text-align:left;">

hsa-mir-335

</td>

<td style="text-align:left;">

MIR-335

</td>

<td style="text-align:left;">

CAAGAGC

</td>

<td style="text-align:left;">

chr7

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">Retinal
Epithelial Cell</span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-451\_5p

</td>

<td style="text-align:right;">

0.87

</td>

<td style="text-align:right;">

0.26

</td>

<td style="text-align:left;">

1.88e-03

</td>

<td style="text-align:right;">

3152

</td>

<td style="text-align:right;">

1565

</td>

<td style="text-align:left;">

hsa-mir-451a

</td>

<td style="text-align:left;">

MIR-451

</td>

<td style="text-align:left;">

AACCGUU

</td>

<td style="text-align:left;">

chr17

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">Red
Blood Cell</span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-24-P1\_3p/P2\_3p

</td>

<td style="text-align:right;">

0.74

</td>

<td style="text-align:right;">

0.14

</td>

<td style="text-align:left;">

7.89e-07

</td>

<td style="text-align:right;">

1289

</td>

<td style="text-align:right;">

626

</td>

<td style="text-align:left;">

hsa-mir-24-2

</td>

<td style="text-align:left;">

MIR-24

</td>

<td style="text-align:left;">

GGCUCAG

</td>

<td style="text-align:left;">

chr19

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">Macrophage</span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-126\_5p

</td>

<td style="text-align:right;">

1.67

</td>

<td style="text-align:right;">

0.15

</td>

<td style="text-align:left;">

6.56e-29

</td>

<td style="text-align:right;">

13328

</td>

<td style="text-align:right;">

3455

</td>

<td style="text-align:left;">

hsa-mir-126

</td>

<td style="text-align:left;">

MIR-126

</td>

<td style="text-align:left;">

AUUAUUA

</td>

<td style="text-align:left;">

chr9

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">c(“Endothelial
Cell”, “Platelet”)</span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-146-P2\_5p

</td>

<td style="text-align:right;">

0.94

</td>

<td style="text-align:right;">

0.23

</td>

<td style="text-align:left;">

1.47e-04

</td>

<td style="text-align:right;">

15403

</td>

<td style="text-align:right;">

6279

</td>

<td style="text-align:left;">

hsa-mir-146b

</td>

<td style="text-align:left;">

MIR-146

</td>

<td style="text-align:left;">

GAGAACU

</td>

<td style="text-align:left;">

chr10

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">c(“Dendritic
Cell”, “Lymphocyte”)</span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-342\_3p

</td>

<td style="text-align:right;">

2.25

</td>

<td style="text-align:right;">

0.18

</td>

<td style="text-align:left;">

6.24e-36

</td>

<td style="text-align:right;">

872

</td>

<td style="text-align:right;">

145

</td>

<td style="text-align:left;">

hsa-mir-342

</td>

<td style="text-align:left;">

MIR-342

</td>

<td style="text-align:left;">

CUCACAC

</td>

<td style="text-align:left;">

chr14

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">c(“Dendritic
Cell”, “Lymphocyte”, “Macrophage”)</span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-34-P2b\_5p

</td>

<td style="text-align:right;">

4.24

</td>

<td style="text-align:right;">

0.24

</td>

<td style="text-align:left;">

2.76e-67

</td>

<td style="text-align:right;">

989

</td>

<td style="text-align:right;">

33

</td>

<td style="text-align:left;">

hsa-mir-34c

</td>

<td style="text-align:left;">

MIR-34

</td>

<td style="text-align:left;">

GGCAGUG

</td>

<td style="text-align:left;">

chr11

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-34-P2a\_5p

</td>

<td style="text-align:right;">

3.97

</td>

<td style="text-align:right;">

0.27

</td>

<td style="text-align:left;">

8.59e-45

</td>

<td style="text-align:right;">

128

</td>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

hsa-mir-34b

</td>

<td style="text-align:left;">

MIR-34

</td>

<td style="text-align:left;">

GGCAGUG

</td>

<td style="text-align:left;">

chr11

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-30-P1a\_5p

</td>

<td style="text-align:right;">

1.74

</td>

<td style="text-align:right;">

0.16

</td>

<td style="text-align:left;">

1.56e-26

</td>

<td style="text-align:right;">

13224

</td>

<td style="text-align:right;">

3165

</td>

<td style="text-align:left;">

hsa-mir-30a

</td>

<td style="text-align:left;">

MIR-30

</td>

<td style="text-align:left;">

GUAAACA

</td>

<td style="text-align:left;">

chr6

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-15-P2c\_5p

</td>

<td style="text-align:right;">

1.64

</td>

<td style="text-align:right;">

0.17

</td>

<td style="text-align:left;">

6.99e-22

</td>

<td style="text-align:right;">

975

</td>

<td style="text-align:right;">

257

</td>

<td style="text-align:left;">

hsa-mir-195

</td>

<td style="text-align:left;">

MIR-15

</td>

<td style="text-align:left;">

AGCAGCA

</td>

<td style="text-align:left;">

chr17

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-218-P1\_5p/P2\_5p

</td>

<td style="text-align:right;">

1.57

</td>

<td style="text-align:right;">

0.19

</td>

<td style="text-align:left;">

3.94e-16

</td>

<td style="text-align:right;">

135

</td>

<td style="text-align:right;">

37

</td>

<td style="text-align:left;">

hsa-mir-218-1

</td>

<td style="text-align:left;">

MIR-218

</td>

<td style="text-align:left;">

UGUGCUU

</td>

<td style="text-align:left;">

chr4

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-374-P2\_5p

</td>

<td style="text-align:right;">

1.14

</td>

<td style="text-align:right;">

0.14

</td>

<td style="text-align:left;">

1.78e-15

</td>

<td style="text-align:right;">

226

</td>

<td style="text-align:right;">

95

</td>

<td style="text-align:left;">

hsa-mir-374b

</td>

<td style="text-align:left;">

MIR-374

</td>

<td style="text-align:left;">

UAUAAUA

</td>

<td style="text-align:left;">

chrX

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-338-P1\_3p

</td>

<td style="text-align:right;">

1.10

</td>

<td style="text-align:right;">

0.22

</td>

<td style="text-align:left;">

1.51e-06

</td>

<td style="text-align:right;">

236

</td>

<td style="text-align:right;">

98

</td>

<td style="text-align:left;">

hsa-mir-338

</td>

<td style="text-align:left;">

MIR-338

</td>

<td style="text-align:left;">

CCAGCAU

</td>

<td style="text-align:left;">

chr17

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-26-P1\_5p/P2\_5p

</td>

<td style="text-align:right;">

1.07

</td>

<td style="text-align:right;">

0.10

</td>

<td style="text-align:left;">

7.47e-28

</td>

<td style="text-align:right;">

84008

</td>

<td style="text-align:right;">

35584

</td>

<td style="text-align:left;">

hsa-mir-26b

</td>

<td style="text-align:left;">

MIR-26

</td>

<td style="text-align:left;">

UCAAGUA

</td>

<td style="text-align:left;">

chr2

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-10-P2c\_5p

</td>

<td style="text-align:right;">

0.63

</td>

<td style="text-align:right;">

0.23

</td>

<td style="text-align:left;">

8.68e-03

</td>

<td style="text-align:right;">

289

</td>

<td style="text-align:right;">

142

</td>

<td style="text-align:left;">

hsa-mir-99a

</td>

<td style="text-align:left;">

MIR-10

</td>

<td style="text-align:left;">

ACCCGUA

</td>

<td style="text-align:left;">

chr21

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-130-P1a\_3p

</td>

<td style="text-align:right;">

0.62

</td>

<td style="text-align:right;">

0.14

</td>

<td style="text-align:left;">

4.12e-05

</td>

<td style="text-align:right;">

729

</td>

<td style="text-align:right;">

421

</td>

<td style="text-align:left;">

hsa-mir-130a

</td>

<td style="text-align:left;">

MIR-130

</td>

<td style="text-align:left;">

AGUGCAA

</td>

<td style="text-align:left;">

chr11

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

</tr>

</tbody>

</table>

</div>

``` r
# Number of upregulated miRNA
signature_mirnas$number_upregulated
```

    ## [1] 31

``` r
# Print list downregulated miRNA
signature_mirnas$down_mirna
```

<div style="border: 1px solid #ddd; padding: 5px; overflow-x: scroll; width:2000px; ">

<table class="table table-striped table-hover table-condensed" style="margin-left: auto; margin-right: auto;">

<caption>

Downregulated in tissue.type\_metastasis.lung\_vs\_tumor.colorect

</caption>

<thead>

<tr>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

miRNA

</th>

<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">

LFC

</th>

<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">

lfcSE

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

FDR

</th>

<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">

RPM metastasis.lung

</th>

<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">

RPM tumor.colorect

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

miRBase\_ID

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

Family

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

Seed

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

Chr

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

Cell-Type Specific

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

Norm Background

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

pCRC Background

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

Hsa-Mir-423\_5p

</td>

<td style="text-align:right;">

\-1.17

</td>

<td style="text-align:right;">

0.17

</td>

<td style="text-align:left;">

1.35e-10

</td>

<td style="text-align:right;">

202

</td>

<td style="text-align:right;">

433

</td>

<td style="text-align:left;">

hsa-mir-423

</td>

<td style="text-align:left;">

MIR-423

</td>

<td style="text-align:left;">

GAGGGGC

</td>

<td style="text-align:left;">

chr17

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Let-7-P1b\_5p

</td>

<td style="text-align:right;">

\-0.86

</td>

<td style="text-align:right;">

0.16

</td>

<td style="text-align:left;">

2.55e-07

</td>

<td style="text-align:right;">

650

</td>

<td style="text-align:right;">

1057

</td>

<td style="text-align:left;">

hsa-let-7e

</td>

<td style="text-align:left;">

LET-7

</td>

<td style="text-align:left;">

GAGGUAG

</td>

<td style="text-align:left;">

chr19

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-197\_3p

</td>

<td style="text-align:right;">

\-0.74

</td>

<td style="text-align:right;">

0.13

</td>

<td style="text-align:left;">

1.76e-07

</td>

<td style="text-align:right;">

113

</td>

<td style="text-align:right;">

181

</td>

<td style="text-align:left;">

hsa-mir-197

</td>

<td style="text-align:left;">

MIR-197

</td>

<td style="text-align:left;">

UCACCAC

</td>

<td style="text-align:left;">

chr1

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-92-P1c\_3p

</td>

<td style="text-align:right;">

\-0.72

</td>

<td style="text-align:right;">

0.18

</td>

<td style="text-align:left;">

2.22e-04

</td>

<td style="text-align:right;">

981

</td>

<td style="text-align:right;">

1398

</td>

<td style="text-align:left;">

hsa-mir-92b

</td>

<td style="text-align:left;">

MIR-92

</td>

<td style="text-align:left;">

AUUGCAC

</td>

<td style="text-align:left;">

chr1

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-362-P2\_3p/P4\_3p

</td>

<td style="text-align:right;">

\-0.65

</td>

<td style="text-align:right;">

0.14

</td>

<td style="text-align:left;">

2.36e-05

</td>

<td style="text-align:right;">

302

</td>

<td style="text-align:right;">

441

</td>

<td style="text-align:left;">

hsa-mir-500a

</td>

<td style="text-align:left;">

MIR-362

</td>

<td style="text-align:left;">

UGCACCU

</td>

<td style="text-align:left;">

chrX

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-221-P2\_3p

</td>

<td style="text-align:right;">

\-0.60

</td>

<td style="text-align:right;">

0.14

</td>

<td style="text-align:left;">

8.50e-05

</td>

<td style="text-align:right;">

1101

</td>

<td style="text-align:right;">

1633

</td>

<td style="text-align:left;">

hsa-mir-222

</td>

<td style="text-align:left;">

MIR-221

</td>

<td style="text-align:left;">

GCUACAU

</td>

<td style="text-align:left;">

chrX

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-7-P1\_5p/P2\_5p/P3\_5p

</td>

<td style="text-align:right;">

\-0.80

</td>

<td style="text-align:right;">

0.25

</td>

<td style="text-align:left;">

1.88e-03

</td>

<td style="text-align:right;">

106

</td>

<td style="text-align:right;">

199

</td>

<td style="text-align:left;">

hsa-mir-7-1

</td>

<td style="text-align:left;">

MIR-7

</td>

<td style="text-align:left;">

GGAAGAC

</td>

<td style="text-align:left;">

chr9

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">c(“Islet
Cell”, “Neural”)</span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-154-P36\_3p

</td>

<td style="text-align:right;">

\-1.40

</td>

<td style="text-align:right;">

0.14

</td>

<td style="text-align:left;">

3.20e-21

</td>

<td style="text-align:right;">

80

</td>

<td style="text-align:right;">

195

</td>

<td style="text-align:left;">

hsa-mir-409

</td>

<td style="text-align:left;">

MIR-154

</td>

<td style="text-align:left;">

AAUGUUG

</td>

<td style="text-align:left;">

chr14

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-154-P12\_3p

</td>

<td style="text-align:right;">

\-1.38

</td>

<td style="text-align:right;">

0.15

</td>

<td style="text-align:left;">

7.64e-18

</td>

<td style="text-align:right;">

46

</td>

<td style="text-align:right;">

115

</td>

<td style="text-align:left;">

hsa-mir-410

</td>

<td style="text-align:left;">

MIR-154

</td>

<td style="text-align:left;">

AUAUAAC

</td>

<td style="text-align:left;">

chr14

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-31\_5p

</td>

<td style="text-align:right;">

\-1.01

</td>

<td style="text-align:right;">

0.34

</td>

<td style="text-align:left;">

1.99e-03

</td>

<td style="text-align:right;">

242

</td>

<td style="text-align:right;">

627

</td>

<td style="text-align:left;">

hsa-mir-31

</td>

<td style="text-align:left;">

MIR-31

</td>

<td style="text-align:left;">

GGCAAGA

</td>

<td style="text-align:left;">

chr9

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-96-P2\_5p

</td>

<td style="text-align:right;">

\-0.78

</td>

<td style="text-align:right;">

0.17

</td>

<td style="text-align:left;">

1.58e-05

</td>

<td style="text-align:right;">

6399

</td>

<td style="text-align:right;">

10417

</td>

<td style="text-align:left;">

hsa-mir-182

</td>

<td style="text-align:left;">

MIR-96

</td>

<td style="text-align:left;">

UUGGCAA

</td>

<td style="text-align:left;">

chr7

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-199-P1\_3p/P2\_3p/P3\_3p

</td>

<td style="text-align:right;">

\-0.65

</td>

<td style="text-align:right;">

0.14

</td>

<td style="text-align:left;">

1.06e-05

</td>

<td style="text-align:right;">

3904

</td>

<td style="text-align:right;">

5772

</td>

<td style="text-align:left;">

hsa-mir-199b

</td>

<td style="text-align:left;">

MIR-199

</td>

<td style="text-align:left;">

CAGUAGU

</td>

<td style="text-align:left;">

chr9

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-15-P1d\_5p

</td>

<td style="text-align:right;">

\-0.63

</td>

<td style="text-align:right;">

0.21

</td>

<td style="text-align:left;">

5.51e-03

</td>

<td style="text-align:right;">

158

</td>

<td style="text-align:right;">

235

</td>

<td style="text-align:left;">

hsa-mir-424

</td>

<td style="text-align:left;">

MIR-15

</td>

<td style="text-align:left;">

AGCAGCA

</td>

<td style="text-align:left;">

chrX

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-92-P1a\_3p/P1b\_3p

</td>

<td style="text-align:right;">

\-0.63

</td>

<td style="text-align:right;">

0.15

</td>

<td style="text-align:left;">

1.76e-04

</td>

<td style="text-align:right;">

26467

</td>

<td style="text-align:right;">

39320

</td>

<td style="text-align:left;">

hsa-mir-92a-1

</td>

<td style="text-align:left;">

MIR-92

</td>

<td style="text-align:left;">

AUUGCAC

</td>

<td style="text-align:left;">

chr13

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-143\_3p

</td>

<td style="text-align:right;">

\-1.34

</td>

<td style="text-align:right;">

0.19

</td>

<td style="text-align:left;">

1.00e-11

</td>

<td style="text-align:right;">

67942

</td>

<td style="text-align:right;">

148431

</td>

<td style="text-align:left;">

hsa-mir-143

</td>

<td style="text-align:left;">

MIR-143

</td>

<td style="text-align:left;">

GAGAUGA

</td>

<td style="text-align:left;">

chr5

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">Mesenchymal</span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-133-P1\_3p/P2\_3p/P3\_3p

</td>

<td style="text-align:right;">

\-1.73

</td>

<td style="text-align:right;">

0.26

</td>

<td style="text-align:left;">

4.97e-11

</td>

<td style="text-align:right;">

37

</td>

<td style="text-align:right;">

116

</td>

<td style="text-align:left;">

hsa-mir-133a-2

</td>

<td style="text-align:left;">

MIR-133

</td>

<td style="text-align:left;">

UUGGUCC

</td>

<td style="text-align:left;">

chr20

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">c(“Skeletal
Myocyte”, “Stem Cell”)</span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-362-P3\_3p

</td>

<td style="text-align:right;">

\-0.89

</td>

<td style="text-align:right;">

0.25

</td>

<td style="text-align:left;">

1.47e-03

</td>

<td style="text-align:right;">

61

</td>

<td style="text-align:right;">

108

</td>

<td style="text-align:left;">

hsa-mir-501

</td>

<td style="text-align:left;">

MIR-362

</td>

<td style="text-align:left;">

AUGCACC

</td>

<td style="text-align:left;">

chrX

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-8-P2a\_3p

</td>

<td style="text-align:right;">

\-0.60

</td>

<td style="text-align:right;">

0.16

</td>

<td style="text-align:left;">

6.49e-04

</td>

<td style="text-align:right;">

6834

</td>

<td style="text-align:right;">

10330

</td>

<td style="text-align:left;">

hsa-mir-200b

</td>

<td style="text-align:left;">

MIR-8

</td>

<td style="text-align:left;">

AAUACUG

</td>

<td style="text-align:left;">

chr1

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">Epithelial
Cell</span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-127\_3p

</td>

<td style="text-align:right;">

\-1.99

</td>

<td style="text-align:right;">

0.22

</td>

<td style="text-align:left;">

3.22e-17

</td>

<td style="text-align:right;">

617

</td>

<td style="text-align:right;">

2263

</td>

<td style="text-align:left;">

hsa-mir-127

</td>

<td style="text-align:left;">

MIR-127

</td>

<td style="text-align:left;">

CGGAUCC

</td>

<td style="text-align:left;">

chr14

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-154-P13\_5p

</td>

<td style="text-align:right;">

\-1.34

</td>

<td style="text-align:right;">

0.15

</td>

<td style="text-align:left;">

3.38e-17

</td>

<td style="text-align:right;">

142

</td>

<td style="text-align:right;">

328

</td>

<td style="text-align:left;">

hsa-mir-411

</td>

<td style="text-align:left;">

MIR-154

</td>

<td style="text-align:left;">

AGUAGAC

</td>

<td style="text-align:left;">

chr14

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-154-P9\_3p

</td>

<td style="text-align:right;">

\-1.34

</td>

<td style="text-align:right;">

0.14

</td>

<td style="text-align:left;">

1.80e-18

</td>

<td style="text-align:right;">

107

</td>

<td style="text-align:right;">

253

</td>

<td style="text-align:left;">

hsa-mir-381

</td>

<td style="text-align:left;">

MIR-154

</td>

<td style="text-align:left;">

AUACAAG

</td>

<td style="text-align:left;">

chr14

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-190-P1\_5p

</td>

<td style="text-align:right;">

\-1.26

</td>

<td style="text-align:right;">

0.18

</td>

<td style="text-align:left;">

2.77e-11

</td>

<td style="text-align:right;">

58

</td>

<td style="text-align:right;">

140

</td>

<td style="text-align:left;">

hsa-mir-190a

</td>

<td style="text-align:left;">

MIR-190

</td>

<td style="text-align:left;">

GAUAUGU

</td>

<td style="text-align:left;">

chr15

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-154-P23\_3p

</td>

<td style="text-align:right;">

\-1.22

</td>

<td style="text-align:right;">

0.14

</td>

<td style="text-align:left;">

1.46e-16

</td>

<td style="text-align:right;">

69

</td>

<td style="text-align:right;">

148

</td>

<td style="text-align:left;">

hsa-mir-654

</td>

<td style="text-align:left;">

MIR-154

</td>

<td style="text-align:left;">

AUGUCUG

</td>

<td style="text-align:left;">

chr14

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-136\_3p

</td>

<td style="text-align:right;">

\-1.09

</td>

<td style="text-align:right;">

0.15

</td>

<td style="text-align:left;">

5.02e-12

</td>

<td style="text-align:right;">

86

</td>

<td style="text-align:right;">

173

</td>

<td style="text-align:left;">

hsa-mir-136

</td>

<td style="text-align:left;">

MIR-136

</td>

<td style="text-align:left;">

AUCAUCG

</td>

<td style="text-align:left;">

chr14

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-378\_3p

</td>

<td style="text-align:right;">

\-0.69

</td>

<td style="text-align:right;">

0.14

</td>

<td style="text-align:left;">

2.25e-06

</td>

<td style="text-align:right;">

3998

</td>

<td style="text-align:right;">

6426

</td>

<td style="text-align:left;">

hsa-mir-378a

</td>

<td style="text-align:left;">

MIR-378

</td>

<td style="text-align:left;">

CUGGACU

</td>

<td style="text-align:left;">

chr5

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-10-P1b\_5p

</td>

<td style="text-align:right;">

\-0.66

</td>

<td style="text-align:right;">

0.25

</td>

<td style="text-align:left;">

1.93e-02

</td>

<td style="text-align:right;">

30285

</td>

<td style="text-align:right;">

39286

</td>

<td style="text-align:left;">

hsa-mir-10b

</td>

<td style="text-align:left;">

MIR-10

</td>

<td style="text-align:left;">

ACCCUGU

</td>

<td style="text-align:left;">

chr2

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

</tr>

</tbody>

</table>

</div>

``` r
# Number of downregulated miRNA
signature_mirnas$number_downregulated
```

    ## [1] 26

``` r
res_tibble <- res$res
res_tibble$miRNA <- rownames(res_tibble)
res_tibble <- as_tibble(res_tibble)

metslfc <- res_dict$tissue.type_metastasis.lung_vs_tumor.colorect$log2FoldChange
normlfc <- res_dict$tissue.type_normal.lung_vs_normal.colorect$log2FoldChange

res_tibble$LFC_adj_background <- mapply(SubtractLFC, metslfc, normlfc)

metsP <- res_dict$tissue.type_metastasis.lung_vs_tumor.colorect$padj
normP <- res_dict$tissue.type_normal.lung_vs_normal.colorect$padj

res_tibble$padj_subt_normal <- mapply(SubtractAdjP, metslfc, normlfc,metsP, normP)

res_tibble %>% select(miRNA, log2FoldChange, lfcSE, LFC_adj_background, padj_subt_normal, baseMean, stat, pvalue, padj) %>% write_csv(path = 
                      '/Users/eirikhoy/Dropbox/projects/comet_analysis/data/Deseq_result_mlu_vs_pcrc.csv')
```

## pCRC vs PM

``` r
#pCRC versus PC metastasis, union of mLi and mLu normal control

column='tissue.type'
tissue_type_A <- 'metastasis.pc'
tissue_type_B <- 'tumor.colorect'
norm_adj_up       = union(dict_sig_mirna$tissue.type_normal.liver_vs_normal.colorect_up, dict_sig_mirna$tissue.type_normal.lung_vs_normal.colorect_up)
norm_adj_down     = union(dict_sig_mirna$tissue.type_normal.liver_vs_normal.colorect_down, dict_sig_mirna$tissue.type_normal.lung_vs_normal.colorect_down)
pCRC_adj_up       = union(dict_sig_mirna$tissue.type_normal.liver_vs_tumor.colorect_up, dict_sig_mirna$tissue.type_normal.lung_vs_tumor.colorect_up)
pCRC_adj_down     = union(dict_sig_mirna$tissue.type_normal.liver_vs_tumor.colorect_down, dict_sig_mirna$tissue.type_normal.lung_vs_tumor.colorect_down)

coef <- paste(column, tissue_type_A, 'vs', tissue_type_B, sep='_')
res <- DeseqResult(dds, column, coef, tissue_type_A, tissue_type_B,
                   lfc.Threshold, rpm.Threshold,
                   norm_adj_up,
                   norm_adj_down,
                   pCRC_adj_up,
                   pCRC_adj_down)

dict_sig_mirna[paste(coef, "up",   sep='_')] <- list(res$up_mirna)
dict_sig_mirna[paste(coef, "down", sep='_')] <- list(res$down_mirna)
res_res <- res$res
res_dict[coef] <- res_res
plotMA(res$res, alpha=0.05)
```

![](comet_analysis_report_automated_no_lfcthreshold_08.12.20_res_alpha_0.05_with_LFC_shrinkage_LFC_0.58_mirge3_7.0_files/figure-gfm/unnamed-chunk-30-1.png)<!-- -->

``` r
# Plot volcano plot
VolcanoPlot(res$res, coef, res$sig,
            res$up_mirna, res$down_mirna,
            norm_adj_up, norm_adj_down,
            pCRC_adj_up, pCRC_adj_down)
```

![](comet_analysis_report_automated_no_lfcthreshold_08.12.20_res_alpha_0.05_with_LFC_shrinkage_LFC_0.58_mirge3_7.0_files/figure-gfm/unnamed-chunk-30-2.png)<!-- -->

``` r
ExpressionPlot(res$res, res$rpm, coef, res$sig,
               tissue_type_A, tissue_type_B,
               res$up_mirna, res$down_mirna,
               norm_adj_up, norm_adj_down,
               pCRC_adj_up, pCRC_adj_down)
```

![](comet_analysis_report_automated_no_lfcthreshold_08.12.20_res_alpha_0.05_with_LFC_shrinkage_LFC_0.58_mirge3_7.0_files/figure-gfm/unnamed-chunk-30-3.png)<!-- -->

``` r
signature_mirnas <- SigList(res, dds, tissue_type_A, tissue_type_B, coef,
                            norm_adj_up, norm_adj_down, 
                            pCRC_adj_up, pCRC_adj_down)

# Print list upregulated miRNA
print("as no normal adjacent PC tissue was available, control was union of lung and liver normal adjacent")
```

    ## [1] "as no normal adjacent PC tissue was available, control was union of lung and liver normal adjacent"

``` r
signature_mirnas$up_mirna
```

<div style="border: 1px solid #ddd; padding: 5px; overflow-x: scroll; width:2000px; ">

<table class="table table-striped table-hover table-condensed" style="margin-left: auto; margin-right: auto;">

<caption>

Upregulated in tissue.type\_metastasis.pc\_vs\_tumor.colorect

</caption>

<thead>

<tr>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

miRNA

</th>

<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">

LFC

</th>

<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">

lfcSE

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

FDR

</th>

<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">

RPM metastasis.pc

</th>

<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">

RPM tumor.colorect

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

miRBase\_ID

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

Family

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

Seed

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

Chr

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

Cell-Type Specific

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

Norm Background

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

pCRC Background

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

Hsa-Mir-155\_5p

</td>

<td style="text-align:right;">

0.70

</td>

<td style="text-align:right;">

0.17

</td>

<td style="text-align:left;">

1.59e-04

</td>

<td style="text-align:right;">

654

</td>

<td style="text-align:right;">

387

</td>

<td style="text-align:left;">

hsa-mir-155

</td>

<td style="text-align:left;">

MIR-155

</td>

<td style="text-align:left;">

UAAUGCU

</td>

<td style="text-align:left;">

chr21

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">c(“Lymphocyte”,
“Macrophage”)</span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-506-P4a1\_3p/P4a2\_3p/P4b\_3p

</td>

<td style="text-align:right;">

2.72

</td>

<td style="text-align:right;">

0.29

</td>

<td style="text-align:left;">

5.75e-11

</td>

<td style="text-align:right;">

191

</td>

<td style="text-align:right;">

14

</td>

<td style="text-align:left;">

hsa-mir-509-1

</td>

<td style="text-align:left;">

MIR-506

</td>

<td style="text-align:left;">

GAUUGGU

</td>

<td style="text-align:left;">

chrX

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-506-P3\_3p

</td>

<td style="text-align:right;">

2.69

</td>

<td style="text-align:right;">

0.29

</td>

<td style="text-align:left;">

6.35e-09

</td>

<td style="text-align:right;">

112

</td>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

hsa-mir-508

</td>

<td style="text-align:left;">

MIR-506

</td>

<td style="text-align:left;">

GAUUGUA

</td>

<td style="text-align:left;">

chrX

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-154-P9\_3p

</td>

<td style="text-align:right;">

0.85

</td>

<td style="text-align:right;">

0.14

</td>

<td style="text-align:left;">

1.43e-08

</td>

<td style="text-align:right;">

455

</td>

<td style="text-align:right;">

253

</td>

<td style="text-align:left;">

hsa-mir-381

</td>

<td style="text-align:left;">

MIR-154

</td>

<td style="text-align:left;">

AUACAAG

</td>

<td style="text-align:left;">

chr14

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-154-P36\_3p

</td>

<td style="text-align:right;">

0.74

</td>

<td style="text-align:right;">

0.13

</td>

<td style="text-align:left;">

7.13e-07

</td>

<td style="text-align:right;">

330

</td>

<td style="text-align:right;">

195

</td>

<td style="text-align:left;">

hsa-mir-409

</td>

<td style="text-align:left;">

MIR-154

</td>

<td style="text-align:left;">

AAUGUUG

</td>

<td style="text-align:left;">

chr14

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-210\_3p

</td>

<td style="text-align:right;">

0.71

</td>

<td style="text-align:right;">

0.18

</td>

<td style="text-align:left;">

3.37e-04

</td>

<td style="text-align:right;">

491

</td>

<td style="text-align:right;">

291

</td>

<td style="text-align:left;">

hsa-mir-210

</td>

<td style="text-align:left;">

MIR-210

</td>

<td style="text-align:left;">

UGUGCGU

</td>

<td style="text-align:left;">

chr11

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-1307\_5p

</td>

<td style="text-align:right;">

0.62

</td>

<td style="text-align:right;">

0.17

</td>

<td style="text-align:left;">

8.29e-04

</td>

<td style="text-align:right;">

797

</td>

<td style="text-align:right;">

492

</td>

<td style="text-align:left;">

hsa-mir-1307

</td>

<td style="text-align:left;">

MIR-1307

</td>

<td style="text-align:left;">

CGACCGG

</td>

<td style="text-align:left;">

chr10

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-127\_3p

</td>

<td style="text-align:right;">

0.61

</td>

<td style="text-align:right;">

0.20

</td>

<td style="text-align:left;">

1.42e-02

</td>

<td style="text-align:right;">

2888

</td>

<td style="text-align:right;">

2263

</td>

<td style="text-align:left;">

hsa-mir-127

</td>

<td style="text-align:left;">

MIR-127

</td>

<td style="text-align:left;">

CGGAUCC

</td>

<td style="text-align:left;">

chr14

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-191\_5p

</td>

<td style="text-align:right;">

0.61

</td>

<td style="text-align:right;">

0.15

</td>

<td style="text-align:left;">

2.41e-04

</td>

<td style="text-align:right;">

20429

</td>

<td style="text-align:right;">

13630

</td>

<td style="text-align:left;">

hsa-mir-191

</td>

<td style="text-align:left;">

MIR-191

</td>

<td style="text-align:left;">

AACGGAA

</td>

<td style="text-align:left;">

chr3

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-150\_5p

</td>

<td style="text-align:right;">

0.70

</td>

<td style="text-align:right;">

0.21

</td>

<td style="text-align:left;">

2.45e-03

</td>

<td style="text-align:right;">

528

</td>

<td style="text-align:right;">

280

</td>

<td style="text-align:left;">

hsa-mir-150

</td>

<td style="text-align:left;">

MIR-150

</td>

<td style="text-align:left;">

CUCCCAA

</td>

<td style="text-align:left;">

chr19

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">Lymphocyte</span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-154-P13\_5p

</td>

<td style="text-align:right;">

0.94

</td>

<td style="text-align:right;">

0.14

</td>

<td style="text-align:left;">

2.08e-09

</td>

<td style="text-align:right;">

628

</td>

<td style="text-align:right;">

328

</td>

<td style="text-align:left;">

hsa-mir-411

</td>

<td style="text-align:left;">

MIR-154

</td>

<td style="text-align:left;">

AGUAGAC

</td>

<td style="text-align:left;">

chr14

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-15-P2a\_5p/P2b\_5p

</td>

<td style="text-align:right;">

0.80

</td>

<td style="text-align:right;">

0.11

</td>

<td style="text-align:left;">

7.44e-11

</td>

<td style="text-align:right;">

8204

</td>

<td style="text-align:right;">

4522

</td>

<td style="text-align:left;">

hsa-mir-16-1

</td>

<td style="text-align:left;">

MIR-15

</td>

<td style="text-align:left;">

AGCAGCA

</td>

<td style="text-align:left;">

chr13

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-136\_3p

</td>

<td style="text-align:right;">

0.69

</td>

<td style="text-align:right;">

0.14

</td>

<td style="text-align:left;">

1.18e-05

</td>

<td style="text-align:right;">

274

</td>

<td style="text-align:right;">

173

</td>

<td style="text-align:left;">

hsa-mir-136

</td>

<td style="text-align:left;">

MIR-136

</td>

<td style="text-align:left;">

AUCAUCG

</td>

<td style="text-align:left;">

chr14

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-335\_5p

</td>

<td style="text-align:right;">

0.59

</td>

<td style="text-align:right;">

0.11

</td>

<td style="text-align:left;">

2.33e-06

</td>

<td style="text-align:right;">

232

</td>

<td style="text-align:right;">

150

</td>

<td style="text-align:left;">

hsa-mir-335

</td>

<td style="text-align:left;">

MIR-335

</td>

<td style="text-align:left;">

CAAGAGC

</td>

<td style="text-align:left;">

chr7

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">Retinal
Epithelial Cell</span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-122\_5p

</td>

<td style="text-align:right;">

1.45

</td>

<td style="text-align:right;">

0.29

</td>

<td style="text-align:left;">

9.08e-08

</td>

<td style="text-align:right;">

244

</td>

<td style="text-align:right;">

9

</td>

<td style="text-align:left;">

hsa-mir-122

</td>

<td style="text-align:left;">

MIR-122

</td>

<td style="text-align:left;">

GGAGUGU

</td>

<td style="text-align:left;">

chr18

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">Hepatocyte</span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-486\_5p

</td>

<td style="text-align:right;">

0.76

</td>

<td style="text-align:right;">

0.22

</td>

<td style="text-align:left;">

2.47e-03

</td>

<td style="text-align:right;">

4442

</td>

<td style="text-align:right;">

2013

</td>

<td style="text-align:left;">

hsa-mir-486-1

</td>

<td style="text-align:left;">

MIR-486

</td>

<td style="text-align:left;">

CCUGUAC

</td>

<td style="text-align:left;">

chr8

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">c(“Platelet”,
“Red Blood Cell”)</span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-126\_5p

</td>

<td style="text-align:right;">

0.63

</td>

<td style="text-align:right;">

0.14

</td>

<td style="text-align:left;">

2.84e-05

</td>

<td style="text-align:right;">

5738

</td>

<td style="text-align:right;">

3455

</td>

<td style="text-align:left;">

hsa-mir-126

</td>

<td style="text-align:left;">

MIR-126

</td>

<td style="text-align:left;">

AUUAUUA

</td>

<td style="text-align:left;">

chr9

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">c(“Endothelial
Cell”, “Platelet”)</span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-342\_3p

</td>

<td style="text-align:right;">

0.76

</td>

<td style="text-align:right;">

0.16

</td>

<td style="text-align:left;">

1.39e-05

</td>

<td style="text-align:right;">

275

</td>

<td style="text-align:right;">

145

</td>

<td style="text-align:left;">

hsa-mir-342

</td>

<td style="text-align:left;">

MIR-342

</td>

<td style="text-align:left;">

CUCACAC

</td>

<td style="text-align:left;">

chr14

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">c(“Dendritic
Cell”, “Lymphocyte”, “Macrophage”)</span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-10-P2c\_5p

</td>

<td style="text-align:right;">

1.10

</td>

<td style="text-align:right;">

0.20

</td>

<td style="text-align:left;">

1.07e-06

</td>

<td style="text-align:right;">

349

</td>

<td style="text-align:right;">

142

</td>

<td style="text-align:left;">

hsa-mir-99a

</td>

<td style="text-align:left;">

MIR-10

</td>

<td style="text-align:left;">

ACCCGUA

</td>

<td style="text-align:left;">

chr21

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Let-7-P1c\_5p

</td>

<td style="text-align:right;">

0.95

</td>

<td style="text-align:right;">

0.18

</td>

<td style="text-align:left;">

1.40e-06

</td>

<td style="text-align:right;">

1282

</td>

<td style="text-align:right;">

606

</td>

<td style="text-align:left;">

hsa-let-7c

</td>

<td style="text-align:left;">

LET-7

</td>

<td style="text-align:left;">

GAGGUAG

</td>

<td style="text-align:left;">

chr21

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-10-P3c\_5p

</td>

<td style="text-align:right;">

0.80

</td>

<td style="text-align:right;">

0.19

</td>

<td style="text-align:left;">

1.02e-04

</td>

<td style="text-align:right;">

497

</td>

<td style="text-align:right;">

239

</td>

<td style="text-align:left;">

hsa-mir-125b-2

</td>

<td style="text-align:left;">

MIR-10

</td>

<td style="text-align:left;">

CCCUGAG

</td>

<td style="text-align:left;">

chr21

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-154-P23\_3p

</td>

<td style="text-align:right;">

0.80

</td>

<td style="text-align:right;">

0.13

</td>

<td style="text-align:left;">

5.21e-08

</td>

<td style="text-align:right;">

256

</td>

<td style="text-align:right;">

148

</td>

<td style="text-align:left;">

hsa-mir-654

</td>

<td style="text-align:left;">

MIR-154

</td>

<td style="text-align:left;">

AUGUCUG

</td>

<td style="text-align:left;">

chr14

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-130-P1a\_3p

</td>

<td style="text-align:right;">

0.70

</td>

<td style="text-align:right;">

0.13

</td>

<td style="text-align:left;">

1.40e-06

</td>

<td style="text-align:right;">

728

</td>

<td style="text-align:right;">

421

</td>

<td style="text-align:left;">

hsa-mir-130a

</td>

<td style="text-align:left;">

MIR-130

</td>

<td style="text-align:left;">

AGUGCAA

</td>

<td style="text-align:left;">

chr11

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-15-P2c\_5p

</td>

<td style="text-align:right;">

0.60

</td>

<td style="text-align:right;">

0.15

</td>

<td style="text-align:left;">

3.85e-04

</td>

<td style="text-align:right;">

413

</td>

<td style="text-align:right;">

257

</td>

<td style="text-align:left;">

hsa-mir-195

</td>

<td style="text-align:left;">

MIR-15

</td>

<td style="text-align:left;">

AGCAGCA

</td>

<td style="text-align:left;">

chr17

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-26-P1\_5p/P2\_5p

</td>

<td style="text-align:right;">

0.59

</td>

<td style="text-align:right;">

0.09

</td>

<td style="text-align:left;">

3.22e-09

</td>

<td style="text-align:right;">

52635

</td>

<td style="text-align:right;">

35584

</td>

<td style="text-align:left;">

hsa-mir-26b

</td>

<td style="text-align:left;">

MIR-26

</td>

<td style="text-align:left;">

UCAAGUA

</td>

<td style="text-align:left;">

chr2

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

</tr>

</tbody>

</table>

</div>

``` r
# Number of upregulated miRNA
signature_mirnas$number_upregulated
```

    ## [1] 25

``` r
# Print list downregulated miRNA
print("as no normal adjacent PC tissue was available, control was union of lung and liver normal adjacent")
```

    ## [1] "as no normal adjacent PC tissue was available, control was union of lung and liver normal adjacent"

``` r
signature_mirnas$down_mirna
```

<div style="border: 1px solid #ddd; padding: 5px; overflow-x: scroll; width:2000px; ">

<table class="table table-striped table-hover table-condensed" style="margin-left: auto; margin-right: auto;">

<caption>

Downregulated in tissue.type\_metastasis.pc\_vs\_tumor.colorect

</caption>

<thead>

<tr>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

miRNA

</th>

<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">

LFC

</th>

<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">

lfcSE

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

FDR

</th>

<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">

RPM metastasis.pc

</th>

<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">

RPM tumor.colorect

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

miRBase\_ID

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

Family

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

Seed

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

Chr

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

Cell-Type Specific

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

Norm Background

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

pCRC Background

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

Hsa-Mir-223\_3p

</td>

<td style="text-align:right;">

\-0.94

</td>

<td style="text-align:right;">

0.21

</td>

<td style="text-align:left;">

4.62e-05

</td>

<td style="text-align:right;">

280

</td>

<td style="text-align:right;">

619

</td>

<td style="text-align:left;">

hsa-mir-223

</td>

<td style="text-align:left;">

MIR-223

</td>

<td style="text-align:left;">

GUCAGUU

</td>

<td style="text-align:left;">

chrX

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">c(“Dendritic
Cell”, “Macrophage”)</span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-143\_3p

</td>

<td style="text-align:right;">

\-0.86

</td>

<td style="text-align:right;">

0.17

</td>

<td style="text-align:left;">

6.37e-06

</td>

<td style="text-align:right;">

77753

</td>

<td style="text-align:right;">

148431

</td>

<td style="text-align:left;">

hsa-mir-143

</td>

<td style="text-align:left;">

MIR-143

</td>

<td style="text-align:left;">

GAGAUGA

</td>

<td style="text-align:left;">

chr5

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">Mesenchymal</span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-133-P1\_3p/P2\_3p/P3\_3p

</td>

<td style="text-align:right;">

\-1.82

</td>

<td style="text-align:right;">

0.23

</td>

<td style="text-align:left;">

3.20e-14

</td>

<td style="text-align:right;">

25

</td>

<td style="text-align:right;">

116

</td>

<td style="text-align:left;">

hsa-mir-133a-2

</td>

<td style="text-align:left;">

MIR-133

</td>

<td style="text-align:left;">

UUGGUCC

</td>

<td style="text-align:left;">

chr20

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">c(“Skeletal
Myocyte”, “Stem Cell”)</span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

</tr>

</tbody>

</table>

</div>

``` r
# Number of downregulated miRNA
signature_mirnas$number_downregulated
```

    ## [1] 3

``` r
res_tibble <- res$res
res_tibble$miRNA <- rownames(res_tibble)
res_tibble <- as_tibble(res_tibble)

metslfc <- res_dict$tissue.type_metastasis.pc_vs_tumor.colorect$log2FoldChange
normlfc <- ( (res_dict$tissue.type_normal.lung_vs_normal.colorect$log2FoldChange + res_dict$tissue.type_normal.liver_vs_normal.colorect$log2FoldChange) / 2)

res_tibble$LFC_adj_background <- mapply(SubtractLFC, metslfc, normlfc)

metsP <- res_dict$tissue.type_metastasis.pc_vs_tumor.colorect$padj
normP <- ( (res_dict$tissue.type_normal.lung_vs_normal.colorect$padj - res_dict$tissue.type_normal.liver_vs_normal.colorect$padj) / 2 )

res_tibble$padj_subt_normal <- mapply(SubtractAdjP, metslfc, normlfc,metsP, normP)


res_tibble %>% select(miRNA, log2FoldChange, lfcSE, LFC_adj_background, padj_subt_normal, baseMean, stat, pvalue, padj) %>% write_csv(path = '/Users/eirikhoy/Dropbox/projects/comet_analysis/data/Deseq_result_pc_vs_pcrc.csv')
```

### All mCRC vs pCRC

``` r
DeseqObject <- function(DESIGN, countdata, coldata, consensus="None", sample_type="None", Ref) {
  "
  Function to create DESeq2 object
  "
  
  dds <- DESeqDataSetFromMatrix(countData = countdata,
                                colData = coldata,
                                design = as.formula(paste("~", DESIGN)))
    # Kick out non-consensus samples
  if (!(consensus == "None")) {
    dds <- dds[, dds$paper %in% consensus]
  }

  # Kick out samples that are not bulk tissue
  if (!(sample_type == "None")) {
    dds <- dds[, dds$sample_type == sample_type]
  }

  dds$type <- relevel(dds$type, ref=ref)
  dds$type <- droplevels(dds$type)
  
  dds <- DESeq(dds,
               parallel=TRUE,
               BPPARAM=MulticoreParam(3)
               )

  return(dds)
}
```

``` r
SigList <- function(res, dds, tissue_type_A, tissue_type_B, coef,
                            norm_adj_up, norm_adj_down, 
                            pCRC_adj_up, pCRC_adj_down){
  "
  Function to create annotated lists of signature miRNA
  Return will print upregulated or downregulated miRNA, 
  by printing <signature_list>$up_mirna
  or          <signature_list>$down_mirna
  "
  group_A_rpm <- rowMeans(res$rpm[res$sig, dds$type == tissue_type_A])
  group_A_rpm_std <- rowSds(res$rpm[res$sig, dds$type == tissue_type_A])
  group_B_rpm <- rowMeans(res$rpm[res$sig, dds$type == tissue_type_B])
  group_B_rpm_std <- rowSds(res$rpm[res$sig, dds$type == tissue_type_B])
  lfc.deseq2  <- res$res[res$sig, ]$log2FoldChange
  lfcSE.deseq2<- res$res[res$sig, ]$lfcSE
  neg.log.10.adj.p <- -log10(res$res[res$sig, ]$padj)
  signature_mirna <- res$sig
  sig_list <- dplyr::tibble(signature_mirna, lfc.deseq2, lfcSE.deseq2,
                            group_A_rpm, #group_A_rpm_std,
                            group_B_rpm, #group_B_rpm_std,
                            neg.log.10.adj.p)
  
  # create list of upregulated mirna
  up_mirna <- sig_list %>%
    filter(lfc.deseq2 > lfc.Threshold) %>% 
    
    # Annotate which miRNA are cell markers
    mutate(
      cell_marker = ifelse(signature_mirna %in% names(cell_spec_dict_inv), cell_spec_dict_inv[signature_mirna], '')) %>%
    mutate(
      cell_marker = cell_spec(cell_marker, color = ifelse(cell_marker != '', 'white', 'black'),
                              background = ifelse(cell_marker != '', 'blue', 'white'),
                              bold = ifelse(cell_marker != '', F, F))) %>%
  # Annotate which miRNA are present in blood cells
    mutate(
      blood_cell = ifelse(signature_mirna %in% blood.cell.mirna, 'yes', '')) %>%
    mutate(
      blood_cell = cell_spec(blood_cell, color = ifelse(blood_cell == 'yes', 'white', 'black'),
                              background = ifelse(blood_cell == 'yes', 'red', 'white'),
                              bold = ifelse(blood_cell == 'yes', T, F)))

    # Annotate which miRNA are in normal_adjacent
    if (norm_adj_up != "None") {
      up_mirna <- up_mirna %>%
        mutate(
          norm_adj = ifelse(signature_mirna %in% norm_adj_up, 'yes', '')) %>%
        mutate(
          norm_adj = cell_spec(norm_adj, color = ifelse(norm_adj == 'yes', 'white', 'black'),
                               background = ifelse(norm_adj == 'yes', 'black', 'white'),
                               bold = ifelse(norm_adj == 'yes', T, F))
        )
    }
    else up_mirna$norm_adj <- "na"
    
  # Annotate which miRNA are in pCRC_adjacent
    if (pCRC_adj_up != "None") {
      up_mirna <- up_mirna %>%
        mutate(
          pCRC_adj = ifelse(signature_mirna %in% pCRC_adj_up, 'yes', '')) %>%
        mutate(
          pCRC_adj = cell_spec(pCRC_adj, color = ifelse(pCRC_adj == 'yes', 'white', 'black'),
                               background = ifelse(pCRC_adj == 'yes', 'black', 'white'),
                               bold = ifelse(pCRC_adj == 'yes', T, F))
        )
    }
    else up_mirna$pCRC_adj <- "na"

    # number of upregulated miRNA
    number_upregulated <- dim(up_mirna)[1]

  # Create kable list with annotations
    up_mirna <- up_mirna %>%
      arrange(-lfc.deseq2) %>%
      arrange(desc(cell_marker)) %>%
      arrange(pCRC_adj) %>%
      arrange(norm_adj) %>%
      kable(col.names = c("miRNA", "LFC", "lfcSE",
                          paste('RPM', tissue_type_A), #paste('std', tissue_type_A), 
                          paste('RPM', tissue_type_B), #paste('std', tissue_type_B), 
                          "-log10(adj p-value)", "cell_marker", "blood_cell", 'norm_adj', 'pCRC_adj'),
            escape = F, booktabs = F, caption = paste("Upregulated in ", coef),
            digits = c(0, 2, 2, 0, 0, 3, 2, 3, 0, 0, 0, 0)) %>%
      kable_styling(bootstrap_options = c("striped", "hover", "condensed"), full_width = T, 
                  fixed_thead = list(enabled = T)) %>%
      column_spec(2, bold = T)
    
  # create list of downregulated miRNA
  down_mirna <- sig_list %>%
    filter(lfc.deseq2 < -lfc.Threshold) %>% 
    
    # Annotate which miRNA are cell markers
    mutate(
      cell_marker = ifelse(signature_mirna %in% names(cell_spec_dict_inv), cell_spec_dict_inv[signature_mirna], '')) %>%
    mutate(
      cell_marker = cell_spec(cell_marker, color = ifelse(cell_marker != '', 'white', 'black'),
                              background = ifelse(cell_marker != '', 'blue', 'white'),
                              bold = ifelse(cell_marker != '', F, F))) %>%
  # Annotate which miRNA are present in blood cells
    mutate(
      blood_cell = ifelse(signature_mirna %in% blood.cell.mirna, 'yes', '')) %>%
    mutate(
      blood_cell = cell_spec(blood_cell, color = ifelse(blood_cell == 'yes', 'white', 'black'),
                              background = ifelse(blood_cell == 'yes', 'red', 'white'),
                              bold = ifelse(blood_cell == 'yes', T, F)))

    # Annotate which miRNA are in normal_adjacent
    if (norm_adj_down != "None") {
      down_mirna <- down_mirna %>%
        mutate(
          norm_adj = ifelse(signature_mirna %in% norm_adj_down, 'yes', '')) %>%
        mutate(
          norm_adj = cell_spec(norm_adj, color = ifelse(norm_adj == 'yes', 'white', 'black'),
                               background = ifelse(norm_adj == 'yes', 'black', 'white'),
                               bold = ifelse(norm_adj == 'yes', T, F))
        )
    }
    else down_mirna$norm_adj <- "na"

    # Annotate which miRNA are in pCRC_adjacent
    if (pCRC_adj_down != "None") {
      down_mirna <- down_mirna %>%
        mutate(
          pCRC_adj = ifelse(signature_mirna %in% pCRC_adj_down, 'yes', '')) %>%
        mutate(
          pCRC_adj = cell_spec(pCRC_adj, color = ifelse(pCRC_adj == 'yes', 'white', 'black'),
                               background = ifelse(pCRC_adj == 'yes', 'black', 'white'),
                               bold = ifelse(pCRC_adj == 'yes', T, F))
        )
    }
    else down_mirna$pCRC_adj <- "na"
  
    # number of upregulated miRNA
    number_downregulated <- dim(down_mirna)[1]

  # Create kable list with annotations    
    down_mirna <- down_mirna %>%
      arrange(lfc.deseq2) %>%
      arrange(desc(cell_marker)) %>%
      arrange(pCRC_adj) %>%
      arrange(norm_adj) %>%
      kable(col.names = c("miRNA", "LFC", "lfcSE",
                          paste('RPM', tissue_type_A), #paste('std', tissue_type_A), 
                          paste('RPM', tissue_type_B), #paste('std', tissue_type_B), 
                          "-log10(adj p-value)", "cell_marker", "blood_cell", 'norm_adj', 'pCRC_adj'),
            escape = F, booktabs = F, caption = paste("Downregulated in ", coef),
            digits = c(0, 2, 2, 0, 0, 3, 2, 3, 0, 0, 0, 0)) %>%
      kable_styling(bootstrap_options = c("striped", "hover", "condensed"), full_width = T, 
                  fixed_thead = list(enabled = T)) %>%
      column_spec(2, bold = T)
    

  # Function return is to print kable, either upregulated or downregulated miRNA
  return_list = list("up_mirna" = up_mirna, "down_mirna" = down_mirna, 
                     "number_upregulated" = number_upregulated, 
                     "number_downregulated" = number_downregulated)
  return(return_list)
}
```

``` r
design <- as.formula(~ type)
ref <- 'tumor'
dds <- DeseqObject(design, countdata, sampleinfo, "None", "None", ref)
#
```

``` r
column='type'
tissue_type_A <- 'metastasis'
tissue_type_B <- 'tumor'
norm_adj_up       = union(dict_sig_mirna$tissue.type_normal.liver_vs_normal.colorect_up, dict_sig_mirna$tissue.type_normal.lung_vs_normal.colorect_up)
norm_adj_down     = union(dict_sig_mirna$tissue.type_normal.liver_vs_normal.colorect_down, dict_sig_mirna$tissue.type_normal.lung_vs_normal.colorect_down)
pCRC_adj_up       = union(dict_sig_mirna$tissue.type_normal.liver_vs_tumor.colorect_up, dict_sig_mirna$tissue.type_normal.lung_vs_tumor.colorect_up)
pCRC_adj_down     = union(dict_sig_mirna$tissue.type_normal.liver_vs_tumor.colorect_down, dict_sig_mirna$tissue.type_normal.lung_vs_tumor.colorect_down)

coef <- paste(column, tissue_type_A, 'vs', tissue_type_B, sep='_')
res <- DeseqResult(dds, column, coef, tissue_type_A, tissue_type_B,
                   lfc.Threshold, rpm.Threshold,
                   norm_adj_up,
                   norm_adj_down,
                   pCRC_adj_up,
                   pCRC_adj_down)
```

    ## using 'normal' for LFC shrinkage, the Normal prior from Love et al (2014).
    ## 
    ## Note that type='apeglm' and type='ashr' have shown to have less bias than type='normal'.
    ## See ?lfcShrink for more details on shrinkage type, and the DESeq2 vignette.
    ## Reference: https://doi.org/10.1093/bioinformatics/bty895

    ## Warning in if (!(norm_adj_up == "None")) {: the condition has length > 1 and
    ## only the first element will be used

    ## Warning in if (!(norm_adj_down == "None")) {: the condition has length > 1 and
    ## only the first element will be used

    ## Warning in if (!(pCRC_adj_up == "None")) {: the condition has length > 1 and
    ## only the first element will be used

    ## Warning in if (!(pCRC_adj_down == "None")) {: the condition has length > 1 and
    ## only the first element will be used

``` r
dict_sig_mirna[paste(coef, "up",   sep='_')] <- list(res$up_mirna)
dict_sig_mirna[paste(coef, "down", sep='_')] <- list(res$down_mirna)
res_res <- res$res
res_dict[coef] <- res_res
```

    ## Warning in `[<-`(`*tmp*`, coef, value = new("DESeqResults", priorInfo = list(:
    ## implicit list embedding of S4 objects is deprecated

``` r
plotMA(res$res, alpha=0.05)
```

![](comet_analysis_report_automated_no_lfcthreshold_08.12.20_res_alpha_0.05_with_LFC_shrinkage_LFC_0.58_mirge3_7.0_files/figure-gfm/unnamed-chunk-34-1.png)<!-- -->

``` r
# Plot volcano plot
VolcanoPlot(res$res, coef, res$sig,
            res$up_mirna, res$down_mirna,
            norm_adj_up, norm_adj_down,
            pCRC_adj_up, pCRC_adj_down)
```

    ## Warning in if (norm_adj_up != "None" & length(up_mirna > 0)) {: the condition
    ## has length > 1 and only the first element will be used

    ## Warning in if (norm_adj_down != "None" & length(down_mirna > 0)) {: the
    ## condition has length > 1 and only the first element will be used

    ## Warning in if (pCRC_adj_up != "None" & length(up_mirna > 0)) {: the condition
    ## has length > 1 and only the first element will be used

    ## Warning in if (pCRC_adj_down != "None" & length(down_mirna > 0)) {: the
    ## condition has length > 1 and only the first element will be used

![](comet_analysis_report_automated_no_lfcthreshold_08.12.20_res_alpha_0.05_with_LFC_shrinkage_LFC_0.58_mirge3_7.0_files/figure-gfm/unnamed-chunk-34-2.png)<!-- -->

``` r
#ExpressionPlot(res$res, res$rpm, coef, res$sig,
#               tissue_type_A, tissue_type_B,
#               res$up_mirna, res$down_mirna,
#               norm_adj_up, norm_adj_down,
#               pCRC_adj_up, pCRC_adj_down)

signature_mirnas <- SigList(res, dds, tissue_type_A, tissue_type_B, coef,
                            norm_adj_up, norm_adj_down, 
                            pCRC_adj_up, pCRC_adj_down)
```

    ## Warning in if (norm_adj_up != "None") {: the condition has length > 1 and only
    ## the first element will be used

    ## Warning in if (pCRC_adj_up != "None") {: the condition has length > 1 and only
    ## the first element will be used

    ## Warning in if (norm_adj_down != "None") {: the condition has length > 1 and only
    ## the first element will be used

    ## Warning in if (pCRC_adj_down != "None") {: the condition has length > 1 and only
    ## the first element will be used

``` r
# Print list upregulated miRNA
print("control was union of lung and liver normal adjacent")
```

    ## [1] "control was union of lung and liver normal adjacent"

``` r
signature_mirnas$up_mirna
```

<table class="table table-striped table-hover table-condensed" style="margin-left: auto; margin-right: auto;">

<caption>

Upregulated in type\_metastasis\_vs\_tumor

</caption>

<thead>

<tr>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

miRNA

</th>

<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">

LFC

</th>

<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">

lfcSE

</th>

<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">

RPM metastasis

</th>

<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">

RPM tumor

</th>

<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">

\-log10(adj p-value)

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

cell\_marker

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

blood\_cell

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

norm\_adj

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

pCRC\_adj

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

Hsa-Mir-210\_3p

</td>

<td style="text-align:right;font-weight: bold;">

1.10

</td>

<td style="text-align:right;">

0.13

</td>

<td style="text-align:right;">

625

</td>

<td style="text-align:right;">

291

</td>

<td style="text-align:right;">

14.551

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-1307\_5p

</td>

<td style="text-align:right;font-weight: bold;">

0.72

</td>

<td style="text-align:right;">

0.12

</td>

<td style="text-align:right;">

837

</td>

<td style="text-align:right;">

492

</td>

<td style="text-align:right;">

7.519

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-191\_5p

</td>

<td style="text-align:right;font-weight: bold;">

0.72

</td>

<td style="text-align:right;">

0.11

</td>

<td style="text-align:right;">

23467

</td>

<td style="text-align:right;">

13630

</td>

<td style="text-align:right;">

9.458

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-8-P1b\_3p

</td>

<td style="text-align:right;font-weight: bold;">

0.66

</td>

<td style="text-align:right;">

0.17

</td>

<td style="text-align:right;">

12715

</td>

<td style="text-align:right;">

7518

</td>

<td style="text-align:right;">

3.068

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-10-P1a\_5p

</td>

<td style="text-align:right;font-weight: bold;">

0.59

</td>

<td style="text-align:right;">

0.15

</td>

<td style="text-align:right;">

140769

</td>

<td style="text-align:right;">

97123

</td>

<td style="text-align:right;">

3.469

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-150\_5p

</td>

<td style="text-align:right;font-weight: bold;">

1.05

</td>

<td style="text-align:right;">

0.16

</td>

<td style="text-align:right;">

664

</td>

<td style="text-align:right;">

280

</td>

<td style="text-align:right;">

9.458

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">Lymphocyte</span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-15-P2a\_5p/P2b\_5p

</td>

<td style="text-align:right;font-weight: bold;">

0.74

</td>

<td style="text-align:right;">

0.08

</td>

<td style="text-align:right;">

7988

</td>

<td style="text-align:right;">

4522

</td>

<td style="text-align:right;">

18.164

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-331\_3p

</td>

<td style="text-align:right;font-weight: bold;">

0.70

</td>

<td style="text-align:right;">

0.12

</td>

<td style="text-align:right;">

104

</td>

<td style="text-align:right;">

60

</td>

<td style="text-align:right;">

7.519

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-335\_5p

</td>

<td style="text-align:right;font-weight: bold;">

0.76

</td>

<td style="text-align:right;">

0.09

</td>

<td style="text-align:right;">

262

</td>

<td style="text-align:right;">

150

</td>

<td style="text-align:right;">

15.670

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">Retinal
Epithelial Cell</span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

<td style="text-align:left;">

<span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-122\_5p

</td>

<td style="text-align:right;font-weight: bold;">

1.56

</td>

<td style="text-align:right;">

0.31

</td>

<td style="text-align:right;">

1459

</td>

<td style="text-align:right;">

9

</td>

<td style="text-align:right;">

24.422

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">Hepatocyte</span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

<td style="text-align:left;">

<span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-126\_5p

</td>

<td style="text-align:right;font-weight: bold;">

0.92

</td>

<td style="text-align:right;">

0.12

</td>

<td style="text-align:right;">

7434

</td>

<td style="text-align:right;">

3455

</td>

<td style="text-align:right;">

12.397

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">c(“Endothelial
Cell”, “Platelet”)</span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

<td style="text-align:left;">

<span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-342\_3p

</td>

<td style="text-align:right;font-weight: bold;">

1.46

</td>

<td style="text-align:right;">

0.13

</td>

<td style="text-align:right;">

472

</td>

<td style="text-align:right;">

145

</td>

<td style="text-align:right;">

27.005

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">c(“Dendritic
Cell”, “Lymphocyte”, “Macrophage”)</span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

<td style="text-align:left;">

<span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-34-P2b\_5p

</td>

<td style="text-align:right;font-weight: bold;">

2.61

</td>

<td style="text-align:right;">

0.21

</td>

<td style="text-align:right;">

332

</td>

<td style="text-align:right;">

33

</td>

<td style="text-align:right;">

19.564

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

<td style="text-align:left;">

<span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-30-P1a\_5p

</td>

<td style="text-align:right;font-weight: bold;">

0.93

</td>

<td style="text-align:right;">

0.13

</td>

<td style="text-align:right;">

7100

</td>

<td style="text-align:right;">

3165

</td>

<td style="text-align:right;">

10.616

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

<td style="text-align:left;">

<span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-15-P2c\_5p

</td>

<td style="text-align:right;font-weight: bold;">

0.87

</td>

<td style="text-align:right;">

0.12

</td>

<td style="text-align:right;">

537

</td>

<td style="text-align:right;">

257

</td>

<td style="text-align:right;">

10.661

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

<td style="text-align:left;">

<span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-10-P3c\_5p

</td>

<td style="text-align:right;font-weight: bold;">

0.82

</td>

<td style="text-align:right;">

0.15

</td>

<td style="text-align:right;">

514

</td>

<td style="text-align:right;">

239

</td>

<td style="text-align:right;">

6.775

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

<td style="text-align:left;">

<span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-26-P1\_5p/P2\_5p

</td>

<td style="text-align:right;font-weight: bold;">

0.71

</td>

<td style="text-align:right;">

0.07

</td>

<td style="text-align:right;">

60843

</td>

<td style="text-align:right;">

35584

</td>

<td style="text-align:right;">

23.244

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

<td style="text-align:left;">

<span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-338-P1\_3p

</td>

<td style="text-align:right;font-weight: bold;">

0.70

</td>

<td style="text-align:right;">

0.16

</td>

<td style="text-align:right;">

169

</td>

<td style="text-align:right;">

98

</td>

<td style="text-align:right;">

4.263

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

<td style="text-align:left;">

<span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-10-P2c\_5p

</td>

<td style="text-align:right;font-weight: bold;">

0.70

</td>

<td style="text-align:right;">

0.17

</td>

<td style="text-align:right;">

264

</td>

<td style="text-align:right;">

142

</td>

<td style="text-align:right;">

3.839

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

<td style="text-align:left;">

<span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

</tr>

</tbody>

</table>

``` r
# Number of upregulated miRNA
signature_mirnas$number_upregulated
```

    ## [1] 19

``` r
# Print list downregulated miRNA
print("control was union of lung and liver normal adjacent")
```

    ## [1] "control was union of lung and liver normal adjacent"

``` r
signature_mirnas$down_mirna
```

<table class="table table-striped table-hover table-condensed" style="margin-left: auto; margin-right: auto;">

<caption>

Downregulated in type\_metastasis\_vs\_tumor

</caption>

<thead>

<tr>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

miRNA

</th>

<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">

LFC

</th>

<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">

lfcSE

</th>

<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">

RPM metastasis

</th>

<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">

RPM tumor

</th>

<th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;">

\-log10(adj p-value)

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

cell\_marker

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

blood\_cell

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

norm\_adj

</th>

<th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">

pCRC\_adj

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

Hsa-Mir-7-P1\_5p/P2\_5p/P3\_5p

</td>

<td style="text-align:right;font-weight: bold;">

\-0.60

</td>

<td style="text-align:right;">

0.18

</td>

<td style="text-align:right;">

112

</td>

<td style="text-align:right;">

199

</td>

<td style="text-align:right;">

2.657

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">c(“Islet
Cell”, “Neural”)</span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-31\_5p

</td>

<td style="text-align:right;font-weight: bold;">

\-0.68

</td>

<td style="text-align:right;">

0.23

</td>

<td style="text-align:right;">

339

</td>

<td style="text-align:right;">

627

</td>

<td style="text-align:right;">

2.235

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-143\_3p

</td>

<td style="text-align:right;font-weight: bold;">

\-1.08

</td>

<td style="text-align:right;">

0.14

</td>

<td style="text-align:right;">

72990

</td>

<td style="text-align:right;">

148431

</td>

<td style="text-align:right;">

13.050

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">Mesenchymal</span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

<td style="text-align:left;">

<span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

</tr>

<tr>

<td style="text-align:left;">

Hsa-Mir-133-P1\_3p/P2\_3p/P3\_3p

</td>

<td style="text-align:right;font-weight: bold;">

\-1.80

</td>

<td style="text-align:right;">

0.19

</td>

<td style="text-align:right;">

33

</td>

<td style="text-align:right;">

116

</td>

<td style="text-align:right;">

18.097

</td>

<td style="text-align:left;">

<span style="     color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: blue !important;">c(“Skeletal
Myocyte”, “Stem Cell”)</span>

</td>

<td style="text-align:left;">

<span style="     color: black !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: white !important;"></span>

</td>

<td style="text-align:left;">

<span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

<td style="text-align:left;">

<span style=" font-weight: bold;    color: white !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: black !important;">yes</span>

</td>

</tr>

</tbody>

</table>

``` r
# Number of downregulated miRNA
signature_mirnas$number_downregulated
```

    ## [1] 4

``` r
res_dict
```

    ## $tissue.type_normal.liver_vs_normal.colorect
    ## log2 fold change (MAP): tissue.type normal.liver vs normal.colorect 
    ## Wald test p-value: tissue.type normal.liver vs normal.colorect 
    ## DataFrame with 389 rows and 6 columns
    ##                                          baseMean     log2FoldChange
    ##                                         <numeric>          <numeric>
    ## Hsa-Let-7-P1a_5p/P2a1_5p/P2a2_5p 89174.7542403281 -0.114702438165718
    ## Hsa-Let-7-P1b_5p                 3359.93816553835 -0.320114053507846
    ## Hsa-Let-7-P1c_5p                 3992.69006418194   3.18371450404926
    ## Hsa-Let-7-P2a3_5p                52000.6252125236  0.135005279299899
    ## Hsa-Let-7-P2b1_5p                9751.60360980849  0.321769922951614
    ## ...                                           ...                ...
    ## Hsa-Mir-95-P2_3p                 390.216333625945 -0.259423173236495
    ## Hsa-Mir-95-P3_5p                 3.95328293251241  0.578777309350614
    ## Hsa-Mir-96-P1_5p                 202.656949570096  -2.86706387684201
    ## Hsa-Mir-96-P2_5p                 27209.4573430933  -2.64828835135873
    ## Hsa-Mir-96-P3_5p                 3344.18853507094  -3.56720648489265
    ##                                              lfcSE               stat
    ##                                          <numeric>          <numeric>
    ## Hsa-Let-7-P1a_5p/P2a1_5p/P2a2_5p 0.137100418774597 -0.876702398466493
    ## Hsa-Let-7-P1b_5p                 0.223907107455248   -1.4689025908534
    ## Hsa-Let-7-P1c_5p                 0.279081137786442   11.1516794079232
    ## Hsa-Let-7-P2a3_5p                0.120493686572101   1.11711255877344
    ## Hsa-Let-7-P2b1_5p                0.122630483027464   2.60055162322238
    ## ...                                            ...                ...
    ## Hsa-Mir-95-P2_3p                  0.19641408635223  -1.08552839444029
    ## Hsa-Mir-95-P3_5p                 0.384387609595463   1.62826495421394
    ## Hsa-Mir-96-P1_5p                 0.333638316705081  -7.80833454343775
    ## Hsa-Mir-96-P2_5p                 0.243241112732279  -10.2148298811257
    ## Hsa-Mir-96-P3_5p                 0.324276388539775  -10.0394868651728
    ##                                                pvalue                 padj
    ##                                             <numeric>            <numeric>
    ## Hsa-Let-7-P1a_5p/P2a1_5p/P2a2_5p    0.380648303687109    0.543582632940632
    ## Hsa-Let-7-P1b_5p                    0.141859211826177    0.240787346389169
    ## Hsa-Let-7-P1c_5p                 7.02674679240312e-29 1.59961824038824e-27
    ## Hsa-Let-7-P2a3_5p                   0.263946201391012    0.404220147427516
    ## Hsa-Let-7-P2b1_5p                  0.0093074014267008   0.0220979408106332
    ## ...                                               ...                  ...
    ## Hsa-Mir-95-P2_3p                    0.277687694717835    0.421431913160009
    ## Hsa-Mir-95-P3_5p                    0.103468717163438     0.18284197964498
    ## Hsa-Mir-96-P1_5p                 5.79486025768406e-15 5.33954980886603e-14
    ## Hsa-Mir-96-P2_5p                  1.7017739490561e-24 2.99357508311232e-23
    ## Hsa-Mir-96-P3_5p                 1.02204432533755e-23 1.71970066915492e-22
    ## 
    ## $tissue.type_normal.lung_vs_normal.colorect
    ## log2 fold change (MAP): tissue.type normal.lung vs normal.colorect 
    ## Wald test p-value: tissue.type normal.lung vs normal.colorect 
    ## DataFrame with 389 rows and 6 columns
    ##                                          baseMean      log2FoldChange
    ##                                         <numeric>           <numeric>
    ## Hsa-Let-7-P1a_5p/P2a1_5p/P2a2_5p 89174.7542403281   0.159488915677012
    ## Hsa-Let-7-P1b_5p                 3359.93816553835   0.197555867999725
    ## Hsa-Let-7-P1c_5p                 3992.69006418194    1.96610542233694
    ## Hsa-Let-7-P2a3_5p                52000.6252125236   0.289175178805092
    ## Hsa-Let-7-P2b1_5p                9751.60360980849    1.10414120876771
    ## ...                                           ...                 ...
    ## Hsa-Mir-95-P2_3p                 390.216333625945   0.165907582077109
    ## Hsa-Mir-95-P3_5p                 3.95328293251241 -0.0273303033157357
    ## Hsa-Mir-96-P1_5p                 202.656949570096  -0.415869888949985
    ## Hsa-Mir-96-P2_5p                 27209.4573430933  -0.107488420248239
    ## Hsa-Mir-96-P3_5p                 3344.18853507094  -0.456955953013122
    ##                                              lfcSE                stat
    ##                                          <numeric>           <numeric>
    ## Hsa-Let-7-P1a_5p/P2a1_5p/P2a2_5p 0.171252727891764   0.886346895351987
    ## Hsa-Let-7-P1b_5p                 0.280263432235887   0.631339921419183
    ## Hsa-Let-7-P1c_5p                 0.350196583213482    5.61986491346816
    ## Hsa-Let-7-P2a3_5p                 0.15044933221432    1.91622249772382
    ## Hsa-Let-7-P2b1_5p                0.153065517284977    7.17242652538056
    ## ...                                            ...                 ...
    ## Hsa-Mir-95-P2_3p                 0.243046643416424   0.843104626057063
    ## Hsa-Mir-95-P3_5p                 0.458035265321474   0.143432621264712
    ## Hsa-Mir-96-P1_5p                 0.388545838443054  -0.690134766070179
    ## Hsa-Mir-96-P2_5p                 0.304517199045329 -0.0704700708063357
    ## Hsa-Mir-96-P3_5p                 0.404724536391786  -0.772763885607241
    ##                                                pvalue                 padj
    ##                                             <numeric>            <numeric>
    ## Hsa-Let-7-P1a_5p/P2a1_5p/P2a2_5p    0.375430626009251    0.530261504618906
    ## Hsa-Let-7-P1b_5p                     0.52781828944235    0.684422033196025
    ## Hsa-Let-7-P1c_5p                 1.91106824923059e-08 1.57358172862179e-07
    ## Hsa-Let-7-P2a3_5p                  0.0553367808976523    0.104465044914105
    ## Hsa-Let-7-P2b1_5p                7.36799615073784e-13 1.05607944827242e-11
    ## ...                                               ...                  ...
    ## Hsa-Mir-95-P2_3p                    0.399169931746439    0.555679005704575
    ## Hsa-Mir-95-P3_5p                    0.885948521274591    0.960397976843884
    ## Hsa-Mir-96-P1_5p                    0.490109441851604    0.651795030916051
    ## Hsa-Mir-96-P2_5p                    0.943819521347161     0.98187676011116
    ## Hsa-Mir-96-P3_5p                    0.439662130177162    0.603366114817595
    ## 
    ## $tissue.type_tumor.colorect_vs_normal.colorect
    ## log2 fold change (MAP): tissue.type tumor.colorect vs normal.colorect 
    ## Wald test p-value: tissue.type tumor.colorect vs normal.colorect 
    ## DataFrame with 389 rows and 6 columns
    ##                                          baseMean      log2FoldChange
    ##                                         <numeric>           <numeric>
    ## Hsa-Let-7-P1a_5p/P2a1_5p/P2a2_5p 89174.7542403281  -0.248058505837001
    ## Hsa-Let-7-P1b_5p                 3359.93816553835  -0.129712872859404
    ## Hsa-Let-7-P1c_5p                 3992.69006418194 -0.0894387113372064
    ## Hsa-Let-7-P2a3_5p                52000.6252125236   0.161274363696562
    ## Hsa-Let-7-P2b1_5p                9751.60360980849  -0.129885530327206
    ## ...                                           ...                 ...
    ## Hsa-Mir-95-P2_3p                 390.216333625945    1.16639047991208
    ## Hsa-Mir-95-P3_5p                 3.95328293251241   0.128041251542255
    ## Hsa-Mir-96-P1_5p                 202.656949570096    1.99605999628903
    ## Hsa-Mir-96-P2_5p                 27209.4573430933     1.8371528804308
    ## Hsa-Mir-96-P3_5p                 3344.18853507094     1.5305097819346
    ##                                               lfcSE               stat
    ##                                           <numeric>          <numeric>
    ## Hsa-Let-7-P1a_5p/P2a1_5p/P2a2_5p   0.09974589771869  -2.50880006172293
    ## Hsa-Let-7-P1b_5p                  0.160828091638605 -0.872616006654668
    ## Hsa-Let-7-P1c_5p                  0.198255228124062  -0.14744126618236
    ## Hsa-Let-7-P2a3_5p                0.0878077644317027   1.81863066456019
    ## Hsa-Let-7-P2b1_5p                0.0893412951279261  -1.43703312496472
    ## ...                                             ...                ...
    ## Hsa-Mir-95-P2_3p                  0.140554605668826   8.28848119000958
    ## Hsa-Mir-95-P3_5p                  0.264177931179274  0.737283525185868
    ## Hsa-Mir-96-P1_5p                  0.217577914511228   8.93516022919441
    ## Hsa-Mir-96-P2_5p                  0.173955103251303   10.4458859814953
    ## Hsa-Mir-96-P3_5p                  0.226074853962932   6.65702796429079
    ##                                                pvalue                 padj
    ##                                             <numeric>            <numeric>
    ## Hsa-Let-7-P1a_5p/P2a1_5p/P2a2_5p   0.0121142030924821   0.0250889779893481
    ## Hsa-Let-7-P1b_5p                     0.38287241295549    0.480492200364671
    ## Hsa-Let-7-P1c_5p                    0.882783735719796    0.926800653290313
    ## Hsa-Let-7-P2a3_5p                  0.0689677962587844    0.122143900850838
    ## Hsa-Let-7-P2b1_5p                   0.150708581526998    0.225175837320617
    ## ...                                               ...                  ...
    ## Hsa-Mir-95-P2_3p                  1.1470369823237e-16 1.40234521387317e-15
    ## Hsa-Mir-95-P3_5p                    0.460949948804283    0.561736432787213
    ## Hsa-Mir-96-P1_5p                 4.06588568948011e-19 8.56094820173867e-18
    ## Hsa-Mir-96-P2_5p                 1.53019828709604e-25 6.44383500899331e-24
    ## Hsa-Mir-96-P3_5p                 2.79420040427442e-11 1.96111472818519e-10
    ## 
    ## $tissue.type_normal.liver_vs_tumor.colorect
    ## log2 fold change (MAP): tissue.type normal.liver vs tumor.colorect 
    ## Wald test p-value: tissue.type normal.liver vs tumor.colorect 
    ## DataFrame with 389 rows and 6 columns
    ##                                          baseMean      log2FoldChange
    ##                                         <numeric>           <numeric>
    ## Hsa-Let-7-P1a_5p/P2a1_5p/P2a2_5p 89174.7542403281   0.134516855928531
    ## Hsa-Let-7-P1b_5p                 3359.93816553835   -0.18431192252938
    ## Hsa-Let-7-P1c_5p                 3992.69006418194    3.24635633273589
    ## Hsa-Let-7-P2a3_5p                52000.6252125236 -0.0245108495144203
    ## Hsa-Let-7-P2b1_5p                9751.60360980849   0.451258084717704
    ## ...                                           ...                 ...
    ## Hsa-Mir-95-P2_3p                 390.216333625945   -1.41840794320829
    ## Hsa-Mir-95-P3_5p                 3.95328293251241   0.414413322947411
    ## Hsa-Mir-96-P1_5p                 202.656949570096   -4.81237877419974
    ## Hsa-Mir-96-P2_5p                 27209.4573430933   -4.46595236546167
    ## Hsa-Mir-96-P3_5p                 3344.18853507094   -5.07557947770261
    ##                                               lfcSE               stat
    ##                                           <numeric>          <numeric>
    ## Hsa-Let-7-P1a_5p/P2a1_5p/P2a2_5p   0.11131421654252   1.19461982343763
    ## Hsa-Let-7-P1b_5p                  0.184409102642062  -1.03024741017946
    ## Hsa-Let-7-P1c_5p                  0.232675048404373    13.987286383856
    ## Hsa-Let-7-P2a3_5p                0.0976442188709745 -0.267757085486139
    ## Hsa-Let-7-P2b1_5p                0.0994097267910815   4.53746602977892
    ## ...                                             ...                ...
    ## Hsa-Mir-95-P2_3p                  0.161409886777642    -8.797157721989
    ## Hsa-Mir-95-P3_5p                  0.332664401974054   1.34224159185026
    ## Hsa-Mir-96-P1_5p                  0.288868128705996  -16.6454924674176
    ## Hsa-Mir-96-P2_5p                  0.201189454623787  -22.1766896561084
    ## Hsa-Mir-96-P3_5p                  0.274372198803725  -18.4172382155403
    ##                                                 pvalue                  padj
    ##                                              <numeric>             <numeric>
    ## Hsa-Let-7-P1a_5p/P2a1_5p/P2a2_5p     0.232235600538294     0.301594555061476
    ## Hsa-Let-7-P1b_5p                     0.302893879225991     0.381823880327226
    ## Hsa-Let-7-P1c_5p                  1.86389083593418e-44  3.13619892828925e-43
    ## Hsa-Let-7-P2a3_5p                    0.788886305673403     0.867326705385247
    ## Hsa-Let-7-P2b1_5p                 5.69341978961745e-06   1.4887523368797e-05
    ## ...                                                ...                   ...
    ## Hsa-Mir-95-P2_3p                  1.40325015866435e-18  7.75796873433006e-18
    ## Hsa-Mir-95-P3_5p                     0.179517674775297     0.240392180408443
    ## Hsa-Mir-96-P1_5p                  3.26270971024334e-62  7.89167911165108e-61
    ## Hsa-Mir-96-P2_5p                 5.76691295122432e-109 4.46359062424763e-107
    ## Hsa-Mir-96-P3_5p                  9.55562651928938e-76  2.84463650997307e-74
    ## 
    ## $tissue.type_normal.lung_vs_tumor.colorect
    ## log2 fold change (MAP): tissue.type normal.lung vs tumor.colorect 
    ## Wald test p-value: tissue.type normal.lung vs tumor.colorect 
    ## DataFrame with 389 rows and 6 columns
    ##                                          baseMean     log2FoldChange
    ##                                         <numeric>          <numeric>
    ## Hsa-Let-7-P1a_5p/P2a1_5p/P2a2_5p 89174.7542403281  0.408125846033093
    ## Hsa-Let-7-P1b_5p                 3359.93816553835  0.332599900506941
    ## Hsa-Let-7-P1c_5p                 3992.69006418194     2.020544361604
    ## Hsa-Let-7-P2a3_5p                52000.6252125236  0.129817388788356
    ## Hsa-Let-7-P2b1_5p                9751.60360980849   1.23298675588116
    ## ...                                           ...                ...
    ## Hsa-Mir-95-P2_3p                 390.216333625945 -0.987860472334413
    ## Hsa-Mir-95-P3_5p                 3.95328293251241 -0.191008452152824
    ## Hsa-Mir-96-P1_5p                 202.656949570096   -2.3341261742075
    ## Hsa-Mir-96-P2_5p                 27209.4573430933  -1.90878481829015
    ## Hsa-Mir-96-P3_5p                 3344.18853507094  -1.93676271625054
    ##                                              lfcSE               stat
    ##                                          <numeric>          <numeric>
    ## Hsa-Let-7-P1a_5p/P2a1_5p/P2a2_5p 0.151440987364688   2.68337077689517
    ## Hsa-Let-7-P1b_5p                 0.249992618790543   1.30066611269162
    ## Hsa-Let-7-P1c_5p                 0.314717268191971   6.48763117983787
    ## Hsa-Let-7-P2a3_5p                0.132884439212234  0.963563940574538
    ## Hsa-Let-7-P2b1_5p                0.135208009802693   9.11432613077411
    ## ...                                            ...                ...
    ## Hsa-Mir-95-P2_3p                 0.215830901283507  -4.59605232940141
    ## Hsa-Mir-95-P3_5p                 0.415244927482672 -0.364000003055949
    ## Hsa-Mir-96-P1_5p                 0.350565916495824  -6.75367079090678
    ## Hsa-Mir-96-P2_5p                 0.272289955724401  -7.05756424076837
    ## Hsa-Mir-96-P3_5p                 0.366207615281681  -5.32493458906083
    ##                                                pvalue                 padj
    ##                                             <numeric>            <numeric>
    ## Hsa-Let-7-P1a_5p/P2a1_5p/P2a2_5p  0.00728841357441722   0.0138946603610811
    ## Hsa-Let-7-P1b_5p                    0.193372766419832    0.281397733924404
    ## Hsa-Let-7-P1c_5p                 8.71964116059556e-11 4.62260428650751e-10
    ## Hsa-Let-7-P2a3_5p                   0.335264592721268    0.444340401997023
    ## Hsa-Let-7-P2b1_5p                7.91608806736923e-20 9.88234220023191e-19
    ## ...                                               ...                  ...
    ## Hsa-Mir-95-P2_3p                  4.3057060677156e-06 1.35472215301296e-05
    ## Hsa-Mir-95-P3_5p                    0.715858007098838    0.805340257986192
    ## Hsa-Mir-96-P1_5p                 1.44150633672609e-11 8.45246897443934e-11
    ## Hsa-Mir-96-P2_5p                 1.69446438026879e-12 1.11145375451529e-11
    ## Hsa-Mir-96-P3_5p                 1.00989374164518e-07   3.947766444613e-07
    ## 
    ## $tissue.type_metastasis.liver_vs_tumor.colorect
    ## log2 fold change (MAP): tissue.type metastasis.liver vs tumor.colorect 
    ## Wald test p-value: tissue.type metastasis.liver vs tumor.colorect 
    ## DataFrame with 389 rows and 6 columns
    ##                                          baseMean      log2FoldChange
    ##                                         <numeric>           <numeric>
    ## Hsa-Let-7-P1a_5p/P2a1_5p/P2a2_5p 89174.7542403281 -0.0291372578148132
    ## Hsa-Let-7-P1b_5p                 3359.93816553835 -0.0384102612949004
    ## Hsa-Let-7-P1c_5p                 3992.69006418194   0.165919575225585
    ## Hsa-Let-7-P2a3_5p                52000.6252125236 -0.0667790393731562
    ## Hsa-Let-7-P2b1_5p                9751.60360980849  -0.022055090655077
    ## ...                                           ...                 ...
    ## Hsa-Mir-95-P2_3p                 390.216333625945   0.197986570574988
    ## Hsa-Mir-95-P3_5p                 3.95328293251241    0.87517582162084
    ## Hsa-Mir-96-P1_5p                 202.656949570096  -0.153456519278348
    ## Hsa-Mir-96-P2_5p                 27209.4573430933  0.0668521603991884
    ## Hsa-Mir-96-P3_5p                 3344.18853507094   0.480378406060166
    ##                                               lfcSE               stat
    ##                                           <numeric>          <numeric>
    ## Hsa-Let-7-P1a_5p/P2a1_5p/P2a2_5p  0.087811822186451  -0.34524409879388
    ## Hsa-Let-7-P1b_5p                  0.143352612343951 -0.309464165623311
    ## Hsa-Let-7-P1c_5p                  0.178645206267475   1.08002602868534
    ## Hsa-Let-7-P2a3_5p                0.0771691545429361 -0.884124420143461
    ## Hsa-Let-7-P2b1_5p                0.0785334182092049 -0.274463101822189
    ## ...                                             ...                ...
    ## Hsa-Mir-95-P2_3p                  0.123880822993862   1.53530515265106
    ## Hsa-Mir-95-P3_5p                  0.226325334424462   3.88027612923505
    ## Hsa-Mir-96-P1_5p                  0.196070125088779 -0.994964941991442
    ## Hsa-Mir-96-P2_5p                  0.155589993744074  0.303386567615509
    ## Hsa-Mir-96-P3_5p                  0.205702483629426   2.17135483299889
    ##                                                pvalue                 padj
    ##                                             <numeric>            <numeric>
    ## Hsa-Let-7-P1a_5p/P2a1_5p/P2a2_5p    0.729910868456285    0.886990854282231
    ## Hsa-Let-7-P1b_5p                    0.756968466884636    0.898589554206931
    ## Hsa-Let-7-P1c_5p                    0.280130589584733    0.477579463300844
    ## Hsa-Let-7-P2a3_5p                   0.376629052036082    0.589523977412743
    ## Hsa-Let-7-P2b1_5p                   0.783728755707429     0.90746054828939
    ## ...                                               ...                  ...
    ## Hsa-Mir-95-P2_3p                    0.124708889065693    0.283896118049548
    ## Hsa-Mir-95-P3_5p                 0.000104337942960983 0.000776515075498082
    ## Hsa-Mir-96-P1_5p                     0.31975331550125    0.522757414533685
    ## Hsa-Mir-96-P2_5p                    0.761595281084944    0.898589554206931
    ## Hsa-Mir-96-P3_5p                   0.0299043603592752   0.0972519954541135
    ## 
    ## $tissue.type_metastasis.lung_vs_tumor.colorect
    ## log2 fold change (MAP): tissue.type metastasis.lung vs tumor.colorect 
    ## Wald test p-value: tissue.type metastasis.lung vs tumor.colorect 
    ## DataFrame with 389 rows and 6 columns
    ##                                          baseMean     log2FoldChange
    ##                                         <numeric>          <numeric>
    ## Hsa-Let-7-P1a_5p/P2a1_5p/P2a2_5p 89174.7542403281 -0.371287975266059
    ## Hsa-Let-7-P1b_5p                 3359.93816553835  -0.86478570992646
    ## Hsa-Let-7-P1c_5p                 3992.69006418194  0.135121670461325
    ## Hsa-Let-7-P2a3_5p                52000.6252125236 -0.362077135000897
    ## Hsa-Let-7-P2b1_5p                9751.60360980849 0.0403928595342307
    ## ...                                           ...                ...
    ## Hsa-Mir-95-P2_3p                 390.216333625945  0.244457301347207
    ## Hsa-Mir-95-P3_5p                 3.95328293251241  0.770986123093289
    ## Hsa-Mir-96-P1_5p                 202.656949570096 -0.746540861417062
    ## Hsa-Mir-96-P2_5p                 27209.4573430933  -0.78147625945163
    ## Hsa-Mir-96-P3_5p                 3344.18853507094 -0.326709592212615
    ##                                               lfcSE              stat
    ##                                           <numeric>         <numeric>
    ## Hsa-Let-7-P1a_5p/P2a1_5p/P2a2_5p 0.0959152266130106 -3.87096089411713
    ## Hsa-Let-7-P1b_5p                  0.156536850438983 -5.51447766923507
    ## Hsa-Let-7-P1c_5p                  0.194969507259272 0.835611129624883
    ## Hsa-Let-7-P2a3_5p                0.0842935080279992 -4.30344211748187
    ## Hsa-Let-7-P2b1_5p                0.0857694713451329 0.474626440937984
    ## ...                                             ...               ...
    ## Hsa-Mir-95-P2_3p                  0.135123404458061  1.75056733087186
    ## Hsa-Mir-95-P3_5p                   0.24193629758336  3.22815141458786
    ## Hsa-Mir-96-P1_5p                  0.213973544808363 -3.63856405996494
    ## Hsa-Mir-96-P2_5p                  0.169868922201866 -4.66162915975394
    ## Hsa-Mir-96-P3_5p                  0.224439059665663 -1.53101854305705
    ##                                                pvalue                 padj
    ##                                             <numeric>            <numeric>
    ## Hsa-Let-7-P1a_5p/P2a1_5p/P2a2_5p   0.0001084071837261  0.00040731631166991
    ## Hsa-Let-7-P1b_5p                 3.49817317913999e-08  2.5543264534475e-07
    ## Hsa-Let-7-P1c_5p                    0.403373705499327    0.532783699755084
    ## Hsa-Let-7-P2a3_5p                1.68164797190917e-05 7.65644429563353e-05
    ## Hsa-Let-7-P2b1_5p                   0.635053256632943    0.757086098142951
    ## ...                                               ...                  ...
    ## Hsa-Mir-95-P2_3p                   0.0800204667405637     0.14746628870761
    ## Hsa-Mir-95-P3_5p                  0.00124593006501447  0.00359832041164626
    ## Hsa-Mir-96-P1_5p                 0.000274162436246517  0.00098241539655002
    ## Hsa-Mir-96-P2_5p                 3.13716090303122e-06 1.57672892139361e-05
    ## Hsa-Mir-96-P3_5p                    0.125764809854213    0.211612962667741
    ## 
    ## $tissue.type_metastasis.pc_vs_tumor.colorect
    ## log2 fold change (MAP): tissue.type metastasis.pc vs tumor.colorect 
    ## Wald test p-value: tissue.type metastasis.pc vs tumor.colorect 
    ## DataFrame with 389 rows and 6 columns
    ##                                          baseMean      log2FoldChange
    ##                                         <numeric>           <numeric>
    ## Hsa-Let-7-P1a_5p/P2a1_5p/P2a2_5p 89174.7542403281  -0.194926960102631
    ## Hsa-Let-7-P1b_5p                 3359.93816553835  -0.170006997664929
    ## Hsa-Let-7-P1c_5p                 3992.69006418194   0.950788049188073
    ## Hsa-Let-7-P2a3_5p                52000.6252125236  -0.235739655653004
    ## Hsa-Let-7-P2b1_5p                9751.60360980849 -0.0275610759199226
    ## ...                                           ...                 ...
    ## Hsa-Mir-95-P2_3p                 390.216333625945  -0.391736823339089
    ## Hsa-Mir-95-P3_5p                 3.95328293251241   0.193480385452117
    ## Hsa-Mir-96-P1_5p                 202.656949570096  -0.537251798926674
    ## Hsa-Mir-96-P2_5p                 27209.4573430933  -0.144468039225267
    ## Hsa-Mir-96-P3_5p                 3344.18853507094  -0.149480720489416
    ##                                               lfcSE               stat
    ##                                           <numeric>          <numeric>
    ## Hsa-Let-7-P1a_5p/P2a1_5p/P2a2_5p 0.0919497346955033  -2.12106919715872
    ## Hsa-Let-7-P1b_5p                  0.146429745186075  -1.18446872666545
    ## Hsa-Let-7-P1c_5p                  0.178397588170818   5.34338828103208
    ## Hsa-Let-7-P2a3_5p                  0.08107727951498  -2.91425069775732
    ## Hsa-Let-7-P2b1_5p                0.0824764385946126  -0.32740564020296
    ## ...                                             ...                ...
    ## Hsa-Mir-95-P2_3p                  0.128019210385648  -3.07056591636714
    ## Hsa-Mir-95-P3_5p                  0.222349699850271  0.967689700538171
    ## Hsa-Mir-96-P1_5p                  0.193381569174697   -2.9161157752392
    ## Hsa-Mir-96-P2_5p                  0.157787812752101  -1.01238717153534
    ## Hsa-Mir-96-P3_5p                  0.201177138993338 -0.825123443737066
    ##                                                pvalue                 padj
    ##                                             <numeric>            <numeric>
    ## Hsa-Let-7-P1a_5p/P2a1_5p/P2a2_5p   0.0339159796632756   0.0901196031052751
    ## Hsa-Let-7-P1b_5p                    0.236227568744681    0.370368857833221
    ## Hsa-Let-7-P1c_5p                 9.12250626572605e-08 1.40060447984391e-06
    ## Hsa-Let-7-P2a3_5p                 0.00356543456685069   0.0141511415829511
    ## Hsa-Let-7-P2b1_5p                   0.743361101075885    0.848252544785979
    ## ...                                               ...                  ...
    ## Hsa-Mir-95-P2_3p                  0.00213653519975517    0.009575796316975
    ## Hsa-Mir-95-P3_5p                    0.333199363429196     0.48229635484693
    ## Hsa-Mir-96-P1_5p                  0.00354418957989829   0.0141511415829511
    ## Hsa-Mir-96-P2_5p                    0.311352969649893     0.46861091937278
    ## Hsa-Mir-96-P3_5p                    0.409301511166841    0.561003252050181
    ## 
    ## $type_metastasis_vs_tumor
    ## log2 fold change (MAP): type metastasis vs tumor 
    ## Wald test p-value: type metastasis vs tumor 
    ## DataFrame with 389 rows and 6 columns
    ##                                          baseMean       log2FoldChange
    ##                                         <numeric>            <numeric>
    ## Hsa-Let-7-P1a_5p/P2a1_5p/P2a2_5p 89174.7542403281   -0.180701783114052
    ## Hsa-Let-7-P1b_5p                 3359.93816553835   -0.294623990307531
    ## Hsa-Let-7-P1c_5p                 3992.69006418194     0.48926666782258
    ## Hsa-Let-7-P2a3_5p                52000.6252125236   -0.206937259890988
    ## Hsa-Let-7-P2b1_5p                9751.60360980849 -0.00482780603867569
    ## ...                                           ...                  ...
    ## Hsa-Mir-95-P2_3p                 390.216333625945   0.0426916660211096
    ## Hsa-Mir-95-P3_5p                 3.95328293251241    0.678957110901101
    ## Hsa-Mir-96-P1_5p                 202.656949570096   -0.464482054965902
    ## Hsa-Mir-96-P2_5p                 27209.4573430933   -0.221824399033929
    ## Hsa-Mir-96-P3_5p                 3344.18853507094   0.0672687256455847
    ##                                               lfcSE               stat
    ##                                           <numeric>          <numeric>
    ## Hsa-Let-7-P1a_5p/P2a1_5p/P2a2_5p 0.0639958634974693  -2.81887235438572
    ## Hsa-Let-7-P1b_5p                  0.106490720796663  -2.76160115635425
    ## Hsa-Let-7-P1c_5p                  0.149925843445819   3.33812788609782
    ## Hsa-Let-7-P2a3_5p                0.0563710904824311  -3.67052297626806
    ## Hsa-Let-7-P2b1_5p                0.0617282834747614 -0.070874324259544
    ## ...                                             ...                ...
    ## Hsa-Mir-95-P2_3p                 0.0918217873522948  0.437722915855047
    ## Hsa-Mir-95-P3_5p                  0.170411797990612   3.95560185906774
    ## Hsa-Mir-96-P1_5p                  0.152114644217898  -3.14528253438741
    ## Hsa-Mir-96-P2_5p                   0.12892939045268  -1.79236001772164
    ## Hsa-Mir-96-P3_5p                  0.164329004650402  0.314884459867911
    ##                                                pvalue                 padj
    ##                                             <numeric>            <numeric>
    ## Hsa-Let-7-P1a_5p/P2a1_5p/P2a2_5p  0.00481926786262815   0.0172690431744175
    ## Hsa-Let-7-P1b_5p                  0.00575186955495521   0.0200538154753844
    ## Hsa-Let-7-P1c_5p                  0.00084344919760993  0.00370925953948912
    ## Hsa-Let-7-P2a3_5p                0.000242054709304128  0.00121656068182724
    ## Hsa-Let-7-P2b1_5p                   0.943497778247006     0.97368970715091
    ## ...                                               ...                  ...
    ## Hsa-Mir-95-P2_3p                    0.661587155561431    0.823261187145574
    ## Hsa-Mir-95-P3_5p                 7.63422117073344e-05 0.000440961730309529
    ## Hsa-Mir-96-P1_5p                  0.00165926501866658  0.00690468346477385
    ## Hsa-Mir-96-P2_5p                   0.0730753150258787     0.17139482978797
    ## Hsa-Mir-96-P3_5p                    0.752849381089381     0.86454810231926
