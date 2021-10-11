# import libraries
library("tidyverse")
library("DESeq2")
library("knitr")
library('kableExtra')
library("RColorBrewer")
library("pheatmap")
library("dplyr")
library("BiocParallel")
library("broman")
library("gridExtra")
library("ggpubr")
library("FactoMineR")
library("factoextra")



############################### DeseqObject ############################### 

DeseqObject <- function(DESIGN, column, countdata, coldata, consensus="None", sample_type="None", ref) {
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
  
  dds[[column]] <- relevel(dds[[column]], ref=ref)
  dds[[column]] <- droplevels(dds[[column]])
  
  #dds$tissue.type <- relevel(dds$tissue.type, ref=ref)
  #dds$tissue.type <- droplevels(dds$tissue.type)
  
  
  dds <- DESeq(dds,
               parallel=TRUE,
               BPPARAM=MulticoreParam(3)
  )
  
  return(dds)
}






############################### DeseqResult ############################### 


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






############################### SigList ############################### 

SigList <- function(res, dds, column, tissue_type_A, tissue_type_B, coef,
                    norm_adj_up, norm_adj_down, 
                    pCRC_adj_up, pCRC_adj_down){
  "
  Function to create annotated lists of signature miRNA
  Return will print upregulated or downregulated miRNA, 
  by printing <signature_list>$up_mirna
  or          <signature_list>$down_mirna
  "
  group_A_rpm <- rowMeans(res$rpm[res$sig, dds[[column]] == tissue_type_A])
  group_A_rpm_std <- rowSds(res$rpm[res$sig, dds[[column]] == tissue_type_A])
  group_B_rpm <- rowMeans(res$rpm[res$sig, dds[[column]] == tissue_type_B])
  group_B_rpm_std <- rowSds(res$rpm[res$sig, dds[[column]] == tissue_type_B])
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





############################### VolcanoPlot ############################### 
"
Function to create volcano plots for DESeq2 results,
annotating signature miRNA, cell marker miRNA and 
normal adjacent tissue effect
"
non_sig_col <- brocolors("crayons")["Black"]
sig_col     <- "#47FF63"#brocolors("crayons")["Asparagus"]
cell_col    <- "#6347FF"#brocolors("crayons")["Navy Blue"]
normal_col  <- "#FF6347"#brocolors("crayons")["Orange Red"]
text_col    <- brocolors("crayons")["Black"]
dot_size    <- 0.4
hits_size   <- 0.8
sig_dot_size<- 1.2
transparent <- 0.4

VolcanoPlot <- function(res, coef, sig, column,
                        up_mirna, down_mirna,
                        norm_adj_up, norm_adj_down,
                        pCRC_adj_up, pCRC_adj_down){
  # Plot log2foldchange againt -log10(padj) for all miRNA
  xmin <- min(res$log2FoldChange, na.rm = TRUE) - 2 
  xmax <- max(res$log2FoldChange, na.rm = TRUE) + 2
  ymin <- 0
  ymax <- max(-log(res$padj, 10), na.rm = TRUE) + 4
  
  plot(res$log2FoldChange, -log(res$padj, 10), ylab="-log10(Adjusted P)", xlab="Log2 FoldChange",
       main=coef, pch=19, cex=dot_size, col=alpha(non_sig_col, transparent), bty='n',
       xlim=c( xmin , xmax ),
       ylim=c( ymin , ymax )
  )
  
  
  # color significant miRNA
  points(res[sig, ]$log2FoldChange, -log(res[sig, ]$padj, 10), 
         pch=21, cex=hits_size, col='black', bg=sig_col)
  
  # color normal adjacent marker miRNA red
  points(res[intersect(norm_adj_up, sig), ]$log2FoldChange, 
         -log(res[intersect(norm_adj_up, sig), ]$padj, 10), 
         pch=19, cex=hits_size, col=normal_col)
  points(res[intersect(norm_adj_down, sig), ]$log2FoldChange, 
         -log(res[intersect(norm_adj_down, sig), ]$padj, 10), 
         pch=19, cex=hits_size, col=normal_col)
  points(res[intersect(pCRC_adj_up, sig), ]$log2FoldChange, 
         -log(res[intersect(pCRC_adj_up, sig), ]$padj, 10), 
         pch=19, cex=hits_size, col=normal_col)
  points(res[intersect(pCRC_adj_down, sig), ]$log2FoldChange, 
         -log(res[intersect(pCRC_adj_down, sig), ]$padj, 10), 
         pch=19, cex=hits_size, col=normal_col)
  
  # color cell marker miRNA
  points(res[intersect(sig, names(cell_spec_dict_inv)), ]$log2FoldChange, 
         -log(res[intersect(sig, names(cell_spec_dict_inv)), ]$padj, 10), 
         pch=19, cex=hits_size, col=cell_col)
  
  
  # Plot Text for signature miRNA, remove upregulated miRNA in 
  # normal colon vs normal adjacent metastatic site
  if (norm_adj_up != "None" & length(up_mirna > 0)) {
    points(res[setdiff(up_mirna, norm_adj_up), ]$log2FoldChange, 
           -log(res[setdiff(up_mirna, norm_adj_up), ]$padj, 10), 
           pch=21, cex=sig_dot_size, col='black', bg=sig_col)
    
    text(res[setdiff(up_mirna, norm_adj_up), ]$log2FoldChange, 
         -log(res[setdiff(up_mirna, norm_adj_up), ]$padj, 10), 
         labels=rownames(res[setdiff(up_mirna, norm_adj_up), ]),
         pos=3, cex=0.7, col=text_col)
  }
  
  # Plot Text for signature miRNA, remove downregulated miRNA in 
  # normal colon vs normal adjacent metastatic site
  if (norm_adj_down != "None" & length(down_mirna > 0)) {
    points(res[setdiff(down_mirna, norm_adj_down), ]$log2FoldChange, 
           -log(res[setdiff(down_mirna, norm_adj_down), ]$padj, 10),
           pch=21, cex=sig_dot_size, col='black', bg=sig_col)
    
    text(res[setdiff(down_mirna, norm_adj_down), ]$log2FoldChange, 
         -log(res[setdiff(down_mirna, norm_adj_down), ]$padj, 10),
         labels=rownames(res[setdiff(down_mirna, norm_adj_down), ]),
         pos=3, cex=0.7, col=text_col)
  }
  
  # Plot Text for signature miRNA, remove upregulated miRNA in 
  # pCRC vs normal adjacent metastatic site
  if (pCRC_adj_up != "None" & length(up_mirna > 0)) {
    points(res[setdiff(up_mirna, pCRC_adj_up), ]$log2FoldChange, 
           -log(res[setdiff(up_mirna, pCRC_adj_up), ]$padj, 10), 
           pch=21, cex=sig_dot_size, col='black', bg=sig_col)
    
    text(res[setdiff(up_mirna, pCRC_adj_up), ]$log2FoldChange, 
         -log(res[setdiff(up_mirna, pCRC_adj_up), ]$padj, 10), 
         labels=rownames(res[setdiff(up_mirna, pCRC_adj_up), ]),
         pos=3, cex=0.7, col=text_col)
  }
  
  # Plot Text for signature miRNA, remove downregulated miRNA in 
  # pCRC vs normal adjacent metastatic site
  if (pCRC_adj_down != "None" & length(down_mirna > 0)) {
    points(res[setdiff(down_mirna, pCRC_adj_down), ]$log2FoldChange, 
           -log(res[setdiff(down_mirna, pCRC_adj_down), ]$padj, 10),
           pch=21, cex=sig_dot_size, col='black', bg=sig_col)
    
    text(res[setdiff(down_mirna, pCRC_adj_down), ]$log2FoldChange, 
         -log(res[setdiff(down_mirna, pCRC_adj_down), ]$padj, 10),
         labels=rownames(res[setdiff(down_mirna, pCRC_adj_down), ]),
         pos=3, cex=0.7, col=text_col)
  }
  
  # Add line marking significance threshold
  abline(h=-log(p.Threshold, 10), lty=3)
  # Add lines marking upper and lower log2foldchange threshold
  abline(v=-lfc.Threshold, lty=3)
  abline(v=lfc.Threshold, lty=3)
  legend("topleft", c("Signature miRNA", "Cell Markers", "Normal Background", "Non-Signature"), 
         fill=c(sig_col, cell_col, normal_col, non_sig_col), bg="transparent", bty="n")
  
}







############################### ExpressionPlot ############################### 
"
Function to plot absolute expression levels against fold changes
"
non_sig_col <- brocolors("crayons")["Black"]
sig_col     <- "#47FF63"#brocolors("crayons")["Asparagus"]
cell_col    <- "#6347FF"#brocolors("crayons")["Navy Blue"]
normal_col  <- "#FF6347"#brocolors("crayons")["Orange Red"]
text_col    <- brocolors("crayons")["Black"]
dot_size    <- 0.4
hits_size   <- 0.8
sig_dot_size<- 1.2
transparent <- 0.4
logscale    <- 10

ExpressionPlot <- function(res, rpm, coef, sig, column,
                           tissue_type_A, tissue_type_B,
                           up_mirna, down_mirna,
                           norm_adj_up, norm_adj_down,
                           pCRC_adj_up, pCRC_adj_down){
  plot(res$log2FoldChange, 
       log(rowMeans(rpm[, dds[[column]] == tissue_type_A]), logscale), 
       xlab="Log2 FoldChange", ylab="Log10 RPM",
       main=coef, pch=19, cex=dot_size, col=alpha(non_sig_col, transparent), bty='n'
  )
  points(res$log2FoldChange, 
         log(rowMeans(rpm[, dds[[column]] == tissue_type_B]), logscale), 
         pch=4, cex=dot_size, col=alpha(non_sig_col, transparent))
  
  # color significant miRNA
  points(res[sig, ]$log2FoldChange, 
         log(rowMeans(rpm[sig, dds$tissue.type == tissue_type_A]), logscale), 
         pch=21, cex=hits_size, col='black', bg=sig_col)
  
  # color normal adjacent marker miRNA red
  points(res[intersect(norm_adj_up, sig), ]$log2FoldChange, 
         log(rowMeans(rpm[intersect(norm_adj_up, sig), dds$tissue.type == tissue_type_A]), logscale), 
         pch=19, cex=hits_size, col=normal_col)
  points(res[intersect(norm_adj_down, sig), ]$log2FoldChange, 
         log(rowMeans(rpm[intersect(norm_adj_down, sig), dds$tissue.type == tissue_type_A]), logscale), 
         pch=19, cex=hits_size, col=normal_col)
  points(res[intersect(pCRC_adj_up, sig), ]$log2FoldChange, 
         log(rowMeans(rpm[intersect(pCRC_adj_up, sig), dds$tissue.type == tissue_type_A]), logscale), 
         pch=19, cex=hits_size, col=normal_col)
  points(res[intersect(pCRC_adj_down, sig), ]$log2FoldChange, 
         log(rowMeans(rpm[intersect(pCRC_adj_down, sig), dds$tissue.type == tissue_type_A]), logscale), 
         pch=19, cex=hits_size, col=normal_col)
  
  # color cell marker miRNA
  points(res[intersect(sig, names(cell_spec_dict_inv)), ]$log2FoldChange,
         log(rowMeans(rpm[intersect(sig, names(cell_spec_dict_inv)), dds$tissue.type == tissue_type_A]), logscale),
         pch=19, cex=hits_size, col=cell_col)
  
  # Plot Text for signature miRNA, remove upregulated miRNA in 
  # normal colon vs normal adjacent metastatic site
  if (norm_adj_up != "None" & length(up_mirna > 0)) {
    points(res[setdiff(up_mirna, norm_adj_up), ]$log2FoldChange, 
           log(rowMeans(rpm[setdiff(up_mirna, norm_adj_up), dds$tissue.type == tissue_type_A]), logscale),
           pch=21, cex=sig_dot_size, col="black", bg=sig_col)
    
    text(res[setdiff(up_mirna, norm_adj_up), ]$log2FoldChange, 
         log(rowMeans(rpm[setdiff(up_mirna, norm_adj_up), dds$tissue.type == tissue_type_A]), logscale), 
         labels=rownames(res[setdiff(up_mirna, norm_adj_up), ]),
         pos=3, cex=0.7, col=text_col)
  }
  # Plot Text for signature miRNA, remove upregulated miRNA in 
  # pCRC vs normal adjacent metastatic site
  if (pCRC_adj_up != "None" & length(up_mirna > 0)) {
    points(res[setdiff(up_mirna, pCRC_adj_up), ]$log2FoldChange, 
           log(rowMeans(rpm[setdiff(up_mirna, pCRC_adj_up), dds$tissue.type == tissue_type_A]), logscale), 
           pch=21, cex=sig_dot_size, col="black", bg=sig_col)
    
    text(res[setdiff(up_mirna, pCRC_adj_up), ]$log2FoldChange, 
         log(rowMeans(rpm[setdiff(up_mirna, pCRC_adj_up), dds$tissue.type == tissue_type_A]), logscale), 
         labels=rownames(res[setdiff(up_mirna, pCRC_adj_up), ]),
         pos=3, cex=0.7, col=text_col)
  }
  abline(v=-lfc.Threshold, lty=3)
  abline(v=lfc.Threshold, lty=3)
  abline(h=log(rpm.Threshold, 10), lty=3)
}



GetResults <- function(dds, column, coef,
                       tissue_type_A, tissue_type_B, 
                       norm_adj_up, norm_adj_down,
                       pCRC_adj_up, pCRC_adj_down,
                       lfc.Threshold, rpm.Threshold,
                       dict_sig_mirna){

    
    res <- DeseqResult(dds, column, coef, tissue_type_A, tissue_type_B,
                       lfc.Threshold, rpm.Threshold,
                       norm_adj_up, norm_adj_down,
                       pCRC_adj_up, pCRC_adj_down)
    VolcanoPlot(res$res, coef, res$sig, column,
                res$up_mirna, res$down_mirna,
                norm_adj_up, norm_adj_down,
                pCRC_adj_up, pCRC_adj_down)
    signature_mirnas <- SigList(res, dds, column, tissue_type_A, tissue_type_B, coef,
                                norm_adj_up, norm_adj_down,
                                pCRC_adj_up, pCRC_adj_down)
    return(list('res'=res, 'signature_mirnas'=signature_mirnas))
    #knitr::kable(signature_mirnas$up_mirna)
    #knitr::kable(signature_mirnas$number_upregulated)
    #knitr::kable(signature_mirnas$down_mirna)
    #knitr::kable(signature_mirnas$number_downregulated)

}

"
Functions to adjust LFC and FDR depending on expression in normal tissue 
"
SubtractLFC <- function(x, y){
  "
  LFC is reduced to zero by value in control group if value in control 
  group has the same sign.
  "
  z = x
  if ( is.na(x) | is.na(y) ){ return( z ) }
  else if (sign(x) == sign(y)){ z = x - y }
  if (sign(z) != sign(x)) { z = 0 }
  return( z )
}

SubtractAdjP <- function(x , y, xP, yP){
  "
  adjP value is increase to 1 by value in control group if 
  LFC in control group has the same sign.
  "
  z = xP
  if ( is.na(xP) | is.na(yP) ){ return( z ) }
  if ( sign(x) == sign(y) ){ 
    z = (xP + ( 1 - yP )) }
  if (z > 1) {z = 1}
  return(z)
}












































