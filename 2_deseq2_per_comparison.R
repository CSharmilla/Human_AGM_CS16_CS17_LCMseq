library(here)
library(patchwork)
library(DESeq2)
library(ashr)
library(sva)
library(edgeR)
library(ggrepel)
library(pheatmap)


txi_cs <- readRDS(here("output", "txi_cs"))
meta_cs <- readRDS(here("output", "meta_cs"))
genes_list <- readRDS(here("output", "genes_list"))
subsets <- readRDS(here("output", "subsets"))
txi_rename <- readRDS(here("output", "txi_rename"))

# removing embryo 2
#for (i in c("V_D"     ,    "V_VL"    ,    "VL_DL"    ,   "DL_D"   ,     "D_DL")){
#  subsets[[i]] <- subsets[[i]][!subsets[[i]]$embryo == 2,]
#}

## pre-processing
rowdata <- list()
txi_subsets <- list()

for (i in names(subsets)){
  
  # subsetting txi per comparison
  txi_subsets[[i]] <- lapply(txi_rename[1:3], function(a){
    a <- a[,rownames(subsets[[i]])]
    return(a)
  })
  
  txi_subsets[[i]]$countsFromAbundance <- "no"
  
  # coldata rownames and txi column names matching
  subsets[[i]] <- subsets[[i]][match(colnames(txi_subsets[[i]]$counts), rownames(subsets[[i]])), ] 
  
  # check if coldata rownames and txi column names match and are in order
  print(all(rownames(subsets[[i]]) == colnames(txi_subsets[[i]]$counts)))
  print(all(rownames(subsets[[i]]) %in% colnames(txi_subsets[[i]]$counts)))
  
  # getting gene list - make sure it is in the same order as txi rows
  rowdata[[i]] <- genes_list[which(genes_list$ensembl_gene_id %in% rownames(txi_subsets[[i]]$counts)),]
  rowdata[[i]] <- rowdata[[i]][match(rownames(txi_subsets[[i]]$counts), rowdata[[i]]$ensembl_gene_id),] 
  rownames(rowdata[[i]]) <- 1:nrow(rowdata[[i]]) 
  
  # check if rowdata and txi rownames  match and are in order
  print(all(rowdata[[i]]$ensembl_gene_id == rownames(txi_subsets[[i]]$counts)))
  print(all(rowdata[[i]]$ensembl_gene_id %in% rownames(txi_subsets[[i]]$counts)))
  
  # factor levels
  subsets[[i]]$embryo <- factor(subsets[[i]]$embryo)
  subsets[[i]]$subdomain <- factor(subsets[[i]]$subdomain)
  subsets[[i]]$title <- factor(subsets[[i]]$title)
}


# creating DESeqDataSet objects

if (!(file.exists(here("output", "dds")))) {
  dds <- list()
  dds_unfil <- list()
  for (i in names(txi_subsets)){
    # converting to deseq2 object
    dds[[i]] <- dds_unfil[[i]] <- DESeqDataSetFromTximport(txi_subsets[[i]], colData = subsets[[i]], 
                                         design = ~embryo + subdomain, 
                                         rowData = rowdata[[i]])
    
    # changing this to make merge easier
    colnames(rowData(dds[[i]]))[1] <- colnames(rowData(dds_unfil[[i]]))[1] <- "row"
    
    # basic pre-filtering
    smallestGroupSize <- 3
    #keep <- rowSums(counts(dds[[i]]) >= 10) >= smallestGroupSize
    keep <- filterByExpr(dds[[i]], group = dds[[i]]$subdomain)
    dds[[i]] <- dds[[i]][keep,]
    
    # estimate size factors - another assay added with normalisation factors
    dds[[i]] <- estimateSizeFactors(dds[[i]])
    dds_unfil[[i]] <- estimateSizeFactors(dds_unfil[[i]])
    
    
  }
  saveRDS(dds, file = here("output", "dds"))
  saveRDS(dds_unfil, file = here("output", "dds_unfil"))
} else{
  dds <- readRDS(here("output", "dds"))
  dds_unfil <- readRDS(here("output", "dds_unfil"))
}




# transforming and plotting to see the data
vst_blind <- list()
pca_plot <- list()
for (i in names(dds)){
  #i <- "V_D"
  # blind=TRUE should be used for comparing samples in a manner unbiased by prior 
  # information on samples, for example to perform sample QA (quality assurance)
  vst_blind[[i]] <- vst(dds[[i]], blind=TRUE)
  
  # plotting pca
  plot <- plotPCA(vst_blind[[i]], intgroup=c("title"), returnData = TRUE) #+ geom_text(label = vst_blind[[i]]$title) 
  
  p <- ggplot(plot, aes(x = PC1, y = PC2, color = title)) + geom_point() +
  geom_text_repel(aes(label = rownames(plot)), max.overlaps = Inf) + ggtitle(paste0(i))
  #pca_plot[[i]] <- p
  print(p)
}
patchwork::wrap_plots(pca_plot)


# running deseq

if (!(file.exists(here("output", "deseq")))) {
  deseq <- list()
  for (i in names(dds)){
    deseq[[i]] <- DESeq(dds[[i]])
  }
  saveRDS(deseq, file = here("output", "deseq"))
} else{
  deseq <- readRDS(here("output", "deseq"))
}



# separate loop since deseq takes longer
for (i in names(deseq)){
  print(resultsNames(deseq[[i]]))
}


# likehood ratio test is better for time series data

# the name of a factor in the design formula, 
# the name of the numerator level for the fold change, and the name of the denominator level for the fold change
#  log2 fold change of 1.5 for a specific gene in the “WT vs KO comparison” means that the 
# expression of that gene is increased in WT relative to KO by a multiplicative factor of 2^1.5 ≈ 2.82

# getting results and transforming
# creating contrasts
contrasts <- list(V_D = c("subdomain", "Ventral_V", "Dorsal_D"), 
                  V_VL = c("subdomain", "Ventral_V", "Ventro-lateral_VL"),
                  D_DL = c("subdomain", "Dorsal_D", "Dorsal-lateral_DL"),
                  DL_D = c("subdomain", "Dorsal-lateral_DL", "Dorsal_D"),
                  VL_DL = c("subdomain", "Ventro-lateral_VL", "Dorsal-lateral_DL"),
                  inner_V_D = c("subdomain", "Ventral_Inner_V_Inner", "Dorsal_Inner_D_Inner"),
                  mid_V_D = c("subdomain", "Ventral_Mid_V_Mid", "Dorsal_Mid_D_Mid"),
                  outer_V_D = c("subdomain", "Ventral_Outer_V_Outer", "Dorsal_Outer_D_Outer"),
                  inner_mid_V = c("subdomain","Ventral_Inner_V_Inner", "Ventral_Mid_V_Mid"),
                  mid_outer_V = c("subdomain","Ventral_Mid_V_Mid", "Ventral_Outer_V_Outer"),
                  inner_mid_D = c("subdomain","Dorsal_Inner_D_Inner", "Dorsal_Mid_D_Mid"),
                  mid_outer_D = c("subdomain","Dorsal_Mid_D_Mid", "Dorsal_Outer_D_Outer"))

results_list <- list() # dataframes with shrunken lfc
results_wald <- list() # wald test dataframes
reslfc <- list() # list of results for plots

if (!file.exists(here("output", "results_list"))){
  for (i in names(deseq)){
    
    #i <- "CS16"
    name1 <- paste0(i, "_apeglm_results")
    name2 <- paste0(i, "_results")
    
    res <- results(deseq[[i]], contrast = contrasts[[i]], alpha = 0.05) # cant use tidy = TRUE here since using lfc shrink
    
    # LFC shrinkage
    reslfc[[i]] <- lfcShrink(deseq[[i]], res = res, type = "ashr")
    apeglm <- data.frame(reslfc[[i]])
    apeglm$row <- rownames(apeglm)
    apeglm <- as.data.frame(merge(apeglm, rowData(dds[[i]]), all.x=TRUE))
    
    assign(name1, apeglm)
    results_list[[i]] <- apeglm
    
    results_wald[[i]] <- results(deseq[[i]], contrast = contrasts[[i]], alpha = 0.05, tidy = TRUE)
    results_wald[[i]] <- merge(results_wald[[i]], rowData(dds[[i]]), all.x=TRUE)
    
    assign(name2, data.frame(results_wald[[i]]))
      
    
  }
  saveRDS(results_list, file = here("output", "results_list"))
  saveRDS(results_wald, file = here("output", "results_wald"))
  saveRDS(reslfc, file = here("output", "reslfc"))
  
}else{
  results_list <- readRDS(here("output", "results_list"))
  results_wald <- readRDS(here("output", "results_wald"))
  reslfc <- readRDS(here("output", "reslfc"))
}

# plots

# transforming and plotting PCA to see the data - not BLIND
vst <- list()
pca_plot_notblind <- list()
for (i in names(deseq)){
  # blind=FALSE should be used for transforming data for downstream analysis, where the full use of the design information should be made. 
  # blind=FALSE will skip re-estimation of the dispersion trend, if this has already been calculated. If many of genes have large 
  # differences in counts due to the experimental design, it is important to set blind=FALSE for downstream analysis
  vst[[i]] <- vst(deseq[[i]], blind=FALSE)
  
  # plotting pca
  pca_plot_notblind[[i]] <- plotPCA(vst[[i]], intgroup=c("title")) #+ geom_text(aes(label = vst[[i]]$title))
  #print(pca_plot_notblind[[i]])
}
patchwork::wrap_plots(pca_plot_notblind)

for (i in names(vst)){
  #i <- "V_D"
  select <- order(results_list[[i]]$padj, decreasing=FALSE)[1:100]
  df <- as.data.frame(colData(vst[[i]])[,c("embryo","subdomain")])
  print(pheatmap(assay(vst[[i]])[select,], cluster_rows=TRUE, show_rownames=FALSE,
                 cluster_cols=TRUE, annotation_col=df, main = i))
}


# # filtering
# results_list_fil <- list()
# results_wald_fil <- list()
# for (i in names(results_list)){
#   for (j in names(results_list[[i]])){
#     results_list_fil[[i]][[j]] <- results_list[[i]][[j]][which(results_list[[i]][[j]]$padj < 0.05),]
#     results_wald_fil[[i]][[j]] <- results_wald[[i]][[j]][which(results_wald[[i]][[j]]$padj < 0.05),]
#     
#   }
# }


# MA plots
# MA plots display a log ratio (M) vs an average (A) in order to visualize the differences between two groups. 
# In general we would expect the expression of genes to remain consistent between conditions and so the MA 
# plot should be similar to the shape of a trumpet with most points residing on a y intercept of 0

for (i in names(reslfc)){
  DESeq2::plotMA(reslfc[[i]], main = i)
}

#subset_V_D <- subsets[["V_D"]]
# plotting counts to look see if genes look ok
# checking KITLG
plotCounts(dds[["VL_DL"]], gene="ENSG00000107779", intgroup="title", normalized = TRUE)


plotDispEsts(deseq[["mid_V_D"]])

boxplot(log10(assays(deseq[["CS16"]])[["cooks"]]), range=0, las=2)
a <- txi_cs$CS17$counts

# plotting kitlg
plot_list <- list()
#dev.new()
for (i in names(dds_unfil)){
  gene <- c("ENSG00000049130")
  #if (gene %in% rownames(dds_unfil[[i]])){
    #i <- "V_D"
  plot <- plotCounts(dds_unfil[[i]], gene=gene, intgroup=c("embryo", "subdomain"), normalized = TRUE, , returnData=TRUE) #main = paste0("KITLG expression in ", i))
  
  print(ggplot(plot, aes(x = subdomain, y = count, color = embryo)) +  geom_point(position=position_jitter(w = 0.1,h = 0)) +
                    geom_text_repel(aes(label = rownames(plot)), max.overlaps = 5) + theme_bw() + ggtitle(paste0("KITLG expression in ", i)) + 
                    theme(plot.title = element_text(hjust = 0.5)))
  #}
}
#patchwork::wrap_plots(plot_list)


plotCounts(dds_unfil[["CS17"]], gene="ENSG00000120253", intgroup=c("embryo", "subdomain"), normalized = TRUE)


