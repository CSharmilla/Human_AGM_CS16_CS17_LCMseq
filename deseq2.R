## not using this script anymore


library(here)
library(patchwork)
library(DESeq2)
library(ashr)
library(sva)
library(edgeR)


txi_cs <- readRDS(here("output", "txi_cs"))
meta_cs <- readRDS(here("output", "meta_cs"))
genes_list <- readRDS(here("output", "genes_list"))
subsets <- readRDS(here("output", "subsets"))


## pre-processing
rowdata <- list()

for (i in names(txi_cs)){
  
  
  # coldata rownames and txi column names matching
  meta_cs[[i]] <- meta_cs[[i]][match(colnames(txi_cs[[i]]$counts), rownames(meta_cs[[i]])), ] 
  
  # check if coldata rownames and txi column names match and are in order
  print(all(rownames(meta_cs[[i]]) == colnames(txi_cs[[i]]$counts)))
  print(all(rownames(meta_cs[[i]]) %in% colnames(txi_cs[[i]]$counts)))
  
  # getting gene list - make sure it is in the same order as txi rows
  rowdata[[i]] <- genes_list[which(genes_list$ensembl_gene_id %in% rownames(txi_cs[[i]]$counts)),]
  rowdata[[i]] <- rowdata[[i]][match(rownames(txi_cs[[i]]$counts), rowdata[[i]]$ensembl_gene_id, ),] 
  rownames(rowdata[[i]]) <- 1:nrow(rowdata[[i]]) 
  
  # check if rowdata and txi rownames  match and are in order
  print(all(rowdata[[i]]$ensembl_gene_id == rownames(txi_cs[[i]]$counts)))
  print(all(rowdata[[i]]$ensembl_gene_id %in% rownames(txi_cs[[i]]$counts)))
  
  # factor levels
  meta_cs[[i]]$embryo <- factor(meta_cs[[i]]$embryo)
  meta_cs[[i]]$subdomain <- factor(meta_cs[[i]]$subdomain)
  meta_cs[[i]]$title <- factor(meta_cs[[i]]$title)
}


# creating DESeqDataSet objects

if (!(file.exists(here("output", "dds_cs")))) {
  dds_cs <- list()
  for (i in names(txi_cs)){
    # converting to deseq2 object
    dds_cs[[i]] <- DESeqDataSetFromTximport(txi_cs[[i]], colData = meta_cs[[i]], 
                                         design = ~embryo + subdomain, 
                                         rowData = rowdata[[i]])
    
    # changing this to make merge easier
    colnames(rowData(dds_cs[[i]]))[1] <- "row"
    
    # basic pre-filtering
    smallestGroupSize <- 3
    #keep <- rowSums(counts(dds_cs[[i]]) >= 10) >= smallestGroupSize
    keep <- filterByExpr(dds_cs[[i]], group = dds_cs[[i]]$subdomain)
    dds_cs[[i]] <- dds_cs[[i]][keep,]
    
    # estimate size factors - another assay added with normalisation factors
    dds_cs[[i]] <- estimateSizeFactors(dds_cs[[i]])
    
  }
  saveRDS(dds_cs, file = here("output", "dds_cs"))
} else{
  dds_cs <- readRDS(here("output", "dds_cs"))
}




# transforming and plotting to see the data
vst_blind <- list()
pca_plot <- list()
for (i in names(dds_cs)){
  # blind=TRUE should be used for comparing samples in a manner unbiased by prior 
  # information on samples, for example to perform sample QA (quality assurance)
  vst_blind[[i]] <- vst(dds_cs[[i]], blind=TRUE)
  
  # plotting pca
  pca_plot[[i]] <- plotPCA(vst_blind[[i]], intgroup=c("subdomain")) #+ geom_text(aes(label = vst_blind[[i]]$title))
  
}
patchwork::wrap_plots(pca_plot)


# removing batch effects

# using sva package which uses surrogate variables. Surrogate variables are covariates constructed directly
# from high-dimensional data that can be used in subsequent analyses to adjust for unknown,
# unmodeled, or latent sources of noise.
# sva package can be used for known and unknown sources of variation. 
# for known sources we used ComBat

# dat2 <- lapply(dds, function(a){
#               dat1 <- counts(a, normalized = TRUE)
#               idx <- rowMeans(dat1) > 1
#               dat1 <- dat1[idx,]
#               return(dat1)
#               })
# 
# modcombat <- lapply(dds, function(b){
#               model.matrix(~ 1, data = colData(b))
# })
# 
# combat <- mapply(function(x,y,z){
#   ComBat(dat=x, batch=colData(y)$embryo, mod=z, par.prior=TRUE, prior.plots=TRUE)},
#   dat, dds, modcombat)




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
contrasts <- list(CS17 = list(V_D = c("subdomain", "Ventral_V", "Dorsal_D"), 
                  V_VL = c("subdomain", "Ventral_V", "Ventro-lateral_VL"),
                  D_DL = c("subdomain", "Dorsal_D", "Dorsal-lateral_DL"),
                  VL_DL = c("subdomain", "Ventro-lateral_VL", "Dorsal-lateral_DL")), 
                  
                  CS16 = list(inner_V_D = c("subdomain", "Ventral_Inner_V_Inner", "Dorsal_Inner_D_Inner"),
                  mid_V_D = c("subdomain", "Ventral_Mid_V_Mid", "Dorsal_Mid_D_Mid"),
                  outer_V_D = c("subdomain", "Ventral_Outer_V_Outer", "Dorsal_Outer_D_Outer")))

results_list <- list() # dataframes with shrunken lfc
results_wald <- list() # wald test dataframes
reslfc <- list() # list of results for plots

if (!file.exists(here("output", "results_list"))){
  for (i in names(deseq)){
     #i <- "CS16"
    for (j in names(contrasts[[i]])){
      #j <- "inner_V_D"
      name1 <- paste0(i, j, "_apeglm_results")
      name2 <- paste0(i, j, "_results")
      
      res <- results(deseq[[i]], contrast = contrasts[[i]][[j]], alpha = 0.05) # cant use tidy = TRUE here since using lfc shrink
      
      # LFC shrinkage
      reslfc[[i]][[j]] <- lfcShrink(deseq[[i]], res = res, type = "ashr")
      apeglm <- data.frame(reslfc[[i]][[j]])
      apeglm$row <- rownames(apeglm)
      apeglm <- as.data.frame(merge(apeglm, rowData(dds[[i]]), all.x=TRUE))
      
      assign(name1, apeglm)
      results_list[[i]][[j]] <- apeglm
      
      results_wald[[i]][[j]] <- results(deseq[[i]], contrast = contrasts[[i]][[j]], alpha = 0.05, tidy = TRUE)
      results_wald[[i]][[j]] <- merge(results_wald[[i]][[j]], rowData(dds[[i]]), all.x=TRUE)
      
      assign(name2, data.frame(results_wald[[i]][[j]]))
      
    }
  }
  saveRDS(results_list, file = here("output", "results_list"))
  saveRDS(results_wald, file = here("output", "results_wald"))
}else{
  results_list <- readRDS(here("output", "results_list"))
  results_wald <- readRDS(here("output", "results_wald"))
  
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

for (i in names(reslfc)){
  for (j in names(reslfc[[i]]))
    plotMA(reslfc[[i]][[j]], main = j)
}


# plotting counts to look at genes with very high logfold change
plotCounts(dds[["CS16"]], gene="ENSG00000260371", intgroup="subdomain", normalized = TRUE)
plotDispEsts(deseq[["CS16"]])

boxplot(log10(assays(deseq[["CS16"]])[["cooks"]]), range=0, las=2)
