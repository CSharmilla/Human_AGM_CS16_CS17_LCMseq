

library(here)
library(patchwork)
library(DESeq2)
library(ashr)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db) # req to run gsea
library(ggplot2)
library(BiocParallel) # needed to update this; trying to fix an error "loadingError in serialize(data, node$con) : error writing to connection"
library(dplyr)

source("functions.R")

# reading all in
deseq <- readRDS(here("output", "deseq"))
dds <- readRDS(here("output", "dds"))
results_list <- readRDS(here("output", "results_list"))
results_wald <- readRDS(here("output", "results_wald"))



# transforming and plotting PCA to see the data - not BLIND

# blind=FALSE should be used for transforming data for downstream analysis, where the full use of the design information should be made. 
# blind=FALSE will skip re-estimation of the dispersion trend, if this has already been calculated. If many of genes have large 
# differences in counts due to the experimental design, it is important to set blind=FALSE for downstream analysis
vst <- lapply(deseq, vst, blind=FALSE)


# filtering
results_list_fil <- list()
results_wald_fil <- list()
for (i in names(results_list)){
  results_list_fil[[i]] <- dplyr::filter(results_list[[i]], padj < 0.05) 
  results_wald_fil[[i]] <- dplyr::filter(data.frame(results_wald[[i]]), padj < 0.05)
  
  # fix rownames
  rownames(results_list_fil[[i]]) <- 1:nrow(results_list_fil[[i]])
  rownames(results_wald_fil[[i]]) <- 1:nrow(results_wald_fil[[i]])
  
  
  # separate data frames so its easier to see
  name <- paste0(i, "_apeglm_results")
  assign(name, results_list[[i]])
}



# plotting top genes - always check with the results lists to make sure nothing funny is happening
selected <- list()
#dev.new()
for (i in names(results_list_fil)){
  plot_top_genes(results_list_fil[[i]], deseq[[i]], dds[[i]], vst[[i]], i)
}


# using clusterprofiler

genelist_ensembl <- list() # all genes used for gsea
genelist_entrez <- list()

fil_genelist_ensembl <- list() # filtered genes used for ORA
fil_genelist_entrez <- list()

# preparing gene lists
for (i in names(results_list)){
  #i <- "VL_DL"
  # extract genes to use
  a <- results_list[[i]]
  #sum(a$padj < 0.05, na.rm=TRUE)
  #a <- a[which(a$padj < 0.05),]
  colnames(a)[1] <- "ENSEMBL"
  
  # getting entrez ids
  entrez <- bitr(a$ENSEMBL, fromType = "ENSEMBL",
                 toType = c("ENTREZID", "SYMBOL"),
                 OrgDb = org.Hs.eg.db)
  # merging the two 
  a <- merge(a, entrez, by = "ENSEMBL")
  #a$ENTREZID <- make.unique(a$ENTREZID, sep = ".") # some genes have duplicated ensembl ids but different entrez ids. not a problem with
  # results tables from deseq
  #a$ENSEMBL <- make.unique(a$ENSEMBL, sep = ".")
  
  # list of ensembl gene lists
  genelist_ensembl[[i]] <- a %>% filter(!is.na(padj)) %>%
    arrange(desc(log2FoldChange)) %>%
    group_by(ENSEMBL) %>% # grouping by emsembl ids
    slice(which.max(abs(log2FoldChange))) %>% # Return a Single Row with the Maximum Value for the Group
    pull(log2FoldChange, name=ENSEMBL) # creating a named vector with logfc
  
  # filtered
  fil_genelist_ensembl[[i]] <- a %>% filter(padj < 0.05) %>%
    group_by(ENSEMBL) %>% # grouping by emsembl ids
    slice(which.max(abs(log2FoldChange))) %>%
    # get ids
    pull(ENSEMBL)
  
  
  # list of entrez gene lists
  genelist_entrez[[i]] <- a %>% filter(!is.na(padj)) %>%
    arrange(desc(log2FoldChange)) %>%
    group_by(ENTREZID) %>% # grouping by emsembl ids
    slice(which.max(abs(log2FoldChange))) %>% # Return a Single Row with the Maximum Value for the Group
    pull(log2FoldChange, name=ENTREZID) # creating anamed vector with logfc
  
  # filtered
  fil_genelist_entrez[[i]] <- a %>% filter(padj < 0.05) %>%
    group_by(ENTREZID) %>% # grouping by emsembl ids
    slice(which.max(abs(log2FoldChange))) %>% # Return a Single Row with the Maximum Value for the Group
    # get ids
    pull(ENTREZID)
  
}  


# Gene ontology GSEA
gse <- list()

for (i in names(genelist_ensembl)){
  
  # running gsego using ensembl ids
  gse[[i]] <- gseGO(geneList=genelist_ensembl[[i]],
                    ont ="ALL",
                    keyType = "ENSEMBL",
                    #nPerm = 10000,
                    #minGSSize = 3,
                    #maxGSSize = 800,
                    pvalueCutoff = 0.05,
                    verbose = TRUE,
                    OrgDb = org.Hs.eg.db,
                    pAdjustMethod = "none")
  
  #print(dotplot(gse[[i]], showCategory=15, split=".sign", font.size=7) + ggtitle(paste0("GSEA dotplot for ", i)) + facet_grid(.~.sign))
  #ridgeplot(gse) + labs(x = "enrichment distribution")
}

# separate loop since gsego takes some time to run
for (i in names(gse)){
  print(dotplot(gse[[i]], showCategory=15, split=".sign", font.size=9) + ggtitle(paste0("GSEA dotplot for ", i)) + facet_grid(.~.sign))
  
}

# GO ORA - nothing
# need entrez ids for this

ggo <- list()
for (i in names(fil_genelist_entrez)){
  ggo[[i]] <- enrichGO(gene          = fil_genelist_entrez[[i]],
                       universe      = fil_genelist_entrez[[i]],
                       OrgDb         = org.Hs.eg.db,
                       ont           = "CC",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.01,
                       qvalueCutoff  = 0.05,
                       readable      = TRUE)
  
}


# KEGG enrichment analysis

# looking at KEGG pathways first

# KEGG pathway ORA analysis - also requires entrez ids.
kk <- list()
for (i in names(fil_genelist_entrez)){
  kk[[i]] <- enrichKEGG(gene  = fil_genelist_entrez[[i]],
                        organism     = 'hsa',
                        pvalueCutoff = 0.05)
  
  # to take a closer look
  name <- paste0("kk_", i)
  assign(name, data.frame(kk[[i]]))
}
head(kk) 

random2 <- kk_VL_DL

# browse KEGG pathways
browseKEGG(kk[["V_VL"]], 'hsa04144')
browseKEGG(kk[["VL_DL"]], 'hsa04110')




# KEGG pathway gsea - nothing here
kk_gsea <- list()

# requires entrez ids
for (i in names(genelist_entrez)){
  #i <- "V_D"
  kk_gsea[[i]] <- gseKEGG(geneList = genelist_entrez[[i]],
                          organism     = 'hsa',
                          minGSSize    = 10,
                          pvalueCutoff = 0.05,
                          verbose      = FALSE)
}


# KEGG module ORA
mkk_ora <- list()
for (i in names(genelist_entrez)){
  mkk_ora[[i]] <- enrichMKEGG(gene = names(genelist_entrez[[i]]),
                              organism = 'hsa',
                              pvalueCutoff = 1,
                              qvalueCutoff = 1)
  
  # to take a closer look
  name <- paste0("mkk_ora", i)
  assign(name, data.frame(mkk_ora[[i]]))
}

# KEGG module GSEA
mkk_gsea <- list()
for (i in names(genelist_entrez)){
  mkk_gsea[[i]] <- gseMKEGG(geneList = genelist_entrez[[i]],
                            organism = 'hsa',
                            pvalueCutoff = 1)
  
  # to take a closer look
  #name <- paste0("mkk_gsea", i)
  #assign(name, data.frame(mkk_gsea[[i]]))
}




