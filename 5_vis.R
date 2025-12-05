library(Glimma)
library(here)
library(dplyr)
library(DESeq2)
library(tibble)
library(pheatmap)
library(EnhancedVolcano)

# reading in data
dds <- readRDS(here("output", "dds"))
deseq <- readRDS(here("output", "deseq"))
results_wald <- readRDS(here("output", "results_wald")) # need stat
results_list <- readRDS(here("output", "results_list"))

# reading in secreted proteins list
sec_prot <-read.csv(here("support", "SEPDB_genes.csv"))
ecm <- read.delim(here("support", "ECM_genes.txt"))

# vst
vst <- lapply(deseq, vst, blind=FALSE)

# trying Glimma
glimmaMDS(deseq[["V_D"]])
glimmaMA(deseq[["V_D"]])


# subsetting secreted factors
sf_list <- list()
sf_assay <- list()

for (i in names(results_list)){
  #i <- "V_D"
  
  # secreted factors from deseq results
  sf_list[[i]] <- results_list[[i]] %>% as.data.frame() %>%
    #filter(padj < 0.05) %>%
    filter(external_gene_name %in% sec_prot$GeneName) %>%
    group_by(external_gene_name) %>%
    dplyr::slice(which.max(abs(log2FoldChange))) %>%
    arrange(desc(log2FoldChange))
  
  # vst assay for secreted proteins
  sf_assay[[i]] <- vst[[i]] %>% assay() %>%
    merge(y = rowData(dds[[i]]), by=0) %>%
    semi_join(y = sf_list[[i]], by=c("row")) %>%
    column_to_rownames(var = "external_gene_name") %>%
    dplyr::select(-c(Row.names, row))
  
  sf_assay[[i]] <- sf_assay[[i]][match(sf_list[[i]]$external_gene_name, rownames(sf_assay[[i]])),]
  # annotation df
  df <- as.data.frame(colData(deseq[[i]])[,c("embryo","subdomain")])
  
  plot_assay <- na.omit(sf_assay[[i]][1:70,])
  # plotting heatmaps
  #print(pheatmap(plot_assay, cluster_rows=FALSE, show_rownames=TRUE,
                 #cluster_cols=TRUE, annotation_col=df, main = i, fontsize_row = 10, fontsize_col = 10, angle_col = 90))
  
  print(EnhancedVolcano(sf_list[[i]],
                        lab = sf_list[[i]]$external_gene_name,
                        x = 'log2FoldChange',
                        y = 'padj',
                        pCutoff = 0.05,
                        FCcutoff = 1, 
                        title = paste0("Differentially expressed secreted factors for ", i),
                        drawConnectors = TRUE,
                        caption = bquote(~Log[2]~ "fold change cutoff, 1; adjp-value cutoff, 0.05"),
                        max.overlaps = 20,
                        labSize = 4,
                        subtitle = ""))
}

EnhancedVolcano(sf_list[[i]],
                lab = sf_list[[i]]$external_gene_name,
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 1, 
                main = paste0("Differentially expressed secreted factors for ", i)) #+
  #ggrepel::geom_label_repel(mapping = ggplot2::aes(label = external_gene_name),max.overlaps = Inf)


# ECM genes

ecm_list <- list()
ecm_assay <- list()

for (i in names(results_list)){
  #i <- "V_D"
  
  # matrisome factors from deseq results
  ecm_list[[i]] <- results_list[[i]] %>% as.data.frame() %>%
    #filter(padj < 0.05) %>%
    filter(external_gene_name %in% ecm$Gene.Name) %>%
    group_by(external_gene_name) %>%
    dplyr::slice(which.max(abs(log2FoldChange))) %>%
    arrange(desc(log2FoldChange))
  
  # vst assay for matrisome factors
  ecm_assay[[i]] <- vst[[i]] %>% assay() %>%
    merge(y = rowData(dds[[i]]), by=0) %>%
    semi_join(y = ecm_list[[i]], by=c("row")) %>%
    column_to_rownames(var = "external_gene_name") %>%
    dplyr::select(-c(Row.names, row))
  
  ecm_assay[[i]] <- ecm_assay[[i]][match(ecm_list[[i]]$external_gene_name, rownames(ecm_assay[[i]])),]
  # annotation df
  df <- as.data.frame(colData(deseq[[i]])[,c("embryo","subdomain")])
  
  plot_assay <- na.omit(ecm_assay[[i]][1:70,])
  # plotting heatmaps
  #print(pheatmap(plot_assay, cluster_rows=FALSE, show_rownames=TRUE,
  #cluster_cols=TRUE, annotation_col=df, main = i, fontsize_row = 10, fontsize_col = 10, angle_col = 90))
  
  print(EnhancedVolcano(ecm_list[[i]],
                        lab = ecm_list[[i]]$external_gene_name,
                        x = 'log2FoldChange',
                        y = 'padj',
                        pCutoff = 0.05,
                        FCcutoff = 1, 
                        title = paste0("Differentially expressed matrisome genes for ", i),
                        caption = bquote(~Log[2]~ "fold change cutoff, 1; adjp-value cutoff, 0.05"),
                        drawConnectors = TRUE,
                        max.overlaps = 20,
                        labSize = 4,
                        subtitle = ""))
}
 vd <- data.frame(results_list[["V_D"]])
