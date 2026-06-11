library(Glimma)
library(here)
library(dplyr)
library(DESeq2)
library(tibble)
library(pheatmap)
library(EnhancedVolcano)
library(decoupleR)

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



a <- results_list[["V_D"]]

for (i in names(results_list[6:8])){
  #i <- "V_D"
  
  
  
  print(EnhancedVolcano(results_list[[i]],
                        lab = results_list[[i]]$external_gene_name,
                        x = 'log2FoldChange',
                        y = 'padj',
                        pCutoff = 0.05,
                        FCcutoff = 1,
                        title = paste0("Differentially expressed genes for ", i),
                        drawConnectors = TRUE,
                        caption = bquote(~Log[2]~ "fold change cutoff, 1; adjp-value cutoff, 0.05"),
                        max.overlaps = 20,
                        labSize = 4,
                        subtitle = ""))
  
}




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
  
  # print(EnhancedVolcano(sf_list[[i]],
  #                       lab = sf_list[[i]]$external_gene_name,
  #                       x = 'log2FoldChange',
  #                       y = 'padj',
  #                       pCutoff = 0.05,
  #                       FCcutoff = 1, 
  #                       title = paste0("Differentially expressed secreted factors for ", i),
  #                       drawConnectors = TRUE,
  #                       caption = bquote(~Log[2]~ "fold change cutoff, 1; adjp-value cutoff, 0.05"),
  #                       max.overlaps = 20,
  #                       labSize = 4,
  #                       subtitle = ""))
  
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
  
  # print(EnhancedVolcano(ecm_list[[i]],
  #                       lab = ecm_list[[i]]$external_gene_name,
  #                       x = 'log2FoldChange',
  #                       y = 'padj',
  #                       pCutoff = 0.05,
  #                       FCcutoff = 1, 
  #                       title = paste0("Differentially expressed matrisome genes for ", i),
  #                       caption = bquote(~Log[2]~ "fold change cutoff, 1; adjp-value cutoff, 0.05"),
  #                       drawConnectors = TRUE,
  #                       max.overlaps = 20,
  #                       labSize = 4,
  #                       subtitle = ""))
}
vd <- data.frame(results_list[["V_D"]])

# transcription factors
net <- get_collectri(organism='human', split_complexes=FALSE) 
 
 
# CS16 mid V vs D (secreted factors and ECM genes) -> CS17 V vs D interactions
sf_mid_V_D <- sf_list[["mid_V_D"]] %>% filter(padj < 0.05, log2FoldChange > 0)
res_vd <- results_list[["V_D"]] %>% as.data.frame() %>%
  filter(padj < 0.05, log2FoldChange > 0)
# now for ECM genes
ecm_mid_V_D <- ecm_list[["mid_V_D"]] %>% filter(padj < 0.05, log2FoldChange > 0)


# creating csv to give to STRING
# Extract gene lists
genes_sf <- data.frame(genes = sf_mid_V_D$external_gene_name, colour = "#F5C827") 
genes_res <- data.frame(genes = res_vd$external_gene_name, colour = "#F54927") 
genes_ecm <- data.frame(genes = ecm_mid_V_D$external_gene_name, colour = "#27ADF5") 

# combining
genes_df <- rbind(genes_ecm, genes_res, genes_sf)
genes_df$colour[genes_df$genes %in% net$source] <- "#a64dff"


# Combine 
all_genes <- genes_df[!duplicated(genes_df$genes) & genes_df$genes != "", ]

write.csv(all_genes$genes,
          file = here("output", "genes_for_STRING.csv"),
          row.names = FALSE)
write.csv(all_genes,
          file = here("output", "labelled_genes_for_STRING.csv"),
          row.names = FALSE)


# CS17 V vs D (secreted factors and ECM genes) -> CS17 mid V vs D interactions
# sf_V_D <- sf_list[["V_D"]] %>% filter(padj < 0.05, log2FoldChange > 0)
# res_mid_vd <- results_list[["mid_V_D"]] %>% as.data.frame() %>%
#   filter(padj < 0.05, log2FoldChange > 0)
# # now for ECM genes
# ecm_V_D <- ecm_list[["V_D"]] %>% filter(padj < 0.05, log2FoldChange > 0)
# 
# 
# # creating csv to give to STRING
# # Extract gene lists
# genes_sf_vd <- data.frame(genes = sf_V_D$external_gene_name, colour = "#F5C827") 
# genes_res_mid_vd <- data.frame(genes = res_mid_vd$external_gene_name, colour = "#F54927") 
# genes_ecm_vd <- data.frame(genes = ecm_V_D$external_gene_name, colour = "#27ADF5") 
# 
# # combining
# genes_df_2 <- rbind(genes_ecm_vd, genes_res_mid_vd, genes_sf_vd)
# # Combine 
# all_genes_2 <- genes_df_2[!duplicated(genes_df_2$genes) & genes_df_2$genes != "", ]
# 
# write.csv(all_genes_2$genes,
#           file = here("output", "genes_for_STRING_2.csv"),
#           row.names = FALSE)
# write.csv(all_genes_2,
#           file = here("output", "labelled_genes_for_STRING_2.csv"),
#           row.names = FALSE)

sf_V_D <- sf_list[["V_D"]] %>% filter(padj < 0.05, log2FoldChange > 0)

# now for ECM genes
ecm_V_D <- ecm_list[["V_D"]] %>% filter(padj < 0.05, log2FoldChange > 0)


# creating csv to give to STRING
# Extract gene lists
genes_res <- data.frame(genes = res_vd$external_gene_name, colour = "#F54927") 

genes_res$colour[genes_res$genes %in% sf_V_D$external_gene_name] <- "#F5C827"
genes_res$colour[genes_res$genes %in% ecm_V_D$external_gene_name] <- "#27ADF5"
genes_res$colour[genes_res$genes %in% net$source] <- "#a64dff"

# Combine 
all_genes <- genes_res[!duplicated(genes_res$genes) & genes_res$genes != "", ]

write.csv(all_genes$genes,
          file = here("output", "vd_genes_for_STRING.csv"),
          row.names = FALSE)
write.csv(all_genes,
          file = here("output", "vd_labelled_genes_for_STRING.csv"),
          row.names = FALSE)


# V_VL
res_vvl <- results_list[["V_VL"]] %>% as.data.frame() %>%
  filter(padj < 0.05, log2FoldChange > 0)

sf_V_VL <- sf_list[["V_VL"]] %>% filter(padj < 0.05, log2FoldChange > 0)

# now for ECM genes
ecm_V_VL <- ecm_list[["V_VL"]] %>% filter(padj < 0.05, log2FoldChange > 0)


# creating csv to give to STRING
# Extract gene lists
genes_res <- data.frame(genes = res_vvl$external_gene_name, colour = "#F54927") 

genes_res$colour[genes_res$genes %in% sf_V_VL$external_gene_name] <- "#F5C827"
genes_res$colour[genes_res$genes %in% ecm_V_VL$external_gene_name] <- "#27ADF5"
genes_res$colour[genes_res$genes %in% net$source] <- "#a64dff"

# Combine 
all_genes <- genes_res[!duplicated(genes_res$genes) & genes_res$genes != "", ]

write.csv(all_genes$genes,
          file = here("output", "vvl_genes_for_STRING.csv"),
          row.names = FALSE)
write.csv(all_genes,
          file = here("output", "vvl_labelled_genes_for_STRING.csv"),
          row.names = FALSE)

# VL to V
res_vvl <- results_list[["V_VL"]] %>% as.data.frame() %>%
  filter(padj < 0.05, log2FoldChange > 0)

sf_V_VL <- sf_list[["V_VL"]] %>% filter(padj < 0.05, log2FoldChange > 0)

# now for ECM genes
ecm_V_VL <- ecm_list[["V_VL"]] %>% filter(padj < 0.05, log2FoldChange > 0)


# creating csv to give to STRING
# Extract gene lists
genes_res <- data.frame(genes = res_vvl$external_gene_name, colour = "#F54927") 

genes_res$colour[genes_res$genes %in% sf_V_VL$external_gene_name] <- "#F5C827"
genes_res$colour[genes_res$genes %in% ecm_V_VL$external_gene_name] <- "#27ADF5"
genes_res$colour[genes_res$genes %in% net$source] <- "#a64dff"

# Combine 
all_genes <- genes_res[!duplicated(genes_res$genes) & genes_res$genes != "", ]

write.csv(all_genes$genes,
          file = here("output", "vvl_genes_for_STRING.csv"),
          row.names = FALSE)
write.csv(all_genes,
          file = here("output", "vvl_labelled_genes_for_STRING.csv"),
          row.names = FALSE)
