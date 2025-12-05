
library(DESeq2)
library(tibble)
library(decoupleR)
library(OmnipathR)
library(tidyr)
library(Glimma)
library(ggplot2)

source("functions.R")

# reading in data
deseq <- readRDS(here("output", "deseq"))
dds <- readRDS(here("output", "dds"))
results_wald <- readRDS(here("output", "results_wald")) # need stat
results_list <- readRDS(here("output", "results_list"))

# vst
vst <- lapply(deseq, vst, blind=FALSE)

# log transform 
rlog <- list()
for (i in names(deseq)){
  rlog[[i]] <- rlog(deseq[[i]])
}


# cleaning up input data
counts <- list()
design <- list()
deg <- list()
deg_all <- list()

# 
# genelist_ensembl[[i]] <- a %>% filter(!is.na(padj)) %>%
#   group_by(ENSEMBL) %>% # grouping by emsembl ids
#   slice(which.max(abs(log2FoldChange))) %>% # Return a Single Row with the Maximum Value for the Group
#   arrange(desc(log2FoldChange)) %>%
#   pull(log2FoldChange, name=ENSEMBL) 

for (i in names(rlog)){
  
  # extract stat from DEG
  deg_all[[i]] <- results_wald[[i]] %>% as.data.frame() %>%
    filter(padj < 0.05) %>%
    filter(!external_gene_name == "") %>%
    group_by(external_gene_name) %>%
    dplyr::slice(which.max(abs(log2FoldChange))) %>%
    dplyr::select(log2FoldChange, stat, pvalue, row, external_gene_name) %>%
    column_to_rownames(var = "external_gene_name")
      
  # just stat  
  deg[[i]] <- deg_all[[i]] %>% dplyr::select(stat)
  
  # create count assays
  counts[[i]] <- rlog[[i]] %>% assay() %>%
    merge(y = rowData(dds[[i]]), by=0) %>%
    semi_join(y = deg_all[[i]], by=c("row")) %>%
    column_to_rownames(var = "external_gene_name") %>%
    dplyr::select(-c(Row.names, row))
    
  # create design
  design[[i]] <- as.data.frame(colData(deseq[[i]])[,c("embryo","subdomain")])  
    
}    
    
   
  
net <- get_collectri(organism='human', split_complexes=FALSE)
net

# jasper <- read.table("M:\\TF_targets\\gene_attribute_edges_jasper.txt", header = TRUE)
# chea <- read.table("M:\\TF_targets\\gene_attribute_edges_chea.txt", header = TRUE)
# encode <- read.table("M:\\TF_targets\\gene_attribute_edges_encode.txt", header = TRUE)
# 
# # targets of pou6f2
# df_pou6f2.1 <- rbind(jasper[jasper$source=="POU6F2",c("source", "target")], chea[chea$source=="POU6F2",c("source", "target")], 
#                      encode[encode$source=="POU6F2",c("source", "target")], net[net$source=="POU6F2",c("source", "target")])
# 
# df_pou6f2 <- unique(df_pou6f2)
# rownames(df_pou6f2) <- 1:nrow(df_pou6f2)
# write.csv(df_pou6f2, "M:\\TF_targets\\POU6F2_targets.csv")
 
# plotting transcription factor heatmaps
assay <- list()
dev.new()
for (i in names(deg_all)){
  
  # vst assay of TFs from net
  assay[[i]] <- vst[[i]] %>% assay() %>%
    merge(y = rowData(dds[[i]]), by=0) %>%
    semi_join(y = deg_all[[i]], by=c("row")) %>%
    filter(external_gene_name %in% net$source) %>%
    column_to_rownames(var = "external_gene_name") %>%
    select(-c(Row.names, row))
  
  # annotation df
  df <- as.data.frame(colData(deseq[[i]])[,c("embryo","subdomain")])
  
  plot_assay <- na.omit(assay[[i]][1:50,])
  # plotting heatmaps
  print(pheatmap(plot_assay, cluster_rows=TRUE, show_rownames=TRUE,
                 cluster_cols=TRUE, annotation_col=df, main = i, fontsize_row = 10, fontsize_col = 10, angle_col = 90))
}

# TF volcano plots
tfs_calvanese <- c("TAL1", "SCL", "SOX17", "HOXA5", "HOXA7", "HOXA9", "RUNX1", "MYB", "GFI1", "GFI1B", "MLLT3", "HLF", "MECOM", "MSI2")



for (i in names(results_list)){
  
  
  #i <- "V_D"
  
  #file_name = paste("DE_volcano_", i, ".tiff", sep="")
  #tiff(here("figs",file_name), width = 1500, height = 1500)
  tfs_label <- results_list[[i]] %>% filter(external_gene_name %in% net$source) #%>% pull(external_gene_name)
  
  # keyvals.shape <- ifelse(tfs_label$external_gene_name %in% tfs_calvanese, 17, 19)
  # names(keyvals.shape)[keyvals.shape == 17] <- 'Key transcription factors functionally implicated in human HSC development'
  # names(keyvals.shape)[keyvals.shape == 19] <- 'Other TFs'
  # 
  # keyvals.colour <- ifelse(tfs_label$external_gene_name %in% tfs_calvanese, "pink", 19)
  # names(keyvals.colour)[keyvals.colour == 17] <- 'Key transcription factors functionally implicated in human HSC development'
  # names(keyvals.colour)[keyvals.colour == 19] <- 'Other TFs'
  
  # print(EnhancedVolcano(tfs_label,
  #                       lab = tfs_label$external_gene_name,
  #                       selectLab = tfs_label$external_gene_name,
  #                       x = 'log2FoldChange',
  #                       y = 'padj',
  #                       pCutoff = 0.05,
  #                       FCcutoff = 1,
  #                       title = paste0("Differentially expressed transcription factors for ", i),
  #                       pointSize = 3,
  #                       #shapeCustom = keyvals.shape,
  #                       #drawConnectors = TRUE,
  #                       #widthConnectors = 0.75,
  #                       lengthConnectors = unit(0.006, 'npc'),
  #                       max.overlaps = 30,
  #                       labSize = 3.5,
  #                       caption = bquote(~Log[2]~ "fold change cutoff, 1; adjp-value cutoff, 0.05"),
  #                       subtitle = ""))
  
  print(EnhancedVolcano(tfs_label,
                        lab = tfs_label$external_gene_name,
                        selectLab = tfs_calvanese,
                        x = 'log2FoldChange',
                        y = 'padj',
                        pCutoff = 0.05,
                        FCcutoff = 1,
                        title = paste0("Key transcription factors functionally implicated in human HSC development for ", i),
                        titleLabSize = 10,
                        drawConnectors = TRUE,
                        #widthConnectors = 0.75,
                        lengthConnectors = unit(0.006, 'npc'),
                        max.overlaps = 30,
                        labSize = 3.5,
                        caption = bquote(~Log[2]~ "fold change cutoff, 1; adjp-value cutoff, 0.05"),
                        subtitle = ""))
  #dev.off()
}


# Run ulm - with input matrix
sample_acts <- lapply(counts, run_ulm, net=net, .source='source', .target='target', .mor='mor', minsize = 5)
sample_acts

for (i in names(sample_acts)){
  a <- sample_acts[[i]]
  name <- paste0("sample_acts_", i)
  assign(name, a)
}

n_tfs <- 100

for (i in names(sample_acts)) {
  
  # Transform to wide matrix
  sample_acts_mat <- sample_acts[[i]] %>%
    pivot_wider(id_cols = 'condition', names_from = 'source',
                values_from = 'score') %>%
    column_to_rownames('condition') %>%
    as.matrix()
  
  # Get top tfs with more variable means across clusters
  tfs <- sample_acts[[i]] %>%
    group_by(source) %>%
    summarise(std = sd(score)) %>%
    arrange(-abs(std)) %>%
    head(n_tfs) %>%
    pull(source)
  sample_acts_mat <- sample_acts_mat[,tfs]
  
  # Scale per sample
  sample_acts_mat <- scale(sample_acts_mat)
  
  # Choose color palette
  palette_length = 100
  my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)
  
  my_breaks <- c(seq(-3, 0, length.out=ceiling(palette_length/2) + 1),
                 seq(0.05, 3, length.out=floor(palette_length/2)))
  
  # Plot
  pheatmap(sample_acts_mat, border_color = NA, color=my_color, breaks = my_breaks,
           main = paste0("Infered TF activity for ", i), angle_col = 45) 
  
}

# having a look at this
test <- V_D_apeglm_results[V_D_apeglm_results$external_gene_name %in% net$target[net$source == "ATF2"],]


# Run ulm - infer TF activities from the t-values of the DEGs
contrast_acts <- lapply(deg, run_ulm, net=net, .source='source', .target='target',
                        .mor='mor', minsize = 5)

contrast_acts_V_D <- contrast_acts$V_D
contrast_acts_V_D <- contrast_acts$V_D
cal_contrast_acts_V_D <- contrast_acts_V_D[contrast_acts_V_D$source %in% tfs_calvanese,]


for (i in names(contrast_acts)){
  a <- contrast_acts[[i]]
  name <- paste0("contrast_acts_", i)
  assign(name, a)
}


# Filter top TFs in both signs
n_tfs <- 100

for (i in names(contrast_acts)){
  #i <- "V_D"
  f_contrast_acts <- contrast_acts[[i]] %>%
    mutate(rnk = NA)
  msk <- f_contrast_acts$score > 0
  f_contrast_acts[msk, 'rnk'] <- rank(-f_contrast_acts[msk, 'score'])
  f_contrast_acts[!msk, 'rnk'] <- rank(-abs(f_contrast_acts[!msk, 'score']))
  tfs <- f_contrast_acts %>%
    arrange(rnk) %>%
    head(n_tfs) %>%
    pull(source)
  f_contrast_acts <- f_contrast_acts %>%
    filter(source %in% tfs)
  
  # Plot
  plot <- ggplot(f_contrast_acts, aes(x = reorder(source, score), y = score)) + 
    geom_bar(aes(fill = score), stat = "identity") +
    scale_fill_gradient2(low = "darkblue", high = "indianred", 
                         mid = "whitesmoke", midpoint = 0) + 
    theme_minimal() +
    theme(axis.title = element_text(face = "bold", size = 12),
          axis.text.x = 
            element_text(angle = 45, hjust = 1, size =10, face= "bold"),
          axis.text.y = element_text(size =10, face= "bold"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) +
    xlab("TFs") + ggtitle(paste0("Infered TF activity for ", i))
  print(plot)
  
  
  
  
}



# looking at specific TFs
TF_targets('MYB', "DL_D")

TF_targets('YAP1', "VL_DL")

TF_targets('MYB', "V_D")


## progeny
pro <- get_progeny(organism = 'human')#, top = 500)
pro
pro <- rbind(pro, list("JAK-STAT", "STAT5A", 3.633489, 5.155198e-17))

sample_acts_pro <- lapply(counts, run_mlm, net=pro, .source='source', .target='target',
                          .mor='weight', minsize = 5)
sample_acts_pro

for (i in names(sample_acts_pro)){
  # Transform to wide matrix
  sample_acts_mat <- sample_acts_pro[[i]] %>%
    pivot_wider(id_cols = 'condition', names_from = 'source',
                values_from = 'score') %>%
    column_to_rownames('condition') %>%
    as.matrix()
  
  # Scale per feature
  sample_acts_mat <- scale(sample_acts_mat)
  
  # Choose color palette
  palette_length = 100
  my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)
  
  my_breaks <- c(seq(-3, 0, length.out=ceiling(palette_length/2) + 1),
                 seq(0.05, 3, length.out=floor(palette_length/2)))
  
  # Plot
  #print(pheatmap(sample_acts_mat, border_color = NA, color=my_color, breaks = my_breaks, main = paste0("Pathway inference for ", i))) 
  
}
vd <- sample_acts_pro[["V_D"]]
a <- deg_all[["V_VL"]]
# looking at specific pathway genes
pathway_genes("MAPK", "V_D")
pathway_genes("MAPK", "V_VL")
pathway_genes("MAPK", "VL_DL")
pathway_genes("MAPK", "DL_D")


pathway_genes("TNFa", "V_D")
pathway_genes("MAPK", "V_VL")


# trying Glimma
glimmaMDS(dds)
