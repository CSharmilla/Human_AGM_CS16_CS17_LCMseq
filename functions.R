

plot_top_genes <- function(res, deseq, dds, vst, title){
  # troubleshooting
  # i <- "V_D"
  # res <- results_list_fil_CS17[[i]]
  # deseq <- deseq_CS17
  # dds <- dds_CS17
  # vst <- vst_CS17

  selected <- res[order(res$padj, decreasing = FALSE)[1:50], ] 
  df <- as.data.frame(colData(deseq)[,c("embryo","subdomain")])
  assay <- assay(vst)[selected$row,]
  assay <- merge(assay, rowData(dds), by=0) # neater without all the rows from running deseq
  assay <- assay[!duplicated(assay$external_gene_name), ]
  rownames(assay) <- ifelse(assay$external_gene_name == "", assay$row, assay$external_gene_name)
  assay <- assay[,colnames(deseq)]
  print(pheatmap(assay, cluster_rows=TRUE, show_rownames=TRUE,
                 cluster_cols=TRUE, annotation_col=df, main = title))
}

plot_top_genes_up <- function(res, deseq, dds, vst, title){
  # troubleshooting
  # i <- "V_D"
  # res <- results_list_fil_CS17[[i]]
  # deseq <- deseq_CS17
  # dds <- dds_CS17
  # vst <- vst_CS17
  
  selected <- res[order(res$log2FoldChange, decreasing = TRUE)[1:50], ] 
  df <- as.data.frame(colData(deseq)[,c("embryo","subdomain")])
  assay <- assay(vst)[selected$row,]
  assay <- merge(assay, rowData(dds), by=0) # neater without all the rows from running deseq
  assay <- assay[!duplicated(assay$external_gene_name), ]
  rownames(assay) <- ifelse(assay$external_gene_name == "", assay$row, assay$external_gene_name)
  assay <- assay[,colnames(deseq)]
  print(pheatmap(assay, cluster_rows=TRUE, show_rownames=TRUE,
                 cluster_cols=TRUE, annotation_col=df, main = paste0("upregulated for ", title)))
}

plot_top_genes_down <- function(res, deseq, dds, vst, title){
  # troubleshooting
  # i <- "V_D"
  # res <- results_list_fil_CS17[[i]]
  # deseq <- deseq_CS17
  # dds <- dds_CS17
  # vst <- vst_CS17
  
  selected <- res[order(res$log2FoldChange, decreasing = FALSE)[1:50], ] 
  df <- as.data.frame(colData(deseq)[,c("embryo","subdomain")])
  assay <- assay(vst)[selected$row,]
  assay <- merge(assay, rowData(dds), by=0) # neater without all the rows from running deseq
  assay <- assay[!duplicated(assay$external_gene_name), ]
  rownames(assay) <- ifelse(assay$external_gene_name == "", assay$row, assay$external_gene_name)
  assay <- assay[,colnames(deseq)]
  print(pheatmap(assay, cluster_rows=TRUE, show_rownames=TRUE,
                 cluster_cols=TRUE, annotation_col=df, main = paste0("downregulated for ", title)))
}


pathview_mod <- function(region, pathway){
  pathview(gene.data  = fil_genelist_entrez[[region]],
           pathway.id = pathway,
           species    = "hsa",
           limit      = list(gene=round(max(abs(fil_genelist_entrez[[region]]))), cpd=1),
           key.pos="topright",
           low=list(gene="red"),
           high=list(gene="green"),
           out.suffix = region)
}


# function to look at specific TFs

#Here blue means that the sign of multiplying the mor and t-value is negative, 
# meaning that these genes are “deactivating” the TF, and red means that the sign is positive, 
# meaning that these genes are “activating” the TF.

TF_targets <- function(TF, region){
  
  df <- net %>%
    filter(source == TF) %>%
    arrange(target) %>%
    mutate(ID = target, color = "3") %>%
    column_to_rownames('target')
  
  inter <- sort(intersect(rownames(deg_all[[region]]),rownames(df)))
  df <- df[inter, ]
  df <- df %>% merge(y=deg_all[[region]][inter, ], by = 0) %>% column_to_rownames(var = "Row.names")
  
  df <- df %>%
    mutate(color = if_else(mor > 0, '1', color)) %>%
    mutate(color = if_else(mor < 0, '2', color)) 
          #%>%
    #mutate(color = if_else(mor > 0 & stat > 0, '1', color)) %>%
    #mutate(color = if_else(mor > 0 & stat < 0, '2', color)) %>%
    #mutate(color = if_else(mor < 0 & stat > 0, '2', color)) %>%
    #mutate(color = if_else(mor < 0 & stat < 0, '1', color))
  
  ggplot(df, aes(x = log2FoldChange, y = -log10(pvalue), color = color, size=abs(mor))) +
    geom_point() +
    scale_colour_manual(values = c("red","royalblue3","grey")) +
    geom_label_repel(aes(label = ID, size=1), max.overlaps = Inf) + 
    theme_minimal() +
    theme(legend.position = "none") +
    geom_vline(xintercept = 0, linetype = 'dotted') +
    geom_hline(yintercept = 0, linetype = 'dotted') +
    ggtitle(TF)
}


# what genes are in specific pathways

pathway_genes <- function(pathway, region){
  df <- pro %>%
    dplyr::filter(source == pathway) %>%
    dplyr::arrange(target) %>%
    dplyr::mutate(ID = target, 
                  color = "3") %>%
    tibble::column_to_rownames('target')
  
  inter <- sort(dplyr::intersect(rownames(deg_all[[region]]), rownames(df)))
  
  df <- df[inter, ]
  
  df <- df %>% merge(y=deg_all[[region]][inter, ], by = 0) %>% column_to_rownames(var = "Row.names")
  
  df <- df %>%
    dplyr::mutate(color = dplyr::if_else(weight > 0 & stat > 0, '1', color)) %>%
    dplyr::mutate(color = dplyr::if_else(weight > 0 & stat < 0, '2', color)) %>%
    dplyr::mutate(color = dplyr::if_else(weight < 0 & stat > 0, '2', color)) %>%
    dplyr::mutate(color = dplyr::if_else(weight < 0 & stat < 0, '1', color))
  
  
    #dplyr::mutate(color = dplyr::if_else(weight > 0, '1', color)) %>%
    #dplyr::mutate(color = dplyr::if_else(weight < 0, '2', color)) 
    
    
    
  
  colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")[c(2, 10)])
  
  p <- ggplot2::ggplot(data = df, 
                       mapping = ggplot2::aes(x = weight, 
                                              y = log2FoldChange, 
                                              color = color)) + 
    ggplot2::geom_point(size = 2.5, 
                        color = "black") + 
    ggplot2::geom_point(size = 1.5) +
    ggplot2::scale_colour_manual(values = c(colors[2], colors[1], "grey")) +
    ggrepel::geom_label_repel(mapping = ggplot2::aes(label = ID),max.overlaps = Inf) + 
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none") +
    ggplot2::geom_vline(xintercept = 0, linetype = 'dotted') +
    ggplot2::geom_hline(yintercept = 0, linetype = 'dotted') +
    ggplot2::ggtitle(paste0(pathway, " genes for ", region))
  
  print(p)
  return(inter)
}


