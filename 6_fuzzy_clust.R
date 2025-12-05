library(here)
library(Mfuzz)
library(marray)
library(dplyr)
library(tibble)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ComplexHeatmap)
library(circlize)
library(decoupleR)
library(OmnipathR)
library(tidyr)
library(DESeq2)


# reading in data separated by CS
txi_cs <- readRDS(here("output", "txi_cs"))
meta_cs <- readRDS(here("output", "meta_cs"))
genes_list <- readRDS(here("output", "genes_list"))

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
  rowdata[[i]] <- rowdata[[i]][match(rownames(txi_cs[[i]]$counts), rowdata[[i]]$ensembl_gene_id),] 
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

if (!(file.exists(here("output", "dds_full")))) {
  dds_full <- list()
  for (i in names(txi_cs)){
    # converting to deseq2 object
    dds_full[[i]] <- DESeqDataSetFromTximport(txi_cs[[i]], colData = meta_cs[[i]], 
                                                           design = ~embryo + subdomain, 
                                                           rowData = rowdata[[i]])
    
    # changing this to make merge easier
    colnames(rowData(dds_full[[i]]))[1] <-  "row"
    
    # basic pre-filtering
    smallestGroupSize <- 3
    #keep <- rowSums(counts(dds[[i]]) >= 10) >= smallestGroupSize
    keep <- filterByExpr(dds_full[[i]], group = dds_full[[i]]$subdomain)
    dds_full[[i]] <- dds_full[[i]][keep,]
    
    # estimate size factors - another assay added with normalisation factors
    dds_full[[i]] <- estimateSizeFactors(dds_full[[i]])
    
    # add normalised counts to assay
    assays(dds_full[[i]])[['sf_norm']] <- counts(dds_full[[i]], normalized = TRUE)
    
    
  }
  saveRDS(dds_full, file = here("output", "dds_full"))
  
} else{
  dds_full <- readRDS(here("output", "dds_full"))
  
}

# averaging replicates
averaged <- list()
for (i in names(dds_full)){
  
  rs <- rowsum(t(assay(dds_full[[i]], "sf_norm")), group = dds_full[[i]]$subdomain)
  rn <- as.numeric(table(dds_full[[i]]$subdomain))
  averaged[[i]] <- t(rs/rn)
}
# changing order of this
averaged[["CS17"]] <- averaged[["CS17"]][ ,c("Ventral_V", "Ventro-lateral_VL", "Dorsal-lateral_DL", "Dorsal_D")]


# sanity check
#a <- averaged$CS17
#b <- assay(dds_full[["CS17"]], "sf_norm")

# converting to expression set 
exprset <- lapply(averaged, ExpressionSet)

# we exclude genes with more than 25% of the measurements missing
exprset_fil <- lapply(exprset, filter.NA, thres=0.25)

# Since  the  clustering  is  performed  in  Euclidian  space,  the  expression  values  of  genes  were
# standardised to have a mean value of zero and a standard deviation of one.  This step ensures
# that vectors of genes with similar changes in expression are close in Euclidean space
# Importantly, Mfuzz assumes that the given expression data are fully preprocessed including  any  data  normalisation.
# The  function standardise does  not  replace  the  normalisation step (eg RPKN normalization).
exprset_std <- lapply(exprset_fil, standardise)
#a <- t(data.frame(exprset_std[["CS17"]]))

# For fuzzy c-means clustering, the fuzzifier m and the number of clusters c has to be chosen in
# advance. For fuzzifier m, we would like to choose a value which prevents clustering of random
# data. Note, that fuzzy clustering can be tuned in such manner, that random data is not
# clustered. This is a clear advantage to hard clustering (such as k-means), which commonly
# detects clusters even in random data. To achieve this, different options exists: Either the
# function partcoef can be used to test, whether random data is clustered for a particular
# setting of m (see example of partcoef) or a direct estimate can be achieved using a relation
# proposed by Schwaemmle and Jensen 

m1 <- lapply(exprset_std, mestimate)

# clustering - separately for cs

## CS17
CS17_expr_std <- t(data.frame(exprset_std[["CS17"]]))

if (!(file.exists(here("output", "c1_CS17")))) {
  c1_CS17 <- mfuzz(exprset_std[["CS17"]], c = 8, m = m1[["CS17"]])
  saveRDS(c1_CS17, file = here("output", "c1_CS17"))
  
} else{
  c1_CS17 <- readRDS(here("output", "c1_CS17"))
  
}


mfuzz.plot2(exprset_std[["CS17"]],
            cl=c1_CS17,
            mfrow=c(4,4), 
            time.labels = c("V", "VL", "DL", "D"),
            xlab = "Subdomain",
            colo = "fancy",
            min.mem = 0)

mfuzzColorBar(col="fancy",main="Membership value", horizontal = TRUE)
#CS17_mem <- as.data.frame(c1_CS17$membership)
#CS17_mem_fil <- CS17_mem %>% filter(if_any(everything(), ~. > 0.5)) %>% setNames(paste0('clust_', names(.))) %>% rownames_to_column(var = "row")

# getting all genes for each cluster with membership above 0.2
CS17_clus_list <- acore(exprset_std[["CS17"]], c1_CS17, min.acore = 0.2)
# adding names
names(CS17_clus_list) <- paste0("clust_",seq_along(CS17_clus_list))

# adding gene names and creating separate vars
for (i in names(CS17_clus_list)){
  
  # adding gene names
  colnames(CS17_clus_list[[i]])[1] <- "row"
  CS17_clus_list[[i]] <- merge(x = CS17_clus_list[[i]], y = rowData(dds_full[["CS17"]]), by=0)
  
  # # easier to look at
  # name <- paste0("CS17_clus__", i)
  # assign(name, data.frame(CS17_clus_list[[i]]))
}


# cheking individual genes - just curious
plotCounts(dds_full[["CS17"]], gene = "ENSG00000204859", intgroup=c("subdomain"))



# preparing gene lists for ORA

CS_17_fuzzy_entrez <- list()
for (i in names(CS17_clus_list)){
  
  # extract genes to use
  a <- CS17_clus_list[[i]]
  colnames(a)[1] <- "ENSEMBL"
  
  # getting entrez ids
  #b <- mapIds(org.Hs.eg.db, keys = a$ENSEMBL, keytype="ENSEMBL", column = c("ENTREZID", "SYMBOL"))
  entrez <- bitr(a$ENSEMBL, fromType = "ENSEMBL",
                 toType = c("ENTREZID"),
                 OrgDb = org.Hs.eg.db)
  # merging the two 
  CS_17_fuzzy_entrez[[i]] <- merge(a, entrez, by = "ENSEMBL")
  
  # easier to look at
  name <- paste0("CS17_", i)
  assign(name, data.frame(CS_17_fuzzy_entrez[[i]]))

}  
#b <- data.frame(CS_17_fuzzy_entrez[[10]])

a <- CS17_clust_8[which(CS17_clust_8$ENTREZID %in% CS17_fuzzy_kegg_ora_8[CS17_fuzzy_kegg_ora_8$Description == "Hippo signalling pathway"])]

# KEGG oRA

# KEGG pathway ORA analysis - also requires entrez ids.
if (!file.exists(here("output", "cs17_fuzzy_kegg_ora"))){
  cs17_fuzzy_kegg_ora <- list()
  for (i in seq_along(CS_17_fuzzy_entrez)){
    cs17_fuzzy_kegg_ora[[i]] <- enrichKEGG(gene  = CS_17_fuzzy_entrez[[i]]$ENTREZID,
                                organism     = 'hsa',
                                pvalueCutoff = 1,
                                qvalueCutoff = 1,
                                pAdjustMethod = "none",
                                minGSSize = 5)
    
    # to take a closer look
    name <- paste0("CS17_fuzzy_kegg_ora_", i)
    assign(name, data.frame(cs17_fuzzy_kegg_ora[[i]]))
  }
  saveRDS(cs17_fuzzy_kegg_ora, file = here("output", "cs17_fuzzy_kegg_ora"))
  head(cs17_fuzzy_kegg_ora) 
} else{
  cs17_fuzzy_kegg_ora <- readRDS(here("output", "cs17_fuzzy_kegg_ora"))
}

CS17_fuzzy_kegg_mod <-  cs17_fuzzy_kegg_ora
CS17_fuzzy_kegg_mod <- lapply(CS17_fuzzy_kegg_mod, function(a){subset(a, !a$category %in% c("Human Diseases"))})
head(CS17_fuzzy_kegg_mod)

# plotting this
dev.new()
for (i in c(1,2,3,4,5,6,8)){
  print(dotplot(CS17_fuzzy_kegg_mod[[i]], showCategory=42, font.size = 8) + ggtitle(paste0("KEGG ORA dotplot for ", i))) + scale_y_discrete(labels=function(x) str_wrap(x, width=40))
  # + theme(text=element_text(size=20), #change font size of all text
  #goplot(go_ora[[i]]) #+ ggtitle(paste0("GO ORA dotplot for ", i)) + facet_grid(.~.sign)
} 

# gene ontology
if (!file.exists(here("output", "CS17_fuzzy_GO_ora"))){
  CS17_fuzzy_GO_ora <- list()
  for (i in seq_along(CS_17_fuzzy_entrez)){
    CS17_fuzzy_GO_ora[[i]] <- enrichGO(gene  = CS_17_fuzzy_entrez[[i]]$ENTREZID,
                                       OrgDb         = org.Hs.eg.db,
                                       ont           = "ALL",
                                       pAdjustMethod = "BH",
                                       pvalueCutoff  = 0.05,
                                       qvalueCutoff  = 0.05,
                                       readable      = TRUE)
    
    # to take a closer look
    name <- paste0("CS17_fuzzy_GO_ora_", i)
    assign(name, data.frame(CS17_fuzzy_GO_ora[[i]]))
  }
  saveRDS(CS17_fuzzy_GO_ora, file = here("output", "CS17_fuzzy_GO_ora"))
  head(CS17_fuzzy_GO_ora) 
} else{
  CS17_fuzzy_GO_ora <- readRDS(here("output", "CS17_fuzzy_GO_ora"))
}

# plotting this
for (i in c("V_D"   ,    "V_VL"  ,    "D_DL"   ,   "VL_DL"  , "outer_V_D")){
  print(dotplot(CS17_fuzzy_GO_ora[[i]], showCategory=30, font.size = 8) + ggtitle(paste0("GO ORA dotplot for ", i))) + scale_y_discrete(labels=function(x) str_wrap(x, width=40))
  # + theme(text=element_text(size=20), #change font size of all text
  #goplot(go_ora[[i]]) #+ ggtitle(paste0("GO ORA dotplot for ", i)) + facet_grid(.~.sign)
}

# # trying progeny
# 
# # create mats
# counts_CS17 <- list()
# for (i in seq_along(CS17_clus_list)){
#   #i <- 1
#   # create count assays
#   counts_CS17[[i]] <- averaged[["CS17"]] %>% 
#     merge(y = rowData(dds_full[["CS17"]]), by=0) %>%
#     merge(y = CS17_clus_list[[i]], by=c("Row.names")) %>%
#     as.data.frame() %>%
#     group_by(external_gene_name.x) %>%
#     dplyr::slice(which.max(MEM.SHIP)) %>%
#     column_to_rownames(var = "external_gene_name.x") %>%
#     dplyr::select(c("Ventral_V"   ,      "Ventro.lateral_VL", "Dorsal.lateral_DL", "Dorsal_D"))
#   
#   # create design
#   #design[[i]] <- as.data.frame(colData(deseq[[i]])[,c("embryo","subdomain")])  
#   
# }    
# 
# #a <- averaged[["CS17"]]
# pro <- get_progeny(organism = 'human')#, top = 500)
# pro
# 
# sample_acts_pro <- lapply(counts_CS17, run_mlm, net=pro, .source='source', .target='target',
#                           .mor='weight', minsize = 5)
#   
# sample_acts_pro
# a <- sample_acts_pro[[1]]
# 
# for (i in seq_along(sample_acts_pro)){
#   # Transform to wide matrix
#   sample_acts_mat <- sample_acts_pro[[i]] %>%
#     pivot_wider(id_cols = 'condition', names_from = 'source',
#                 values_from = 'score') %>%
#     column_to_rownames('condition') %>%
#     as.matrix()
#   
#   # Scale per feature
#   sample_acts_mat <- scale(sample_acts_mat)
#   
#   # Choose color palette
#   palette_length = 100
#   my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)
#   
#   my_breaks <- c(seq(-3, 0, length.out=ceiling(palette_length/2) + 1),
#                  seq(0.05, 3, length.out=floor(palette_length/2)))
#   
#   # Plot
#   print(pheatmap(sample_acts_mat, 
#                  border_color = NA, 
#                  color=my_color, 
#                  breaks = my_breaks, 
#                  main = paste0("Pathway inference for ", i), 
#                  cluster_rows = FALSE, 
#                  cluster_cols = FALSE)) 
#   
# }

# plot heatmaps
CS17_fuzzy_assay <- list()
col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
h <- list()
# adding row annotation
CS17_kegg_text_list <- lapply(CS17_fuzzy_kegg_mod, function(a){return(a$Description)})
names(CS17_kegg_text_list) <- paste0("clust_",seq_along(CS17_kegg_text_list))

for (i in c("clust_1" ,"clust_2" ,"clust_3" , "clust_5" ,"clust_6" , "clust_8")){
  dev.new()
  i <- "clust_1"
  std <- t(data.frame(exprset_std[["CS17"]]))
  # genes from cluster
  #CS17_fuzzy_assay[[i]] <- averaged[["CS17"]][which(rownames(averaged[["CS17"]]) %in% CS17_clus_list[[i]]$Row.names),]
  CS17_fuzzy_assay[[i]] <- std[which(rownames(std) %in% CS17_clus_list[[i]]$Row.names),]
  # annotation df
  #df <- as.data.frame(colData(deseq[[i]])[,c("embryo","subdomain")])
  
  # plotting heatmaps
  # print(pheatmap(CS17_fuzzy_assay[[i]], cluster_rows=FALSE, show_rownames=FALSE,
  # cluster_cols=FALSE, 
  # #annotation_col=df, 
  # main = i, fontsize_row = 10, fontsize_col = 10, angle_col = 90, color=colorRampPalette(c("navy","lightblue", "white","pink", "red"))(50)))
  if (length(CS17_kegg_text_list[[i]]) == 0) {
    ha <- NULL
  } else {
    #ha <- rowAnnotation(foo = anno_empty(border = FALSE, width = max_text_width(CS17_kegg_text_list[[i]]) + unit(2, "mm")))
    ha <- rowAnnotation(textbox = anno_textbox(
      list(i = 1:nrow(CS17_fuzzy_assay[[i]])), 
      list(i = CS17_kegg_text_list[[i]]), 
      word_wrap = TRUE,
      add_new_line = TRUE,
      gp = gpar(fontsize = 19),
      padding = unit(1, "mm")))
  }
  
  #h[[i]] <- 
    print(Heatmap(CS17_fuzzy_assay[[i]], name = paste0("Cluster ",i), col = col_fun, 
                  cluster_rows = FALSE, cluster_columns = FALSE, show_row_names = FALSE,
                  row_title=paste0("Cluster ",i),heatmap_legend_param = list(direction = "horizontal"),
                  row_title_gp = grid::gpar(fontsize = 15), 
                  column_title_gp = grid::gpar(fontsize = 20),
                  column_names_gp = gpar(fontsize = 20),
                  right_annotation = ha,
                  width = unit(30, "cm"),
                  height = unit(20,"cm"),
                  column_names_rot = 45))
    
}

# adding row annotation
#names(CS17_kegg_text_list) <- names(h) <- paste0("clust_",seq_along(CS17_kegg_text_list))
#names(h)
png(file=here("figs", "heat.png"),width = 2000, height = 2000)
draw(h[[1]] %v% h[[2]] %v% h[[3]] %v% h[[4]] %v% h[[5]] %v% h[[6]] %v% h[[7]] %v% h[[8]], merge_legend = TRUE,  heatmap_legend_side = "bottom")
dev.off()

for (i in seq_along(CS17_kegg_text_list)){
  if (!length(CS17_kegg_text_list[[i]]) == 0){
    dev.new()
    
    print(grid.draw(textbox_grob(CS17_kegg_text_list[[i]])))
  }
}


# trying to add annotation - didnt work
for(i in seq_along(CS17_kegg_text_list)) {
  decorate_annotation("foo", slice = i, {
    grid.rect(x = 0, width = unit(2, "mm"), gp = gpar(fill = i, col = NA), just = "left")
    grid.text(paste(CS17_kegg_text_list, collapse = "\n"), x = unit(4, "mm"), just = "left")
  })
}




## CS16
c1_CS16 <- mfuzz(exprset_std[["CS16"]], c = 16, m = m1[["CS16"]])
mfuzz.plot2(exprset_std[["CS16"]],
            cl=c1_CS16,
            mfrow=c(4,5), 
            time.labels = c("D_inner", "D_mid", "D_outer", 
                            "V_inner", "V_mid", "V_outer"),
            xlab = "Subdomain",
            colo = "fancy",
            min.mem = 0.3, 
            cex.axis = 0.65)

mfuzzColorBar(col="fancy",main="Membership value", horizontal = TRUE)
