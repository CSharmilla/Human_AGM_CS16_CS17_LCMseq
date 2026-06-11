library(here)
library(DESeq2)
library(tibble)
library(decoupleR)
library(OmnipathR)
library(tidyr)
library(Glimma)
library(ggplot2)
library(dplyr)
library(BulkSignalR)

library(CellChat)
library(igraph)
library(ggraph)

source("functions.R")

# reading in data
deseq <- readRDS(here("output", "deseq"))
dds <- readRDS(here("output", "dds"))
results_wald <- readRDS(here("output", "results_wald")) # need stat
results_list <- readRDS(here("output", "results_list"))



# # cleaning up input data
# dds_norm_bulksignalr <- dds
# 
# 
# for (i in names(dds_norm_bulksignalr)){
#   # add normalised counts to assay
#   assays(dds_norm_bulksignalr[[i]])[['sf_norm']] <- counts(dds_norm_bulksignalr[[i]], normalized = TRUE)
# }    
# 
# # averaging replicates
# averaged <- list()
# counts_bulksignalr <- list()
# deg_all <- list()
# 
# for (i in names(dds_norm_bulksignalr)){
#   
#   deg_all[[i]] <- results_wald[[i]] %>% as.data.frame() %>%
#     #filter(padj < 0.05) %>%
#     filter(!external_gene_name == "") %>%
#     group_by(external_gene_name) %>%
#     dplyr::slice(which.max(abs(log2FoldChange))) %>%
#     dplyr::select(log2FoldChange, stat, pvalue, row, external_gene_name) %>%
#     column_to_rownames(var = "external_gene_name")
#   
#   
#   # create count assays
#   counts_bulksignalr[[i]] <- assay(dds_norm_bulksignalr[[i]], "sf_norm")
#   
#   counts_bulksignalr[[i]] <- counts_bulksignalr[[i]] %>%
#     merge(y = rowData(dds[[i]]), by=0) %>%
#     semi_join(y = deg_all[[i]], by=c("row")) %>%
#     column_to_rownames(var = "external_gene_name") %>%
#     dplyr::select(-c(Row.names, row))
#   
#   
#   
#   rs <- rowsum(t(counts_bulksignalr[[i]]), group = dds_norm_bulksignalr[[i]]$subdomain)
#   rn <- as.numeric(table(dds_norm_bulksignalr[[i]]$subdomain))
#   averaged[[i]] <- t(rs/rn)
# }

# # ligand receptor interaction using bulksignalr
# 
# # Building a BSRDataModel object
# bsrdm <- lapply(averaged, function(x) {
#   BSRDataModel(counts = x)
# })
# 
# bsrdm
# 
# bsrdm <- lapply(bsrdm, function(y) {
#   learnParameters(y, quick=TRUE)
# })
# 
# ls(BulkSignalR:::.SignalR)
# 
# bsrinf <- lapply(bsrdm[c("inner_mid_V" ,"mid_outer_V")], function(z) {
#   BSRInference(z, min.cor = 0.3,
#                reference="REACTOME-GOBP")
# })
# saveRDS(bsrinf, file = here("output", "bsrinf"))
# 
# # reducing to best pathways
# bsrinf.redBP <- lapply(bsrinf, function(a){
#   reduceToBestPathway(a)
# })
# 
# # reducing to ligands or receptors
# bsrinf.L <- lapply(bsrinf, function(b){
#   reduceToLigand(b)
# })
# 
# bsrinf.R <- lapply(bsrinf, function(c){
#   reduceToLigand(c)
# })
# 
# LRinter.dataframe <- lapply(bsrinf, function(k){
#   LRinter(k)
# })
# 
# LR_inner_mid_V <- LRinter.dataframe$inner_mid_V
# 
# expr_inner_V  <- averaged[["inner_mid_V"]]
# expr_mid_V <- averaged[["inner_mid_V"]][, 1]
# 
# 
# # Building a BSRSignature object
# # Scoring by ligand-receptor
# 
# bsrsig.redBP <- lapply(bsrinf.redBP, function(d){
#   BSRSignature(d, qval.thres=0.001)
# })
# 
# #scoresLR <- mapply(function(x, y){scoreLRGeneSignatures(bsrdm[c("V_D","V_VL")], bsrsig.redBP, name.by.pathway=FALSE)}
# #                   , bsrdm, bsrsig.redBP, SIMPLIFY = FALSE)
#   
#   
# scoresLR_vd <- scoreLRGeneSignatures(bsrdm[["V_D"]], bsrsig.redBP[["V_D"]],
#                                   name.by.pathway=FALSE
# )  
# 
# scoresLR_vvl <- scoreLRGeneSignatures(bsrdm[["V_VL"]], bsrsig.redBP[["V_VL"]],
#                                      name.by.pathway=FALSE
# )
# 
# simpleHeatmap(scoresLR_vd[1:50, ],
#               hcl.palette="Cividis",
#               pointsize=8)
# 
# simpleHeatmap(scoresLR_vvl[1:50, ],
#               hcl.palette="Cividis",
#               pointsize=8)
# 
# scoresPathway_vd <- scoreLRGeneSignatures(bsrdm[["V_D"]], bsrsig.redBP[["V_D"]],
#                                        name.by.pathway=TRUE
# )
# 
# simpleHeatmap(scoresPathway_vd[1:10, ],
#               hcl.palette="Blue-Red 2",
#               pointsize=8)
# 
# scoresPathway_vvl <- scoreLRGeneSignatures(bsrdm[["V_VL"]], bsrsig.redBP[["V_VL"]],
#                                           name.by.pathway=TRUE
# )
# 
# simpleHeatmap(scoresPathway_vvl[1:10, ],
#               hcl.palette="Blue-Red 2",
#               pointsize=8)

# trying with cellchat
counts_cellchat <- list()
deg_all <- list()

for (i in names(dds)){
  
  # to remove duplicates
  deg_all[[i]] <- results_wald[[i]] %>% as.data.frame() %>%
    #filter(padj < 0.05) %>%
    filter(!external_gene_name == "") %>%
    group_by(external_gene_name) %>%
    dplyr::slice(which.max(abs(log2FoldChange))) %>%
    dplyr::select(log2FoldChange, stat, pvalue, row, external_gene_name) %>%
    column_to_rownames(var = "external_gene_name")
  
  
  # create count assays
  counts_cellchat[[i]] <- assay(dds[[i]], "counts")
  
  counts_cellchat[[i]] <- counts_cellchat[[i]] %>%
    merge(y = rowData(dds[[i]]), by=0) %>%
    semi_join(y = deg_all[[i]], by=c("row")) %>%
    column_to_rownames(var = "external_gene_name") %>%
    dplyr::select(-c(Row.names, row))
  
  # log-normalise with pseudocount
  counts_cellchat[[i]] <- as.matrix(log2(counts_cellchat[[i]] + 1))
  
  
}

# inner mid v 
counts_inner_mid_v <- counts_cellchat$inner_mid_V

meta_inner_mid_v <- data.frame(
  samples = as.factor(colnames(counts_inner_mid_v)),
  sampleType = ifelse(grepl("Inner", colnames(counts_inner_mid_v)),
                      "V_inner", "V_mid"),
  stringsAsFactors = FALSE
)
rownames(meta_inner_mid_v) <- colnames(counts_inner_mid_v)

#Creating CellChat object
cellchat_inner_mid_v <- createCellChat(object = as.matrix(counts_inner_mid_v),
                               meta = meta_inner_mid_v,
                               group.by = "sampleType")

cellchat_inner_mid_v <- setIdent(cellchat_inner_mid_v, ident.use = "sampleType")


# Loading in the mouse ligand-receptor database
showDatabaseCategory(CellChatDB.human)
cellchat_inner_mid_v@DB <- subsetDB(CellChatDB.human, search = c("Secreted Signaling", "ECM_receptor", "Cell-Cell Contact"), key = "annotation")
  #CellChatDB.human
#a <- cellchat_inner_mid_v@DB$interaction

# Subsetting to relevant ligand/receptor genes only
cellchat_inner_mid_v <- subsetData(cellchat_inner_mid_v)

#Checking subset results
cellchat_inner_mid_v@data.signaling

#Identifying over expressed genes/interactions
cellchat_inner_mid_v <- identifyOverExpressedGenes(cellchat_inner_mid_v, do.fast = FALSE, thresh.p = 1)
cellchat_inner_mid_v@var.features
cellchat_inner_mid_v <- identifyOverExpressedInteractions(cellchat_inner_mid_v)

#nrow(cellchat_inner_mid_v@LR$LRsig)
#dim(cellchat_inner_mid_v@data.signaling)


# Computing strength of ligand-receptor interactions
if (!(file.exists(here("output", "cellchat_inner_mid_v")))) {
  cellchat_inner_mid_v <- computeCommunProb(cellchat_inner_mid_v)
} else {
  cellchat_inner_mid_v <- readRDS(here("output", "cellchat_inner_mid_v"))
}

df.net <- subsetCommunication(cellchat_inner_mid_v)
df.net_mid_inner_v <- subsetCommunication(cellchat_inner_mid_v, sources.use = c("V_mid"), targets.use = c("V_inner"))


netVisual_bubble(
  cellchat_inner_mid_v,
  sources.use = "V_mid",
  targets.use = "V_inner",
  signaling = c(df.net_mid_inner_v$pathway_name)
)




# ggplot(df.net_mid_inner_v,
#        aes(x = interaction_name_2,
#            y = pathway_name,
#            fill = prob)) +
#   geom_tile() +
#   scale_fill_gradient(low = "white", high = "red") +
#   theme_bw() +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1),
#     panel.grid = element_blank()
#   ) +
#   labs(
#     x = "Ligand–Receptor Pair",
#     y = "Pathway",
#     fill = "Interaction strength"
#   )


# df <- df.net_mid_inner_v %>%
#   arrange(pathway_name, desc(prob)) %>%
#   mutate(interaction_name_2 = factor(
#     interaction_name_2,
#     levels = unique(interaction_name_2)
#   ))
# df$pathway_name <- factor(df$pathway_name,
#                           levels = unique(df$pathway_name))
# 
# ggplot(df,
#        aes(x = "V_mid → V_inner",
#            y = interaction_name_2,
#            fill = prob)) +
#   geom_tile() +
#   facet_grid(pathway_name ~ ., scales = "free_y", space = "free_y") +
#   scale_fill_gradient(low = "white", high = "red") +
#   theme_bw() +
#   theme(
#     strip.text.y = element_text(angle = 0, face = "bold"),
#     axis.text.x = element_blank(),
#     axis.ticks.x = element_blank()
#   )

df <- df.net_mid_inner_v %>%
  arrange(pathway_name, desc(prob)) %>%
  group_by(pathway_name) %>%
  mutate(interaction_name_2 = factor(
    interaction_name_2,
    levels = unique(interaction_name_2)
  )) %>%
  ungroup()

ggplot(df,
       aes(x = "V_mid → V_inner",
           y = interaction_name_2)) +
  geom_point(aes(size = prob, fill = prob),
             shape = 21,
             color = "black",
             stroke = 0.5,
             position = position_jitter(width = 0.1, height = 0)) +
  scale_fill_gradient(low = "lightblue", high = "red") +
  scale_size(range = c(2, 6)) +
  facet_grid(pathway_name ~ ., scales = "free_y", space = "free_y") +
  theme_bw() +
  theme(
    strip.text.y = element_text(angle = 0, face = "bold"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank()
  ) +
  labs(
    title = "Ligand–receptor communication",
    subtitle = "Ligands: Ventral_mid → Receptors: Ventral_inner",
    y = "Ligand–Receptor Pair",
    size = "Interaction strength",
    fill = "Interaction strength",
    x = ""
  )





# cellchat_inner_mid_v <- computeCommunProbPathway(cellchat_inner_mid_v)
# 
# cellchat_inner_mid_v <- aggregateNet(cellchat_inner_mid_v)
# 
# cellchat_inner_mid_v <- netAnalysis_computeCentrality(cellchat_inner_mid_v, slot.name = "netP")
# 
# netAnalysis_signalingRole_heatmap(cellchat_inner_mid_v, signaling = c(df.net_mid_inner_v$pathway_name))

run_cellchat_pipeline <- function(counts_mat, name_label, label1, label2) {
  
  # counts_mat <- counts_cellchat$inner_mid_D 
  # name_label <- "inner_mid_D"
  # label1 <- "D_Inner" 
  # label2 <- "D_Mid"
  
  meta <- data.frame(
    samples = colnames(counts_mat),
    sampleType = ifelse(grepl(label1, colnames(counts_mat)), label1, label2),
    stringsAsFactors = FALSE
  )
  rownames(meta) <- colnames(counts_mat)
  
  cellchat <- createCellChat(object = as.matrix(counts_mat),
                             meta = meta,
                             group.by = "sampleType")
  
  cellchat <- setIdent(cellchat, ident.use = "sampleType")
  
  cellchat@DB <- subsetDB(CellChatDB.human,
                          search = c("Secreted Signaling",
                                     "ECM_receptor",
                                     "Cell-Cell Contact"),
                          key = "annotation")
  
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat, do.fast = FALSE, thresh.p = 1)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(cellchat)
  
  df <- subsetCommunication(cellchat,
                            sources.use = label1, #unique(meta$sampleType)[1],
                            targets.use = label2) #unique(meta$sampleType)[2])
  
  df_plot <- df %>%
    arrange(pathway_name, desc(prob)) %>%
    group_by(pathway_name) %>%
    mutate(interaction_name_2 = factor(interaction_name_2,
                                       levels = unique(interaction_name_2))) %>%
    ungroup()
  
  p <- ggplot(df_plot,
              aes(x = 1,
                  y = interaction_name_2)) +
    geom_point(aes(fill = prob),
               shape = 21,
               size = 4,            # fixed size
               color = "black",
               stroke = 0.5) +
    scale_fill_gradient(low = "lightblue", high = "red") +
    scale_x_continuous(breaks = NULL) +
    facet_grid(pathway_name ~ ., scales = "free_y", space = "free_y") +
    theme_bw() +
    theme(
      strip.text.y = element_text(angle = 0, face = "bold"),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid = element_blank()
    ) +
    labs(
      title = name_label,
      subtitle = paste(label1, "→", label2),
      y = "Ligand–Receptor Pair",
      fill = "Interaction strength",
      x = ""
    )
    
    
    # ggplot(df_plot,
    #           aes(x = paste0(name_label),
    #               y = interaction_name_2)) +
    # geom_point(aes(size = prob, 
    #                fill = prob),
    #            shape = 21,
    #            color = "black",
    #            stroke = 0.5,
    #            position = position_jitter(width = 0.1, height = 0)) +
    # scale_fill_gradient(low = "lightblue", high = "red") +
    # scale_size(range = c(2, 6)) +
    # facet_grid(pathway_name ~ ., scales = "free_y", space = "free_y") +
    # theme_bw() +
    # theme(
    #   strip.text.y = element_text(angle = 0, face = "bold"),
    #   axis.text.x = element_blank(),
    #   axis.ticks.x = element_blank(),
    #   panel.grid = element_blank()
    # ) +
    # labs(
    #   title = name_label,
    #   subtitle = paste(label1, "→", label2),
    #   y = "Ligand–Receptor Pair",
    #   fill = "Interaction strength",
    #   size = "Interaction strength",
    #   x = ""
    # )
    # 
  print(p)
  
  return(list(cellchat = cellchat, plot = p))
}

res_mid_outer_V <- run_cellchat_pipeline(counts_cellchat$mid_outer_V, "mid_outer_V", "V_Mid", "V_Outer")
res_mid_inner_D <- run_cellchat_pipeline(counts_cellchat$inner_mid_D, "inner_mid_D", "D_Mid", "D_Inner")
res_inner_mid_V <- run_cellchat_pipeline(counts_cellchat$inner_mid_V, "inner_mid_V", "V_Mid", "V_Inner")

