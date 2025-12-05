

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
library(pathview)
library(ReactomePA) # reactome pathway analysis
library(SBGNview)
library(EnhancedVolcano)

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


# filtering and plotting
results_list_fil <- list()
results_wald_fil <- list()
selected <- list()
for (i in names(results_list)){
  i <- "V_VL"
  # just in case i need this at some point 
  results_list_fil[[i]] <- dplyr::filter(results_list[[i]], padj < 0.05) 
  results_wald_fil[[i]] <- dplyr::filter(data.frame(results_wald[[i]]), padj < 0.05)
  
  # plotting top genes - always check with the results lists to make sure nothing funny is happening
  #plot_top_genes(results_list_fil[[i]], deseq[[i]], dds[[i]], vst[[i]], i)
  
  #plot_top_genes_up(results_list_fil[[i]], deseq[[i]], dds[[i]], vst[[i]], i)
  #plot_top_genes_down(results_list_fil[[i]], deseq[[i]], dds[[i]], vst[[i]], i)
  
  #file_name = paste("DE_volcano_", i, ".tiff", sep="")
  #tiff(here("figs",file_name), width = 1500, height = 1500)
  print(EnhancedVolcano(results_list[[i]],
                        lab = results_list[[i]]$external_gene_name,
                        x = 'log2FoldChange',
                        y = 'padj',
                        pCutoff = 0.05,
                        FCcutoff = 1, 
                        title = paste0("Differentially expressed genes for ", i),
                        drawConnectors = TRUE,
                        #widthConnectors = 0.75,
                        lengthConnectors = unit(0.006, 'npc'),
                        max.overlaps = 25,
                        labSize = 3.5,
                        caption = bquote(~Log[2]~ "fold change cutoff, 1; adjp-value cutoff, 0.05"),
                        subtitle = ""))
  #dev.off()
  # separate data frames so its easier to see
  name <- paste0(i, "_apeglm_results")
  assign(name, results_list_fil[[i]])
}



# using clusterprofiler

genelist_ensembl <- list() # all genes used for gsea
genelist_entrez <- list()

fil_genelist_ensembl <- list() # filtered genes used for ORA
fil_genelist_entrez <- list()

# preparing gene lists
for (i in names(results_list)){
  
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
    group_by(ENSEMBL) %>% # grouping by emsembl ids
    dplyr::slice(which.max(abs(log2FoldChange))) %>% # Return a Single Row with the Maximum Value for the Group
    arrange(desc(log2FoldChange)) %>%
    pull(log2FoldChange, name=ENSEMBL) # creating a named vector with logfc
  
  # filtered
  fil_genelist_ensembl[[i]] <- a %>% filter(padj < 0.05) %>%
    group_by(ENSEMBL) %>% # grouping by emsembl ids
    dplyr::slice(which.max(abs(log2FoldChange))) %>%
    arrange(desc(log2FoldChange)) %>%
    pull(log2FoldChange, name=ENSEMBL) # creating a named vector with logfc
  
  
  # list of entrez gene lists
  genelist_entrez[[i]] <- a %>% filter(!is.na(padj)) %>%
    group_by(ENTREZID) %>% # grouping by emsembl ids
    dplyr::slice(which.max(abs(log2FoldChange))) %>% # Return a Single Row with the Maximum Value for the Group
    arrange(desc(log2FoldChange)) %>%
    pull(log2FoldChange, name=ENTREZID) # creating anamed vector with logfc
  
  # filtered
  fil_genelist_entrez[[i]] <- a %>% filter(padj < 0.05) %>%
    group_by(ENTREZID) %>% # grouping by emsembl ids
    dplyr::slice(which.max(abs(log2FoldChange))) %>% # Return a Single Row with the Maximum Value for the Group
    arrange(desc(log2FoldChange)) %>%
    pull(log2FoldChange, name=ENTREZID) # creating anamed vector with logfc
  
}  


# Gene ontology GSEA

if (!file.exists(here("output", "gse"))){
  gse <- list()
  for (i in names(genelist_ensembl)){
    
    # running gsego using ensembl ids
    gse[[i]] <- gseGO(geneList=genelist_ensembl[[i]],
                      ont ="ALL",
                      keyType = "ENSEMBL",
                      #nPerm = 10000,
                      #minGSSize = 3,
                      #maxGSSize = 800,
                      pvalueCutoff = 0.01,
                      verbose = TRUE,
                      OrgDb = org.Hs.eg.db,
                      pAdjustMethod = "none")
    
    # to take a closer look
    name <- paste0("go_gsea_", i)
    assign(name, data.frame(gse[[i]]))
    
    #print(dotplot(gse[[i]], showCategory=15, split=".sign", font.size=7) + ggtitle(paste0("GSEA dotplot for ", i)) + facet_grid(.~.sign))
    #ridgeplot(gse) + labs(x = "enrichment distribution")
  }
  saveRDS(gse, file = here("output", "gse"))
  
}else{
  gse <- readRDS(here("output", "gse"))
  
}


# separate loop since gsego takes some time to run
for (i in names(gse)){
  print(dotplot(gse[[i]], showCategory=15, split=".sign", font.size=8) + ggtitle(paste0("GO GSEA dotplot for ", i)) + facet_grid(.~.sign))
  
}

# GO ORA 
# need entrez ids for this
if (!file.exists(here("output", "go_ora"))){
 
  go_ora <- list()
  for (i in names(fil_genelist_entrez)){
    go_ora[[i]] <- enrichGO(gene        = names(fil_genelist_entrez[[i]]),
                         #universe      = names(fil_genelist_entrez[[i]]),
                         OrgDb         = org.Hs.eg.db,
                         ont           = "ALL",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05,
                         qvalueCutoff  = 0.05,
                         readable      = TRUE)
    
    # to take a closer look
    name <- paste0("go_ora_", i)
    assign(name, data.frame(go_ora[[i]]))
  
  }
  saveRDS(go_ora, file = here("output", "go_ora"))
} else{
  go_ora <- readRDS(here("output", "go_ora"))
}

# plotting this
for (i in c("V_D"   ,    "V_VL"  ,    "D_DL"   ,   "VL_DL"  , "outer_V_D")){
  print(dotplot(go_ora[[i]], showCategory=30, font.size = 8) + ggtitle(paste0("GO ORA dotplot for ", i))) + scale_y_discrete(labels=function(x) str_wrap(x, width=40))
  # + theme(text=element_text(size=20), #change font size of all text
  #goplot(go_ora[[i]]) #+ ggtitle(paste0("GO ORA dotplot for ", i)) + facet_grid(.~.sign)
}

# this doesnt work
#goplot(go_ora[["V_D"]])



# KEGG enrichment analysis

# looking at KEGG pathways first

# KEGG pathway ORA analysis - also requires entrez ids.

if (!file.exists(here("output", "kegg_ora"))){
  kegg_ora <- list()
  for (i in names(fil_genelist_entrez)){
    kegg_ora[[i]] <- enrichKEGG(gene  = names(fil_genelist_entrez[[i]]),
                          organism     = 'hsa',
                          pvalueCutoff = 0.05)
    
    # to take a closer look
    name <- paste0("kegg_ora_", i)
    assign(name, data.frame(kegg_ora[[i]]))
  }
  saveRDS(kegg_ora, file = here("output", "kegg_ora"))
  head(kegg_ora) 
} else{
  kegg_ora <- readRDS(here("output", "kegg_ora"))
}

# plotting this
for (i in c("V_VL"  ,    "D_DL"   ,   "VL_DL")){
  print(dotplot(kegg_ora[[i]], showCategory=30, font.size = 8) + ggtitle(paste0("KEGG ORA dotplot for ", i))) + scale_y_discrete(labels=function(x) str_wrap(x, width=40))
  # + theme(text=element_text(size=20), #change font size of all text
  #goplot(go_ora[[i]]) #+ ggtitle(paste0("GO ORA dotplot for ", i)) + facet_grid(.~.sign)
}


# browse KEGG pathways
browseKEGG(kegg_ora[["V_D"]], 'hsa04630')
browseKEGG(kegg_ora[["VL_DL"]], 'hsa04110')

# looking at different pathways


# V_VL
pathview_mod("V_VL", "hsa04144")

# VL_DL
pathview_mod("VL_DL", "hsa04152")
pathview_mod("VL_DL", "hsa03082")
pathview_mod("VL_DL", "hsa05221")


# D_DL
pathview_mod("D_DL", "hsa05020")
pathview_mod("D_DL", "hsa04140")
pathview_mod("D_DL", "hsa04141")


# KEGG pathway gsea 

if (!file.exists(here("output", "kegg_gsea"))){
  kegg_gsea <- list()
  
  # requires entrez ids
  for (i in names(genelist_entrez)){
    #i <- "V_D"
    kegg_gsea[[i]] <- gseKEGG(geneList = genelist_entrez[[i]],
                            organism     = 'hsa',
                            minGSSize    = 20,
                            pvalueCutoff = 0.05,
                            pAdjustMethod = "none",
                            verbose      = FALSE)
    
    # to take a closer look
    name <- paste0("kegg_gsea_", i)
    assign(name, data.frame(kegg_gsea[[i]]))
  }
  saveRDS(kegg_gsea, file = here("output", "kegg_gsea"))
} else{
  kegg_gsea <- readRDS(here("output", "kegg_gsea"))
  
}

# plotting this
for (i in names(kegg_gsea)){
  print(dotplot(kegg_gsea[[i]], showCategory=30, split=".sign", font.size=8) + ggtitle(paste0("KEGG GSEA dotplot for ", i)) + facet_grid(.~.sign))

}

# hippo pathway
pathview_mod("V_VL", "hsa04520")
pathview_mod("VL_DL", "hsa04520")
pathview_mod("DL_D", "hsa04520")

pathview_mod("V_D", "hsa04520")

# NFKB pathway
pathview_mod("V_VL", "hsa04064")
pathview_mod("VL_DL", "hsa04064")
pathview_mod("DL_D", "hsa04064")
pathview_mod("V_D", "hsa04064")

# TNF signalling pathway
pathview_mod("V_VL", "hsa04668")
pathview_mod("VL_DL", "hsa04668")
pathview_mod("DL_D", "hsa04668")
pathview_mod("V_D", "hsa04668")

# Toll-like receptor signaling pathway
pathview_mod("V_VL", "hsa04620")
pathview_mod("VL_DL", "hsa04620")
pathview_mod("DL_D", "hsa04620")

# Signaling pathways regulating pluripotency of stem cells
pathview_mod("V_VL", "hsa04550")
pathview_mod("VL_DL", "hsa04550")
pathview_mod("DL_D", "hsa04550")

# Regulation of actin cytoskeleton
pathview_mod("V_VL", "hsa04510")
pathview_mod("VL_DL", "hsa04510")
pathview_mod("DL_D", "hsa04510")

pathview_mod("V_D", "hsa04810")


# Cell adhesion molecules
pathview_mod("inner_mid_V", "hsa04514")
pathview_mod("mid_outer_V", "hsa04514")

pathview_mod("inner_mid_D", "hsa04514")
pathview_mod("mid_outer_D", "hsa04514")

#cAMP pathway
pathview_mod("inner_mid_V", "hsa04024")
pathview_mod("mid_outer_V", "hsa04024")

pathview_mod("inner_mid_D", "hsa04024")
pathview_mod("mid_outer_D", "hsa04024")





# KEGG module ORA - not sure about this...
if (!file.exists(here("output", "modkegg_ora"))){
  modkegg_ora <- list()
  for (i in names(fil_genelist_entrez)){
    modkegg_ora[[i]] <- enrichMKEGG(gene = names(fil_genelist_entrez[[i]]),
                                organism = 'hsa',
                                minGSSize = 10,
                                pvalueCutoff = 1,
                                qvalueCutoff = 1)
    
    # to take a closer look
    name <- paste0("modkegg_ora_", i)
    assign(name, data.frame(modkegg_ora[[i]]))
  }
  saveRDS(modkegg_ora, file = here("output", "modkegg_ora"))
} else{
  modkegg_ora <- readRDS(here("output", "modkegg_ora"))
}

# KEGG module GSEA
if (!file.exists(here("output", "modkegg_gsea"))){
  modkegg_gsea <- list()
  for (i in names(genelist_entrez)){
    modkegg_gsea[[i]] <- gseMKEGG(geneList = genelist_entrez[[i]],
                              organism = 'hsa',
                              pvalueCutoff = 1)
    
    # to take a closer look
    name <- paste0("modkegg_gsea_", i)
    assign(name, data.frame(modkegg_gsea[[i]]))
  }
  saveRDS(modkegg_gsea, file = here("output", "modkegg_gsea"))
} else{
  modkegg_gsea <- readRDS(here("output", "modkegg_gsea"))
}

dotplot(mkk_gsea[[1]])


# wikipathways

# gsea
if (!file.exists(here("output", "wikip_gsea"))){
  wikip_gsea <- list()
  
  for (i in names(genelist_entrez)){
    wikip_gsea[[i]] <- gseWP(genelist_entrez[[i]], 
                          pvalueCutoff = 0.05, 
                          organism = "Homo sapiens", 
                          pAdjustMethod = "none")
    
    # to take a closer look
    name <- paste0("wikip_gsea_", i)
    assign(name, data.frame(wikip_gsea[[i]]))
  
  }
  saveRDS(wikip_gsea, file = here("output", "wikip_gsea"))
} else{
  wikip_gsea <- readRDS(here("output", "wikip_gsea"))
}

# cant be visualised using pathview

# ora
if (!file.exists(here("output", "wikip_ora"))){
  wikip_ora <- list()
  
  for (i in names(fil_genelist_entrez)){
    wikip_ora[[i]] <- enrichWP(names(fil_genelist_entrez[[i]]), 
                             pvalueCutoff = 0.05, 
                             organism = "Homo sapiens", 
                             pAdjustMethod = "none")
    
    # to take a closer look
    name <- paste0("wikip_ora_", i)
    assign(name, data.frame(wikip_ora[[i]]))
    
  }
  saveRDS(wikip_ora, file = here("output", "wikip_ora"))
} else{
  wikip_ora <- readRDS(here("output", "wikip_ora"))
}


# reactome pathway analysis


# ORA
if (!file.exists(here("output", "react_ora"))){
  react_ora <- list()
  
  for (i in names(fil_genelist_entrez)){
    react_ora[[i]] <- enrichPathway(names(fil_genelist_entrez[[i]]), 
                                    #universe = names(fil_genelist_entrez[[i]]),
                               pvalueCutoff = 0.5, 
                               #qvalueCutoff = 1,
                               organism = "human",
                               readable=TRUE)
    
    # to take a closer look
    name <- paste0("react_ora_", i)
    assign(name, data.frame(react_ora[[i]]))
    
  }
  saveRDS(react_ora, file = here("output", "react_ora"))
} else{
  react_ora <- readRDS(here("output", "react_ora"))
}

# plotting this
for (i in c("V_VL"  ,    "D_DL"   ,   "VL_DL", "inner_V_D")){
  print(dotplot(react_ora[[i]], showCategory=30, font.size=7) + ggtitle(paste0("Reactome ORA dotplot for ", i)))
  
}

# visualising - this is not very useful, hard to read
viewPathway("Signaling by NOTCH", 
            readable = TRUE, 
            foldChange = fil_genelist_entrez[["VL_DL"]])


# SBGNview visualisation
data("sbgn.xmls")
data("pathways.info", "pathways.stats")

is.reactome <- pathways.info[,"sub.database"]== "reactome"

reactome.ids <- pathways.info[is.reactome ,"pathway.id"]

#input.pathways <- findPathways("R-HSA-69620")
SBGNview.obj <- SBGNview(
  gene.data = fil_genelist_entrez[["VL_DL"]], 
  gene.id.type = "entrez",
  input.sbgn = "R-HSA-1912408",
  output.file = "VL_DL", 
  output.formats =  c("svg"),
  key.pos = "topright"#,
  #min.gene.value = round(max(abs(fil_genelist_entrez[["V_VL"]]))),
  #max.gene.value = round(max(abs(fil_genelist_entrez[["V_VL"]])))
) 
print(SBGNview.obj)

edox <- setReadable(react_ora[["D_DL"]], 'org.Hs.eg.db', 'ENTREZID')
cnetplot(edox, foldChange=fil_genelist_entrez[["VL_DL"]], showCategory = c("TRAF6 mediated NF-kB activation",
                                                                           "RIP-mediated NFkB activation via ZBP1",
                                                                           "TNF signaling",
                                                                           "NIK-->noncanonical NF-kB signaling",
                                                                           "Signaling by NOTCH",
                                                                           "Pre-NOTCH Transcription and Translation"))

cnetplot(edox, foldChange=fil_genelist_entrez[["VL_DL"]], showCategory = c("Signaling by NOTCH",
                                                                           "Pre-NOTCH Transcription and Translation",
                                                                           "Pre-NOTCH Expression and Processing",
                                                                           "Negative regulation of NOTCH4 signaling",
                                                                           "NOTCH1 Intracellular Domain Regulates Transcription"))
heatplot(edox, foldChange=fil_genelist_entrez[["VL_DL"]], showCategory=5)

cnetplot(edox, foldChange=fil_genelist_entrez[["V_VL"]], showCategory = c("SMAD2/SMAD3:SMAD4 heterotrimer regulates transcription",
                                                                           "Signaling by TGFB family members",
                                                                           "Signaling by TGF-beta Receptor Complex"))

cnetplot(edox, foldChange=fil_genelist_entrez[["D_DL"]], showCategory = c("M Phase",
                                                                          "Transcriptional Regulation by TP53",
                                                                          "Cell Cycle Checkpoints"))

# gsea

if (!file.exists(here("output", "react_gsea"))){
  react_gsea <- list()
  
  for (i in names(genelist_entrez)){
    react_gsea[[i]] <- gsePathway(genelist_entrez[[i]], 
                             pvalueCutoff = 0.6, 
                             pAdjustMethod = "BH",
                             organism = "human", 
                             verbose = FALSE)
    
    # to take a closer look
    name <- paste0("react_gsea_", i)
    assign(name, data.frame(react_gsea[[i]]))
    
  }
  saveRDS(react_gsea, file = here("output", "react_gsea"))
} else{
  react_gsea <- readRDS(here("output", "react_gsea"))
}

# dotplots
for (i in c("mid_V_D")){
  print(dotplot(react_gsea[[i]], showCategory=17, split=".sign", font.size=8) + ggtitle(paste0("Reactome GSEA dotplot for ", i)) + facet_grid(.~.sign))
  
}

# visualising 
viewPathway("TCF dependent signaling in response to WNT", 
            readable = TRUE, 
            foldChange = fil_genelist_entrez[["VL_DL"]])

viewPathway("E2F mediated regulation of DNA replication", 
            readable = TRUE, 
            foldChange = fil_genelist_entrez[["VL_DL"]])



# SBGNview visualisation
data("sbgn.xmls")
data("pathways.info", "pathways.stats")

is.reactome <- pathways.info[,"sub.database"]== "reactome"

reactome.ids <- pathways.info[is.reactome ,"pathway.id"]

#input.pathways <- findPathways("R-HSA-69620")
SBGNview.obj <- SBGNview(
  gene.data = fil_genelist_entrez[["VL_DL"]], 
  gene.id.type = "entrez",
  input.sbgn = "R-HSA-1912408",
  output.file = "TEST", 
  output.formats =  c("png")
) 
print(SBGNview.obj)


edox <- setReadable(react_gsea[["mid_V_D"]], 'org.Hs.eg.db', 'ENTREZID')
cnetplot(edox, foldChange=fil_genelist_entrez[["mid_V_D"]], showCategory = c("Mitochondrial translation termination",
                                                                          "Mitochondrial translation elongation",
                                                                          "Mitochondrial translation initiation",
                                                                          "Mitochondrial translation",
                                                                          "Receptor Mediated Mitophagy"))


