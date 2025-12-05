library(here)
library(biomaRt)
library(tximport)
library(readr)
library(stringr)

meta <- read.table(here("support", "LCMseq_CS16_17_metadata.txt"), header = TRUE, sep = "\t")
meta$matname <- gsub(x = meta$raw.file, pattern = "_1.fastq.gz", replacement = "")
meta$title <- make.unique(meta$title, sep = "_")

# reading in quant files
quant_files <- list.files(path = 'salmon/nodecoy_k21_quant', pattern = "_cut", full.names = TRUE)
quant_files
quant_files <- paste0(quant_files, "/quant.sf")
quant_files
names(quant_files) <- basename(gsub("/quant.sf", "", quant_files))
names(quant_files) <- gsub("nodecoy_k21_quant_", "", names(quant_files))
names(quant_files) <- gsub("_cut", "", names(quant_files))
quant_files


# creating tx2gene to use tximport
## change version to see which can annotate all genes; mirror will not work with version option.
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", version = 113) 
filters <- data.frame(listFilters(ensembl))
# trying different versions. 
tx2gene.111 <- getBM(attributes= c("ensembl_transcript_id_version", "ensembl_gene_id"), mart= ensembl)
tx2gene.113 <- getBM(attributes= c("ensembl_transcript_id_version", "ensembl_gene_id"), mart= ensembl) # has more

txi <- tximport(files = quant_files, type = "salmon", txIn = TRUE, txOut = FALSE, countsFromAbundance = "no", 
                tx2gene = tx2gene.113) #, ignoreTxVersion = TRUE)

txi_rename <- txi # just to avoid running tximport again

txi_rename <- lapply(txi_rename[1:3], function(a){
  colnames(a) <- meta[match(colnames(a), meta$matname), "title"]
  return(a)
  })
txi_rename$countsFromAbundance <- "no" # adding this back because it is required later for deseq2

# saving txi
saveRDS(txi_rename, file = here("output", "txi_rename"))

# saving txi for each CS
txi_CS16 <- lapply(txi_rename[1:3], function(a){
  a <- a[,grep("CS16_", colnames(a))]
  return(a)
})
txi_CS16$countsFromAbundance <- "no" # adding this back because it is required later for deseq2

txi_CS17 <- lapply(txi_rename[1:3], function(a){
  a <- a[,grep("CS17_", colnames(a))]
  return(a)
})
txi_CS17$countsFromAbundance <- "no" # adding this back because it is required later for deseq2

txi_cs <- list(CS16 = txi_CS16, CS17 = txi_CS17)
saveRDS(txi_cs, file = here("output", "txi_cs"))

# saving txi for each subset of comparison
test <- txi_rename$counts


# also get list of genes to save for later
genes_list <- getBM(attributes= c("ensembl_gene_id", "external_gene_name"), mart= ensembl)
saveRDS(genes_list, file = here("output", "genes_list"))




# cleaning up the meta file also then saving
# title to row names
rownames(meta) <- meta$title

# removing square brackets from meta
meta$characteristics..subdomain <- gsub("[][]","", meta$characteristics..subdomain)
# removing underscore 
meta$characteristics..subdomain <- gsub(" ","_", meta$characteristics..subdomain)
# removing brakcets
meta$characteristics..subdomain <- gsub("\\*|\\(|\\)","", meta$characteristics..subdomain)

# changing names of some columns
colnames(meta)[7] <- "subdomain"
colnames(meta)[5] <- "embryo"
colnames(meta)[6] <- "stage"

# saving meta
saveRDS(meta, file = here("output", "meta"))

# meta for CS16
meta_CS16 <- meta[meta$stage == "CS16", ]

# meta for CS17
meta_CS17 <- meta[meta$stage == "CS17", ]

meta_cs <- list(CS16 = meta_CS16, CS17 = meta_CS17)
saveRDS(meta_cs, file = here("output", "meta_cs"))

# meta for every subset - because running deseq2 separately for each one due to how dispersed samples are
subsets <- list(V_D = meta[which(meta$subdomain == "Ventral_V" | meta$subdomain == "Dorsal_D"),], 
                V_VL = meta[which(meta$subdomain == "Ventral_V" | meta$subdomain == "Ventro-lateral_VL"),],
                VL_DL = meta[which(meta$subdomain == "Ventro-lateral_VL" | meta$subdomain == "Dorsal-lateral_DL"),],
                DL_D = meta[which(meta$subdomain == "Dorsal-lateral_DL" | meta$subdomain == "Dorsal_D"),],
                D_DL = meta[which(meta$subdomain == "Dorsal_D" | meta$subdomain == "Dorsal-lateral_DL"),],
                
                inner_V_D = meta[which(meta$subdomain == "Ventral_Inner_V_Inner" | meta$subdomain == "Dorsal_Inner_D_Inner"),],
                mid_V_D = meta[which(meta$subdomain == "Ventral_Mid_V_Mid" | meta$subdomain == "Dorsal_Mid_D_Mid"),],
                outer_V_D = meta[which(meta$subdomain == "Ventral_Outer_V_Outer" | meta$subdomain == "Dorsal_Outer_D_Outer"),],
                inner_mid_V = meta[which(meta$subdomain == "Ventral_Inner_V_Inner" | meta$subdomain == "Ventral_Mid_V_Mid"),],
                mid_outer_V = meta[which(meta$subdomain == "Ventral_Mid_V_Mid" | meta$subdomain == "Ventral_Outer_V_Outer"),],
                inner_mid_D = meta[which(meta$subdomain == "Dorsal_Inner_D_Inner" | meta$subdomain == "Dorsal_Mid_D_Mid"),],
                mid_outer_D = meta[which(meta$subdomain == "Dorsal_Mid_D_Mid" | meta$subdomain == "Dorsal_Outer_D_Outer"),])

# saving meta for every subset
saveRDS(subsets, file = here("output", "subsets"))

# sanity check
# V_D = meta[which(meta$subdomain == "Ventral_V" | meta$subdomain == "Dorsal_D"),]
# V_VL = meta[which(meta$subdomain == "Ventral_V" | meta$subdomain == "Ventro-lateral_VL"),]
# D_DL = meta[which(meta$subdomain == "Dorsal_D" | meta$subdomain == "Dorsal-lateral_DL"),]
# VL_DL = meta[which(meta$subdomain == "Ventro-lateral_VL" | meta$subdomain == "Dorsal-lateral_DL"),]
# inner_V_D = meta[which(meta$subdomain == "Ventral_Inner_V_Inner" | meta$subdomain == "Dorsal_Inner_D_Inner"),]
# mid_V_D = meta[which(meta$subdomain == "Ventral_Mid_V_Mid" | meta$subdomain == "Dorsal_Mid_D_Mid"),]
# outer_V_D = meta[which(meta$subdomain == "Ventral_Outer_V_Outer" | meta$subdomain == "Dorsal_Outer_D_Outer"),]

