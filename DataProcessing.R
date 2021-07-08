# Data processing

# 1. Transcriptome data

# Transcriptomics data can be downloaded from GDC data portal (https://portal.gdc.cancer.gov/)

library(openxlsx)

# Convertin FPKM values to TPM

#Convert FPKM to TPM

fpkmToTpm <- function(fpkm){exp(log(fpkm) - log(sum(fpkm)) + log(1e6))}

fpkm_tables = list.files(path = "/path/to/all/FPKM/files/downloaded/from/GDC",
                         pattern = "*.txt", recursive = FALSE, full.names = TRUE)
c.types = c()
for (i in 1:length(fpkm_tables)) {
  fpkm_table = read.delim(fpkm_tables[i], header = TRUE)
  rownames(fpkm_table) = as.character(fpkm_table$ensembl_id)
  fpkm_table = fpkm_table[,-1]
  tpm_mm = c()
  c = sub("TCGA-","",strsplit(basename(fpkm_tables[i]),"_")[[1]][1])
  c.types = c(c.types,c)
  for (coln in 1:length(colnames(fpkm_table))) {
    column = fpkmToTpm(fpkm_table[,coln])
    tpm_mm = cbind(tpm_mm, column)
    colnames(tpm_mm)[coln] = colnames(fpkm_table)[coln]}
  rownames(tpm_mm) = rownames(fpkm_table)
  tpm_mm = as.data.frame(tpm_mm)
  save(tpm_mm, file = paste0("~/",c,"_tpm.RData"))}

#TPM table for unique sample names
for (c in c.types) {
  load(paste0("~/",c,"_tpm.RData"))
  sample_names = colnames(tpm_mm)
  sample_names = sapply(sample_names, function(x) substr(x,1,15))
  sample_names =unique(sample_names)
  tpm_samples = c()
  for (cname in 1:length(sample_names)) {
    pattern = paste0("^",sample_names[cname])
    data = as.data.frame(tpm_mm[, grep(pattern, colnames(tpm_mm))])
    if(dim(data)[2] > 1){
      RM = rowMeans(data[])
      tpm_samples = cbind(tpm_samples, RM)
      a = which(colnames(tpm_samples) == "RM")
      colnames(tpm_samples)[a] = sample_names[cname]}
    else if (dim(data)[2] == 1){
      RM = data[,1]
      tpm_samples = cbind(tpm_samples, RM)
      a = which(colnames(tpm_samples) == "RM")
      colnames(tpm_samples)[a] = sample_names[cname]}}
  gene_ids = c(as.character(rownames(tpm_mm)))
  gene_ids_w_f = sapply(gene_ids, function(x) substr(x,1,15))
  rownames(tpm_samples) = gene_ids_w_f
  sample_ids = as.character(colnames(tpm_samples))
  sample_ids_new = sapply(sample_ids, function(x) gsub("\\.","-",x))
  colnames(tpm_samples) = sample_ids_new
  save(tpm_samples, file = paste0("~/",c,"_Samples_tpm.RData"))}

#Separate different tissue types within the data
human_genes = read.delim("~/human_genes.txt", header = TRUE) # human_genes.txt file can de downloaded from the repository
for (c in c.types) {
  load(paste0("~/",c,"_Samples_tpm.RData"))
  mt_genes = c(as.character(human_genes[human_genes$Chromosome.scaffold.name == "MT",]$Gene.stable.ID)) #mit. genes were removed
  tpm_samples = tpm_samples[!(rownames(tpm_samples) %in% mt_genes),]
  tpm_samples = as.data.frame(tpm_samples)
  primary_samples = c()
  for(col in colnames(tpm_samples)){
    if(strsplit(col,"-")[[1]][4] == "01" | strsplit(col,"-")[[1]][4] == "03" |
       strsplit(col,"-")[[1]][4] == "09"){primary_samples = c(primary_samples,col)}}
  tpm_samples = tpm_samples[,colnames(tpm_samples) %in% primary_samples]
  save(tpm_samples, file = paste0("~/",c,"_tpm_primary.RData"))}

#Convert GeneIDs to Gene names (if there are more than one genes mapped to one ID, take the mean)
#Filter low expressing genes, genes having 0 TPM for all samples were removed
human_genes = read.delim("~/human_genes.txt", header = TRUE) # human_genes.txt file can de downloaded from the repository
c.types = c()
#read the first lines
for (i in 1:length(fpkm_tables)) {
  c = sub("TCGA-","",strsplit(basename(fpkm_tables[i]),"_")[[1]][1])
  c.types = c(c.types,c)}
library(dplyr)
for (type in c.types) {
  assign("tpm_primary",get(load(paste0("~/",type,"_tpm_primary.RData"))))
  tpm_primary = as.data.frame(tpm_primary)
  tpm_primary$Gene.stable.ID = rownames(tpm_primary)
  tpm_primary = merge(tpm_primary, human_genes[,c(1,3)], by = "Gene.stable.ID")
  tpm_primary = tpm_primary[,-1]
  tpm_primary <- tpm_primary %>% group_by(Gene.name) %>% summarise_all(mean, na.rm = TRUE)
  tpm_primary = as.data.frame(tpm_primary)
  rownames(tpm_primary) = as.character(tpm_primary$Gene.name)
  tpm_primary = tpm_primary[,-1]
  tpm_primary = tpm_primary[which(rowSums(tpm_primary) > 0),] #filtered out genes
  save(tpm_primary, file = paste0("~/",type,"_tpm_primary_genenames.RData"))}

# 2. Proteomic data

# Proteomics data can be downloaded from CPTAC (https://cptac-data-portal.georgetown.edu/)

# Preparation of proteomic data after downloading from CPTAC

# COREAD
require(openxlsx)
COREAD_proteome = read.delim("~/TCGA_Colon_VU_Proteome_CDAP_Protein_Report.r2/TCGA_Colon_VU_Proteome_CDAP.r2.summary.tsv", sep = "\t", header = TRUE)
COREAD_pcounts = COREAD_proteome[, grepl("Spectral.Counts",names(COREAD_proteome))]
COREAD_pcounts = COREAD_pcounts[,-c(96)]

#Converting CPTAC barcodes to TCGA sample names
ConvertBarcodes <- function(x){
  name = gsub("\\.","-",x)
  name = strsplit(name, "-")
  final_name = paste("TCGA",name[[1]][1],name[[1]][2],name[[1]][3], sep = "-")
  return(final_name)}

cnames = c()
for (i in colnames(COREAD_pcounts)) {
  new_name = ConvertBarcodes(i)
  cnames = c(cnames,new_name)}
cnames = sapply(cnames, function(x) substr(x,1,15))
colnames(COREAD_pcounts) = cnames
rownames(COREAD_pcounts) = c(as.character(COREAD_proteome$Gene))

#Some column names are same! For those, the mean was taken
nms = unique(names(COREAD_pcounts))
COREAD_pcounts = sapply(nms, function(x) rowMeans(COREAD_pcounts[names(COREAD_pcounts) %in% x]))
COREAD_pcounts = as.data.frame(COREAD_pcounts) #5561 genes x 90 samples
save(COREAD_pcounts, file = "~/COREAD_pcounts.RData")

# OV
OV_proteome_JHU = read.delim("~/TCGA_Ovarian_JHU_Proteome_CDAP_Protein_Report.r3/TCGA_Ovarian_JHU_Proteome.itraq.tsv", sep = "\t", header = TRUE)
OV_proteome_PNNL = read.delim("~/TCGA_Ovarian_PNNL_Proteome_CDAP_Protein_Report.r3/TCGA_Ovarian_PNNL_Proteome.itraq.tsv", sep = "\t", header = TRUE)
OV_proteome = merge(OV_proteome_JHU[-c(1,2,3),-c(266,267,268,269,270,271)], OV_proteome_PNNL[-c(1,2,3),-c(170,171,172,173,174,175)], by = "Gene")
OV_proteome = OV_proteome[, !grepl("OVARIAN.CONTROL", names(OV_proteome))]
OV_pcounts = OV_proteome[, !grepl("Unshared.Log.Ratio",names(OV_proteome))] #Take the "log ratio" samples -relative abundances with respect to the pooled reference
rownames(OV_pcounts) = as.character(OV_pcounts$Gene)
OV_pcounts = OV_pcounts[,-c(1)]

name.converter <- function(x){
  name = gsub("X", "", x)
  name = gsub("\\.","-",name)
  name = substr(name,1,10)
  name = paste("TCGA",name,sep = "-")
  return(name)}

new_colnames = as.vector(as.character(sapply(colnames(OV_pcounts), name.converter)))
colnames(OV_pcounts) = new_colnames
#Some column names are same! For those, the mean was taken.
nms = unique(names(OV_pcounts))
OV_pcounts = sapply(nms, function(x) rowMeans(OV_pcounts[names(OV_pcounts) %in% x]))
OV_pcounts = as.data.frame(OV_pcounts) #7169 genes x 174 samples
save(OV_pcounts, file = "~/OV_pcounts.RData")

# BRCA
BRCA_proteome = read.delim("~/TCGA_Breast_BI_Proteome_CDAP_Protein_Report.r3/TCGA_Breast_BI_Proteome.itraq.tsv", sep = "\t", header = TRUE)
BRCA_proteome = BRCA_proteome[-c(1,2,3),-c(224,225,226,227,228,229)]
BRCA_pcounts = BRCA_proteome[, !grepl("Unshared.Log.Ratio",names(BRCA_proteome))]
rownames(BRCA_pcounts) = as.character(BRCA_pcounts$Gene)
BRCA_pcounts = BRCA_pcounts[,-c(1,2,3,4)]

name.converter <- function(x){
  name = gsub("\\.","-",x)
  name = substr(name,1,10)
  name = paste("TCGA",name,sep = "-")
  return(name)}

new_colnames = as.vector(as.character(sapply(colnames(BRCA_pcounts), name.converter)))
colnames(BRCA_pcounts) = new_colnames
#Some column names are same! For those, the mean was taken. 
nms = unique(names(BRCA_pcounts))
BRCA_pcounts = sapply(nms, function(x) rowMeans(BRCA_pcounts[names(BRCA_pcounts) %in% x]))
BRCA_pcounts = as.data.frame(BRCA_pcounts) #10625 genes x 105 samples
save(BRCA_pcounts, file = "~/BRCA_pcounts.RData")

# Take the overlap samples between CPTAC and TCGA*

# Take the overlap genes between CPTAC and TCGA*

# Normalization - quantile normalization followed by log2 transformation (for COREAD)

c.types = c("OV","BRCA","COREAD")
for (type in c.types){
  assign("tpm_primary",get(load(paste0("~/",type,"_tpm_primary_genenames.RData")))) # Output files produced by transcriptome data processing code
  tpm_primary = as.data.frame(tpm_primary)
  assign("pcounts_primaryI",get(load(paste0("~/",type,"_pcounts.RData")))) # Folder storing *_pcounts.RData files for COREAD, OV, and BRCA cancer types
  pcounts_primaryI = as.data.frame(pcounts_primaryI)
  overlap_samples = intersect(colnames(tpm_primary),colnames(pcounts_primaryI))
  overlap_genes = intersect(rownames(tpm_primary),rownames(pcounts_primaryI))
  pcounts_primaryI = pcounts_primaryI[rownames(pcounts_primaryI) %in% overlap_genes, colnames(pcounts_primaryI) %in% overlap_samples]
  if(type == "COREAD"){
    #Normalization for COREAD data
    pcounts_primary.m = as.matrix(pcounts_primaryI)
    pcounts_primary.m = normalize.quantiles(pcounts_primary.m, copy = TRUE)
    pcounts_primary.m = log2(pcounts_primary.m + 1)
    colnames(pcounts_primary.m) = colnames(pcounts_primaryI)
    rownames(pcounts_primary.m) = rownames(pcounts_primaryI)
    assign("pcounts_primary",pcounts_primary.m)}
  else{
    assign("pcounts_primary",pcounts_primaryI)}
  save(pcounts_primary, file = paste0("~/",type,"_pcounts.RData"))}