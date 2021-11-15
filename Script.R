
# Set a working directory

setwd(dir = "./Github") # Path to where you download the files from the repository

# Required libraries

library(openxlsx)
library(dplyr)
library(Hmisc)
library(gdata)
library(tidyr)
library(survival)
library(survminer)

## CALCULATING WHOLE-CHROMOSOME LEVEL ANEUPLOIDY SCORES

url = "https://ars.els-cdn.com/content/image/1-s2.0-S1535610818301119-mmc2.xlsx" # Arm-level aneuploidy scores from Taylor et al., 2018
download.file(url, "./Aneuploidyscores.xlsx") # It can be found in the "Data" folder

supp_mmc2 = read.xlsx("./Aneuploidyscores.xlsx", 
                      rowNames = FALSE, colNames = TRUE, startRow = 2)
Aneuploidy_table = data.frame("Sample" = supp_mmc2$Sample,
                              "Type" = supp_mmc2$Type)
Aneuploidy_table[,1:2] = sapply(Aneuploidy_table[,1:2], as.character)

chromosomes = seq(1,22)
arms = c("p","q")
j = 3 #Aneuploidy scores start from column 3

for (i in chromosomes) {
  if(i %in% c(13, 14, 15, 21, 22)){
    columns = paste0(i,"q")
    chr.df = as.data.frame(supp_mmc2[,columns])
    colnames(chr.df) = "q"
    Aneuploidy_table = cbind(Aneuploidy_table,chr.df$q)
    colnames(Aneuploidy_table)[j] = paste0("chr",i)
    j = j + 1}
  else{
    columns = paste0(i,arms)
    chr.df = supp_mmc2[,columns]
    colnames(chr.df) = arms
    chr.df$score = ifelse(chr.df$p == 0 & chr.df$q == 0,0,
                          ifelse(chr.df$p == 1 & chr.df$q == 1,1,
                                 ifelse(chr.df$p == -1 & chr.df$q == -1,-1,NA)))
    Aneuploidy_table = cbind(Aneuploidy_table,chr.df$score)
    colnames(Aneuploidy_table)[j] = paste0("chr",i)
    j = j + 1}}

Aneuploidy_table$Type = ifelse(Aneuploidy_table$Type %in% c("COAD","READ"),"COREAD",Aneuploidy_table$Type) #Because we considered colon and rectal tumors together as COREAD
# Aneuploidy_table can be found in the "Data" folder as "AneuploidyTable.RData"

## DETECTING CANCER TYPE-SPECIFIC WHOLE CHROMOSOME-LEVEL ANEUPLOIDIES

load("./Data/AneuploidyTable.RData")
chromosomes = paste0("chr",as.character(seq(1,22)))
c.types = unique(Aneuploidy_table$Type)
Summary.df = c()
for (t in c.types) {
  mm = Aneuploidy_table[Aneuploidy_table$Type == t,]
  mm.Amp.sum = sum(mm == 1, na.rm = TRUE)
  mm.Del.sum = sum(mm == -1, na.rm = TRUE)
  mm.Diploid.sum = sum(mm == 0, na.rm = TRUE)
  for (c in chromosomes) {
    mm.chrI.Amp = sum(mm[,c] == 1, na.rm = TRUE)
    mm.chrI.Del = sum(mm[,c] == -1, na.rm = TRUE)
    mm.chrI.Diploid = sum(mm[,c] == 0, na.rm = TRUE)
    # Chi-Square for Chromosome of interest (ChrI). ChrO stands for other chromosomes
    df.amp = matrix(c(mm.chrI.Amp,(mm.Amp.sum - mm.chrI.Amp),(mm.chrI.Del+mm.chrI.Diploid),(mm.Del.sum - mm.chrI.Del) + (mm.Diploid.sum - mm.chrI.Diploid)), ncol = 2)
    colnames(df.amp) = c("Amp","Del.Diploid")
    rownames(df.amp) = c("ChrI","ChrO")
    df.del = matrix(c(mm.chrI.Del,(mm.Del.sum - mm.chrI.Del),(mm.chrI.Amp+mm.chrI.Diploid),(mm.Amp.sum - mm.chrI.Amp) + (mm.Diploid.sum - mm.chrI.Diploid)), ncol = 2)
    colnames(df.del) = c("Del","Amp.Diploid")
    rownames(df.del) = c("ChrI","ChrO")
    test.res.amp = chisq.test(df.amp)
    OR.amp = (mm.chrI.Amp * ((mm.Del.sum - mm.chrI.Del) + (mm.Diploid.sum - mm.chrI.Diploid))) /
      ((mm.chrI.Del+mm.chrI.Diploid) * (mm.Amp.sum - mm.chrI.Amp))
    rw.amp = c(t,c,"Amp", test.res.amp$statistic, test.res.amp$p.value, test.res.amp$stdres[1],OR.amp)
    Summary.df = rbind(Summary.df,rw.amp)
    test.res.del = chisq.test(df.del)
    OR.del = mm.chrI.Del * ((mm.Amp.sum - mm.chrI.Amp) + (mm.Diploid.sum - mm.chrI.Diploid)) /
      ((mm.chrI.Amp+mm.chrI.Diploid) * (mm.Del.sum - mm.chrI.Del))
    rw.del = c(t,c,"Del", test.res.del$statistic, test.res.del$p.value, test.res.del$stdres[1],OR.del)
    Summary.df = rbind(Summary.df,rw.del)}}
colnames(Summary.df) = c("Cancer.Type","Chromosomal.Aneuploidy","TestType","X-squared","p.value","StdResidual","OddsRatio")
Summary.df = as.data.frame(Summary.df)

# Multiple Testing Correction after chi-square test

p = as.numeric(as.character(Summary.df$p.value))
adjusted_pvalues = p.adjust(p, method = "holm", n = length(p))
Summary.df$p.adj = adjusted_pvalues
Summary.df[,1:7] = sapply(Summary.df[,1:7], as.character)
Summary.df[,4:7] = sapply(Summary.df[,4:7], as.numeric)

# Selecting most frequently occured aneuploidies in each cancer type
Summary.df = Summary.df[Summary.df$p.adj <= 0.05 & Summary.df$StdResidual >= 2,] # Summary.df can be found in the "Data" folder as "FrequentAneuploidies.RData"

## DATA PROCESSING

# 1. Transcriptome data

# Transcriptomics data needs to be downloaded from GDC data portal (https://portal.gdc.cancer.gov/)

# Converting FPKM values to TPM

# Convert FPKM to TPM

fpkmToTpm <- function(fpkm){exp(log(fpkm) - log(sum(fpkm)) + log(1e6))}

fpkm_tables = list.files(path = "path/to/where/fpkm tables/are",
                         pattern = "*.txt", recursive = FALSE, full.names = TRUE) # FPKM tables downloaded from GDC
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
  save(tpm_mm, file = paste0("./",c,"_tpm.RData"))}

# TPM table for unique sample names

for (c in c.types) {
  load(paste0("./",c,"_tpm.RData"))
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
  save(tpm_samples, file = paste0("./",c,"_Samples_tpm.RData"))}

# Separate different tissue types within the data

human_genes = read.delim("./Data/human_genes.txt", header = TRUE) 

for (c in c.types) {
  load(paste0("./",c,"_Samples_tpm.RData"))
  mt_genes = c(as.character(human_genes[human_genes$Chromosome.scaffold.name == "MT",]$Gene.stable.ID)) #mitochondrial genes were removed
  tpm_samples = tpm_samples[!(rownames(tpm_samples) %in% mt_genes),]
  tpm_samples = as.data.frame(tpm_samples)
  primary_samples = c()
  for(col in colnames(tpm_samples)){
    if(strsplit(col,"-")[[1]][4] == "01" | strsplit(col,"-")[[1]][4] == "03" |
       strsplit(col,"-")[[1]][4] == "09"){primary_samples = c(primary_samples,col)}}
  tpm_samples = tpm_samples[,colnames(tpm_samples) %in% primary_samples]
  save(tpm_samples, file = paste0("./",c,"_tpm_primary.RData"))}

# Convert GeneIDs to Gene names (if there are more than one genes mapped to one ID, take the mean)
# Filter low expressing genes, genes having 0 TPM for all samples were removed

human_genes = read.delim("./Data/human_genes.txt", header = TRUE) 
c.types = c()

# Read the first lines

for (i in 1:length(fpkm_tables)) {
  c = sub("TCGA-","",strsplit(basename(fpkm_tables[i]),"_")[[1]][1])
  c.types = c(c.types,c)}

for (type in c.types) {
  assign("tpm_primary",get(load(paste0("./",type,"_tpm_primary.RData"))))
  tpm_primary = as.data.frame(tpm_primary)
  tpm_primary$Gene.stable.ID = rownames(tpm_primary)
  tpm_primary = merge(tpm_primary, human_genes[,c(1,3)], by = "Gene.stable.ID")
  tpm_primary = tpm_primary[,-1]
  tpm_primary <- tpm_primary %>% group_by(Gene.name) %>% summarise_all(mean, na.rm = TRUE)
  tpm_primary = as.data.frame(tpm_primary)
  rownames(tpm_primary) = as.character(tpm_primary$Gene.name)
  tpm_primary = tpm_primary[,-1]
  tpm_primary = tpm_primary[which(rowSums(tpm_primary) > 0),] #filtered out genes
  save(tpm_primary, file = paste0("./",type,"_tpm_primary_genenames.RData"))} 

# 2. Proteomic data

# Proteomics data needs to be downloaded from CPTAC (https://cptac-data-portal.georgetown.edu/). 

# Preparation of proteomic data after downloading from CPTAC

# COREAD

COREAD_proteome = read.delim("path/to/where/CPTAC data/is/TCGA_Colon_VU_Proteome_CDAP_Protein_Report.r2/TCGA_Colon_VU_Proteome_CDAP.r2.summary.tsv", sep = "\t", header = TRUE)
COREAD_pcounts = COREAD_proteome[, grepl("Spectral.Counts",names(COREAD_proteome))]
COREAD_pcounts = COREAD_pcounts[,-c(96)]

# Converting CPTAC barcodes to TCGA sample names

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

# Some column names are same! For those, the mean was taken

nms = unique(names(COREAD_pcounts))
COREAD_pcounts = sapply(nms, function(x) rowMeans(COREAD_pcounts[names(COREAD_pcounts) %in% x]))
COREAD_pcounts = as.data.frame(COREAD_pcounts) # 5561 genes x 90 samples - It can be found in ./Data/CPTAC_Data

# OV

OV_proteome_JHU = read.delim("path/to/where/CPTAC data/is//TCGA_Ovarian_JHU_Proteome_CDAP_Protein_Report.r3/TCGA_Ovarian_JHU_Proteome.itraq.tsv", sep = "\t", header = TRUE)
OV_proteome_PNNL = read.delim("path/to/where/CPTAC data/is//TCGA_Ovarian_PNNL_Proteome_CDAP_Protein_Report.r3/TCGA_Ovarian_PNNL_Proteome.itraq.tsv", sep = "\t", header = TRUE)
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
OV_pcounts = as.data.frame(OV_pcounts) #7169 genes x 174 samples - It can be found in ./Data/CPTAC_Data

# BRCA

BRCA_proteome = read.delim("path/to/where/CPTAC data/is//TCGA_Breast_BI_Proteome_CDAP_Protein_Report.r3/TCGA_Breast_BI_Proteome.itraq.tsv", sep = "\t", header = TRUE)
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

# Some column names are same! For those, the mean was taken. 

nms = unique(names(BRCA_pcounts))
BRCA_pcounts = sapply(nms, function(x) rowMeans(BRCA_pcounts[names(BRCA_pcounts) %in% x]))
BRCA_pcounts = as.data.frame(BRCA_pcounts) # 10625 genes x 105 samples - It can be found in ./Data/CPTAC_Data

# To run the below code either save the *_tpm_primary_genenames.RData and *_pcounts.RData that can be obtained by running the above code or use the data in the "Data" folder.

# Take the overlap samples/genes between CPTAC and TCGA*

# Normalization - quantile normalization followed by log2 transformation (for COREAD)

c.types = c("OV","BRCA","COREAD")

for (type in c.types){
  assign("tpm_primary",get(load(paste0("./",type,"_tpm_primary_genenames.RData")))) # Output files produced by transcriptome data processing code
  tpm_primary = as.data.frame(tpm_primary)
  assign("pcounts_primaryI",get(load(paste0("./Data/CPTAC_Data/",type,"_pcounts.RData")))) # Output files produced by proteome data processing code
  pcounts_primaryI = as.data.frame(pcounts_primaryI)
  overlap_samples = intersect(colnames(tpm_primary),colnames(pcounts_primaryI))
  overlap_genes = intersect(rownames(tpm_primary),rownames(pcounts_primaryI))
  pcounts_primaryI = pcounts_primaryI[rownames(pcounts_primaryI) %in% overlap_genes, colnames(pcounts_primaryI) %in% overlap_samples]
  if(type == "COREAD"){
    # Normalization for COREAD data
    pcounts_primary.m = as.matrix(pcounts_primaryI)
    pcounts_primary.m = normalize.quantiles(pcounts_primary.m, copy = TRUE)
    pcounts_primary.m = log2(pcounts_primary.m + 1)
    colnames(pcounts_primary.m) = colnames(pcounts_primaryI)
    rownames(pcounts_primary.m) = rownames(pcounts_primaryI)
    assign("pcounts_primary",pcounts_primary.m)}
  else{
    assign("pcounts_primary",pcounts_primaryI)}
  save(pcounts_primary, file = paste0("./",type,"_pcounts.RData"))} # All the files can be found in the "Data/CPTAC_TCGA_Intersection" folder.

## DETECTING TRANSCRIPTOMIC AND PROTEOMIC CHANGES

load("./Data/AneuploidyTable.RData") # This file can be produced by running the AneuploidyScores.R or can be downloaded from the repository
load("./Data/FrequentAneuploidies.RData")

# 1. Whole chromosome-level amplifications

cases_slctd = Summary.df[Summary.df$TestType == "Amp",]

# Detecting transcriptomic changes induced by amplifications 

c.types = unique(cases_slctd$Cancer.Type)
Stat = c()

for (type in c.types) {
  assign("tpm_primary",get(load(paste0("./",type,"_tpm_primary_genenames.RData")))) # Output files produced by transcriptome data processing code
  cases = as.character(cases_slctd[cases_slctd$Cancer.Type == type,"Chromosomal.Aneuploidy"])
  for (c in cases) {
    df = Aneuploidy_table[Aneuploidy_table$Type == type,colnames(Aneuploidy_table) %in% c("Sample","Type",c)]
    df = df[complete.cases(df),]
    colnames(df)[3] = "chr"
    case_samples = df[df$chr == 1,"Sample"] # Samples with amplification of corresponding chromosome
    noch_samples = df[df$chr == 0,"Sample"] # Samples that are diploid for corresponding chromosome
    wilcox.df = c()
    noch.df = tpm_primary[, colnames(tpm_primary) %in% noch_samples]
    colnames(noch.df) = paste0("noch",seq(1,length(colnames(noch.df))))
    case.df = tpm_primary[, colnames(tpm_primary) %in% case_samples]
    colnames(case.df) = paste0("case",seq(1,length(colnames(case.df))))
    wilcox.df = cbind(case.df,noch.df)
    wilcox.df = wilcox.df[which(rowSums(wilcox.df) > 0),] #filtering low expressed genes
    wilcox.df = as.data.frame(wilcox.df)
    rw = c(type, c, dim(case.df)[2], dim(noch.df)[2], dim(wilcox.df)[2])
    Stat = rbind(Stat,rw)
    data_stat = c()
    for (i in 1:length(rownames(wilcox.df))) {
      case = as.numeric(as.vector(wilcox.df[i,grepl("case", names(wilcox.df))]))
      noch = as.numeric(as.vector(wilcox.df[i,grepl("noch", names(wilcox.df))]))
      res = wilcox.test(case,noch)
      roww = c(mean(case, na.rm = TRUE), mean(noch, na.rm = TRUE),median(case), median(noch), res$p.value)
      data_stat = rbind(data_stat, roww)}
    rownames(data_stat) = c(as.character(rownames(wilcox.df)))
    colnames(data_stat) = c("Mean.case","Mean.noch", "Median.case","Median.noch","p.value")
    data_stat = as.data.frame(data_stat)
    p = data_stat$p.value
    adjusted_pvalues = p.adjust(p, method = "BH")
    data_stat$p.adj = adjusted_pvalues
    data_stat = data_stat[complete.cases(data_stat),]
    name = paste("wilcox",type,c,"amp","vs_noch",sep = "_")
    assign(name,data_stat)}}

df.tobesaved = mget(ls(pattern = "^wilcox_")) # The output file can be found in ./Data/DE_Transcriptome.RData

# Detecting proteomic changes induced by amplifications

cases_slctd = cases_slctd[cases_slctd$Cancer.Type %in% c("BRCA","OV","COREAD"),] #13 amplification cases for which proteomic data is available
c.types = c("COREAD", "BRCA","OV") 

sample.info = c()
final.dataset = list()
for (type in c.types) {
  assign("pcounts_primary",get(load(paste0("./Data/CPTAC_TCGA_Intersection/",type,"_pcounts.RData")))) # Output files produced by proteome data processing code
  pcounts_primary = as.data.frame(pcounts_primary)
  print(dim(pcounts_primary))
  cases = as.character(cases_slctd[cases_slctd$Cancer.Type == type,"Chromosomal.Aneuploidy"])
  for (c in cases) {
    df = Aneuploidy_table[Aneuploidy_table$Type == type,colnames(Aneuploidy_table) %in% c("Sample","Type",c)]
    df = df[complete.cases(df),]
    colnames(df)[3] = "chr"
    case_samples = df[df$chr == 1,"Sample"] # Samples with amplification of corresponding chromosome
    noch_samples = df[df$chr == 0,"Sample"] # Samples that are diploid for corresponding chromosome
    wilcox.df = c()
    noch.df = as.data.frame(pcounts_primary[, colnames(pcounts_primary) %in% noch_samples])
    colnames(noch.df) = paste0("noch",seq(1,length(colnames(noch.df))))
    case.df = as.data.frame(pcounts_primary[, colnames(pcounts_primary) %in% case_samples])
    colnames(case.df) = paste0("case",seq(1,length(colnames(case.df))))
    wilcox.df = cbind(case.df,noch.df)
    wilcox.df = as.data.frame(wilcox.df)
    wilcox.df = wilcox.df[which(rowMeans(!is.na(wilcox.df)) > 0.5), ] #Removing genes which are not detected in the half of the samples
    rw = c("CPTAC-TCGA",type, c, dim(case.df)[2], dim(noch.df)[2], dim(wilcox.df)[2])
    sample.info = rbind(sample.info,rw)
    data_stat = c()
    for (i in 1:length(rownames(wilcox.df))) {
      case = as.numeric(as.vector(wilcox.df[i,grepl("case", names(wilcox.df))]))
      noch = as.numeric(as.vector(wilcox.df[i,grepl("noch", names(wilcox.df))]))
      res = wilcox.test(case,noch)
      roww = c(mean(case, na.rm = TRUE), mean(noch, na.rm = TRUE),median(case), median(noch), res$p.value)
      data_stat = rbind(data_stat, roww)}
    rownames(data_stat) = c(as.character(rownames(wilcox.df)))
    colnames(data_stat) = c("Mean.case","Mean.noch", "Median.case","Median.noch","p.value")
    data_stat = as.data.frame(data_stat)
    p = data_stat$p.value
    adjusted_pvalues = p.adjust(p, method = "BH")
    data_stat$p.adj = adjusted_pvalues
    name = paste("wilcox",type,c,"amp","vs_noch",sep = "_")
    final.dataset[[name]] = data_stat}} # The output file can be found in ./Data/DE_Proteome.RData

# 2. Whole chromosome-level deletions

load("./Data/FrequentAneuploidies.RData")
cases_slctd = Summary.df[Summary.df$TestType == "Del",]

# Detecting transcriptomic changes induced by deletions

c.types = unique(cases_slctd$Cancer.Type)
Stat = c()

for (type in c.types) {
  assign("tpm_primary",get(load(paste0("./",type,"_tpm_primary_genenames.RData"))))
  cases = as.character(cases_slctd[cases_slctd$Cancer.Type == type,"Chromosomal.Aneuploidy"])
  for (c in cases) {
    # Select samples which have chr ane (chr1 amp for example) and samples which do not, then prepare wilcox matrix for test
    df = Aneuploidy_table[Aneuploidy_table$Type == type,colnames(Aneuploidy_table) %in% c("Sample","Type",c)]
    df = df[complete.cases(df),]
    colnames(df)[3] = "chr"
    case_samples = df[df$chr == -1,"Sample"]
    noch_samples = df[df$chr == 0,"Sample"]
    wilcox.df = c()
    noch.df = tpm_primary[, colnames(tpm_primary) %in% noch_samples]
    colnames(noch.df) = paste0("noch",seq(1,length(colnames(noch.df))))
    case.df = tpm_primary[, colnames(tpm_primary) %in% case_samples]
    colnames(case.df) = paste0("case",seq(1,length(colnames(case.df))))
    wilcox.df = cbind(case.df,noch.df)
    wilcox.df = wilcox.df[which(rowSums(wilcox.df) > 0),] #filtering low expressed genes
    wilcox.df = as.matrix(wilcox.df)
    rw = c(type, c, dim(case.df)[2], dim(noch.df)[2], dim(wilcox.df)[2])
    Stat = rbind(Stat,rw)
    data_stat = c()
    for (i in 1:length(rownames(wilcox.df))) {
      case = as.numeric(as.vector(wilcox.df[i,grepl("case", colnames(wilcox.df))])) 
      noch = as.numeric(as.vector(wilcox.df[i,grepl("noch", colnames(wilcox.df))]))
      res = wilcox.test(case,noch, exact = FALSE)
      roww = c(mean(case, na.rm = TRUE), mean(noch, na.rm = TRUE),median(case), median(noch), res$p.value)
      data_stat = rbind(data_stat, roww)}
    rownames(data_stat) = c(as.character(rownames(wilcox.df)))
    colnames(data_stat) = c("Mean.case","Mean.noch", "Median.case","Median.noch","p.value")
    data_stat = as.data.frame(data_stat)
    p = data_stat$p.value
    adjusted_pvalues = p.adjust(p, method = "BH")
    data_stat$p.adj = adjusted_pvalues
    data_stat = data_stat[complete.cases(data_stat),]
    name = paste("wilcox",type,c,"del","vs_noch",sep = "_")
    assign(name,data_stat)}}

df.tobesaved = mget(ls(pattern = "^wilcox_"))  # The output file can be found in ./Data/DE_Transcriptome_deletions.RData

# Detecting proteomic changes induced by deletions

cases_slctd = cases_slctd[cases_slctd$Cancer.Type %in% c("BRCA","OV","COREAD"),] #20 deletion cases for which proteomic data is available
c.types = c("COREAD", "BRCA","OV") 
sample.info = c()
final.dataset = list()

for (type in c.types) {
  assign("pcounts_primary",get(load(paste0("./Data/CPTAC_TCGA_Intersection/",type,"_pcounts.RData"))))
  pcounts_primary = as.data.frame(pcounts_primary)
  print(dim(pcounts_primary))
  cases = as.character(cases_slctd[cases_slctd$Cancer.Type == type,"Chromosomal.Aneuploidy"])
  for (c in cases) {
    df = Aneuploidy_table[Aneuploidy_table$Type == type,colnames(Aneuploidy_table) %in% c("Sample","Type",c)]
    df = df[complete.cases(df),]
    colnames(df)[3] = "chr"
    case_samples = df[df$chr == -1,"Sample"]
    noch_samples = df[df$chr == 0,"Sample"]
    wilcox.df = c()
    noch.df = as.data.frame(pcounts_primary[, colnames(pcounts_primary) %in% noch_samples])
    colnames(noch.df) = paste0("noch",seq(1,length(colnames(noch.df))))
    case.df = as.data.frame(pcounts_primary[, colnames(pcounts_primary) %in% case_samples])
    colnames(case.df) = paste0("case",seq(1,length(colnames(case.df))))
    wilcox.df = cbind(case.df,noch.df)
    wilcox.df = as.data.frame(wilcox.df)
    wilcox.df = wilcox.df[which(rowMeans(!is.na(wilcox.df)) > 0.5), ] #It is needed to put this filtering as wilcoxon cannot be performed for some genes having a lot of NAs
    rw = c("CPTAC-TCGA",type, c, dim(case.df)[2], dim(noch.df)[2], dim(wilcox.df)[2])
    sample.info = rbind(sample.info,rw)
    data_stat = c()
    for (i in 1:length(rownames(wilcox.df))) {
      case = as.numeric(as.vector(wilcox.df[i,grepl("case", names(wilcox.df))]))
      noch = as.numeric(as.vector(wilcox.df[i,grepl("noch", names(wilcox.df))]))
      res = wilcox.test(case,noch)
      roww = c(mean(case, na.rm = TRUE), mean(noch, na.rm = TRUE),median(case, na.rm = TRUE), median(noch, na.rm = TRUE), res$p.value)
      data_stat = rbind(data_stat, roww)}
    rownames(data_stat) = c(as.character(rownames(wilcox.df)))
    colnames(data_stat) = c("Mean.case","Mean.noch", "Median.case","Median.noch","p.value")
    data_stat = as.data.frame(data_stat)
    p = data_stat$p.value
    adjusted_pvalues = p.adjust(p, method = "BH")
    data_stat$p.adj = adjusted_pvalues
    name = paste("wilcox",type,c,"del","vs_noch",sep = "_")
    final.dataset[[name]] = data_stat}} # The output file can be found in ./Data/DE_Proteome_deletions.RData

# Detecting significantly changed genes in different cut-offs

human_genes = read.delim("./Data/human_genes.txt", header = TRUE)
human_genes = human_genes[!duplicated(human_genes$Gene.name),]

assign("df.lists", get(load("./Data/DE_Transcriptome.RData"))) # To run the same code for deletions and the proteome data, change the input data!

sig.levels = c(0.1,0.05,0.01)
Stat = c() # To check number of genes after the cut-off

for (i in 1:length(df.lists)) {
  df = df.lists[[i]]
  chr = gsub("chr","",strsplit(names(df.lists)[i],"_")[[1]][3])
  df$Gene.name = rownames(df)
  df = merge(df, human_genes[,c(2,3,4)], by = "Gene.name")
  df$Chr.Info = ifelse(as.character(df$Chromosome.scaffold.name) == chr, "ChrI","ChrO")
  for (s in sig.levels) {
    df$DE.info.pvalue = ifelse(df$p.value < s, "DE","NotDE")
    df$DE.info.padj = ifelse(df$p.adj < s, "DE","NotDE")
    rw = c(names(df.lists)[i], s, dim(df[df$DE.info.pvalue == "DE",])[1],dim(df[df$DE.info.padj == "DE",])[1],dim(df[df$DE.info.pvalue == "DE" & df$Gene.type == "protein_coding",])[1], dim(df[df$DE.info.padj == "DE" & df$Gene.type == "protein_coding",])[1] )
    Stat = rbind(Stat,rw)
    assign(paste(names(df.lists)[i],s,sep = "_"),df)}}

colnames(Stat) = c("Case","Significance","Number of genes.pvalue","Number of genes.padj","Number of pcgenes.pvalue","Number of pcgenes.padj")
df.lists.sig = mget(ls(pattern = "^wilcox_")) # Outputs can be found in the "Data" folder: DE_Transcriptome_Sig.RData, DE_Proteome_Sig.RData, DE_Transcriptome_Sig_deletions.RData, DE_Proteome_Sig_deletions.RData

## STATISTICAL ANALYSES

# 1. Chi-Square tests

load("./Data/HumanPCgenes_Subunits_allchr.RData")

# 1.1 Transcriptome Level Analyses

# For amplifications

load("./Data/DE_Transcriptome_Sig.RData")
Info.table = read.xlsx("./Data/Table_cutoff.xlsx",colNames = TRUE,rowNames = FALSE, sheet = 1) # File for which cut-off were used to detect DEGs in each aneuploidy case

# For deletions

#load("./Data/DE_Transcriptome_Sig_deletions.RData")
#Info.table = read.xlsx("./Data/Table_cutoff.xlsx",colNames = TRUE,rowNames = FALSE, sheet = 2) 

padj_0.1 = Info.table[Info.table$Sig == "padj<0.1",] 
pvalue_0.05 = Info.table[Info.table$Sig == "pvalue<0.05",] 

cases = Info.table$Case
Summary = c()

# 1.1.1 Overlap with complex subunits

for (case in cases) {
  if(case %in% padj_0.1$Case){
    sig = "0.1"
    val = "p.adj"}
  else{
    sig = "0.05"
    val = "p.value"}
  name = paste("wilcox",case,"amp_vs_noch",sig,sep = "_")
  #name = paste("wilcox",case,"del_vs_noch",sig,sep = "_") # For deletions
  df = df.lists.sig[[which(names(df.lists.sig)==name)]]
  df = df[complete.cases(df),]
  df = df[df$Gene.type == "protein_coding",]
  df$LogFC = log2(df$Median.case / df$Median.noch)
  df = df[complete.cases(df),]
  df = df[df$Chr.Info == "ChrO",]
  if(val == "p.adj"){
    ChrO_DE = df[df$Chr.Info == "ChrO" & df$DE.info.padj == "DE" & abs(df$LogFC) >= 0,]
    if(dim(ChrO_DE)[1] == 0)next
    ChrO_DE_subunit = ChrO_DE[as.character(ChrO_DE$Gene.name) %in% as.character(Subunits_Info_pc_genes$Gene.name),]
    ChrO_NotDE = df[!(df$Gene.name %in% ChrO_DE$Gene.name),]
    ChrO_NotDE_subunit = ChrO_NotDE[as.character(ChrO_NotDE$Gene.name) %in% as.character(Subunits_Info_pc_genes$Gene.name),]}
  else{
    ChrO_DE = df[df$Chr.Info == "ChrO" & df$DE.info.pvalue == "DE" & abs(df$LogFC) >= 0,]
    if(dim(ChrO_DE)[1] == 0)next
    ChrO_DE_subunit = ChrO_DE[as.character(ChrO_DE$Gene.name) %in% as.character(Subunits_Info_pc_genes$Gene.name),]
    ChrO_NotDE = df[!(df$Gene.name %in% ChrO_DE$Gene.name),]
    ChrO_NotDE_subunit = ChrO_NotDE[as.character(ChrO_NotDE$Gene.name) %in% as.character(Subunits_Info_pc_genes$Gene.name),]}
  mm_ChrO = matrix(c(dim(ChrO_DE_subunit)[1],dim(ChrO_NotDE_subunit)[1],
                     (dim(ChrO_DE)[1] - dim(ChrO_DE_subunit)[1]),(dim(ChrO_NotDE)[1]) - dim(ChrO_NotDE_subunit)[1]), ncol = 2)
  colnames(mm_ChrO) = c("Complex","NotComplex")
  rownames(mm_ChrO) = c("DE","NotDE")
  mm = mm_ChrO
  res = chisq.test(mm)
  mm.res = cbind("Status" = c("DE.sub","NotDE.sub"),"Res/StdRes" = c(res$stdres[1],res$stdres[2]))
  final.mm = cbind(mm, mm.res,
                   data.frame("p.value" = c(res$p.value,""),
                              "Per.DE.sub" = c((dim(ChrO_DE_subunit)[1]/dim(ChrO_DE)[1])*100,""),
                              "Per.NotDE.sub" = c((dim(ChrO_NotDE_subunit)[1]/dim(ChrO_NotDE)[1])*100,""),
                              "Case" = c(case,case),
                              "Cutoff" = c(paste(val,sig,sep = "_"),paste(val,sig,sep = "_"))))
  Summary = rbind(Summary,final.mm)}

# Multiple testing correction
Summary = Summary[grep(pattern = "^DE", rownames(Summary)),]
p = as.numeric(as.character(Summary$p.value))
adjusted_pvalues = p.adjust(p, method = "holm", n = length(p))
Summary$p.adj = adjusted_pvalues

## 1.1.2 Overlap with partners

for (case in cases) {
  if(case %in% padj_0.1$Case){
    sig = "0.1"
    val = "p.adj"}
  else{
    sig = "0.05"
    val = "p.value"}
  name = paste("wilcox",case,"amp_vs_noch",sig,sep = "_")
  #name = paste("wilcox",case,"del_vs_noch",sig,sep = "_") # For deletions
  df = df.lists.sig[[which(names(df.lists.sig)==name)]]
  df = df[complete.cases(df),]
  df = df[df$Gene.type == "protein_coding",]
  df$LogFC = log2(df$Median.case / df$Median.noch)
  df = df[complete.cases(df),]
  if(val == "p.adj"){df_ChrI_DE = df[as.character(df$Chr.Info) == "ChrI" & as.character(df$DE.info.padj) == "DE" & abs(df$LogFC) >= 0,]} #DEGs on aneuploid chromosomes
  else{df_ChrI_DE = df[as.character(df$Chr.Info) == "ChrI" & as.character(df$DE.info.pvalue) == "DE" & abs(df$LogFC) >= 0,]}
  #Check point
  if(dim(df_ChrI_DE)[1] == 0)next
  all.partners = Subunits.Partner.Info[Subunits.Partner.Info$Subunit %in% df_ChrI_DE$Gene.name,]
  all.partners = unlist(strsplit(all.partners$Partners, split=";"))
  #Check point
  if(length(all.partners) == 0)next
  #Overlap between Subunits and proteins encoded on other chromosomes
  data = df[as.character(df$Chr.Info) == "ChrO",]
  data = data[as.character(data$Gene.name) %in% as.character(Subunits_Info_pc_genes$Gene.name),]
  #check point
  if(dim(data)[1] < 10)next
  data$Partner = ifelse(as.character(data$Gene.name) %in% all.partners, "Partner","Not.Partner")
  partner.df = data[as.character(data$Partner) == "Partner",]
  notpartner.df = data[as.character(data$Partner) == "Not.Partner",]
  if(val == "p.adj"){
    #Check point
    if(dim(data[as.character(data$DE.info.padj) == "DE" & abs(data$LogFC) >= 0,])[1] < 10)next
    mm = matrix(c(dim(partner.df[as.character(partner.df$DE.info.padj) == "DE" & abs(partner.df$LogFC) >= 0,])[1],
                  dim(partner.df[!(abs(partner.df$LogFC) >= 0) | as.character(partner.df$DE.info.padj) == "NotDE",])[1],
                  dim(notpartner.df[as.character(notpartner.df$DE.info.padj) == "DE" & abs(notpartner.df$LogFC) >= 0,])[1],
                  dim(notpartner.df[!(abs(notpartner.df$LogFC) >= 0) | as.character(notpartner.df$DE.info.padj) == "NotDE",])[1]), ncol = 2)}
  else{
    #Check point
    if(dim(data[as.character(data$DE.info.pvalue) == "DE" & abs(data$LogFC) >= 0,])[1] < 10)next
    mm = matrix(c(dim(partner.df[as.character(partner.df$DE.info.pvalue) == "DE" & abs(partner.df$LogFC) >= 0,])[1],
                  dim(partner.df[!(abs(partner.df$LogFC) >= 0) | as.character(partner.df$DE.info.pvalue) == "NotDE",])[1],
                  dim(notpartner.df[as.character(notpartner.df$DE.info.pvalue) == "DE" & abs(notpartner.df$LogFC) >= 0,])[1],
                  dim(notpartner.df[!(abs(notpartner.df$LogFC) >= 0) | as.character(notpartner.df$DE.info.pvalue) == "NotDE",])[1]), ncol = 2)}
  colnames(mm) = c("Partner","NotPartner")
  rownames(mm) = c("DE","NotDE")
  res = chisq.test(mm)
  mm.res = cbind("Status" = c("Partner.DE","Partner.NotDE"),"Res/StdRes" = c(res$stdres[1],res$stdres[2]))
  final.mm = cbind(mm, mm.res,
                   data.frame("p.value" = c(res$p.value,""),
                              "Per.DE.Partner" = c((mm[1,1]/(mm[1,1]+mm[1,2]))*100,""),
                              "Per.NotDE.Partner" = c((mm[2,1]/(mm[2,1]+mm[2,2]))*100,""),
                              "Case" = c(case,case),
                              "Cutoff" = c(paste(val,sig,sep = "_"),paste(val,sig,sep = "_"))))
  Summary = rbind(Summary,final.mm)}

# Multiple testing correction
Summary = Summary[grep(pattern = "^DE", rownames(Summary)),]
p = as.numeric(as.character(Summary$p.value))
adjusted_pvalues = p.adjust(p, method = "holm", n = length(p))
Summary$p.adj = adjusted_pvalues

# 1.2 Proteome Level Analyses

# For amplifications

load("./Data/DE_Proteome_Sig.RData")
Info.table = read.xlsx("./Data/Table_cutoff.xlsx",colNames = TRUE,rowNames = FALSE, sheet = 1)

# For deletions

#load("./Data/DE_Proteome_Sig_deletions.RData")
#Info.table = read.xlsx("./Data/Table_cutoff.xlsx",colNames = TRUE,rowNames = FALSE, sheet = 2)

Info.table$Type = unlist(lapply(Info.table$Case, function(x) strsplit(x, "_")[[1]][1]))
Info.table = Info.table[Info.table$Type %in% c("BRCA","OV","COREAD"),]

cases = Info.table$Case
Summary = c()

# 1.2.1 Overlap with complex subunits

for (case in cases) {
  sig = "0.1"
  val = "p.value"
  name = paste("wilcox",case,"amp_vs_noch",sig,sep = "_")
  #name = paste("wilcox",case,"del_vs_noch",sig,sep = "_")
  df = df.lists.sig[[which(names(df.lists.sig)==name)]]
  df = df[complete.cases(df),]
  df = df[df$Gene.type == "protein_coding",]
  df = df[df$Chr.Info == "ChrO",]
  ChrO_DE = df[df$Chr.Info == "ChrO" & df$DE.info.pvalue == "DE",]
  if(dim(ChrO_DE)[1] == 0)next
  ChrO_DE_subunit = ChrO_DE[as.character(ChrO_DE$Gene.name) %in% as.character(Subunits_Info_pc_genes$Gene.name),]
  ChrO_NotDE = df[!(df$Gene.name %in% ChrO_DE$Gene.name),]
  ChrO_NotDE_subunit = ChrO_NotDE[as.character(ChrO_NotDE$Gene.name) %in% as.character(Subunits_Info_pc_genes$Gene.name),]
  mm_ChrO = matrix(c(dim(ChrO_DE_subunit)[1],dim(ChrO_NotDE_subunit)[1],
                     (dim(ChrO_DE)[1] - dim(ChrO_DE_subunit)[1]),(dim(ChrO_NotDE)[1]) - dim(ChrO_NotDE_subunit)[1]), ncol = 2)
  colnames(mm_ChrO) = c("Complex","NotComplex")
  rownames(mm_ChrO) = c("DE","NotDE")
  mm = mm_ChrO
  res = chisq.test(mm)
  mm.res = cbind("Status" = c("DE.sub","NotDE.sub"),"Res/StdRes" = c(res$stdres[1],res$stdres[2]))
  final.mm = cbind(mm, mm.res,
                   data.frame("p.value" = c(res$p.value,""),
                              "Per.DE.sub" = c((dim(ChrO_DE_subunit)[1]/dim(ChrO_DE)[1])*100,""),
                              "Per.NotDE.sub" = c((dim(ChrO_NotDE_subunit)[1]/dim(ChrO_NotDE)[1])*100,""),
                              "Case" = c(case,case),
                              "Cutoff" = c(paste(val,sig,sep = "_"),paste(val,sig,sep = "_"))))
  Summary = rbind(Summary,final.mm)}

## 1.2.2 Overlap with partners

for (case in cases) {
  sig = "0.1"
  val = "p.value"
  name = paste("wilcox",case,"amp_vs_noch",sig,sep = "_")
  #name = paste("wilcox",case,"del_vs_noch",sig,sep = "_")
  df = df.lists.sig[[which(names(df.lists.sig)==name)]]
  df = df[complete.cases(df),]
  df = df[df$Gene.type == "protein_coding",]
  df_ChrI_DE = df[as.character(df$Chr.Info) == "ChrI" & as.character(df$DE.info.pvalue) == "DE",]
  #Check point
  if(dim(df_ChrI_DE)[1] == 0)next
  all.partners = Subunits.Partner.Info[Subunits.Partner.Info$Subunit %in% df_ChrI_DE$Gene.name,]
  all.partners = unlist(strsplit(all.partners$Partners, split=";"))
  #Check point
  if(length(all.partners) == 0)next
  #Overlap between Subunits and proteins encoded on other chromosomes
  data = df[as.character(df$Chr.Info) == "ChrO",]
  data = data[as.character(data$Gene.name) %in% as.character(Subunits_Info_pc_genes$Gene.name),] #intersect between corum subunits
  #check point
  if(dim(data)[1] < 10)next
  data$Partner = ifelse(as.character(data$Gene.name) %in% all.partners, "Partner","Not.Partner")
  partner.df = data[as.character(data$Partner) == "Partner",]
  notpartner.df = data[as.character(data$Partner) == "Not.Partner",]
  #Check point
  if(dim(data[as.character(data$DE.info.pvalue) == "DE",])[1] < 10)next
  mm = matrix(c(dim(partner.df[as.character(partner.df$DE.info.pvalue) == "DE",])[1],
                dim(partner.df[as.character(partner.df$DE.info.pvalue) == "NotDE",])[1],
                dim(notpartner.df[as.character(notpartner.df$DE.info.pvalue) == "DE",])[1],
                dim(notpartner.df[as.character(notpartner.df$DE.info.pvalue) == "NotDE",])[1]), ncol = 2)
  colnames(mm) = c("Partner","NotPartner")
  rownames(mm) = c("DE","NotDE")
  res = chisq.test(mm)
  mm.res = cbind("Status" = c("Partner.DE","Partner.NotDE"),"Res/StdRes" = c(res$stdres[1],res$stdres[2]))
  final.mm = cbind(mm, mm.res,
                   data.frame("p.value" = c(res$p.value,""),
                              "Per.DE.Partner" = c((mm[1,1]/(mm[1,1]+mm[1,2]))*100,""),
                              "Per.NotDE.Partner" = c((mm[2,1]/(mm[2,1]+mm[2,2]))*100,""),
                              "Case" = c(case,case),
                              "Cutoff" = c(paste(val,sig,sep = "_"),paste(val,sig,sep = "_"))))
  Summary = rbind(Summary,final.mm)}


# 2. Correlations between protein pairs

# For COREAD, OV, and BRCA primary tumour samples and genes (covered by both transcriptomic and proteomic data), all possible pairwise correlations at both mRNA 
# and proteomic levels were calculated. Then we focused on the correlations between proteins that encoded by genes found on amplified chromosome and 
# show significant changes in abundances and their partners encoded by genes found on other chromosomes.

# Function to make df for correlation results

flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut])}

c.types = c("COREAD","OV","BRCA")

for (type in c.types) {
  # Transcriptomic measurements
  assign("tpm_primary",get(load(paste0("./Data/TCGA_Data/",type,"_tpm_primary_genenames.RData")))) # Output files can be produced by transcriptome data processing code
  # Proteomic measurements
  assign("pcounts_primary",get(load(paste0("./Data/CPTAC_TCGA_Intersection/",type,"_pcounts.RData")))) # Output files can be produced by proteome data processing code
  pcounts_primary = as.data.frame(pcounts_primary)
  # Make sure transcriptome data cover the same samples and genes
  tpm_primary = tpm_primary[rownames(pcounts_primary),colnames(pcounts_primary)]
  # Calculating correlations
  # mRNA
  mm.t = t(tpm_primary)
  res.t <- rcorr(as.matrix(mm.t), type = "spearman")
  Corr.mat.t = flattenCorrMatrix(res.t$r, res.t$P)
  # Protein
  mm.p = t(pcounts_primary)
  res.p <- rcorr(as.matrix(mm.p), type = "spearman")
  Corr.mat.p = flattenCorrMatrix(res.p$r, res.p$P)
  print(type)
  save(Corr.mat.p,Corr.mat.t, 
       file = paste0("./SpearmanCorr_allpairs_mrna&protein_",type,".RData"))} # The outputs can be found in the Data folder

# 3. Making unique correlation tables for amplification and deletion cases (Unique set of aneuploid proteins and their partners)

# To obtain differentially abundant aneuploid proteins, load proteomic changes data

# For amplifications
load("./Data/DE_Proteome_Sig.RData")

# For deletions
#load("./Data/DE_Proteome_Sig_deletions.RData")

cases = unlist(lapply(names(df.lists.sig), function(x) paste(strsplit(x,"_")[[1]][2],strsplit(x,"_")[[1]][3],sep = "_")))
cases = unique(cases)

Cor_DAP_Partners = c()
types = c("COREAD","OV","BRCA")

for(type in types){
  load(paste0("./Data/SpearmanCorr_allpairs_mrna&protein_",type,".RData"))
  type.cases = cases[grep(pattern = type, cases)]
  for(case in type.cases){
    # Find differentiall up(down)-regulated proteins
    name = paste("wilcox",case,"amp_vs_noch","0.1",sep = "_")
    #name = paste("wilcox",case,"del_vs_noch","0.1",sep = "_") # If the code is running for deletion cases
    df = df.lists.sig[[which(names(df.lists.sig)==name)]]
    df = df[complete.cases(df),]
    df = df[df$Gene.type == "protein_coding",]
    df = df[df$Chr.Info == "ChrI" & df$DE.info.pvalue == "DE",] #Differentially abundant proteins on the aneuploid chromosome
    df$Diff = df$Median.case - df$Median.noch
    df = df[df$Diff > 0,]
    #df = df[df$Diff < 0,] # If the code is running for deletion cases
    # Take correlations for those proteins and their pairs
    # protein correlations
    protein.corr.row = Corr.mat.p[Corr.mat.p$row %in% df$Gene.name,]
    colnames(protein.corr.row) = c("DAP","Partner","p.cor","p.p")
    protein.corr.col = Corr.mat.p[Corr.mat.p$column %in% df$Gene.name,]
    colnames(protein.corr.col) = c("Partner","DAP","p.cor","p.p")
    protein.corr = rbind(protein.corr.row,protein.corr.col)
    protein.corr$Pair = paste(protein.corr$DAP,protein.corr$Partner,sep = "_")
    # mrna correlations
    mrna.corr.row = Corr.mat.t[Corr.mat.t$row %in% df$Gene.name,]
    colnames(mrna.corr.row) = c("DAP","Partner","t.cor","t.p")
    mrna.corr.col = Corr.mat.t[Corr.mat.t$column %in% df$Gene.name,]
    colnames(mrna.corr.col) = c("Partner","DAP","t.cor","t.p")
    mrna.corr = rbind(mrna.corr.row,mrna.corr.col)
    mrna.corr$Pair = paste(mrna.corr$DAP,mrna.corr$Partner,sep = "_")
    # Merging data from each case
    corr.df = merge(mrna.corr[,-c(1,2)],protein.corr[,-c(1,2)], by = "Pair")
    corr.df$Case = case
    # Final data
    Cor_DAP_Partners = rbind(Cor_DAP_Partners,corr.df)
    print(case)}
  keep(cases,Cor_DAP_Partners,types,df.lists.sig, sure = TRUE)}

# Correlations between DAP and their partners of other chromosomes

load("./Data/HumanPCgenes_Subunits_allchr.RData")

Cor_DAP_Partners_new = c() #the correlated proteins that are covered by the CORUM database

for(case in cases){
  case.cor = Cor_DAP_Partners[Cor_DAP_Partners$Case == case,]
  case.cor = case.cor %>% separate(Pair, c("DAP", "Partner"), sep = "_", remove = FALSE)
  #For further analyses, continue with the proteins covered by the CORUM database
  case.cor = case.cor[case.cor$DAP %in% Subunits_Info_pc_genes$Gene.name &
                        case.cor$Partner %in% Subunits_Info_pc_genes$Gene.name,]
  case.cor = merge(case.cor, 
                   Subunits_Info_pc_genes[Subunits_Info_pc_genes$Chromosome.scaffold.name %in% c(as.character(seq(1,22)),"X","Y"),c(2,3)], 
                   by.x = "Partner", by.y = "Gene.name")
  case.cor$Chr.info = ifelse(case.cor$Chromosome.scaffold.name != gsub("chr","",strsplit(case,"_")[[1]][2]),"Other","Aneuploid")
  case.cor = case.cor[!(duplicated(case.cor$Pair)),]
  DAPs = unique(case.cor$DAP)
  case.cor.new = c()
  for(p in DAPs){
    all.partners = unlist(strsplit(Subunits.Partner.Info[Subunits.Partner.Info$Subunit == p,"Partners"], split = ";"))
    protein.cor = case.cor[case.cor$DAP == p,]
    protein.cor$Cocomplex.member = ifelse(protein.cor$Partner %in% all.partners, "YES","NO")
    case.cor.new = rbind(case.cor.new,protein.cor)}
  Cor_DAP_Partners_new = rbind(Cor_DAP_Partners_new,case.cor.new)}

# Correlation between co-complex members only

Final.df = Cor_DAP_Partners_new[Cor_DAP_Partners_new$Cocomplex.member == "YES" &
                                  Cor_DAP_Partners_new$Chr.info == "Other",]

# Remove duplicated pairs by condering the one with the maximum correlation in absolute

Final.df = Final.df %>% group_by(Pair) %>% mutate(duplicate.flag = n() > 1)
repeated.pairs = Final.df %>% filter(duplicate.flag)
unique.pairs = Final.df %>% filter(!duplicate.flag)

repeated.pairs = repeated.pairs[order(abs(repeated.pairs$p.cor), decreasing = TRUE),]
repeated.pairs = repeated.pairs[!duplicated(repeated.pairs$Pair),]
Final.df = rbind(unique.pairs,repeated.pairs)

# Check for interacting pairs like A-B and B-A, and keep only one

Final.df.intbased = c()

for(i in 1:length(rownames(Final.df))){
  pair = Final.df$Pair[i]
  pair.r = paste(Final.df$Partner[i],Final.df$DAP[i],sep = "_")
  df = Final.df[Final.df$Pair %in% c(pair,pair.r),]
  print(paste(i, dim(df)[1],sep = " "))
  df = df[order(df$p.cor, decreasing = TRUE),]
  Final.df.intbased = rbind(Final.df.intbased,df[1,])}

Final.df.intbased = Final.df.intbased[!duplicated(Final.df.intbased$Pair),] # Correlations between differentially abundant proteins of aneuploid chromosome and their co-complex members of other chromosomes: A unique set

# The output files can be found in ./Data folder: Correlations_cocomplexmembers_amplifications.RData and Correlations_cocomplexmembers_deletions.RData

## GROUPING PROTEINS AND PROTEIN PAIRS

load("./Data/HumanPCgenes_Subunits_allchr.RData")

# Correlations between aneuploid proteins and their partners 

load("./Data/Correlations_cocomplexmembers_amplifications.RData") # For amplification cases
#load("./Data/Correlations_cocomplexmembers_deletions.RData") # For deletion cases

Correlations = Final.df.intbased

# The number of complexes a protein involves and number of partners

for(i in 1:length(rownames(Correlations))){
  p = as.character(unlist(Correlations$DAP))[i]
  p.comp = human_proteincomplex[grepl(paste0('\\b',p,'\\b'),human_proteincomplex$subunits.Gene.name.),]
  p.partners = strsplit(Subunits.Partner.Info[Subunits.Partner.Info$Subunit == p,"Partners"],";")[[1]]
  Correlations$Number.complex[i] = length(rownames(p.comp))
  Correlations$Number.partner[i] = length(p.partners)}
Correlations$Promiscuity = ifelse(Correlations$Number.complex > 5,"Promiscuous","Non-promiscuous")

# Aggregation propensity

aggregators = read.delim("./Data/aggregators.csv")
Correlations$Aggregation = ifelse(Correlations$DAP %in% aggregators$aggregators, "Aggregation prone","Non-aggregation prone")

# Co-occurrence frequency

DAP_proteins = unique(Correlations$DAP)
Final.df = c()

for(p in DAP_proteins){
  df.p = Correlations[Correlations$DAP == p,]
  p.comp = human_proteincomplex[grepl(paste0('\\b',p,'\\b'),human_proteincomplex$subunits.Gene.name.),]
  for (i in 1:length(rownames(df.p))) {
    B.comp = human_proteincomplex[grepl(paste0('\\b',df.p$Partner[i],'\\b'),human_proteincomplex$subunits.Gene.name.),]
    together = intersect(p.comp$ComplexID,B.comp$ComplexID)
    atleastone = union(p.comp$ComplexID,B.comp$ComplexID)
    df.p$number_of_complexes_with_protein_A_and_B[i] = length(together)
    df.p$number_of_complexes_with_protein_A_or_B[i] = length(atleastone)}
  Final.df = rbind(Final.df,df.p)}

Final.df$JaccardIndex = Final.df$number_of_complexes_with_protein_A_and_B / Final.df$number_of_complexes_with_protein_A_or_B

## FUNCTIONAL ANNOTATION OF PROTEIN COMPLEXES

# CORUM protein complexes and their GO terms
allcomplexes = read.delim("./Data/allComplexes.txt")
allcomplexes = allcomplexes[allcomplexes$Organism == "Human",]

# Making a table with GO IDs and their description
GO.terms = allcomplexes[,9:10]
GO.terms[,1:2] <- sapply(GO.terms[,1:2], as.character)
GO.terms = separate_rows(GO.terms,1:2,sep = ";")
GO.terms = GO.terms[!duplicated(GO.terms$GO.ID),]
GO.terms = GO.terms[GO.terms$GO.ID != "None",] # Table with GO ID and GO description

# Correlations between aneuploid proteins and their partners 

# Correlation data

load("./Data/SpearmanCorr_DAP&Partners_mrna&protein_13amplificationcases.RData") # Amplifications
#load("./Data/SpearmanCorr_DAP&Partners_mrna&protein_20deletioncases.RData") # Deletions

cases = unique(Cor_DAP_Partners_new$Case)

# Step1: Grouping correlations and retrieveing complexes

Info.df = c() # How many pair after grouping as top, middle and down and how many complexes 
for(case in cases){
  case.df = Cor_DAP_Partners_new[Cor_DAP_Partners_new$Case == case,]
  # Decision point 1: proteins on other chromosomes
  case.df = case.df[case.df$Chr.info == "Other",]
  # Decision point 2: only partner proteins so for a pair/cor we can be sure they are in the same complex
  case.df = case.df[case.df$Cocomplex.member == "YES",]
  case.df = case.df[order(case.df$p.cor, decreasing = TRUE),]
  case.df$Group[case.df$p.cor < 0.2 & case.df$p.cor > -0.2] <- "Middle"
  case.df.top = case.df[case.df$p.cor > 0,]
  case.df.down = case.df[case.df$p.cor < 0,]
  if(nrow(case.df.down) < 20 | nrow(case.df.top) < 20) next
  if(min(case.df.top$p.cor[1:20]) < 0.2 | max(case.df.down$p.cor[(nrow(case.df.down) - 19):nrow(case.df.down)]) > -0.2) next
  case.df$Group[1:20] <- "Top"
  case.df$Group[(nrow(case.df) - 19):nrow(case.df)] <- "Down"
  # Top 20 and down 20
  case.df = case.df[complete.cases(case.df),]
  groups = unique(case.df$Group)
  for(g in groups){
    df = case.df[case.df$Group == g,]
    proteins = union(df$Partner,df$DAP)
    complexes = c()
    for(p in proteins){complexes = rbind(complexes,allcomplexes[grep(p,as.character(allcomplexes$subunits.Gene.name.)),])}
    complexes = complexes[!duplicated(complexes$ComplexID),]
    assign(g,complexes)}
  df.lists = list("Top" = Top, "Middle" = Middle, "Down" = Down)
  assign(case, df.lists)
  line = c(case, as.numeric(table(case.df$Group)), nrow(Down), nrow(Middle), nrow(Top))
  Info.df = rbind(Info.df, line)}
colnames(Info.df) = c("Case","Number.of.pair.Down","Number.of.pair.Middle","Number.of.pair.Top",
                      "Number.of.complex.Down","Number.of.complex.Middle","Number.of.complex.Top")
Info.df = as.data.frame(Info.df)
Info.df[,1:7] <- sapply(Info.df[,1:7], as.character)
Info.df[,2:7] <- sapply(Info.df[,2:7], as.numeric)

# Step 2: Comparing enriched GO terms in top 20 + down 20 correlations pooled to those in the background correlations

for(l in 1:length(case.lists)){
  case.dfs = case.lists[[l]]
  df1 = rbind(case.dfs[["Top"]], case.dfs[["Down"]])
  df1_go = unlist(lapply(df1$GO.ID, function(x) strsplit(as.character(x),";")[[1]]))
  df1_go = unique(df1_go)
  df2 = case.dfs[["Middle"]]
  df2_go = unlist(lapply(df2$GO.ID, function(x) strsplit(as.character(x),";")[[1]]))
  df2_go = unique(df2_go)
  common_goterms = intersect(df1_go,df2_go)
  common_goterms = common_goterms[common_goterms != "None"]
  for(j in 1:length(common_goterms)){
    GO = common_goterms[j]
    GO.m = matrix(c(nrow(df1[grep(GO,as.character(df1$GO.ID)),]),
                    (nrow(df1) - nrow(df1[grep(GO,as.character(df1$GO.ID)),])),
                    nrow(df2[grep(GO,as.character(df2$GO.ID)),]),
                    (nrow(df2) - nrow(df2[grep(GO,as.character(df2$GO.ID)),]))), nrow = 2, ncol = 2)
    colnames(GO.m) = comparisons[[i]]
    rownames(GO.m) = c("wGO","woGO")
    res = chisq.test(GO.m)
    line = c(names(case.lists)[l],"TopDown;Middle", GO, res$p.value, res$stdres[1], res$stdres[3]) #stdres for first row
    Stat = rbind(Stat, line)}}

colnames(Stat) = c("Case","Comparison","GO.ID","pvalue","Stdres.w.TopDown","Stdres.w.Middle")
Stat = as.data.frame(Stat)
Stat[,1:6] <- sapply(Stat[,1:6], as.character)
Stat[,4:6] <- sapply(Stat[,4:6], as.numeric)
Stat = merge(Stat, GO.terms, by = "GO.ID", all.x = TRUE)

## CALCULATION OF THE STOICHIOMETRY SCORE AND SURVIVAL ANALYSIS

# PART I: Calculating the mean residuals for the samples (The stoichiometry deviation score)

# Correlation data (Amplification cases)

load("./Data/SpearmanCorr_DAP&Partners_mrna&protein_13amplificationcases.RData") 
c.types = unique(gsub("_.*","",Cor_DAP_Partners_new$Case))
Cor_DAP_Partners_new$Type = gsub("_.*","",Cor_DAP_Partners_new$Case)

all.data = list()
for(type in c.types){
  #proteomic measurements
  assign("pcounts_primary",get(load(paste0("./Data/CPTAC_TCGA_Intersection/",type,"_pcounts.RData"))))
  pcounts_primary = as.data.frame(pcounts_primary)
  # Correlation data
  df.type = Cor_DAP_Partners_new[Cor_DAP_Partners_new$Type == type,]
  df.type = df.type[df.type$Chr.info == "Other" & df.type$Cocomplex.member == "YES",]
  df.type.ordered = df.type[order(abs(df.type$p.cor), decreasing = TRUE),]
  df.top = df.type.ordered[1:30,] #The top 30 strongest correlation
  # Linear regression model
  allres.df = as.data.frame(colnames(pcounts_primary))
  colnames(allres.df) = "samples"
  for(r in 1:30){
    dap = df.top$DAP[r]
    partner = df.top$Partner[r]
    m = t(pcounts_primary[c(partner,dap),])
    colnames(m) = c("y","x")
    m = as.data.frame(m)
    # Partner is dependent variable and aneuploid protein is independent variable
    model = lm(y~x,data = m)
    res.df = as.data.frame(model$residuals)
    colnames(res.df) = df.top$Pair[r]
    res.df$samples = rownames(res.df)
    allres.df = merge(allres.df,res.df, by = "samples", all.x = TRUE)}
  name = paste(type,paste0("top","30"),sep = ".")
  all.data[[name]] = allres.df}

Final.df = c()

for(i in c("BRCA","COREAD","OV")){
  sub.data = all.data[grep(i,names(all.data))]
  type.df = c()
  cnames = c()
  for(j in 1:length(sub.data)){
    df = sub.data[[j]]
    df.abs = abs(df[,-1])
    avg = rowMeans(df.abs, na.rm = TRUE)
    type.df = cbind(type.df,avg)
    cnames = c(cnames, strsplit(names(sub.data)[j],"\\.")[[1]][2])}
  type.df = as.data.frame(type.df)
  colnames(type.df) = cnames
  type.df$Samples = df$samples
  type.df$Tissue = i
  Final.df = rbind(Final.df,type.df)}

# Clinical data from cBioPortal (TCGA Pancancer)
clinicaldata.brca = read.delim("./Data/Clinical_data/brca_tcga_pan_can_atlas_2018_clinical_data.tsv")
clinicaldata.brca = clinicaldata.brca[,c(3,10,16,17,35,36)]

clinicaldata.ov = read.delim("./Data/Clinical_data/ov_tcga_pan_can_atlas_2018_clinical_data.tsv")
clinicaldata.ov = clinicaldata.ov[,c(3,10,16,17,35,36)]

clinicaldata.coread = read.delim("./Data/Clinical_data/coadread_tcga_pan_can_atlas_2018_clinical_data.tsv")
clinicaldata.coread = clinicaldata.coread[,c(3,10,16,17,35,36)]

clinical.data = rbind(clinicaldata.brca, clinicaldata.coread,
                      clinicaldata.ov)

clinical.data$status = as.numeric(unlist(lapply(as.character(clinical.data$Overall.Survival.Status), function(x) strsplit(x,":")[[1]][1])))
clinical.data$status.diseasefree = as.numeric(unlist(lapply(as.character(clinical.data$Disease.Free.Status), function(x) strsplit(x,":")[[1]][1])))

m = merge(Final.df, clinical.data, by.x = "Samples", by.y = "Sample.ID")

# Grouping samples - By using Surv_cut package

tissues = unique(m$Tissue)

for(tissue in tissues){
  df = m[m$Tissue == tissue,]
  
  # Overall Survival
  res.cut = surv_cutpoint(df,time = "Overall.Survival..Months.", event = "status", variables = "top30")
  res.cat <- surv_categorize(res.cut)

  fit = survfit(Surv(Overall.Survival..Months., status) ~ top30, data = res.cat)
  png(filename = paste0("./","Survplot_survcut_os","_",tissue,".png"))
  print(ggsurvplot(fit, data = res.cat, risk.table = TRUE, conf.int = TRUE,pval = TRUE))
  dev.off()
  
  # Disease free survival
  res.cut = surv_cutpoint(df,time = "Disease.Free..Months.", event = "status.diseasefree", variables = "top30")
  res.cat <- surv_categorize(res.cut)

  fit = survfit(Surv(Disease.Free..Months., status.diseasefree) ~ top30, data = res.cat)
  png(filename = paste0("./","Survplot_survcut_dfs","_",tissue,".png"))
  print(ggsurvplot(fit, data = res.cat, risk.table = TRUE, conf.int = TRUE,pval = TRUE))
  dev.off()}
