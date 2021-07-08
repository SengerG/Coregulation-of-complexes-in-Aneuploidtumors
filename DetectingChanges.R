
library(openxlsx)

load("~/AneuploidyTable.RData") # This file can be produced by running the AneuploidyScores.R or can be downloaded from the repository
cases_slctd = read.xlsx("~/CasesSelected_alltumours_amp.xlsx") # 86 selected amplifications after multiple testing correction

# Detecting transcriptomic changes induced by amplifications 

c.types = unique(cases_slctd$Cancer.Type)
Stat = c()

for (type in c.types) {
  assign("tpm_primary",get(load(paste0("~/",type,"_tpm_primary_genenames.RData")))) # Output files produced by transcriptome data processing code
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

df.tobesaved = mget(ls(pattern = "^wilcox_"))
save(df.tobesaved, file = "~/DE_Transcriptome.RData")

# Detecting proteomic changes induced by amplifications

cases_slctd = cases_slctd[cases_slctd$Cancer.Type %in% c("BRCA","OV","COREAD"),] #13 amplification cases for which proteomic data is available
c.types = c("COREAD", "BRCA","OV") 
path = "~/" # Path output files produced by proteome data processing code
sample.info = c()
final.dataset = list()
for (type in c.types) {
  assign("pcounts_primary",get(load(paste0(path,type,"_pcounts.RData")))) # Output files produced by proteome data processing code
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
    final.dataset[[name]] = data_stat}}
save(final.dataset, file = "~/DE_Proteome.RData")