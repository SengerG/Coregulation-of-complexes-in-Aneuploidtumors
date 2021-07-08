
# For COREAD, OV, and BRCA primary tumour samples and genes (covered by both transcriptomic and proteomic data), all possible pairwise correlations at both mRNA 
# and proteomic levels were calculated. Then we focused on the correlations between proteins that encoded by genes found on amplified chromosome and 
# show significant increases in abundances and their partners encoded by genes found on other chromosomes.

#Function to make df for correlation results
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut])}


#Required Libraries
library(Hmisc)
library(openxlsx)

#Analysis
c.types = c("COREAD","OV","BRCA")
for (type in c.types) {
  #transcriptomic measurements
  assign("tpm_primary",get(load(paste0("~/",type,"_tpm_primary_genenames.RData")))) # Output files produced by transcriptome data processing code
  #proteomic measurements
  assign("pcounts_primary",get(load(paste0("~/",type,"_pcounts.RData")))) # Output files produced by proteome data processing code
  pcounts_primary = as.data.frame(pcounts_primary)
  #Make sure transcriptome data cover the same samples and genes
  tpm_primary = tpm_primary[rownames(pcounts_primary),colnames(pcounts_primary)]
  #Calculating correlations
  #mRNA
  mm.t = t(tpm_primary)
  res.t <- rcorr(as.matrix(mm.t), type = "spearman")
  Corr.mat.t = flattenCorrMatrix(res.t$r, res.t$P)
  #Protein
  mm.p = t(pcounts_primary)
  res.p <- rcorr(as.matrix(mm.p), type = "spearman")
  Corr.mat.p = flattenCorrMatrix(res.p$r, res.p$P)
  print(type)
  save(Corr.mat.p,Corr.mat.t, 
       file = paste0("~/SpearmanCorr_allpairs_mrna&protein_",type,".RData"))}