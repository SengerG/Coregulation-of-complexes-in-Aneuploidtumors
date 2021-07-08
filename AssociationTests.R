
library(openxlsx)

load("~/HumanPCgenes_Subunits_allchr.RData")

# 1. Transcriptome Level Analyses

## 1.1. Overlap with complex subunits

load("~/DE_Transcriptome_Sig.RData")

Info.table = read.xlsx("~/Table_cutoff.xlsx",colNames = TRUE,rowNames = FALSE) # File for which cut-off were used to detect DEGs in each aneuploidy case

padj_0.1 = Info.table[Info.table$Sig == "padj<0.1",] # 64
pvalue_0.05 = Info.table[Info.table$Sig == "pvalue<0.05",] # 22

cases = Info.table$Case
Summary = c()

for (case in cases) {
  if(case %in% padj_0.1$Case){
    sig = "0.1"
    val = "p.adj"}
  else{
    sig = "0.05"
    val = "p.value"}
  name = paste("wilcox",case,"amp_vs_noch",sig,sep = "_")
  df = df.tobesaved.sig[[which(names(df.tobesaved.sig)==name)]]
  df = df[complete.cases(df),]
  df = df[df$Gene.type == "protein_coding",]
  df = df[df$Chr.Info == "ChrO",]
  if(val == "p.adj"){
    ChrO_DE = df[df$Chr.Info == "ChrO" & df$DE.info.padj == "DE",]
    if(dim(ChrO_DE)[1] == 0)next
    ChrO_DE_subunit = ChrO_DE[as.character(ChrO_DE$Gene.name) %in% as.character(Subunits_Info_pc_genes$Gene.name),]
    ChrO_NotDE = df[!(df$Gene.name %in% ChrO_DE$Gene.name),]
    ChrO_NotDE_subunit = ChrO_NotDE[as.character(ChrO_NotDE$Gene.name) %in% as.character(Subunits_Info_pc_genes$Gene.name),]}
  else{
    ChrO_DE = df[df$Chr.Info == "ChrO" & df$DE.info.pvalue == "DE",]
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

save(Summary, file = "~/Transcriptome_ComplexSubunitOverlap.RData")

## 1.2. Overlap with partners

Summary = c()

for (case in cases) {
  if(case %in% padj_0.1$Case){
    sig = "0.1"
    val = "p.adj"}
  else{
    sig = "0.05"
    val = "p.value"}
  name = paste("wilcox",case,"amp_vs_noch",sig,sep = "_")
  df = df.tobesaved.sig[[which(names(df.tobesaved.sig)==name)]]
  df = df[complete.cases(df),]
  df = df[df$Gene.type == "protein_coding",]
  if(val == "p.adj"){df_ChrI_DE = df[as.character(df$Chr.Info) == "ChrI" & as.character(df$DE.info.padj) == "DE",]} #DEGs on aneuploid chromosomes
  else{df_ChrI_DE = df[as.character(df$Chr.Info) == "ChrI" & as.character(df$DE.info.pvalue) == "DE",]}
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
    if(dim(data[as.character(data$DE.info.padj) == "DE",])[1] < 10)next
    mm = matrix(c(dim(partner.df[as.character(partner.df$DE.info.padj) == "DE",])[1],
                  dim(partner.df[as.character(partner.df$DE.info.padj) == "NotDE",])[1],
                  dim(notpartner.df[as.character(notpartner.df$DE.info.padj) == "DE",])[1],
                  dim(notpartner.df[as.character(notpartner.df$DE.info.padj) == "NotDE",])[1]), ncol = 2)}
  else{
    #Check point
    if(dim(data[as.character(data$DE.info.pvalue) == "DE",])[1] < 10)next
    mm = matrix(c(dim(partner.df[as.character(partner.df$DE.info.pvalue) == "DE",])[1],
                  dim(partner.df[as.character(partner.df$DE.info.pvalue) == "NotDE",])[1],
                  dim(notpartner.df[as.character(notpartner.df$DE.info.pvalue) == "DE",])[1],
                  dim(notpartner.df[as.character(notpartner.df$DE.info.pvalue) == "NotDE",])[1]), ncol = 2)}
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
save(Summary, file = "~/Transcriptome_PartnersOverlap.RData")

# 2. Proteome Level Analyses

## 2.1. Overlap with complex subunits

load("~/DE_Proteome_Sig.RData")
cases_slctd = read.xlsx("~/CasesSelected_alltumours_amp.xlsx")

cases_slctd = cases_slctd[cases_slctd$Cancer.Type %in% c("BRCA","OV","COREAD"),] #13 amplification cases for which proteomic data is available

cases = paste(cases_slctd$Cancer.Type, cases_slctd$Chromosomal.Aneuploidy, sep = "_")
Summary = c()

for (case in cases) {
  sig = "0.1"
  val = "p.value"
  name = paste("wilcox",case,"amp_vs_noch",sig,sep = "_")
  df = final.dataset.sig[[which(names(final.dataset.sig)==name)]]
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

save(Summary, file = "~/Proteome_ComplexSubunitOverlap.RData")

## 2.2. Overlap with partners

Summary = c()

for (case in cases) {
  sig = "0.1"
  val = "p.value"
  name = paste("wilcox",case,"amp_vs_noch",sig,sep = "_")
  df = final.dataset.sig[[which(names(final.dataset.sig)==name)]]
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
save(Summary, file = "~/Proteome_PartnersOverlap.RData")