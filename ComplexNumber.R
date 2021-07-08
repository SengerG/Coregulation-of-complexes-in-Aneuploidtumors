# Find the number of complexes a protein involves and number of partners (Promiscuous vs non-promiscuous proteins)

# Correlation Data

load("~/SpearmanCorr_mrna&protein.RData")

#Human protein complexes - CORUM Database

load("~/HumanPCgenes_Subunits_allchr.RData")

#Find the number of complexes a protein involves and number of partners

for(i in 1:length(rownames(Correlations))){
  p = as.character(unlist(Correlations$DAP))[i]
  p.comp = human_proteincomplex[grepl(paste0('\\b',p,'\\b'),human_proteincomplex$subunits.Gene.name.),]
  p.partners = strsplit(Subunits.Partner.Info[Subunits.Partner.Info$Subunit == p,"Partners"],";")[[1]]
  Correlations$Number.complex[i] = length(rownames(p.comp))
  Correlations$Number.partner[i] = length(p.partners)}

Correlations$Promiscuity = ifelse(Correlations$Number.complex > 5,"Promiscuous","Non-promiscuous")