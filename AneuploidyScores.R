# Calculation whole-chromosome level aneuploidy scores

library(openxlsx)

url = "https://ars.els-cdn.com/content/image/1-s2.0-S1535610818301119-mmc2.xlsx" # Arm-level aneuploidy scores from Taylor et al., 2018
destfile = "~/Aneuploidyscores.xlsx"
download.file(url, destfile = destfile)

supp_mmc2 = read.xlsx("~/Aneuploidyscores.xlsx", 
                      rowNames = FALSE, colNames = TRUE, startRow = 2)
Aneuploidy_table = data.frame("Sample" = supp_mmc2$Sample,
                              "Type" = supp_mmc2$Type)
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
save(Aneuploidy_table,file = "~/AneuploidyTable.RData")

# Detecting cancer type-specific whole chromosome-level aneuploidies
load("~/AneuploidyTable.RData")
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
    #Chi-Square for Chromosome of interest (ChrI). ChrO stands for other chromosomes
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