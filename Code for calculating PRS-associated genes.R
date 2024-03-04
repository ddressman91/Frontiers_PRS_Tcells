##Run association of PRS with gene expression

library(ggplot2)
library(stringr)
library(SummarizedExperiment)
library(edgeR)
library(SCeQTL)
library(dplyr)
library(data.table)

#Read in key file for FID to projid conversion for ROSMAP_n381
keyfile<-fread("/chop_1k_iid_projid.txt",fill=TRUE)
keyfile<-keyfile[,1:2]
colnames(keyfile)<-c("FID","projid")

#Read in all_score files, convert all IDs to projid, merge ROSMAP data into one table
AD.Kunkle<-fread("ROSMAP.all.noAPOE.finalQC.AD.Kunkle.PRS.all_score",fill=TRUE)
AD.Kunkle.381<-AD.Kunkle %>% filter(nchar(FID)==9)
AD.Kunkle.other<-AD.Kunkle %>% filter(nchar(FID)!=9)
AD.Kunkle.381<-inner_join(AD.Kunkle.381,keyfile,by="FID")
AD.Kunkle.other$projid<-str_replace(AD.Kunkle.other$FID,"MAP","")
AD.Kunkle.other$projid<-str_replace(AD.Kunkle.other$projid,"ROS","")
AD.Kunkle.other$projid<-as.character(AD.Kunkle.other$projid)
AD.Kunkle.381$projid<-as.character(AD.Kunkle.381$projid)
AD.Kunkle.all<-bind_rows(AD.Kunkle.other,AD.Kunkle.381)
AD.Kunkle.all$projid<-str_pad(AD.Kunkle.all$projid,8,pad="0")

#Standardize PRS distribution to be used for correlation with gene expression
AD.Kunkle.all$PRS_5e8_standard <- scale(AD.Kunkle.all$`Pt_5e-08`,center=T,scale=T)
AD.Kunkle.all$PRS_1e7_standard <- scale(AD.Kunkle.all$`Pt_1e-07`,center=T,scale=T)
AD.Kunkle.all$PRS_1e6_standard <- scale(AD.Kunkle.all$`Pt_1e-06`,center=T,scale=T)
AD.Kunkle.all$PRS_1e5_standard <- scale(AD.Kunkle.all$`Pt_1e-05`,center=T,scale=T)
AD.Kunkle.all$PRS_1e4_standard <- scale(AD.Kunkle.all$`Pt_0.0001`,center=T,scale=T)
AD.Kunkle.all$PRS_1e3_standard <- scale(AD.Kunkle.all$`Pt_0.001`,center=T,scale=T)
AD.Kunkle.all$PRS_0.05_standard <- scale(AD.Kunkle.all$`Pt_0.05`,center=T,scale=T)
AD.Kunkle.all$PRS_0.25_standard <- scale(AD.Kunkle.all$`Pt_0.25`,center=T,scale=T)
AD.Kunkle.all$PRS_0.5_standard <- scale(AD.Kunkle.all$`Pt_0.5`,center=T,scale=T)
AD.Kunkle.all$PRS_0.75_standard <- scale(AD.Kunkle.all$`Pt_0.75`,center=T,scale=T)
AD.Kunkle.all$PRS_1_standard <- scale(AD.Kunkle.all$`Pt_1`,center=T,scale=T)

#Check shape of PRS distribution
ggplot(AD.Kunkle.all,aes(AD.Kunkle.all$PRS_5e8_standard))+geom_histogram()
ggplot(AD.Kunkle.all,aes(AD.Kunkle.all$PRS_1e7_standard))+geom_histogram()
ggplot(AD.Kunkle.all,aes(AD.Kunkle.all$PRS_1e6_standard))+geom_histogram()
ggplot(AD.Kunkle.all,aes(AD.Kunkle.all$PRS_1e5_standard))+geom_histogram()
ggplot(AD.Kunkle.all,aes(AD.Kunkle.all$PRS_1e4_standard))+geom_histogram()
ggplot(AD.Kunkle.all,aes(AD.Kunkle.all$PRS_1e3_standard))+geom_histogram()
ggplot(AD.Kunkle.all,aes(AD.Kunkle.all$PRS_0.05_standard))+geom_histogram()
ggplot(AD.Kunkle.all,aes(AD.Kunkle.all$PRS_0.25_standard))+geom_histogram()
ggplot(AD.Kunkle.all,aes(AD.Kunkle.all$PRS_0.5_standard))+geom_histogram()
ggplot(AD.Kunkle.all,aes(AD.Kunkle.all$PRS_0.75_standard))+geom_histogram()
ggplot(AD.Kunkle.all,aes(AD.Kunkle.all$PRS_1_standard))+geom_histogram()

#Join PRS file with ID's from plate map
platemap<-fread("/ROSMAP_Tcell_DGE_PlateMap_20170613.csv",fill=TRUE)
platemap<-unique(platemap,by="Sample Identifier")
platemap$projid<-platemap$`Sample Identifier`
platemap$projid<-str_pad(platemap$projid,8,pad="0")
platemaptrim<-platemap %>% select('Sample Identifier', projid)
AD.Kunklernaseq<-inner_join(AD.Kunkle.all,platemaptrim,by="projid")
AD.Kunklernaseq<-AD.Kunklernaseq[,1:(ncol(AD.Kunklernaseq)-1)]

#Check shape of PRS distribution for sequencing samples
ggplot(AD.Kunklernaseq,aes(AD.Kunklernaseq$PRS_5e8_standard))+geom_histogram()
ggplot(AD.Kunklernaseq,aes(AD.Kunklernaseq$PRS_1e7_standard))+geom_histogram()
ggplot(AD.Kunklernaseq,aes(AD.Kunklernaseq$PRS_1e6_standard))+geom_histogram()
ggplot(AD.Kunklernaseq,aes(AD.Kunklernaseq$PRS_1e5_standard))+geom_histogram()
ggplot(AD.Kunklernaseq,aes(AD.Kunklernaseq$PRS_1e4_standard))+geom_histogram()
ggplot(AD.Kunklernaseq,aes(AD.Kunklernaseq$PRS_1e3_standard))+geom_histogram()
ggplot(AD.Kunklernaseq,aes(AD.Kunklernaseq$PRS_0.05_standard))+geom_histogram()
ggplot(AD.Kunklernaseq,aes(AD.Kunklernaseq$PRS_0.25_standard))+geom_histogram()
ggplot(AD.Kunklernaseq,aes(AD.Kunklernaseq$PRS_0.5_standard))+geom_histogram()
ggplot(AD.Kunklernaseq,aes(AD.Kunklernaseq$PRS_0.75_standard))+geom_histogram()
ggplot(AD.Kunklernaseq,aes(AD.Kunklernaseq$PRS_1_standard))+geom_histogram()

#Read expression counts file, filter to retain genes with at least 20% expression and max at least 3
scrnaseq<-read.table("SSF-11688.HYLTJBGX2.unq.refseq.umi.dat",header=T,as.is=T)
colnames(scrnaseq)<-gsub("T384s1_","",colnames(scrnaseq))
scrnaseq <- scrnaseq[rowSums(scrnaseq == 0) <= 307, ]
scrnaseq$Max <- apply(X = scrnaseq, MARGIN = 1, FUN = max)
scrnaseq <- scrnaseq %>% filter(Max >= 3)
scrnaseq <- scrnaseq[,1:384]

map<-read.csv("/ROSMAP_Tcell_DGE_PlateMap_20170613.csv",header=T,as.is=T)
map[,"ID"]<-paste(map$Sample.Identifier,map$Cell.or.Tissue.Type,sep="_")

dge <- DGEList(counts=scrnaseq, samples=map)
dge.norm<-calcNormFactors(dge)
dge.norm$samples$Sample.Identifier<-str_pad(dge.norm$samples$Sample.Identifier, 
                                            8, pad = "0")
rownames(dge.norm$counts)<-gsub("-","_",rownames(dge.norm$counts))

dge.norm.c1<-dge.norm[,dge.norm$samples$Cell.or.Tissue.Type=="CD4+CD45RO-"]
colnames(dge.norm.c1$counts)<-dge.norm.c1$samples$Sample.Identifier

#############get phenotypes for covariates
full.pheno<-read.table("/rosmap_xsect_slopes_n3328_012618.csv",
                       header=T,as.is=T,sep=",")
full.pheno[,"ID"]<-paste(str_pad(full.pheno$projid, 8 ,pad="0"),
                         str_pad(full.pheno$projid, 8 ,pad="0"),sep="_")
full.pheno[,"clinicalAD"]<-ifelse(full.pheno$dcfdx_lv %in% c(4,5,6),1,0)

full.pheno$projid<-str_pad(full.pheno$projid, 8, pad = "0")

covar.pheno<-subset(full.pheno,select=c(projid,msex,age_bl,clinicalAD,pathoAD,EV1,EV2,EV3,EV4,EV5,EV6,EV7,EV8,EV9,EV10))
AD.Kunkle.prs<-merge(AD.Kunklernaseq,covar.pheno,by="projid")

##################################################################

##run PRS association with each cell type
##c1.prs contains expression (samples in rows, genes in columns)
##followed by PRS and phenotypes
##use the PRS for all variants (P<1), column X1 to run QTL analysis
##if needed the PRS QTL analysis can be run for PRS computed using variants at different p-value cutoffs

c1.expr<-as.data.frame(t(dge.norm.c1$counts))
c1.expr[,"Sample"]<-dge.norm.c1$samples$Sample.Identifier
c1.prs<-merge(c1.cases,AD.Kunkle.prs.cases,by.x="Sample",by.y="projid")
tail(colnames(c1.prs),50)

c1.qtl<-c()
for(i in 2:(ncol(c1.expr)))
{
  t<-summary(lm(as.formula(paste(colnames(c1.prs)[i],"~", colnames(c1.prs)[ncol(c1.prs)-14], "+",
                                 colnames(c1.prs)[ncol(c1.prs)-13],"+",
                                 colnames(c1.prs)[ncol(c1.prs)-12],"+",
                                 colnames(c1.prs)[ncol(c1.prs)-11],"+",
                                 colnames(c1.prs)[ncol(c1.prs)-10],"+",
                                 colnames(c1.prs)[ncol(c1.prs)-9],"+",
                                 colnames(c1.prs)[ncol(c1.prs)-8],"+",
                                 colnames(c1.prs)[ncol(c1.prs)-7],"+",
                                 colnames(c1.prs)[ncol(c1.prs)-6],"+",
                                 colnames(c1.prs)[ncol(c1.prs)-5],"+",
                                 colnames(c1.prs)[ncol(c1.prs)-4],"+",
                                 colnames(c1.prs)[ncol(c1.prs)-3],"+",
                                 colnames(c1.prs)[ncol(c1.prs)-2],"+",
                                 colnames(c1.prs)[ncol(c1.prs)-1],"+",
                                 colnames(c1.prs)[ncol(c1.prs)],
                                 sep="")), data =c1.prs ))$coef
  c1.qtl<-rbind(c1.qtl,c(colnames(c1.prs)[i],t[2,]))
}

#Convert to data frame, change column names and add columns for
#Bonferroni significance and FDR q-values
c1.qtl <- as.data.frame(c1.qtl)
colnames(c1.qtl)<-c("Gene","Effect.size","SE","t.value","p.value")
c1.qtl$p.value <- as.numeric(c1.qtl$p.value)
c1.qtl$Bonferroni.qval <- p.adjust(c1.qtl$p.value, method = "bonferroni")
c1.qtl$FDR.qval <- p.adjust(c1.qtl$p.value, method = "fdr")

################C2
dge.norm.c2<-dge.norm[,dge.norm$samples$Cell.or.Tissue.Type=="CD4+CD45RO+"]
colnames(dge.norm.c2$counts)<-dge.norm.c2$samples$Sample.Identifier

c2.expr<-as.data.frame(t(dge.norm.c2$counts))
c2.expr[,"Sample"]<-dge.norm.c2$samples$Sample.Identifier
c2.prs<-merge(c2.cases,AD.Kunkle.prs.cases,by.x="Sample",by.y="projid")

c2.qtl<-c()
for(i in 2:(ncol(c2.expr)))
{
  t<-summary(lm(as.formula(paste(colnames(c2.prs)[i],"~", colnames(c2.prs)[ncol(c2.prs)-14], "+",
                                 colnames(c2.prs)[ncol(c2.prs)-13],"+",
                                 colnames(c2.prs)[ncol(c2.prs)-12],"+",
                                 colnames(c2.prs)[ncol(c2.prs)-11],"+",
                                 colnames(c2.prs)[ncol(c2.prs)-10],"+",
                                 colnames(c2.prs)[ncol(c2.prs)-9],"+",
                                 colnames(c2.prs)[ncol(c2.prs)-8],"+",
                                 colnames(c2.prs)[ncol(c2.prs)-7],"+",
                                 colnames(c2.prs)[ncol(c2.prs)-6],"+",
                                 colnames(c2.prs)[ncol(c2.prs)-5],"+",
                                 colnames(c2.prs)[ncol(c2.prs)-4],"+",
                                 colnames(c2.prs)[ncol(c2.prs)-3],"+",
                                 colnames(c2.prs)[ncol(c2.prs)-2],"+",
                                 colnames(c2.prs)[ncol(c2.prs)-1],"+",
                                 colnames(c2.prs)[ncol(c2.prs)],
                                 sep="")), data =c2.prs ))$coef
  c2.qtl<-rbind(c2.qtl,c(colnames(c2.prs)[i],t[2,]))
}

#Convert to data frame, change column names and add columns for
#Bonferroni significance and FDR q-values
c2.qtl <- as.data.frame(c2.qtl)
colnames(c2.qtl)<-c("Gene","Effect.size","SE","t.value","p.value")
c2.qtl$p.value <- as.numeric(c2.qtl$p.value)
c2.qtl$Bonferroni.qval <- p.adjust(c2.qtl$p.value, method = "bonferroni")
c2.qtl$FDR.qval <- p.adjust(c2.qtl$p.value, method = "fdr")

################C3
dge.norm.c3<-dge.norm[,dge.norm$samples$Cell.or.Tissue.Type=="CD8+CD45RO-"]
colnames(dge.norm.c3$counts)<-dge.norm.c3$samples$Sample.Identifier

c3.expr<-as.data.frame(t(dge.norm.c3$counts))
c3.expr[,"Sample"]<-dge.norm.c3$samples$Sample.Identifier
c3.prs<-merge(c3.cases,AD.Kunkle.prs.cases,by.x="Sample",by.y="projid")

c3.qtl<-c()
for(i in 2:(ncol(c3.expr)))
{
  t<-summary(lm(as.formula(paste(colnames(c3.prs)[i],"~", colnames(c3.prs)[ncol(c3.prs)-14], "+",
                                 colnames(c3.prs)[ncol(c3.prs)-13],"+",
                                 colnames(c3.prs)[ncol(c3.prs)-12],"+",
                                 colnames(c3.prs)[ncol(c3.prs)-11],"+",
                                 colnames(c3.prs)[ncol(c3.prs)-10],"+",
                                 colnames(c3.prs)[ncol(c3.prs)-9],"+",
                                 colnames(c3.prs)[ncol(c3.prs)-8],"+",
                                 colnames(c3.prs)[ncol(c3.prs)-7],"+",
                                 colnames(c3.prs)[ncol(c3.prs)-6],"+",
                                 colnames(c3.prs)[ncol(c3.prs)-5],"+",
                                 colnames(c3.prs)[ncol(c3.prs)-4],"+",
                                 colnames(c3.prs)[ncol(c3.prs)-3],"+",
                                 colnames(c3.prs)[ncol(c3.prs)-2],"+",
                                 colnames(c3.prs)[ncol(c3.prs)-1],"+",
                                 colnames(c3.prs)[ncol(c3.prs)],
                                 sep="")), data =c3.prs ))$coef
  c3.qtl<-rbind(c3.qtl,c(colnames(c3.prs)[i],t[2,]))
}

#Convert to data frame, change column names and add columns for
#Bonferroni significance and FDR q-values
c3.qtl <- as.data.frame(c3.qtl)
colnames(c3.qtl)<-c("Gene","Effect.size","SE","t.value","p.value")
c3.qtl$p.value <- as.numeric(c3.qtl$p.value)
c3.qtl$Bonferroni.qval <- p.adjust(c3.qtl$p.value, method = "bonferroni")
c3.qtl$FDR.qval <- p.adjust(c3.qtl$p.value, method = "fdr")

###############C4
dge.norm.c4<-dge.norm[,dge.norm$samples$Cell.or.Tissue.Type=="CD8+CD45RO+"]
colnames(dge.norm.c4$counts)<-dge.norm.c4$samples$Sample.Identifier

c4.expr<-as.data.frame(t(dge.norm.c4$counts))
c4.expr[,"Sample"]<-dge.norm.c4$samples$Sample.Identifier
c4.prs<-merge(c4.cases,AD.Kunkle.prs.cases,by.x="Sample",by.y="projid")

c4.qtl<-c()
for(i in 2:(ncol(c4.expr)))
{
  t<-summary(lm(as.formula(paste(colnames(c4.prs)[i],"~", colnames(c4.prs)[ncol(c4.prs)-14], "+",
                                 colnames(c4.prs)[ncol(c4.prs)-13],"+",
                                 colnames(c4.prs)[ncol(c4.prs)-12],"+",
                                 colnames(c4.prs)[ncol(c4.prs)-11],"+",
                                 colnames(c4.prs)[ncol(c4.prs)-10],"+",
                                 colnames(c4.prs)[ncol(c4.prs)-9],"+",
                                 colnames(c4.prs)[ncol(c4.prs)-8],"+",
                                 colnames(c4.prs)[ncol(c4.prs)-7],"+",
                                 colnames(c4.prs)[ncol(c4.prs)-6],"+",
                                 colnames(c4.prs)[ncol(c4.prs)-5],"+",
                                 colnames(c4.prs)[ncol(c4.prs)-4],"+",
                                 colnames(c4.prs)[ncol(c4.prs)-3],"+",
                                 colnames(c4.prs)[ncol(c4.prs)-2],"+",
                                 colnames(c4.prs)[ncol(c4.prs)-1],"+",
                                 colnames(c4.prs)[ncol(c4.prs)],
                                 sep="")), data =c4.prs ))$coef
  c4.qtl<-rbind(c4.qtl,c(colnames(c4.prs)[i],t[2,]))
}

#Convert to data frame, change column names and add columns for
#Bonferroni significance and FDR q-values
c4.qtl <- as.data.frame(c4.qtl)
colnames(c4.qtl)<-c("Gene","Effect.size","SE","t.value","p.value")
c4.qtl$p.value <- as.numeric(c4.qtl$p.value)
c4.qtl$Bonferroni.qval <- p.adjust(c4.qtl$p.value, method = "bonferroni")
c4.qtl$FDR.qval <- p.adjust(c4.qtl$p.value, method = "fdr")

#Save files of PRS-associated genes
setwd("C:/Users/dalli/Documents/Elyaman lab/PRS paper/PRS associations/Associations.after.final.PRS.QC/")
write.table(c1.qtl,"c1.AD.Kunkle.prs.genes.covar.qval",col.names=T,row.names=F,sep="\t",quote=F)
write.table(c2.qtl,"c2.AD.Kunkle.prs.genes.covar.qval",col.names=T,row.names=F,sep="\t",quote=F)
write.table(c3.qtl,"c3.AD.Kunkle.prs.genes.covar.qval",col.names=T,row.names=F,sep="\t",quote=F)
write.table(c4.qtl,"c4.AD.Kunkle.prs.genes.covar.qval",col.names=T,row.names=F,sep="\t",quote=F)