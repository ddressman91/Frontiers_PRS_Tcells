#Load AD PRS file and packages
library(pROC)
AD.Kunklepheno <- read.csv(file="AD new PRS with pheno.csv")

#Create ROC objects for AD PRS at each significance level
clinicalAD.roc.5e8<-AD.Kunklepheno %>% roc(clinicalAD,PRS_5e8_standard)
clinicalAD.roc.1e7<-AD.Kunklepheno %>% roc(clinicalAD,PRS_1e7_standard)
clinicalAD.roc.1e6<-AD.Kunklepheno %>% roc(clinicalAD,PRS_1e6_standard)
clinicalAD.roc.1e5<-AD.Kunklepheno %>% roc(clinicalAD,PRS_1e5_standard)
clinicalAD.roc.1e4<-AD.Kunklepheno %>% roc(clinicalAD,PRS_1e4_standard)
clinicalAD.roc.1e3<-AD.Kunklepheno %>% roc(clinicalAD,PRS_1e3_standard)
clinicalAD.roc.05<-AD.Kunklepheno %>% roc(clinicalAD,PRS_0.05_standard)
clinicalAD.roc.25<-AD.Kunklepheno %>% roc(clinicalAD,PRS_0.25_standard)
clinicalAD.roc.5<-AD.Kunklepheno %>% roc(clinicalAD,PRS_0.5_standard)
clinicalAD.roc.75<-AD.Kunklepheno %>% roc(clinicalAD,PRS_0.75_standard)
clinicalAD.roc.1<-AD.Kunklepheno %>% roc(clinicalAD,PRS_1_standard)

pathoAD.roc.5e8<-AD.Kunklepheno %>% roc(pathoAD,PRS_5e8_standard)
pathoAD.roc.1e7<-AD.Kunklepheno %>% roc(pathoAD,PRS_1e7_standard)
pathoAD.roc.1e6<-AD.Kunklepheno %>% roc(pathoAD,PRS_1e6_standard)
pathoAD.roc.1e5<-AD.Kunklepheno %>% roc(pathoAD,PRS_1e5_standard)
pathoAD.roc.1e4<-AD.Kunklepheno %>% roc(pathoAD,PRS_1e4_standard)
pathoAD.roc.1e3<-AD.Kunklepheno %>% roc(pathoAD,PRS_1e3_standard)
pathoAD.roc.05<-AD.Kunklepheno %>% roc(pathoAD,PRS_0.05_standard)
pathoAD.roc.25<-AD.Kunklepheno %>% roc(pathoAD,PRS_0.25_standard)
pathoAD.roc.5<-AD.Kunklepheno %>% roc(pathoAD,PRS_0.5_standard)
pathoAD.roc.75<-AD.Kunklepheno %>% roc(pathoAD,PRS_0.75_standard)
pathoAD.roc.1<-AD.Kunklepheno %>% roc(pathoAD,PRS_1_standard)

#Plot all significance levels on a ROC curve
plot(clinicalAD.roc.5e8,col="#a6cee3",legacy.axes=TRUE,xlim=c(1,0),
     xlab="False positive rate",ylab="True positive rate")
plot(clinicalAD.roc.1e6,col="#b2df8a",add=TRUE)
plot(clinicalAD.roc.1e5,col="#33a02c",add=TRUE)
plot(clinicalAD.roc.1e4,col="#fb9a99",add=TRUE)
plot(clinicalAD.roc.1e3,col="#e31a1c",add=TRUE)
plot(clinicalAD.roc.05,col="#fdbf6f",add=TRUE)
plot(clinicalAD.roc.25,col="#ff7f00",add=TRUE)
plot(clinicalAD.roc.5,col="#cab2d6", add=TRUE)
plot(clinicalAD.roc.75,col="#6a3d9a", lwd = 4, add=TRUE)
plot(clinicalAD.roc.1,col="#ffff99", add=TRUE)
legend(x = 0.2, y = 0.7, title= "SNP pvalue", legend=c("5e-8",
      "1e-6","1e-5","1e-4","0.001","0.05","0.25","0.5","0.75","1"),
      fill=c("#a6cee3","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f",
            "#ff7f00","#cab2d6","#6a3d9a","#ffff99"))
      title(main = "Predictive value of PRS against clinical AD", line = 3)
      
plot(pathoAD.roc.5e8,col="#a6cee3",legacy.axes=TRUE,xlim=c(1,0),
     xlab="False positive rate",ylab="True positive rate")
plot(pathoAD.roc.1e6,col="#b2df8a",add=TRUE)
plot(pathoAD.roc.1e5,col="#33a02c",add=TRUE)
plot(pathoAD.roc.1e4,col="#fb9a99",add=TRUE)
plot(pathoAD.roc.1e3,col="#e31a1c",add=TRUE)
plot(pathoAD.roc.05,col="#fdbf6f",add=TRUE)
plot(pathoAD.roc.25,col="#ff7f00",add=TRUE)
plot(pathoAD.roc.5,col="#cab2d6", add=TRUE)
plot(pathoAD.roc.75,col="#6a3d9a", lwd = 4, add=TRUE)
plot(pathoAD.roc.1,col="#ffff99", add=TRUE)
legend(x = 0.2, y = 0.7, title= "SNP pvalue", legend=c("5e-8",
      "1e-6","1e-5","1e-4","0.001","0.05","0.25","0.5","0.75","1"),
      fill=c("#a6cee3","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f",
            "#ff7f00","#cab2d6","#6a3d9a","#ffff99"))
      title(main = "Predictive value of PRS against patho AD", line = 3)

#Legend hexadecimal 11-color scheme obtained from colorbrewer2.org