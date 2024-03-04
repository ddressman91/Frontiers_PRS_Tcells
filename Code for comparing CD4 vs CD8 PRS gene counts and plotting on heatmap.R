library(dplyr)
library(ComplexHeatmap)
library(data.table)

#Load table of significant PRS-associated genes by trait
sigprs <-read.csv(file = "sig.prs.no.psych.csv")

#Make data frame with numbers of total, positive t-value, and negative t-value PRS-associated genes from each trait
cd4vscd8<-data.frame(
  Trait = c("Lymphocyte counts","WBC counts","C-reactive protein","Ulcerative colitis",
            "Crohn's disease","Multiple sclerosis","Rheumatoid arthritis","Systemic lupus erythematosus",
            "Type 1 diabetes","Alzheimer's disease", "Parkinson's disease",
            "Amyotrophic lateral sclerosis","Epilepsy","Stroke"),
  Total.CD4 = c(sum(!is.na(sigprs$c1.lymph | sigprs$c2.lymph)),sum(!is.na(sigprs$c1.wbc | sigprs$c2.wbc)),
                sum(!is.na(sigprs$c1.crp | sigprs$c2.crp)),sum(!is.na(sigprs$c1.ulccol | sigprs$c2.ulccol)),
                sum(!is.na(sigprs$c1.crohns | sigprs$c2.crohns)),sum(!is.na(sigprs$c1.MS | sigprs$c2.MS)),
                sum(!is.na(sigprs$c1.RA | sigprs$c2.RA)),sum(!is.na(sigprs$c1.SLE | sigprs$c2.SLE)),
                sum(!is.na(sigprs$c1.T1D | sigprs$c2.T1D)),sum(!is.na(sigprs$c1.AD | sigprs$c2.AD)),
                sum(!is.na(sigprs$c1.PD | sigprs$c2.PD)),sum(!is.na(sigprs$c1.ALS | sigprs$c2.ALS)),
                sum(!is.na(sigprs$c1.epilepsy | sigprs$c2.epilepsy)),sum(!is.na(sigprs$c1.stroke | sigprs$c2.stroke))),
  Total.CD8 = c(sum(!is.na(sigprs$c3.lymph | sigprs$c4.lymph)),sum(!is.na(sigprs$c3.wbc | sigprs$c4.wbc)),
                sum(!is.na(sigprs$c3.crp | sigprs$c4.crp)),sum(!is.na(sigprs$c3.ulccol | sigprs$c4.ulccol)),
                sum(!is.na(sigprs$c3.crohns | sigprs$c4.crohns)),sum(!is.na(sigprs$c3.MS | sigprs$c4.MS)),
                sum(!is.na(sigprs$c3.RA | sigprs$c4.RA)),sum(!is.na(sigprs$c3.SLE | sigprs$c4.SLE)),
                sum(!is.na(sigprs$c3.T1D | sigprs$c4.T1D)),sum(!is.na(sigprs$c3.AD | sigprs$c4.AD)),
                sum(!is.na(sigprs$c3.PD | sigprs$c4.PD)),sum(!is.na(sigprs$c3.ALS | sigprs$c4.ALS)),
                sum(!is.na(sigprs$c3.epilepsy | sigprs$c4.epilepsy)),sum(!is.na(sigprs$c3.stroke | sigprs$c4.stroke))),
  Pos.CD4 = c(sum((sigprs$c1.lymph > 0 | sigprs$c2.lymph > 0), na.rm = T),
              sum((sigprs$c1.wbc > 0 | sigprs$c2.wbc > 0), na.rm = T),
              sum((sigprs$c1.crp > 0 | sigprs$c2.crp > 0), na.rm = T),
              sum((sigprs$c1.ulccol > 0 | sigprs$c2.ulccol > 0), na.rm = T),
              sum((sigprs$c1.crohns > 0 | sigprs$c2.crohns > 0), na.rm = T),
              sum((sigprs$c1.MS > 0 | sigprs$c2.MS > 0), na.rm = T),
              sum((sigprs$c1.RA > 0 | sigprs$c2.RA > 0), na.rm = T),
              sum((sigprs$c1.SLE > 0 | sigprs$c2.SLE > 0), na.rm = T),
              sum((sigprs$c1.T1D > 0 | sigprs$c2.T1D > 0), na.rm = T),
              sum((sigprs$c1.AD > 0 | sigprs$c2.AD > 0), na.rm = T),
              sum((sigprs$c1.PD > 0 | sigprs$c2.PD > 0), na.rm = T),
              sum((sigprs$c1.ALS > 0 | sigprs$c2.ALS > 0), na.rm = T),
              sum((sigprs$c1.epilepsy > 0 | sigprs$c2.epilepsy > 0), na.rm = T),
              sum((sigprs$c1.stroke > 0 | sigprs$c2.stroke > 0), na.rm = T)),
  Pos.CD8 = c(sum((sigprs$c3.lymph > 0 | sigprs$c4.lymph > 0), na.rm = T),
              sum((sigprs$c3.wbc > 0 | sigprs$c4.wbc > 0), na.rm = T),
              sum((sigprs$c3.crp > 0 | sigprs$c4.crp > 0), na.rm = T),
              sum((sigprs$c3.ulccol > 0 | sigprs$c4.ulccol > 0), na.rm = T),
              sum((sigprs$c3.crohns > 0 | sigprs$c4.crohns > 0), na.rm = T),
              sum((sigprs$c3.MS > 0 | sigprs$c4.MS > 0), na.rm = T),
              sum((sigprs$c3.RA > 0 | sigprs$c4.RA > 0), na.rm = T),
              sum((sigprs$c3.SLE > 0 | sigprs$c4.SLE > 0), na.rm = T),
              sum((sigprs$c3.T1D > 0 | sigprs$c4.T1D > 0), na.rm = T),
              sum((sigprs$c3.AD > 0 | sigprs$c4.AD > 0), na.rm = T),
              sum((sigprs$c3.PD > 0 | sigprs$c4.PD > 0), na.rm = T),
              sum((sigprs$c3.ALS > 0 | sigprs$c4.ALS > 0), na.rm = T),
              sum((sigprs$c3.epilepsy > 0 | sigprs$c4.epilepsy > 0), na.rm = T),
              sum((sigprs$c3.stroke > 0 | sigprs$c4.stroke > 0), na.rm = T)),
  Neg.CD4 = c(sum((sigprs$c1.lymph < 0 | sigprs$c2.lymph < 0), na.rm = T),
              sum((sigprs$c1.wbc < 0 | sigprs$c2.wbc < 0), na.rm = T),
              sum((sigprs$c1.crp < 0 | sigprs$c2.crp < 0), na.rm = T),
              sum((sigprs$c1.ulccol < 0 | sigprs$c2.ulccol < 0), na.rm = T),
              sum((sigprs$c1.crohns < 0 | sigprs$c2.crohns < 0), na.rm = T),
              sum((sigprs$c1.MS < 0 | sigprs$c2.MS < 0), na.rm = T),
              sum((sigprs$c1.RA < 0 | sigprs$c2.RA < 0), na.rm = T),
              sum((sigprs$c1.SLE < 0 | sigprs$c2.SLE < 0), na.rm = T),
              sum((sigprs$c1.T1D < 0 | sigprs$c2.T1D < 0), na.rm = T),
              sum((sigprs$c1.AD < 0 | sigprs$c2.AD < 0), na.rm = T),
              sum((sigprs$c1.PD < 0 | sigprs$c2.PD < 0), na.rm = T),
              sum((sigprs$c1.ALS < 0 | sigprs$c2.ALS < 0), na.rm = T),
              sum((sigprs$c1.epilepsy < 0 | sigprs$c2.epilepsy < 0), na.rm = T),
              sum((sigprs$c1.stroke < 0 | sigprs$c2.stroke < 0), na.rm = T)),
  Neg.CD8 = c(sum((sigprs$c3.lymph < 0 | sigprs$c4.lymph < 0), na.rm = T),
              sum((sigprs$c3.wbc < 0 | sigprs$c4.wbc < 0), na.rm = T),
              sum((sigprs$c3.crp < 0 | sigprs$c4.crp < 0), na.rm = T),
              sum((sigprs$c3.ulccol < 0 | sigprs$c4.ulccol < 0), na.rm = T),
              sum((sigprs$c3.crohns < 0 | sigprs$c4.crohns < 0), na.rm = T),
              sum((sigprs$c3.MS < 0 | sigprs$c4.MS < 0), na.rm = T),
              sum((sigprs$c3.RA < 0 | sigprs$c4.RA < 0), na.rm = T),
              sum((sigprs$c3.SLE < 0 | sigprs$c4.SLE < 0), na.rm = T),
              sum((sigprs$c3.T1D < 0 | sigprs$c4.T1D < 0), na.rm = T),
              sum((sigprs$c3.AD < 0 | sigprs$c4.AD < 0), na.rm = T),
              sum((sigprs$c3.PD < 0 | sigprs$c4.PD < 0), na.rm = T),
              sum((sigprs$c3.ALS < 0 | sigprs$c4.ALS < 0), na.rm = T),
              sum((sigprs$c3.epilepsy < 0 | sigprs$c4.epilepsy < 0), na.rm = T),
              sum((sigprs$c3.stroke < 0 | sigprs$c4.stroke < 0), na.rm = T))
  )

#Add columns for percentages of total, positive t-value, and negative t-value genes in CD4 vs CD8
cd4vscd8$Percent.Total.CD4 <- (cd4vscd8$Total.CD4/(cd4vscd8$Total.CD4+cd4vscd8$Total.CD8))*100
cd4vscd8$Percent.Total.CD8 <- (cd4vscd8$Total.CD8/(cd4vscd8$Total.CD4+cd4vscd8$Total.CD8))*100
cd4vscd8$Percent.Pos.CD4 <- (cd4vscd8$Pos.CD4/(cd4vscd8$Pos.CD4+cd4vscd8$Pos.CD8))*100
cd4vscd8$Percent.Pos.CD8 <- (cd4vscd8$Pos.CD8/(cd4vscd8$Pos.CD4+cd4vscd8$Pos.CD8))*100
cd4vscd8$Percent.Neg.CD4 <- (cd4vscd8$Neg.CD4/(cd4vscd8$Neg.CD4+cd4vscd8$Neg.CD8))*100
cd4vscd8$Percent.Neg.CD8 <- (cd4vscd8$Neg.CD8/(cd4vscd8$Neg.CD4+cd4vscd8$Neg.CD8))*100
cd4vscd8$Total.CD4.skew <- cd4vscd8$Percent.Total.CD4 - cd4vscd8$Percent.Total.CD8
cd4vscd8$Pos.CD4.skew <- cd4vscd8$Percent.Pos.CD4 - cd4vscd8$Percent.Pos.CD8
cd4vscd8$Neg.CD4.skew <- cd4vscd8$Percent.Neg.CD4 - cd4vscd8$Percent.Neg.CD8

#Form matrix for heatmap
cd4vscd8matrix <- as.matrix(cd4vscd8[,14:16])
rownames(cd4vscd8matrix) <- cd4vscd8$Trait
colnames(cd4vscd8matrix) <- c("All genes",
                              "+ t-value",
                              "- t-value")

#Plot heatmap
cd4.cd8.heatmap <- Heatmap(cd4vscd8matrix, row_names_side = "left", column_names_side = "top",
                           cluster_rows = F, cluster_columns = F, column_names_rot = 45,
                           border = T, heatmap_legend_param = list(title = "CD4 vs CD8",
                            at = c(-100, 0, 100), labels = c("100% CD8", "50/50", "100% CD4")),
                           column_title = "Abundance of PRS genes by cell type")

#Save heatmap
tiff(filename = "CD4.vs.CD8 heatmap.postQC.tiff", width = 5, height = 6, units = "in", res = 300)
plot(cd4.cd8.heatmap)
dev.off()