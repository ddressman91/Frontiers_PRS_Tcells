library(data.table)
library(dplyr)

#Read in file with all t values for PRS-associated genes, regardless of significance
allprs <- fread("all.prs.genes.filtered.postQC.txt", fill = T)

#Run correlation tests between all combinations of traits for a given
## T cell subtype, where c1 is CD4+CD45RO-, c2 is CD4+CD45RO+, c3 is
### CD8+CD45RO-, and c4 is CD8+CD45RO+, extract p and r values
cor.test(allprs$c3.lymph, allprs$c3.wbc, method = "pearson")[c("p.value", "estimate")]
cor.test(allprs$c3.lymph, allprs$c3.crp, method = "pearson")[c("p.value", "estimate")]
cor.test(allprs$c3.lymph, allprs$c3.ulccol, method = "pearson")[c("p.value", "estimate")]
cor.test(allprs$c3.lymph, allprs$c3.crohns, method = "pearson")[c("p.value", "estimate")]
cor.test(allprs$c3.lymph, allprs$c3.MS, method = "pearson")[c("p.value", "estimate")]
cor.test(allprs$c3.lymph, allprs$c3.RA, method = "pearson")[c("p.value", "estimate")]
cor.test(allprs$c3.lymph, allprs$c3.SLE, method = "pearson")[c("p.value", "estimate")]
cor.test(allprs$c3.lymph, allprs$c3.T1D, method = "pearson")[c("p.value", "estimate")]
cor.test(allprs$c3.lymph, allprs$c3.AD, method = "pearson")[c("p.value", "estimate")]
cor.test(allprs$c3.lymph, allprs$c3.PD, method = "pearson")[c("p.value", "estimate")]
cor.test(allprs$c3.lymph, allprs$c3.ALS, method = "pearson")[c("p.value", "estimate")]
cor.test(allprs$c3.lymph, allprs$c3.epilepsy, method = "pearson")[c("p.value", "estimate")]
cor.test(allprs$c3.lymph, allprs$c3.stroke, method = "pearson")[c("p.value", "estimate")]
cor.test(allprs$c3.wbc, allprs$c3.crp, method = "pearson")[c("p.value", "estimate")]
cor.test(allprs$c3.wbc, allprs$c3.ulccol, method = "pearson")[c("p.value", "estimate")]
cor.test(allprs$c3.wbc, allprs$c3.crohns, method = "pearson")[c("p.value", "estimate")]
cor.test(allprs$c3.wbc, allprs$c3.MS, method = "pearson")[c("p.value", "estimate")]
cor.test(allprs$c3.wbc, allprs$c3.RA, method = "pearson")[c("p.value", "estimate")]
cor.test(allprs$c3.wbc, allprs$c3.SLE, method = "pearson")[c("p.value", "estimate")]
cor.test(allprs$c3.wbc, allprs$c3.T1D, method = "pearson")[c("p.value", "estimate")]
cor.test(allprs$c3.wbc, allprs$c3.AD, method = "pearson")[c("p.value", "estimate")]
cor.test(allprs$c3.wbc, allprs$c3.PD, method = "pearson")[c("p.value", "estimate")]
cor.test(allprs$c3.wbc, allprs$c3.ALS, method = "pearson")[c("p.value", "estimate")]
cor.test(allprs$c3.wbc, allprs$c3.epilepsy, method = "pearson")[c("p.value", "estimate")]
cor.test(allprs$c3.wbc, allprs$c3.stroke, method = "pearson")[c("p.value", "estimate")]
cor.test(allprs$c3.crp, allprs$c3.ulccol, method = "pearson")[c("p.value", "estimate")]
cor.test(allprs$c3.crp, allprs$c3.crohns, method = "pearson")[c("p.value", "estimate")]
cor.test(allprs$c3.crp, allprs$c3.MS, method = "pearson")[c("p.value", "estimate")]
cor.test(allprs$c3.crp, allprs$c3.RA, method = "pearson")[c("p.value", "estimate")]
cor.test(allprs$c3.crp, allprs$c3.SLE, method = "pearson")[c("p.value", "estimate")]
cor.test(allprs$c3.crp, allprs$c3.T1D, method = "pearson")[c("p.value", "estimate")]
cor.test(allprs$c3.crp, allprs$c3.AD, method = "pearson")[c("p.value", "estimate")]
cor.test(allprs$c3.crp, allprs$c3.PD, method = "pearson")[c("p.value", "estimate")]
cor.test(allprs$c3.crp, allprs$c3.ALS, method = "pearson")[c("p.value", "estimate")]
cor.test(allprs$c3.crp, allprs$c3.epilepsy, method = "pearson")[c("p.value", "estimate")]
cor.test(allprs$c3.crp, allprs$c3.stroke, method = "pearson")[c("p.value", "estimate")]
cor.test(allprs$c3.ulccol, allprs$c3.crohns, method = "pearson")[c("p.value", "estimate")]
cor.test(allprs$c3.ulccol, allprs$c3.MS, method = "pearson")[c("p.value", "estimate")]
cor.test(allprs$c3.ulccol, allprs$c3.RA, method = "pearson")[c("p.value", "estimate")]
cor.test(allprs$c3.ulccol, allprs$c3.SLE, method = "pearson")[c("p.value", "estimate")]
cor.test(allprs$c3.ulccol, allprs$c3.T1D, method = "pearson")[c("p.value", "estimate")]
cor.test(allprs$c3.ulccol, allprs$c3.AD, method = "pearson")[c("p.value", "estimate")]
cor.test(allprs$c3.ulccol, allprs$c3.PD, method = "pearson")[c("p.value", "estimate")]
cor.test(allprs$c3.ulccol, allprs$c3.ALS, method = "pearson")[c("p.value", "estimate")]
cor.test(allprs$c3.ulccol, allprs$c3.epilepsy, method = "pearson")[c("p.value", "estimate")]
cor.test(allprs$c3.ulccol, allprs$c3.stroke, method = "pearson")[c("p.value", "estimate")]
cor.test(allprs$c3.crohns, allprs$c3.MS, method = "pearson")[c("p.value", "estimate")]
cor.test(allprs$c3.crohns, allprs$c3.RA, method = "pearson")[c("p.value", "estimate")]
cor.test(allprs$c3.crohns, allprs$c3.SLE, method = "pearson")[c("p.value", "estimate")]
cor.test(allprs$c3.crohns, allprs$c3.T1D, method = "pearson")[c("p.value", "estimate")]
cor.test(allprs$c3.crohns, allprs$c3.AD, method = "pearson")[c("p.value", "estimate")]
cor.test(allprs$c3.crohns, allprs$c3.PD, method = "pearson")[c("p.value", "estimate")]
cor.test(allprs$c3.crohns, allprs$c3.ALS, method = "pearson")[c("p.value", "estimate")]
cor.test(allprs$c3.crohns, allprs$c3.epilepsy, method = "pearson")[c("p.value", "estimate")]
cor.test(allprs$c3.crohns, allprs$c3.stroke, method = "pearson")[c("p.value", "estimate")]
cor.test(allprs$c3.MS, allprs$c3.RA, method = "pearson")[c("p.value", "estimate")]
cor.test(allprs$c3.MS, allprs$c3.SLE, method = "pearson")[c("p.value", "estimate")]
cor.test(allprs$c3.MS, allprs$c3.T1D, method = "pearson")[c("p.value", "estimate")]
cor.test(allprs$c3.MS, allprs$c3.AD, method = "pearson")[c("p.value", "estimate")]
cor.test(allprs$c3.MS, allprs$c3.PD, method = "pearson")[c("p.value", "estimate")]
cor.test(allprs$c3.MS, allprs$c3.ALS, method = "pearson")[c("p.value", "estimate")]
cor.test(allprs$c3.MS, allprs$c3.epilepsy, method = "pearson")[c("p.value", "estimate")]
cor.test(allprs$c3.MS, allprs$c3.stroke, method = "pearson")[c("p.value", "estimate")]
cor.test(allprs$c3.RA, allprs$c3.SLE, method = "pearson")[c("p.value", "estimate")]
cor.test(allprs$c3.RA, allprs$c3.T1D, method = "pearson")[c("p.value", "estimate")]
cor.test(allprs$c3.RA, allprs$c3.AD, method = "pearson")[c("p.value", "estimate")]
cor.test(allprs$c3.RA, allprs$c3.PD, method = "pearson")[c("p.value", "estimate")]
cor.test(allprs$c3.RA, allprs$c3.ALS, method = "pearson")[c("p.value", "estimate")]
cor.test(allprs$c3.RA, allprs$c3.epilepsy, method = "pearson")[c("p.value", "estimate")]
cor.test(allprs$c3.RA, allprs$c3.stroke, method = "pearson")[c("p.value", "estimate")]
cor.test(allprs$c3.SLE, allprs$c3.T1D, method = "pearson")[c("p.value", "estimate")]
cor.test(allprs$c3.SLE, allprs$c3.AD, method = "pearson")[c("p.value", "estimate")]
cor.test(allprs$c3.SLE, allprs$c3.PD, method = "pearson")[c("p.value", "estimate")]
cor.test(allprs$c3.SLE, allprs$c3.ALS, method = "pearson")[c("p.value", "estimate")]
cor.test(allprs$c3.SLE, allprs$c3.epilepsy, method = "pearson")[c("p.value", "estimate")]
cor.test(allprs$c3.SLE, allprs$c3.stroke, method = "pearson")[c("p.value", "estimate")]
cor.test(allprs$c3.T1D, allprs$c3.AD, method = "pearson")[c("p.value", "estimate")]
cor.test(allprs$c3.T1D, allprs$c3.PD, method = "pearson")[c("p.value", "estimate")]
cor.test(allprs$c3.T1D, allprs$c3.ALS, method = "pearson")[c("p.value", "estimate")]
cor.test(allprs$c3.T1D, allprs$c3.epilepsy, method = "pearson")[c("p.value", "estimate")]
cor.test(allprs$c3.T1D, allprs$c3.stroke, method = "pearson")[c("p.value", "estimate")]
cor.test(allprs$c3.AD, allprs$c3.PD, method = "pearson")[c("p.value", "estimate")]
cor.test(allprs$c3.AD, allprs$c3.ALS, method = "pearson")[c("p.value", "estimate")]
cor.test(allprs$c3.AD, allprs$c3.epilepsy, method = "pearson")[c("p.value", "estimate")]
cor.test(allprs$c3.AD, allprs$c3.stroke, method = "pearson")[c("p.value", "estimate")]
cor.test(allprs$c3.PD, allprs$c3.ALS, method = "pearson")[c("p.value", "estimate")]
cor.test(allprs$c3.PD, allprs$c3.epilepsy, method = "pearson")[c("p.value", "estimate")]
cor.test(allprs$c3.PD, allprs$c3.stroke, method = "pearson")[c("p.value", "estimate")]
cor.test(allprs$c3.ALS, allprs$c3.epilepsy, method = "pearson")[c("p.value", "estimate")]
cor.test(allprs$c3.ALS, allprs$c3.stroke, method = "pearson")[c("p.value", "estimate")]
cor.test(allprs$c3.epilepsy, allprs$c3.stroke, method = "pearson")[c("p.value", "estimate")]

#Count absolute numbers of genes shared between any two traits
nrow(sigprs[(!is.na(sigprs$c3.lymph) & !is.na(sigprs$c3.wbc)),])
nrow(sigprs[(!is.na(sigprs$c3.lymph) & !is.na(sigprs$c3.crp)),])
nrow(sigprs[(!is.na(sigprs$c3.lymph) & !is.na(sigprs$c3.ulccol)),])
nrow(sigprs[(!is.na(sigprs$c3.lymph) & !is.na(sigprs$c3.crohns)),])
nrow(sigprs[(!is.na(sigprs$c3.lymph) & !is.na(sigprs$c3.MS)),])
nrow(sigprs[(!is.na(sigprs$c3.lymph) & !is.na(sigprs$c3.RA)),])
nrow(sigprs[(!is.na(sigprs$c3.lymph) & !is.na(sigprs$c3.SLE)),])
nrow(sigprs[(!is.na(sigprs$c3.lymph) & !is.na(sigprs$c3.T1D)),])
nrow(sigprs[(!is.na(sigprs$c3.lymph) & !is.na(sigprs$c3.AD)),])
nrow(sigprs[(!is.na(sigprs$c3.lymph) & !is.na(sigprs$c3.PD)),])
nrow(sigprs[(!is.na(sigprs$c3.lymph) & !is.na(sigprs$c3.ALS)),])
nrow(sigprs[(!is.na(sigprs$c3.lymph) & !is.na(sigprs$c3.epilepsy)),])
nrow(sigprs[(!is.na(sigprs$c3.lymph) & !is.na(sigprs$c3.stroke)),])
nrow(sigprs[(!is.na(sigprs$c3.wbc) & !is.na(sigprs$c3.crp)),])
nrow(sigprs[(!is.na(sigprs$c3.wbc) & !is.na(sigprs$c3.ulccol)),])
nrow(sigprs[(!is.na(sigprs$c3.wbc) & !is.na(sigprs$c3.crohns)),])
nrow(sigprs[(!is.na(sigprs$c3.wbc) & !is.na(sigprs$c3.MS)),])
nrow(sigprs[(!is.na(sigprs$c3.wbc) & !is.na(sigprs$c3.RA)),])
nrow(sigprs[(!is.na(sigprs$c3.wbc) & !is.na(sigprs$c3.SLE)),])
nrow(sigprs[(!is.na(sigprs$c3.wbc) & !is.na(sigprs$c3.T1D)),])
nrow(sigprs[(!is.na(sigprs$c3.wbc) & !is.na(sigprs$c3.AD)),])
nrow(sigprs[(!is.na(sigprs$c3.wbc) & !is.na(sigprs$c3.PD)),])
nrow(sigprs[(!is.na(sigprs$c3.wbc) & !is.na(sigprs$c3.ALS)),])
nrow(sigprs[(!is.na(sigprs$c3.wbc) & !is.na(sigprs$c3.epilepsy)),])
nrow(sigprs[(!is.na(sigprs$c3.wbc) & !is.na(sigprs$c3.stroke)),])
nrow(sigprs[(!is.na(sigprs$c3.crp) & !is.na(sigprs$c3.ulccol)),])
nrow(sigprs[(!is.na(sigprs$c3.crp) & !is.na(sigprs$c3.crohns)),])
nrow(sigprs[(!is.na(sigprs$c3.crp) & !is.na(sigprs$c3.MS)),])
nrow(sigprs[(!is.na(sigprs$c3.crp) & !is.na(sigprs$c3.RA)),])
nrow(sigprs[(!is.na(sigprs$c3.crp) & !is.na(sigprs$c3.SLE)),])
nrow(sigprs[(!is.na(sigprs$c3.crp) & !is.na(sigprs$c3.T1D)),])
nrow(sigprs[(!is.na(sigprs$c3.crp) & !is.na(sigprs$c3.AD)),])
nrow(sigprs[(!is.na(sigprs$c3.crp) & !is.na(sigprs$c3.PD)),])
nrow(sigprs[(!is.na(sigprs$c3.crp) & !is.na(sigprs$c3.ALS)),])
nrow(sigprs[(!is.na(sigprs$c3.crp) & !is.na(sigprs$c3.epilepsy)),])
nrow(sigprs[(!is.na(sigprs$c3.crp) & !is.na(sigprs$c3.stroke)),])
nrow(sigprs[(!is.na(sigprs$c3.ulccol) & !is.na(sigprs$c3.crohns)),])
nrow(sigprs[(!is.na(sigprs$c3.ulccol) & !is.na(sigprs$c3.MS)),])
nrow(sigprs[(!is.na(sigprs$c3.ulccol) & !is.na(sigprs$c3.RA)),])
nrow(sigprs[(!is.na(sigprs$c3.ulccol) & !is.na(sigprs$c3.SLE)),])
nrow(sigprs[(!is.na(sigprs$c3.ulccol) & !is.na(sigprs$c3.T1D)),])
nrow(sigprs[(!is.na(sigprs$c3.ulccol) & !is.na(sigprs$c3.AD)),])
nrow(sigprs[(!is.na(sigprs$c3.ulccol) & !is.na(sigprs$c3.PD)),])
nrow(sigprs[(!is.na(sigprs$c3.ulccol) & !is.na(sigprs$c3.ALS)),])
nrow(sigprs[(!is.na(sigprs$c3.ulccol) & !is.na(sigprs$c3.epilepsy)),])
nrow(sigprs[(!is.na(sigprs$c3.ulccol) & !is.na(sigprs$c3.stroke)),])
nrow(sigprs[(!is.na(sigprs$c3.crohns) & !is.na(sigprs$c3.MS)),])
nrow(sigprs[(!is.na(sigprs$c3.crohns) & !is.na(sigprs$c3.RA)),])
nrow(sigprs[(!is.na(sigprs$c3.crohns) & !is.na(sigprs$c3.SLE)),])
nrow(sigprs[(!is.na(sigprs$c3.crohns) & !is.na(sigprs$c3.T1D)),])
nrow(sigprs[(!is.na(sigprs$c3.crohns) & !is.na(sigprs$c3.AD)),])
nrow(sigprs[(!is.na(sigprs$c3.crohns) & !is.na(sigprs$c3.PD)),])
nrow(sigprs[(!is.na(sigprs$c3.crohns) & !is.na(sigprs$c3.ALS)),])
nrow(sigprs[(!is.na(sigprs$c3.crohns) & !is.na(sigprs$c3.epilepsy)),])
nrow(sigprs[(!is.na(sigprs$c3.crohns) & !is.na(sigprs$c3.stroke)),])
nrow(sigprs[(!is.na(sigprs$c3.MS) & !is.na(sigprs$c3.RA)),])
nrow(sigprs[(!is.na(sigprs$c3.MS) & !is.na(sigprs$c3.SLE)),])
nrow(sigprs[(!is.na(sigprs$c3.MS) & !is.na(sigprs$c3.T1D)),])
nrow(sigprs[(!is.na(sigprs$c3.MS) & !is.na(sigprs$c3.AD)),])
nrow(sigprs[(!is.na(sigprs$c3.MS) & !is.na(sigprs$c3.PD)),])
nrow(sigprs[(!is.na(sigprs$c3.MS) & !is.na(sigprs$c3.ALS)),])
nrow(sigprs[(!is.na(sigprs$c3.MS) & !is.na(sigprs$c3.epilepsy)),])
nrow(sigprs[(!is.na(sigprs$c3.MS) & !is.na(sigprs$c3.stroke)),])
nrow(sigprs[(!is.na(sigprs$c3.RA) & !is.na(sigprs$c3.SLE)),])
nrow(sigprs[(!is.na(sigprs$c3.RA) & !is.na(sigprs$c3.T1D)),])
nrow(sigprs[(!is.na(sigprs$c3.RA) & !is.na(sigprs$c3.AD)),])
nrow(sigprs[(!is.na(sigprs$c3.RA) & !is.na(sigprs$c3.PD)),])
nrow(sigprs[(!is.na(sigprs$c3.RA) & !is.na(sigprs$c3.ALS)),])
nrow(sigprs[(!is.na(sigprs$c3.RA) & !is.na(sigprs$c3.epilepsy)),])
nrow(sigprs[(!is.na(sigprs$c3.RA) & !is.na(sigprs$c3.stroke)),])
nrow(sigprs[(!is.na(sigprs$c3.SLE) & !is.na(sigprs$c3.T1D)),])
nrow(sigprs[(!is.na(sigprs$c3.SLE) & !is.na(sigprs$c3.AD)),])
nrow(sigprs[(!is.na(sigprs$c3.SLE) & !is.na(sigprs$c3.PD)),])
nrow(sigprs[(!is.na(sigprs$c3.SLE) & !is.na(sigprs$c3.ALS)),])
nrow(sigprs[(!is.na(sigprs$c3.SLE) & !is.na(sigprs$c3.epilepsy)),])
nrow(sigprs[(!is.na(sigprs$c3.SLE) & !is.na(sigprs$c3.stroke)),])
nrow(sigprs[(!is.na(sigprs$c3.T1D) & !is.na(sigprs$c3.AD)),])
nrow(sigprs[(!is.na(sigprs$c3.T1D) & !is.na(sigprs$c3.PD)),])
nrow(sigprs[(!is.na(sigprs$c3.T1D) & !is.na(sigprs$c3.ALS)),])
nrow(sigprs[(!is.na(sigprs$c3.T1D) & !is.na(sigprs$c3.epilepsy)),])
nrow(sigprs[(!is.na(sigprs$c3.T1D) & !is.na(sigprs$c3.stroke)),])
nrow(sigprs[(!is.na(sigprs$c3.AD) & !is.na(sigprs$c3.PD)),])
nrow(sigprs[(!is.na(sigprs$c3.AD) & !is.na(sigprs$c3.ALS)),])
nrow(sigprs[(!is.na(sigprs$c3.AD) & !is.na(sigprs$c3.epilepsy)),])
nrow(sigprs[(!is.na(sigprs$c3.AD) & !is.na(sigprs$c3.stroke)),])
nrow(sigprs[(!is.na(sigprs$c3.PD) & !is.na(sigprs$c3.ALS)),])
nrow(sigprs[(!is.na(sigprs$c3.PD) & !is.na(sigprs$c3.epilepsy)),])
nrow(sigprs[(!is.na(sigprs$c3.PD) & !is.na(sigprs$c3.stroke)),])
nrow(sigprs[(!is.na(sigprs$c3.ALS) & !is.na(sigprs$c3.epilepsy)),])
nrow(sigprs[(!is.na(sigprs$c3.ALS) & !is.na(sigprs$c3.stroke)),])
nrow(sigprs[(!is.na(sigprs$c3.epilepsy) & !is.na(sigprs$c3.stroke)),])