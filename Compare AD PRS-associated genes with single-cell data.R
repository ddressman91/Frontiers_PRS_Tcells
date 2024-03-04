library(Seurat)
library(dplyr)
library(ggplot2)
library(data.table)
library(cowplot)

#Read in downloaded datasets from Xu & Jia 2021, create Seurat objects
AD1 <- Read10X(data.dir = "/AD1/")
AD2 <- Read10X(data.dir = "/AD2/")
AD3 <- Read10X(data.dir = "/AD3/")
HC1 <- Read10X(data.dir = "/HC1/")
HC2 <- Read10X(data.dir = "/HC2/")

AD1 <- CreateSeuratObject(counts = AD1, project = "AD1")
AD2 <- CreateSeuratObject(counts = AD2, project = "AD2")
AD3 <- CreateSeuratObject(counts = AD3, project = "AD3")
HC1 <- CreateSeuratObject(counts = HC1, project = "HC1")
HC2 <- CreateSeuratObject(counts = HC2, project = "HC2")

#Run QC on individual Seurat objects and prepare for integration
AD1 <- RenameCells(AD1, add.cell.id = "AD1")
AD2 <- RenameCells(AD2, add.cell.id = "AD2")
AD3 <- RenameCells(AD3, add.cell.id = "AD3")
HC1 <- RenameCells(HC1, add.cell.id = "HC1")
HC2 <- RenameCells(HC2, add.cell.id = "HC2")

#Prepare Seurat objects for integration
object.list <- c(AD1, AD2, AD3, HC1, HC2)
for (i in 1:length(object.list)) {
  object.list[[i]] <- NormalizeData(object.list[[i]], verbose = F)
  object.list[[i]] <- FindVariableFeatures(object.list[[i]], verbose = F)
  object.list[[i]] <- ScaleData(object.list[[i]], verbose = F)
  object.list[[i]] <- RunPCA(object.list[[i]], verbose = F)
  object.list[[i]] <- FindNeighbors(object.list[[i]], dims = 1:30,
                                    reduction = "pca", verbose = F)
  object.list[[i]] <- FindClusters(object.list[[i]], resolution = 2,
                                   cluster.name = "unintegrated_clusters")
  object.list[[i]] <- RunUMAP(object.list[[i]], dims = 1:30, reduction = "pca",
                  reduction.name = "umap.unintegrated")
}

features <- SelectIntegrationFeatures(object.list = object.list, nfeatures = 5000)

#Integrate Seurat objects from individual patients
anchors <- FindIntegrationAnchors(object.list = object.list,
                                  anchor.features = features)
adpbmc <- IntegrateData(anchorset = anchors)
setwd("C:/Users/dalli/Documents/PRS paper/Comparison with single-cell data/")
save(adpbmc, file = "Xu.Jia.2021.AD.PBMCs.integrated.RData")

#Filter Seurat objects and run QC
DefaultAssay(adpbmc) <- "RNA"
adpbmc[["percent.mt"]] <- PercentageFeatureSet(adpbmc, pattern = "^MT-")
adpbmc <- subset(adpbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & 
                   percent.mt < 15 & nCount_RNA < 10000)
adpbmc <- SCTransform(adpbmc)
adpbmc <- RunPCA(adpbmc, verbose = F)
adpbmc <- RunUMAP(adpbmc, dims = 1:30, verbose = F)
adpbmc <- FindNeighbors(adpbmc, reduction = "pca", dims = 1:30, verbose = F)
adpbmc <- FindClusters(adpbmc, verbose = F)

adpbmc@meta.data$Condition[adpbmc@meta.data$orig.ident %in% c("AD1","AD2","AD3")] <- "AD"
adpbmc@meta.data$Condition[adpbmc@meta.data$orig.ident %in% c("HC1","HC2")] <- "HC"
save(adpbmc, file = "Xu.Jia.2021.AD.PBMCs.clustered.res0.8.RData")

#Subset T cells, recluster, and annotate T cell subtypes
adpbmc@meta.data$is_tcell<-"No"
adpbmc@meta.data$is_tcell[adpbmc@meta.data$seurat_clusters %in%
                              c("0","1","2","3","4","5","6","12","13","14")]<-"Yes"
Idents(adpbmc) <- "is_tcell"
adtcells <- subset(adpbmc, subset = is_tcell == "Yes")
adtcells <- SCTransform(adtcells)
adtcells <- RunPCA(adtcells, verbose = F)
adtcells <- RunUMAP(adtcells, dims = 1:30, verbose = F)
adtcells <- FindNeighbors(adtcells, reduction = "pca", dims = 1:30, verbose = F)
adtcells <- FindClusters(adtcells, resolution = 0.8, verbose = F)

save(adtcells, file = "Xu.Jia.2021.AD.Tcells.clustered.res0.8.RData")

#Fine cell type annotation with marker expression
adtcells <- PrepSCTFindMarkers(adtcells)
Tcellclustermarkers <- FindAllMarkers(adtcells, assay = "SCT")

DotPlot(adtcells, features = c("CD3E", "CD4", "IL7R", "CD8A", "CD8B", "LEF1",
                               "SELL", "TCF7", "NOSIP", "CCR7", "LTB", "TRAC",
                               "CCL5", "GZMK", "GZMH", "GZMA", "GZMB", "KLRD1",
                               "NKG7", "CST7", "CXCR3", "CXCR4", "CX3CR1",
                               "FOXP3","TIGIT","IL2RA","CTLA4","TRDC","TRGC1",
                               "NCR3","FCER1G","CD79A","LYZ","CD14","CD74",
                               "CD27","CD28","AQP3")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

adtcells@meta.data$CellType[adtcells@meta.data$seurat_clusters %in% c(0,1,4,8,9)] <- "CD4 naive"
adtcells@meta.data$CellType[adtcells@meta.data$seurat_clusters %in% c(3)] <- "CD4 memory"
adtcells@meta.data$CellType[adtcells@meta.data$seurat_clusters %in% c(11)] <- "CD8 naive"
adtcells@meta.data$CellType[adtcells@meta.data$seurat_clusters %in% c(5,6,7)] <- "CD8 memory"
adtcells@meta.data$CellType[adtcells@meta.data$seurat_clusters %in% c(2,13)] <- "CD4/CD8 TEMRA"
adtcells@meta.data$CellType[adtcells@meta.data$seurat_clusters %in% c(10,12)] <- "CD4 low quality"
adtcells@meta.data$CellType[adtcells@meta.data$seurat_clusters %in% c(14)] <- "Treg"

# Join layers of RNA assay, designate cell identities for expression testing,
## DGE comparisons between control and AD for CD4 naive, CD4 memory, CD8 naive,
### and CD8 memory (excluding TEMRA from memory populations)
adtcells[['RNA']] <- JoinLayers(adtcells[['RNA']])
naivecd4ad <- rownames(adtcells@meta.data %>% filter(Condition == "AD",
                                                     seurat_clusters %in% c("0","1","4","8","9")))
naivecd4hc <- rownames(adtcells@meta.data %>% filter(Condition == "HC",
                                                     seurat_clusters %in% c("0","1","4","8","9")))
memorycd4ad <- rownames(adtcells@meta.data %>% filter(Condition == "AD",
                                                    seurat_clusters == "3"))
memorycd4hc <- rownames(adtcells@meta.data %>% filter(Condition == "HC",
                                                    seurat_clusters == "3"))
naivecd8ad <- rownames(adtcells@meta.data %>% filter(Condition == "AD",
                                                   seurat_clusters == "11"))
naivecd8hc <- rownames(adtcells@meta.data %>% filter(Condition == "HC",
                                                    seurat_clusters == "11"))
memorycd8ad <- rownames(adtcells@meta.data %>% filter(Condition == "AD",
                                                    seurat_clusters %in% c("5", "6", "7")))
memorycd8hc <- rownames(adtcells@meta.data %>% filter(Condition == "HC",
                                                    seurat_clusters %in% c("5", "6", "7")))
naivecd4ADdegs <- FindMarkers(adtcells, assay = "RNA", ident.1 = naivecd4ad, ident.2 = naivecd4hc)
memorycd4ADdegs <- FindMarkers(adtcells, assay = "RNA", ident.1 = memorycd4ad, ident.2 = memorycd4hc)
naivecd8ADdegs <- FindMarkers(adtcells, assay = "RNA", ident.1 = naivecd8ad, ident.2 = naivecd8hc)
memorycd8ADdegs <- FindMarkers(adtcells, assay = "RNA", ident.1 = memorycd8ad, ident.2 = memorycd8hc)

write.csv(naivecd4ADdegs, file = "Xu.Jia.naive.CD4.AD.vs.HC.DEGs.fromRNAcounts.csv", row.names = T)
write.csv(memorycd4ADdegs, file = "Xu.Jia.memory.CD4.AD.vs.HC.DEGs.fromRNAcounts.csv", row.names = T)
write.csv(naivecd8ADdegs, file = "Xu.Jia.naive.CD8.AD.vs.HC.DEGs.fromRNAcounts.csv", row.names = T)
write.csv(memorycd8ADdegs, file = "Xu.Jia.memory.CD8.AD.vs.HC.DEGs.fromRNAcounts.csv", row.names = T)

#Plot genes highlighted in text from comparison dataset
DefaultAssay(adtcells) <- "RNA"
Idents(adtcells) <- "CellType"
umap <- DimPlot(adtcells) + ggtitle("T cell subtypes in comparison dataset")
mapk1 <- VlnPlot(adtcells, features = "MAPK1", idents = "CD8 memory", group.by = "Condition") +
  theme(legend.position = "none")
c1qbp <- VlnPlot(adtcells, features = "C1QBP", idents = "CD4 memory", group.by = "Condition") +
  theme(legend.position = "none")
dbnl <- VlnPlot(adtcells, features = "DBNL", idents = "CD4 naive", group.by = "Condition") +
  theme(legend.position = "none")
itsn2 <- VlnPlot(adtcells, features = "ITSN2", idents = "CD8 memory", group.by = "Condition") +
  theme(legend.position = "none")
nfatc2 <- VlnPlot(adtcells, features = "NFATC2", idents = "CD4 memory", group.by = "Condition") +
  theme(legend.position = "none")
geneplot <- plot_grid(mapk1, c1qbp, dbnl, itsn2, nfatc2,
                      nrow = 1, labels = c("B", "C", "D", "E", "F"))
suppfig2 <- plot_grid(umap, geneplot, nrow = 2, labels = c("A", ""),
                      rel_heights = c(2,1))