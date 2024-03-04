library(dplyr)
library(ComplexHeatmap)
library(ggplot2)
library(circlize)

#Save CSV file, and reload later if needed, add a column to sigprs with # of non-NA values per row

sigprs <-read.csv(file = "sig.prs.no.psych.csv")
sigprs <- sigprs %>% mutate(Num.Traits = rowSums(!is.na(select(., -Gene))))

#Plot results with ComplexHeatmap
ann_colors = list(Cell.Type = c("CD4+CD45RO-" = "#e41a1c", "CD4+CD45RO+" = "#4daf4a",
                                "CD8+CD45RO-" = "#377eb8", "CD8+CD45RO+" = "#984ea3"), 
                  Trait.Type = c("Immunity and Inflammation" = "#66c2a5",
                                 "Autoimmune Disease" = "#fc8d62",
                                 "Neurological Disease" = "#8da0cb"),
                  Trait = c("Lymphocyte counts" = "dodgerblue2",
                            "White blood cell counts"  = "#FB9A99",
                            "C-reactive protein" = "green4",
                            "Ulcerative colitis" = "#6A3D9A",
                            "Crohn's disease" = "#FF7F00",
                            "Multiple sclerosis" = "black",
                            "Rheumatoid arthritis" = "gold1",
                            "Systemic lupus erythematosus" = "skyblue2",
                            "Type 1 diabetes" = "#E31A1C",
                            "Alzheimer's disease" = "palegreen2",
                            "Parkinson's disease" = "blue1",
                            "Amyotrophic lateral sclerosis" = "#FDBF6F",
                            "Epilepsy" = "deeppink1", "Stroke" = "gray70"))

ann_labels = list(
  Cell.Type = list(
    title = "Cell type",
    at = c("CD4+CD45RO-", "CD4+CD45RO+", "CD8+CD45RO-", "CD8+CD45RO+"),
    labels = c("CD4+CD45RO-", "CD4+CD45RO+", "CD8+CD45RO-", "CD8+CD45RO+")),
  Trait.Type = list(
    title = "Trait type",
    at = c("Immunity and Inflammation", "Autoimmune Disease",
           "Neurological Disease"),
    labels = c("Immunity and Inflammation", "Autoimmune Disease",
               "Neurological Disease")),
  Trait = list(
    title = "Trait",
    at = c("Lymphocyte counts", "White blood cell counts", "C-reactive protein",
           "Ulcerative colitis", "Crohn's disease", "Multiple sclerosis",
           "Rheumatoid arthritis", "Systemic lupus erythematosus",
           "Type 1 diabetes", "Alzheimer's disease", "Parkinson's disease",
           "Amyotrophic lateral sclerosis", "Epilepsy", "Stroke"),
    labels = c("Lymphocyte counts", "White blood cell counts", "C-reactive protein",
               "Ulcerative colitis", "Crohn's disease", "Multiple sclerosis",
               "Rheumatoid arthritis", "Systemic lupus erythematosus",
               "Type 1 diabetes", "Alzheimer's disease", "Parkinson's disease",
               "Amyotrophic lateral sclerosis", "Epilepsy", "Stroke")))

col_annotations = read.csv("C:/Users/dalli/Documents/PRS paper/Updated column annotations for heatmap.csv")

ha<-HeatmapAnnotation(df = col_annotations, col = ann_colors, which = "column",
                      annotation_legend_param = ann_labels)

prs.heatmap <- Heatmap(as.matrix(signew[,2:57]), cluster_rows = F, cluster_columns = F,
                       col = colorRamp2(c(-5, 0, 5), c("blue", "lightyellow", "red")),
                       border = T, show_row_names = F, show_column_names = F, top_annotation = ha,
                       heatmap_legend_param = list(title = "t value"))
  
prs.heatmap <- draw(prs.heatmap, merge_legend = TRUE)

#Save the heatmap
tiff(filename="prs.heatmap.all.traits.postQC.tiff",width=180,height=160,units="mm",res=300)
plot(prs.heatmap)
dev.off()