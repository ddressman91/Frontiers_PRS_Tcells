#Generating enrichR results from GO Biological Processes and plotting top pathways

library(dplyr)
library(ggplot2)
library(scales)
library(cowplot)
library(ggpubr)
library(readxl)

#Load supplementary table 4 with significant GSEA pathways
all.sig <- read_xlsx("/Updated supplementary table 4.xlsx")

#Plot top pathways
p1 <- ggplot(all.sig[all.sig$Cell.Type == "CD4+CD45RO-", ], aes(x = FDR.qval,
          y = NAME, color = Trait, shape = t.value.sign)) + geom_point() +
  expand_limits(x=0) + labs(x="FDR q-value", y=NULL, title="Top pathways in CD4+CD45RO-") +
  scale_y_discrete(labels = wrap_format(50))+
  theme(axis.text = element_text(size = 9),
        panel.background = element_rect(fill="white", colour="white"),
        panel.grid.major = element_line(linewidth = 0.25, linetype = 'solid',colour = "lightgray"),
        panel.grid.minor = element_line(linewidth = 0.25, linetype = 'solid',colour = "lightgray"),
        panel.border = element_rect(colour = "black",fill=NA),
        plot.title = element_text(hjust = 1),
        legend.position = "none")+
  scale_color_manual(values = c("Lymphocyte counts" = "dodgerblue2", "WBC counts"  = "#FB9A99",
                                "C-reactive protein" = "green4", "Ulcerative colitis" = "#6A3D9A",
                                "Crohn's disease" = "#FF7F00", "Multiple sclerosis" = "black",
                                "Rheumatoid arthritis" = "gold1", "Lupus" = "skyblue2",
                                "Type 1 diabetes" = "#E31A1C", "Alzheimer's" = "palegreen2",
                                "Parkinson's" = "blue1", "ALS" = "#FDBF6F",
                                "Epilepsy" = "deeppink1", "Stroke" = "gray70")) +
  scale_shape_manual(name = "t value sign", values = c("-" = 16, "+" = 17))

p2 <- ggplot(all.sig[all.sig$Cell.Type == "CD4+CD45RO+",], aes(x = FDR.qval,
          y = NAME, color = Trait, shape = t.value.sign)) + geom_point() +
  expand_limits(x=0) + labs(x="FDR q-value", y=NULL, title="Top pathways in CD4+CD45RO+") +
  scale_y_discrete(labels = wrap_format(40))+
  theme(axis.text = element_text(size = 9),
        panel.background = element_rect(fill="white", colour="white"),
        panel.grid.major = element_line(linewidth = 0.25, linetype = 'solid',colour = "lightgray"),
        panel.grid.minor = element_line(linewidth = 0.25, linetype = 'solid',colour = "lightgray"),
        panel.border = element_rect(colour = "black",fill=NA),
        plot.title = element_text(hjust = 1),
        legend.position = "none")+
  scale_color_manual(values = c("Lymphocyte counts" = "dodgerblue2", "WBC counts"  = "#FB9A99",
                                "C-reactive protein" = "green4", "Ulcerative colitis" = "#6A3D9A",
                                "Crohn's disease" = "#FF7F00", "Multiple sclerosis" = "black",
                                "Rheumatoid arthritis" = "gold1", "Lupus" = "skyblue2",
                                "Type 1 diabetes" = "#E31A1C", "Alzheimer's" = "palegreen2",
                                "Parkinson's" = "blue1", "ALS" = "#FDBF6F",
                                "Epilepsy" = "deeppink1", "Stroke" = "gray70")) +
  scale_shape_manual(name = "t value sign", values = c("-" = 16, "+" = 17))

p3 <- ggplot(all.sig[all.sig$Cell.Type %in% c("CD8+CD45RO-", "CD8+CD45RO+"),],
             aes(x = FDR.qval, y = NAME, color = Trait, shape = t.value.sign)) + geom_point() +
  expand_limits(x=0) + labs(x="FDR q-value", y=NULL, title="Top pathways in CD8+") +
  scale_y_discrete(labels = wrap_format(40))+
  theme(axis.text = element_text(size = 9),
        panel.background = element_rect(fill="white", colour="white"),
        panel.grid.major = element_line(linewidth = 0.25, linetype = 'solid',colour = "lightgray"),
        panel.grid.minor = element_line(linewidth = 0.25, linetype = 'solid',colour = "lightgray"),
        panel.border = element_rect(colour = "black",fill=NA),
        plot.title = element_text(hjust = 1),
        legend.position = "none")+
  scale_color_manual(values = c("Lymphocyte counts" = "dodgerblue2", "WBC counts"  = "#FB9A99",
                                "C-reactive protein" = "green4", "Ulcerative colitis" = "#6A3D9A",
                                "Crohn's disease" = "#FF7F00", "Multiple sclerosis" = "black",
                                "Rheumatoid arthritis" = "gold1", "Lupus" = "skyblue2",
                                "Type 1 diabetes" = "#E31A1C", "Alzheimer's" = "palegreen2",
                                "Parkinson's" = "blue1", "ALS" = "#FDBF6F",
                                "Epilepsy" = "deeppink1", "Stroke" = "gray70")) +
  scale_shape_manual(name = "t value sign", values = c("-" = 16, "+" = 17))

p4 <- ggplot(all.sig,aes(FDR.qval,y=NAME,color=Trait, shape = t.value.sign)) + geom_point() +
  expand_limits(x=0) + labs(x="FDR q-value", y=NULL, title="Top pathways in CD4+CD45RO-") +
  scale_y_discrete(labels = wrap_format(40))+
  theme(axis.text = element_text(size = 9),
        panel.background = element_rect(fill="white", colour="white"),
        panel.grid.major = element_line(linewidth = 0.25, linetype = 'solid',colour = "lightgray"),
        panel.grid.minor = element_line(linewidth = 0.25, linetype = 'solid',colour = "lightgray"),
        panel.border = element_rect(colour = "black",fill=NA),
        plot.title = element_text(hjust = 1)) +
  scale_color_manual(values = c("Lymphocyte counts" = "dodgerblue2", "WBC counts"  = "#FB9A99",
                                "C-reactive protein" = "green4", "Ulcerative colitis" = "#6A3D9A",
                                "Crohn's disease" = "#FF7F00", "Multiple sclerosis" = "black",
                                "Rheumatoid arthritis" = "gold1", "Lupus" = "skyblue2",
                                "Type 1 diabetes" = "#E31A1C", "Alzheimer's" = "palegreen2",
                                "Parkinson's" = "blue1", "ALS" = "#FDBF6F",
                                "Epilepsy" = "deeppink1", "Stroke" = "gray70")) +
  scale_shape_manual(name = "t value sign", values = c("-" = 16, "+" = 17))

leg <- ggpubr::get_legend(p4)
leg <- as_ggplot(leg)

col1 <- plot_grid(p1, p3, labels = c('A', 'C'), label_size = 12, ncol = 1)
composite <- plot_grid(col1, p2, leg, labels = c('', 'B', ''),
                       label_size = 12, ncol = 3, rel_widths = c(2,2,0.8))

#Save plot
ggsave(filename="Figure 3.tiff",
       width = 180, height = 75, plot = composite, device = 'tiff', units = "mm")