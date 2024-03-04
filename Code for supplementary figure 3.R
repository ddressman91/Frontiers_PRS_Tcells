library(dplyr)
library(ggplot2)
library(cowplot)

#Read in file with all nominally significant PRS-associated genes
sigprs <- read.csv("sig.prs.no.psych.csv", header = T)

#Create tables for each cell type, add a column to each table for
## the number of traits that each gene is significant in
c1sig <- sigprs[,c(1:15)]
c2sig <- sigprs[,c(1,16:29)]
c3sig <- sigprs[,c(1,30:43)]
c4sig <- sigprs[,c(1,44:57)]

sigprs <- sigprs %>% mutate(Num.Traits = rowSums(!is.na(select(., -Gene))))
c1sig <- c1sig %>% mutate(Num.Traits = rowSums(!is.na(select(., -Gene))))
c2sig <- c2sig %>% mutate(Num.Traits = rowSums(!is.na(select(., -Gene))))
c3sig <- c3sig %>% mutate(Num.Traits = rowSums(!is.na(select(., -Gene))))
c4sig <- c4sig %>% mutate(Num.Traits = rowSums(!is.na(select(., -Gene))))

#Plot histogram of the number of genes shared by n traits across cell types
genehist <- ggplot(sigprs, aes(Num.Traits)) + geom_histogram(bins = 14) +
  xlab(label = "Number of traits across cell types") +
  ylab(label = "Number of genes shared by n traits") +
  scale_x_continuous(breaks = c(0:13), limits = c(0,13)) +
  scale_y_continuous(breaks = c(0,200,400,600,800,1000,1200,1400)) +
  ggtitle("Shared PRS-associated genes across traits and cell types")

#Plot stacked bar chart of number of genes shared by n traits within each cell type
c1sig$Cell.Type <- "CD4+CD45RO-"
c2sig$Cell.Type <- "CD4+CD45RO+"
c3sig$Cell.Type <- "CD8+CD45RO-"
c4sig$Cell.Type <- "CD8+CD45RO+"
allcells <- bind_rows(c1sig, c2sig, c3sig, c4sig) %>%
  select(Gene, Num.Traits, Cell.Type)
allcells$Num.Traits <- as.factor(allcells$Num.Traits)

cellhist <- ggplot(allcells, aes(x = Num.Traits, fill = Cell.Type)) +
  geom_histogram(position = "dodge", stat = "count") +
  scale_fill_manual(values = c("CD4+CD45RO-" = "#e41a1c",
                                "CD4+CD45RO+" = "#4daf4a",
                                "CD8+CD45RO-" = "#377eb8",
                                "CD8+CD45RO+" = "#984ea3")) +
  labs(fill = "Cell type") + xlab(label = "Number of traits") +
  scale_x_discrete(breaks = c(0:6)) + ylab(label = "Number of genes") +
  scale_y_continuous(breaks = c(0,200,400,600,800,1000,1200,1400,1600,
                                1800,2000,2200,2400,2600,2800,3000)) +
  ggtitle("Shared PRS-associated genes within cell types")

supp2 <- plot_grid(genehist, cellhist)