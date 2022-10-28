https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE163943
```
suppressPackageStartupMessages({
	source("/usr/local/prog/scripts/plotThemes.R")
	source("/usr/local/prog/scripts/readmatrixkits.R")
	source("/usr/local/prog/scripts/enrichmentkits.R")
})
library(tidyverse)

###matrix <- readMatrixFromFile (matrix = "TCGA_COAD.TPM.tab")
matrix <- read.table(file = "GSE163943.protein_coding.txt", row.names = 1, header = T)
colnames(matrix) <- c("Blood_CTL_R1", "Blood_CTL_R2", "Blood_CTL_R3", "Blood_CTL_R4", "Blood_POD_R1", "Blood_POD_R2", "Blood_POD_R3", "Blood_POD_R4")

metasheet <- read.table(file = "metasheet.txt", row.names = NULL, header = T)
metasheet <- metasheet[,c("SID", "TYPE")]
colnames(metasheet) <- c("ID", "GROUP")
rownames(metasheet) <- c("Blood_CTL_R1", "Blood_CTL_R2", "Blood_CTL_R3", "Blood_CTL_R4", "Blood_POD_R1", "Blood_POD_R2", "Blood_POD_R3", "Blood_POD_R4")

sig.matrix <- diffIimma (matrix = matrix, metasheet = metasheet, ref = "CTL", exp = "POD", logFC.cutoff = log(1.5, 2), adj.p.cutoff = 1, output = "GSE163943.DEG")
sig.matrix <- sig.matrix[sig.matrix$P.Value <= 0.05,]
sig.matrix <- arrange(sig.matrix, desc(logFC), P.Value)

library(tidyverse)
library(tidyr)
sig.matrix$NAME <- str_extract(rownames(sig.matrix), "[|].*[|]")
sig.matrix$NAME <- gsub("[|]", "", sig.matrix$NAME)

write.table(x = sig.matrix, file = "GSE163943.DEG.POD_vs_CTL.diff_limma.significant.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
```

heatmap
```

dat <- sig.matrix[,1:8]
colnames(dat) <- c("Blood_CTL_R1", "Blood_CTL_R2", "Blood_CTL_R3", "Blood_CTL_R4", "Blood_POD_R1", "Blood_POD_R2", "Blood_POD_R3", "Blood_POD_R4")
dat <- data.frame(dat)

library(ggplot2)
library(pheatmap)
anno_col <- data.frame(row.names = metasheet$ID, GROUP = metasheet$GROUP)
rownames(anno_col) <- c("Blood_CTL_R1", "Blood_CTL_R2", "Blood_CTL_R3", "Blood_CTL_R4", "Blood_POD_R1", "Blood_POD_R2", "Blood_POD_R3", "Blood_POD_R4")
anno_col$GROUP <- factor(anno_col$GROUP)

p <- pheatmap(mat = dat,
         show_rownames = F,
         show_colnames = T,
         scale = "row",
         cluster_rows = T,
         cluster_cols = F,
         #cutree_col = 2,
         cutree_row = 2,
         #gaps_row = c(4),
         gaps_col = c(4),
         treeheight_row = 40,
         treeheight_col = 40,
         color = colorRampPalette(c("navy", "white", "red"))(100),
         border = FALSE,
         #border_color = "white",
         #annotation_colors = ann_colors,
         #annotation_row = anno_row,
         annotation_col = anno_col,
         #cellwidth = 6,
         #cellheight = 6,
         fontsize = 6)

#w_size <- ncol(matrix) * 0.25 + 10
#h_size <- nrow(matrix) * 0.25 + 1
w_size <- 6
h_size <- 8

ggsave(filename = paste("GSE163943.DEG", "heatmap.pdf", sep = "."), plot = p, width = w_size, height = h_size, units = "cm")
```

TOP20基因
```

matrix <- read.table(file = "GSE163943.protein_coding.txt", row.names = 1, header = T)
metasheet <- read.table(file = "metasheet.txt", row.names = NULL, header = T)
metasheet <- metasheet[,c("GID", "TYPE")]
colnames(metasheet) <- c("ID", "GROUP")

sig.matrix <- diffIimma (matrix = matrix, metasheet = metasheet, ref = "CTL", exp = "POD", logFC.cutoff = log(1.5, 2), adj.p.cutoff = 1, output = "GSE163943.DEG")
sig.matrix <- sig.matrix[sig.matrix$P.Value <= 0.05,]
sig.matrix <- arrange(sig.matrix, desc(logFC), P.Value)


dat <- sig.matrix[,1:8]
colnames(dat) <- c("Blood_CTL_R1", "Blood_CTL_R2", "Blood_CTL_R3", "Blood_CTL_R4", "Blood_POD_R1", "Blood_POD_R2", "Blood_POD_R3", "Blood_POD_R4")
dat <- data.frame(dat)


anno_col <- data.frame(row.names = metasheet$ID, GROUP = metasheet$GROUP)
rownames(anno_col) <- c("Blood_CTL_R1", "Blood_CTL_R2", "Blood_CTL_R3", "Blood_CTL_R4", "Blood_POD_R1", "Blood_POD_R2", "Blood_POD_R3", "Blood_POD_R4")
anno_col$GROUP <- factor(anno_col$GROUP)

dat <- rbind(head(dat, 10), tail(dat, 10))
rownames(dat) <- str_extract(string = rownames(dat), pattern = "[|].*[|]")
rownames(dat) <- gsub(pattern = "[|]", replacement = "", x = rownames(dat))

library(ggplot2)
library(pheatmap)


p <- pheatmap(mat = dat,
         show_rownames = T,
         show_colnames = T,
         scale = "row",
         cluster_rows = T,
         cluster_cols = F,
         #cutree_col = 2,
         cutree_row = 2,
         #gaps_row = c(4),
         gaps_col = c(4),
         treeheight_row = 30,
         treeheight_col = 30,
         color = colorRampPalette(c("navy", "white", "red"))(100),
         border = FALSE,
         #border_color = "white",
         #annotation_colors = ann_colors,
         #annotation_row = anno_row,
         annotation_col = anno_col,
         #cellwidth = 6,
         #cellheight = 6,
         fontsize = 6)

#w_size <- ncol(matrix) * 0.25 + 10
#h_size <- nrow(matrix) * 0.25 + 1
w_size <- 7
h_size <- 8

ggsave(filename = paste("GSE163943.DEG", "TOP20.heatmap.pdf", sep = "."), plot = p, width = w_size, height = h_size, units = "cm")

```

3. 火山图
```
suppressPackageStartupMessages({
	source("/usr/local/prog/scripts/plotThemes.R")
	source("/usr/local/prog/scripts/readmatrixkits.R")
	source("/usr/local/prog/scripts/enrichmentkits.R")
})

matrix <- read.table(file = "GSE163943.protein_coding.txt", row.names = 1, header = T)
metasheet <- read.table(file = "metasheet.txt", row.names = NULL, header = T)
metasheet <- metasheet[,c("GID", "TYPE")]
colnames(metasheet) <- c("ID", "GROUP")

sig.matrix <- diffIimma (matrix = matrix, metasheet = metasheet, ref = "CTL", exp = "POD", logFC.cutoff = log(1.0, 2), adj.p.cutoff = 1, output = "GSE163943.DEG")
executeVolcanoPlot (matrix = sig.matrix, output = "GSE163943.DEG", logFC_cutoff = log(1.5, 2), pvalue_cutoff = -log(0.05, 10), w_size = 5, h_size = 4.5)
```

4. boxplot
```
suppressPackageStartupMessages({
	source("/usr/local/prog/scripts/plotThemes.R")
	source("/usr/local/prog/scripts/readmatrixkits.R")
	source("/usr/local/prog/scripts/enrichmentkits.R")
})
library(tidyverse)

matrix <- read.table(file = "GSE163943.protein_coding.txt", row.names = 1, header = T)
metasheet <- read.table(file = "metasheet.txt", row.names = NULL, header = T)
metasheet <- metasheet[,c("GID", "TYPE")]
colnames(metasheet) <- c("ID", "GROUP")
rownames(metasheet) <- c("Blood_CTL_R1", "Blood_CTL_R2", "Blood_CTL_R3", "Blood_CTL_R4", "Blood_POD_R1", "Blood_POD_R2", "Blood_POD_R3", "Blood_POD_R4")
metasheet$ID <- rownames(metasheet)

sig.matrix <- diffIimma (matrix = matrix, metasheet = metasheet, ref = "CTL", exp = "POD", logFC.cutoff = log(1.5, 2), adj.p.cutoff = 1, output = "GSE163943.DEG")
sig.matrix <- sig.matrix[sig.matrix$P.Value <= 0.05,]
sig.matrix <- arrange(sig.matrix, desc(logFC), P.Value)

dat <- sig.matrix[,1:8]
colnames(dat) <- c("Blood_CTL_R1", "Blood_CTL_R2", "Blood_CTL_R3", "Blood_CTL_R4", "Blood_POD_R1", "Blood_POD_R2", "Blood_POD_R3", "Blood_POD_R4")
dat <- data.frame(dat)

dat <- rbind(head(dat, 10), tail(dat, 10))
rownames(dat) <- str_extract(string = rownames(dat), pattern = "[|].*[|]")
rownames(dat) <- gsub(pattern = "[|]", replacement = "", x = rownames(dat))

executeBoxPlot (matrix = dat, metasheet = metasheet, compaired = list(c("POD", "CTL")), method = "t.test", output = "GSE163943.DEG", col.num = 5) 
```

5. DEG功能分析
```
comparisons.grp <- "POD_vs_CTL"
geneset <- sig.matrix$NAME

enrichmentGeneSet(geneset = geneset, ontology = "BP", species = "human", output = paste("DEG.BP", comparisons.grp, sep = "."))
enrichmentGeneSet(geneset = geneset, ontology = "MF", species = "human", output = paste("DEG.MF", comparisons.grp, sep = "."))
enrichmentGeneSet(geneset = geneset, ontology = "CC", species = "human", output = paste("DEG.CC", comparisons.grp, sep = "."))
enrichmentGeneSet(geneset = geneset, ontology = "KEGG", species = "human", output = paste("DEG.KEGG", comparisons.grp, sep = "."))

plotGSVAHeatmap (matrix = sig.matrix, metasheet = metasheet, ref = "POD", exp = "CTL", output = "DEG.GSVA")
```

6. aging.genes
```
aging.genes <- read.table("/home/wangk/lab/yinyibo/aging.gene.txt")
colnames(aging.genes) <- "GENEID"
aging.genes <- aging.genes$GENEID

aging.genes.sig.matrix <- sig.matrix[sig.matrix$NAME %in% aging.genes,]
write.table(x = aging.genes.sig.matrix, file = "GSE163943.DEG.aging.genes.POD_vs_CTL.diff_limma.significant.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)


############
### heatmap
############
dat <- aging.genes.sig.matrix
rownames(dat) <- dat$NAME
dat <- dat[,1:8]

anno_col <- data.frame(row.names = metasheet$ID, GROUP = metasheet$GROUP)
rownames(anno_col) <- c("Blood_CTL_R1", "Blood_CTL_R2", "Blood_CTL_R3", "Blood_CTL_R4", "Blood_POD_R1", "Blood_POD_R2", "Blood_POD_R3", "Blood_POD_R4")
anno_col$GROUP <- factor(anno_col$GROUP)

library(ggplot2)
library(pheatmap)

p <- pheatmap(mat = dat,
         show_rownames = T,
         show_colnames = T,
         scale = "row",
         cluster_rows = T,
         cluster_cols = F,
         #cutree_col = 2,
         cutree_row = 2,
         #gaps_row = c(4),
         gaps_col = c(4),
         treeheight_row = 30,
         treeheight_col = 30,
         color = colorRampPalette(c("navy", "white", "red"))(100),
         border = FALSE,
         #border_color = "white",
         #annotation_colors = ann_colors,
         #annotation_row = anno_row,
         annotation_col = anno_col,
         #cellwidth = 6,
         #cellheight = 6,
         fontsize = 6)

#w_size <- ncol(matrix) * 0.25 + 10
#h_size <- nrow(matrix) * 0.25 + 1
w_size <- 7
h_size <- 6

ggsave(filename = paste("GSE163943.DEG", "AGING.heatmap.pdf", sep = "."), plot = p, width = w_size, height = h_size, units = "cm")
```

