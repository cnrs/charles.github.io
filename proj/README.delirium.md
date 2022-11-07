https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE163943
```
suppressPackageStartupMessages({
	source("/usr/local/prog/scripts/plotThemes.R")
	source("/usr/local/prog/scripts/readmatrixkits.R")
	source("/usr/local/prog/scripts/enrichmentkits.R")
})
library(tidyverse)

matrix <- readMatrixFromFile (matrix = "GSE163943.protein_coding.txt")
colnames(matrix) <- c("Blood_CTL_R1", "Blood_CTL_R2", "Blood_CTL_R3", "Blood_CTL_R4", "Blood_POD_R1", "Blood_POD_R2", "Blood_POD_R3", "Blood_POD_R4")

metasheet <- read.table(file = "metasheet.txt", row.names = NULL, header = T)
metasheet <- metasheet[,c("SID", "TYPE")]
colnames(metasheet) <- c("ID", "GROUP")
rownames(metasheet) <- metasheet$ID

sig.matrix <- diffIimma (matrix = matrix, metasheet = metasheet, ref = "CTL", exp = "POD", logFC.cutoff = log(1.5, 2), adj.p.cutoff = 1, output = "GSE163943.DEG")
sig.matrix <- sig.matrix[sig.matrix$P.Value <= 0.05,]
sig.matrix <- arrange(sig.matrix, desc(logFC), P.Value)


saveRDS(object = matrix, file = "GSE163943.matrix.rds")
saveRDS(object = metasheet, file = "GSE163943.metasheet.rds")
saveRDS(object = sig.matrix, file = "GSE163943.sig.matrix.rds")

library(tidyverse)
library(tidyr)
sig.matrix$NAME <- str_extract(rownames(sig.matrix), "[|].*[|]")
sig.matrix$NAME <- gsub("[|]", "", sig.matrix$NAME)

write.table(x = sig.matrix, file = "GSE163943.DEG.POD_vs_CTL.diff_limma.significant.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
```

heatmap
```

dat <- sig.matrix[,1:8]
dat <- data.frame(dat)

library(ggplot2)
library(pheatmap)
anno_col <- data.frame(row.names = metasheet$ID, GROUP = metasheet$GROUP)
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
dat <- data.frame(dat)


anno_col <- data.frame(row.names = metasheet$ID, GROUP = metasheet$GROUP)
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

matrix <- readMatrixFromFile (matrix = "GSE163943.protein_coding.txt")

sig.matrix <- diffIimma (matrix = matrix, metasheet = metasheet, ref = "CTL", exp = "POD", logFC.cutoff = log(1.5, 2), adj.p.cutoff = 1, output = "GSE163943.DEG")
sig.matrix <- sig.matrix[sig.matrix$P.Value <= 0.05,]
sig.matrix <- arrange(sig.matrix, desc(logFC), P.Value)

dat <- sig.matrix[,1:8]
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

7. GSEA & GSVA
```
sed -e 's/[[:space:]]/\t/g' GSE163943.protein_coding.txt > GSE163943.protein_coding.txt.1
mv GSE163943.protein_coding.txt.1 GSE163943.protein_coding.txt

suppressPackageStartupMessages({
	source("/usr/local/prog/scripts/plotThemes.R")
	source("/usr/local/prog/scripts/readmatrixkits.R")
	source("/usr/local/prog/scripts/enrichmentkits.R")
})

matrix <- readMatrixFromFile (matrix = "GSE163943.protein_coding.txt")
colnames(matrix) <- c("Blood_CTL_R1", "Blood_CTL_R2", "Blood_CTL_R3", "Blood_CTL_R4", "Blood_POD_R1", "Blood_POD_R2", "Blood_POD_R3", "Blood_POD_R4")

metasheet <- read.table(file = "metasheet.txt", row.names = NULL, header = T)
metasheet <- metasheet[,c("SID", "TYPE")]
colnames(metasheet) <- c("ID", "GROUP")
rownames(metasheet) <- metasheet$ID

sig.matrix <- diffIimma (matrix = matrix, metasheet = metasheet, ref = "CTL", exp = "POD", logFC.cutoff = log(1.0, 2), adj.p.cutoff = 1, output = "GSE163943.DEG")

enrichmentGSEA(matrix = sig.matrix, ontology = "BP", species = "human", pvalue.cutoff = 0.05, output = "GSE163943.ALLGENE", w_size = 7, h_size = 7)
enrichmentGSEA(matrix = sig.matrix, ontology = "KEGG", species = "human", pvalue.cutoff = 0.05, output = "GSE163943.ALLGENE", w_size = 7, h_size = 7)

gs.matrix <- enrichmentGSVA (matrix = matrix[,1:8], method = "gsva", species = "Homo sapiens", ontology = "KEGG", m_s_mnsize = 1, m_x_mxsize = 500, output = "GSE163943.ALLGENE")
sig.matrix <- diffIimma (matrix = gs.matrix, metasheet = metasheet, ref = "CTL", exp = "POD", logFC.cutoff = 0.00, adj.p.cutoff = 1, output = "GSE163943.ALLGENE.GSVA")
plotHistogramForGSVA (matrix = sig.matrix, output = "GSE163943.ALLGENE.GSVA", logFC.cutoff = 0.00, pvalue.cutoff = 0.05, ns = T)



sig.matrix <- diffIimma (matrix = matrix, metasheet = metasheet, ref = "CTL", exp = "POD", logFC.cutoff = log(1.0, 2), adj.p.cutoff = 1, output = "GSE163943.DEG")
sig.matrix <- sig.matrix[abs(sig.matrix$logFC) >= log(1.5, 2) & sig.matrix$P.Value <= 0.05,]
sig.matrix <- arrange(sig.matrix, desc(logFC), P.Value)

enrichmentGSEA(matrix = sig.matrix, ontology = "BP", species = "human", pvalue.cutoff = 0.05, output = "GSE163943.DEG", w_size = 7, h_size = 7)
enrichmentGSEA(matrix = sig.matrix, ontology = "KEGG", species = "human", pvalue.cutoff = 0.05, output = "GSE163943.DEG", w_size = 7, h_size = 7)

gs.matrix <- enrichmentGSVA (matrix = matrix[rownames(sig.matrix), 1:8], method = "gsva", species = "Homo sapiens", ontology = "KEGG", m_s_mnsize = 1, m_x_mxsize = 500, output = "GSE163943.DEG")
sig.matrix <- diffIimma (matrix = gs.matrix, metasheet = metasheet, ref = "CTL", exp = "POD", logFC.cutoff = 0.00, adj.p.cutoff = 1, output = "GSE163943.DEG.GSVA")
plotHistogramForGSVA (matrix = sig.matrix, output = "GSE163943.DEG.GSVA", logFC.cutoff = 0.00, pvalue.cutoff = 0.05, ns = T)

```

9. aging.genes
```
executeCircosPlot (covariates = unique(aging.genes.sig.matrix$NAME), output = "DEG.AGING")
```


10. LASSO鉴定关键衰老基因
```
matrix <- readMatrixFromFile (matrix = "GSE163943.protein_coding.txt")
colnames(matrix) <- c("Blood_CTL_R1", "Blood_CTL_R2", "Blood_CTL_R3", "Blood_CTL_R4", "Blood_POD_R1", "Blood_POD_R2", "Blood_POD_R3", "Blood_POD_R4")

metasheet <- read.table(file = "metasheet.txt", row.names = NULL, header = T)
metasheet <- metasheet[,c("SID", "TYPE")]
colnames(metasheet) <- c("ID", "GROUP")
rownames(metasheet) <- metasheet$ID

### sig.matrix <- diffIimma (matrix = matrix, metasheet = metasheet, ref = "CTL", exp = "POD", logFC.cutoff = log(1.0, 2), adj.p.cutoff = 1, output = "GSE163943.DEG")
### sig.matrix <- sig.matrix[abs(sig.matrix$logFC) >= log(1.5, 2) & sig.matrix$P.Value <= 0.05,]
### sig.matrix <- arrange(sig.matrix, desc(logFC), P.Value)

aging.degs = c("ADCY5", "AKT1", "APP", "BCL6", "COL4A2", "E2F1", "HOXB7", "IGF1", "IL6", "KRT14", "NFE2L2", "RGN", "SERPINF1", "TBX2")
matrix <- mergeMatrixMetasheet (matrix = matrix, metasheet = metasheet)

matrix$GROUP <- ifelse(matrix$GROUP == "POD", 1, 0)
executeLASSOAnalysis(matrix = matrix, phenotype = "GROUP", covariates = aging.degs, model = "binomial", measure = "default", output = "GSE163943.AGING.DEG", w_size = 6.5, h_szie = 8)
```

预测到2个关键基因（GSE163943.AGING.DEG.LASSO.coefficient.txt）：
```
ID	coefficient
COL4A2	0.328030958447826
NFE2L2	3.49659258777546
```

