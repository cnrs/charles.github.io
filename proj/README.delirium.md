https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE163943
```
suppressPackageStartupMessages({
	source("/usr/local/prog/scripts/plotThemes.R")
	source("/usr/local/prog/scripts/readmatrixkits.R")
	source("/usr/local/prog/scripts/enrichmentkits.R")
})

###matrix <- readMatrixFromFile (matrix = "TCGA_COAD.TPM.tab")
matrix <- read.table(file = "GSE163943.protein_coding.txt", row.names = 1, header = T)
metasheet <- readMetasheetFromFile(matrix = "consensusclusterplus.sample.clusters.txt")
metasheet <- readMetasheetFromFile(matrix = "metasheet.txt")

sig.matrix <- diffIimma (matrix = matrix, metasheet = "metasheet.txt", ref = "Normal", exp = "Tumor", logFC.cutoff = log(2, 2), adj.p.cutoff = 0.05, output = "DEG")
gene_list_degs <- unique(rownames(sig.matrix))

gene_list <- c()
cluster.num = 6
for (i in 1:cluster.num){
	flag = 0
	for (j in 1:cluster.num){
		if (i == j) {next}
		comparisons.grp <- paste (i, "_vs_", j, sep = "")
		print(paste(i, ", ", j, comparisons.grp, sep = " "))
		sig.matrix <- diffIimma (matrix = matrix, metasheet = "consensusclusterplus.sample.clusters.txt", ref = i, exp = j, logFC.cutoff = log(1.5, 2), adj.p.cutoff = 0.05, output = "DEG.CLUSTER")
		
		gene_list <- c(gene_list, rownames(sig.matrix))
	}
}
gene_list <- unique(gene_list)
gene_list <- intersect (gene_list, gene_list_degs)
write.table(x = gene_list, file = "m1A_phenotype_related_DEGs.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
```
