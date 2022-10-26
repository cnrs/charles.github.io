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
metasheet <- read.table(file = "metasheet.txt", row.names = NULL, header = T)
metasheet <- metasheet[,c("GID", "TYPE")]
colnames(metasheet) <- c("ID", "GROUP")

sig.matrix <- diffIimma (matrix = matrix, metasheet = metasheet, ref = "CTL", exp = "POD", logFC.cutoff = log(1.5, 2), adj.p.cutoff = 1, output = "GSE163943.DEG")
sig.matrix <- sig.matrix[sig.matrix$P.Value <= 0.05,]
sig.matrix <- arrange(sig.matrix, desc(logFC), P.Value)

write.table(x = sig.matrix, file = "GSE163943.sig.matrix.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

library(tidyverse)
library(tidyr)
sig.matrix$NAME <- str_extract(rownames(sig.matrix), "[|].*[|]")
sig.matrix$NAME <- gsub("[|]", "", sig.matrix$NAME)
```

heatmap
```

```
