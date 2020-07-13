library(scRNAseq)
data(allen)
allen
sce <- as(allen, "SingleCellExperiment")
sce
isSpike(sce, "ERCC") <- grepl("^ERCC-", rownames(sce))
sce
table(isSpike(sce, "ERCC"))
isSpike(sce, "Adam") <- grepl("^Adam[0-9]", rownames(sce))
sce

sizeFactors(sce) <- colSums(assay(sce))
head(sizeFactors(sce))

colData(sce)
rowData(sce)

library(magrittr)
assay(sce) %>% log1p %>% rowVars -> vars
names(vars) <- rownames(sce)
vars <- sort(vars, decreasing = TRUE)

sce_sub <- sce[names(vars[1:100]),]
sce_sub

library(Rtsne)
set.seed(5252)

pca_data <- prcomp(t(log1p(assay(sce_sub))))
tsne_data <- Rtsne(pca_data$x[,1:50], pca = FALSE)

reducedDims(sce_sub) <- SimpleList(PCA=pca_data$x, TSNE=tsne_data$Y)
sce_sub

reducedDims(sce_sub)
reducedDimNames(sce_sub)
head(reducedDim(sce_sub, "PCA")[,1:3])
head(reducedDim(sce_sub, "TSNE")[,1:2])

dim(reducedDim(sce_sub, "PCA"))
dim(reducedDim(sce_sub[,1:10], "PCA"))


##----scater----
library(scater)
data("sc_example_counts")
data("sc_example_cell_info")
example_sce <- SingleCellExperiment(assays = list(counts = sc_example_counts), colData = sc_example_cell_info)
example_sce

str(counts(example_sce))
example_sce$whee <- sample(LETTERS, ncol(example_sce), replace=TRUE)
colData(example_sce)
rowData(example_sce)



































