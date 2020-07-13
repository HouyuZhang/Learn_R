if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("simpleSingleCell", version = "3.8")
library(simpleSingleCell)
#---- Analyzing single-cell RNA-seq data containing read counts----
###1. Setting up the data
#1.1 Loading in the count matrix
lun.sdrf <- read.delim("raw_data/326010d4cda_E-MTAB-5522.sdrf.txt")
plate1 <- read.delim("raw_data/counts_Calero_20160113.tsv", header=TRUE, row.names=1, check.names=FALSE)
plate2 <- read.delim("raw_data/counts_Calero_20160325.tsv", header=TRUE, row.names=1, check.names=FALSE)

gene.lengths <- plate1$Length # First column is the gene length.
plate1 <- as.matrix(plate1[,-1]) # Discarding gene length (as it is not a cell).
plate2 <- as.matrix(plate2[,-1])
rbind(Plate1=dim(plate1), Plate2=dim(plate2))

#verifying that the genes are in the same order between the two matrices
stopifnot(identical(rownames(plate1), rownames(plate2)))
all.counts <- cbind(plate1, plate2)

library(SingleCellExperiment)
sce <- SingleCellExperiment(list(counts=as.matrix(all.counts)))
rowData(sce)$GeneLength <- gene.lengths
sce

#We identify the rows corresponding to ERCC spike-in transcripts from the row names. We store this information in the SingleCellExperiment object for future use.
isSpike(sce, "ERCC") <- grepl("^ERCC", rownames(sce))
summary(isSpike(sce, "ERCC"))

#we must remove the rows corresponding to the SIRV transcripts prior to further analysis
is.sirv <- grepl("^SIRV", rownames(sce))
sce <- sce[!is.sirv,] 
summary(is.sirv)

#1.2 Incorporating cell-based annotation
metadata <- read.delim("raw_data/326010d4cda_E-MTAB-5522.sdrf.txt", check.names=FALSE, header=TRUE)
m <- match(colnames(sce), metadata[["Source Name"]]) # Enforcing identical order.
stopifnot(all(!is.na(m))) # Checking that nothing's missing.
metadata <- metadata[m,]
head(colnames(metadata))

colData(sce)$Plate <- factor(metadata[["Factor Value[block]"]])
pheno <- metadata[["Factor Value[phenotype]"]]
levels(pheno) <- c("induced", "control")
colData(sce)$Oncogene <- pheno
table(colData(sce)$Oncogene, colData(sce)$Plate)

#1.3 Incorporating gene-based annotation
library(org.Mm.eg.db)
symb <- mapIds(org.Mm.eg.db, keys=rownames(sce), keytype="ENSEMBL", column="SYMBOL")
rowData(sce)$ENSEMBL <- rownames(sce)
rowData(sce)$SYMBOL <- symb
head(rowData(sce))

# replace missing symbols with the Ensembl identifier and concatenate duplicated symbols with the (unique) Ensembl identifiers.
new.names <- rowData(sce)$SYMBOL
missing.name <- is.na(new.names)
new.names[missing.name] <- rowData(sce)$ENSEMBL[missing.name]

dup.name <- new.names %in% new.names[duplicated(new.names)]
new.names[dup.name] <- paste0(new.names, "_", rowData(sce)$ENSEMBL)[dup.name]
rownames(sce) <- new.names
head(rownames(sce))

#determine the chromosomal location for each gene
library(TxDb.Mmusculus.UCSC.mm10.ensGene)
location <- mapIds(TxDb.Mmusculus.UCSC.mm10.ensGene, keys=rowData(sce)$ENSEMBL, column="CDSCHROM", keytype="GENEID")
rowData(sce)$CHR <- location
summary(location=="chrM")

###2. Quality control on the cells

#2.1 Defining the quality control metrics
library(scater)
mito <- which(rowData(sce)$CHR=="chrM")
sce <- calculateQCMetrics(sce, feature_controls=list(Mt=mito))
head(colnames(colData(sce)), 10)

#to remove putative low-quality cells that have low library sizes, low numbers of expressed features, and high spike-in (or mitochondrial) proportions.
sce$PlateOnco <- paste0(sce$Oncogene, ".", sce$Plate)
multiplot(
  plotColData(sce, y="total_counts", x="PlateOnco"),
  plotColData(sce, y="total_features_by_counts", x="PlateOnco"),
  plotColData(sce, y="pct_counts_ERCC", x="PlateOnco"),
  plotColData(sce, y="pct_counts_Mt", x="PlateOnco"),
  cols=2)

#Clear discrepancies may correspond to technical differences between batches of cells or genuine biological differences in RNA content.
par(mfrow=c(1,3))
plot(sce$total_features_by_counts, sce$total_counts/1e6, xlab="Number of expressed genes",
     ylab="Library size (millions)")
plot(sce$total_features_by_counts, sce$pct_counts_ERCC, xlab="Number of expressed genes",
     ylab="ERCC proportion (%)")
plot(sce$total_features_by_counts, sce$pct_counts_Mt, xlab="Number of expressed genes",
     ylab="Mitochondrial proportion (%)")

#2.2 Identifying outliers for each metric

libsize.drop <- isOutlier(sce$total_counts, nmads=3, type="lower", log=TRUE, batch=sce$PlateOnco)
feature.drop <- isOutlier(sce$total_features_by_counts, nmads=3, type="lower", log=TRUE, batch=sce$PlateOnco)

spike.drop <- isOutlier(sce$pct_counts_ERCC, nmads=3, type="higher", batch=sce$PlateOnco)

keep <- !(libsize.drop | feature.drop | spike.drop)
data.frame(ByLibSize=sum(libsize.drop), ByFeature=sum(feature.drop),BySpike=sum(spike.drop), Remaining=sum(keep))

sce$PassQC <- keep
saveRDS(sce, file="416B_preQC.rds")
sce <- sce[,keep]
dim(sce)

attr(libsize.drop, "thresholds")
attr(spike.drop, "thresholds")

###3. Classification of cell cycle phase
#Pairs with changes in the sign across cell cycle phases were chosen as markers. 
#Cells in a test dataset can then be classified into the appropriate phase, 
#based on whether the observed sign for each marker pair is consistent with one phase or another.
set.seed(100)
library(scran)
mm.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", package="scran"))
assignments <- cyclone(sce, mm.pairs, gene.names=rowData(sce)$ENSEMBL)
plot(assignments$score$G1, assignments$score$G2M, xlab="G1 score", ylab="G2/M score", pch=16)

sce$phases <- assignments$phases
table(sce$phases)

###4. Examining gene-level expression metrics
#4.1 Inspecting the most highly expressed genes
fontsize <- theme(axis.text=element_text(size=12), axis.title=element_text(size=16))
plotHighestExprs(sce, n=50) + fontsize

#4.2 Filtering out low-abundance genes

#the average count for each gene
ave.counts <- calcAverage(sce, use_size_factors=FALSE)
hist(log10(ave.counts), breaks=100, main="", col="grey80", xlab=expression(Log[10]~"average count"))

demo.keep <- ave.counts >= 1
filtered.sce <- sce[demo.keep,]
summary(demo.keep)

#the number of cells that express each gene
num.cells <- nexprs(sce, byrow=TRUE)
smoothScatter(log10(ave.counts), num.cells, ylab="Number of cells", xlab=expression(Log[10]~"average count"))

#remove genes that are not expressed in any cell
to.keep <- num.cells > 0
sce <- sce[to.keep,]
summary(to.keep)

###5. Normalization of cell-specific biases
#5.1 Using the deconvolution method to deal with zero counts
#we pool counts from many cells to increase the count size for accurate size factor estimation.
#Pool-based size factors are then ¡°deconvolved¡± into cell-based factors for cell-specific normalization.
sce <- computeSumFactors(sce)
summary(sizeFactors(sce))

plot(sce$total_counts/1e6, sizeFactors(sce), log="xy",
     xlab="Library size (millions)", ylab="Size factor",
     col=c("red", "black")[sce$Oncogene], pch=16)
legend("bottomright", col=c("red", "black"), pch=16, cex=1.2,legend=levels(sce$Oncogene))