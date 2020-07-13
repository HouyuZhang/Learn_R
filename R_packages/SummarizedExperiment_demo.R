#这个包有两个类:
#  SummarizedExperiment
#    RangedSummarizedExperiment
  
library(SummarizedExperiment)
data(airway, package="airway")
se <- airway
se
#******1. Anatomy of a SummarizedExperiment******
assays(se)$counts

rowRanges(se)

colData(se)
# subset for only those samples treated with dexamethasone
se[, se$dex == "trt"]

metadata(se)

metadata(se)$formula <- counts ~ dex + albut

metadata(se)


#******2. Constructing a SummarizedExperiment******
nrows <- 200
ncols <- 6
counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
rowRanges <- GRanges(rep(c("chr1", "chr2"), c(50, 150)),
                     IRanges(floor(runif(200, 1e5, 1e6)), width=100),
                     strand=sample(c("+", "-"), 200, TRUE),
                     feature_id=sprintf("ID%03d", 1:200))
colData <- DataFrame(Treatment=rep(c("ChIP", "Input"), 3),
                     row.names=LETTERS[1:6])

SummarizedExperiment(assays=list(counts=counts),
                     rowRanges=rowRanges, colData=colData)

SummarizedExperiment(assays=list(counts=counts), colData=colData)

#******3. Common operations on SummarizedExperiment******







