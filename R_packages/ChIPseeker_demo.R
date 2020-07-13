#Gene Ontology (GO) annotates genes to biological processes, molecular functions, and cellular components in a directed acyclic graph structure, 
#Kyoto Encyclopedia of Genes and Genomes (KEGG) annotates genes to pathways 
#Disease Ontology (DO) annotates genes with human disease association
#Reactome annotates gene to pathways and reactions.

#ChIPseeker was developed for annotating nearest genes and genomic features to peaks.
#ChIPseeker包是个好东西,以后要用

## loading packages
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
library(clusterProfiler)
library(org.Hs.eg.db)
#******1. ChIP profiling******
files <- getSampleFiles()
files

peak <- readPeakFile(files[[4]])
peak

#***ChIP peaks coverage plot
covplot(peak, weightCol="V5")
covplot(peak, weightCol="V5", chrs=c("chr17", "chr18"), xlim=c(4.5e7, 5e7))

#***Profile of ChIP peaks binding to TSS regions
#promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
#tagMatrix <- getTagMatrix(peak, windows=promoter)

## to speed up the compilation of this vignettes, we use a precalculated tagMatrix
data("tagMatrixList")
tagMatrix <- tagMatrixList[[4]]

#Heatmap of ChIP binding to TSS regions
tagHeatmap(tagMatrix, xlim=c(-3000, 3000), color="red")  # ==
#peakHeatmap(files[[4]], TxDb=txdb, upstream=3000, downstream=3000, color="red")

#Average Profile of ChIP peaks binding to TSS region
plotAvgProf(tagMatrix, xlim=c(-3000, 3000),
            xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")

#plotAvgProf2(files[[4]], TxDb=txdb, upstream=3000, downstream=3000,xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")

#Confidence interval estimated by bootstrap method is also supported for characterizing ChIP binding profiles.
plotAvgProf(tagMatrix, xlim=c(-3000, 3000), conf = 0.95, resample = 1000)

#***Profile of ChIP peaks binding to start site of Exon/Intron
getBioRegion(TxDb = txdb)



#******2. Peak Annotation******
peakAnno <- annotatePeak(files[[4]], tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Hs.eg.db")

#***Visualize Genomic Annotation
plotAnnoPie(peakAnno)
plotAnnoBar(peakAnno)
vennpie(peakAnno)
upsetplot(peakAnno)
upsetplot(peakAnno, vennpie=TRUE)

#***Visualize distribution of TF-binding loci relative to TSS
plotDistToTSS(peakAnno,
              title="Distribution of transcription factor-binding loci\nrelative to TSS")

#******3. Functional enrichment analysis******
library(ReactomePA)

pathway1 <- enrichPathway(as.data.frame(peakAnno)$geneId)
head(pathway1, 2)

gene <- seq2gene(peak, tssRegion = c(-1000, 1000), flankDistance = 3000, TxDb=txdb)
pathway2 <- enrichPathway(gene)
head(pathway2, 2)

dotplot(pathway2)

#******4. ChIP peak data set comparison******
#***Profile of several ChIP peak data binding to TSS region
#Average profiles

## promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
## tagMatrixList <- lapply(files, getTagMatrix, windows=promoter)
##
## to speed up the compilation of this vigenette, we load a precaculated tagMatrixList
data("tagMatrixList")
plotAvgProf(tagMatrixList, xlim=c(-3000, 3000))

plotAvgProf(tagMatrixList, xlim=c(-3000, 3000), conf=0.95,resample=500, facet="row")

#Peak heatmaps
tagHeatmap(tagMatrixList, xlim=c(-3000, 3000), color=NULL)

#***ChIP peak annotation comparision
peakAnnoList <- lapply(files, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=FALSE)

plotAnnoBar(peakAnnoList)
plotDistToTSS(peakAnnoList)

#***Functional profiles comparison
genes = lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
names(genes) = sub("_", "\n", names(genes))
compKEGG <- compareCluster(geneCluster   = genes,
                           fun           = "enrichKEGG",
                           pvalueCutoff  = 0.05,
                           pAdjustMethod = "BH")
plot(compKEGG, showCategory = 15, title = "KEGG Pathway Enrichment Analysis")

#***Overlap of peaks and annotated genes
genes= lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
vennplot(genes)

#******5. Statistical testing of ChIP seq overlap******
#***Shuffle genome coordination
p <- GRanges(seqnames=c("chr1", "chr3"),
             ranges=IRanges(start=c(1, 100), end=c(50, 130)))
shuffle(p, TxDb=txdb)

#***Peak overlap enrichment analysis
enrichPeakOverlap(queryPeak     = files[[5]],
                  targetPeak    = unlist(files[1:4]),
                  TxDb          = txdb,
                  pAdjustMethod = "BH",
                  nShuffle      = 50,
                  chainFile     = NULL,
                  verbose       = FALSE)


#******6. Data Mining with ChIP seq data deposited in GEO******
#***GEO data collection
getGEOspecies()
getGEOgenomeVersion()
hg19 <- getGEOInfo(genome="hg19", simplify=TRUE)
head(hg19)

#***Download GEO ChIP data sets
downloadGEObedFiles(genome="hg19", destDir="hg19")
gsm <- hg19$gsm[sample(nrow(hg19), 10)]
downloadGSMbedFiles(gsm, destDir="hg19")

#***Overlap significant testing


sessionInfo()