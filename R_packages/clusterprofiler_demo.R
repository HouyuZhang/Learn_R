#clusterProfiler包主要用来:
#1. GO分析: GO classification -- groupGO()
#           GO over-representation test -- enrichGO()
#           GO Gene Set Enrichment Analysis -- gseGO()
#2. KEGG分析
#           KEGG over-representation test -- enrichKEGG()
#           KEGG Gene Set Enrichment Analysis -- gseKEGG()
#3. DAVID分析 -- enrichDAVID()
#4. 导入GSEA的gmt进行分析 -- read.gmt()
#5. Biological theme comparison -- compareCluster()

#Gene Ontology (GO) annotates genes to biological processes, molecular functions, and cellular components in a directed acyclic graph structure, 
#Kyoto Encyclopedia of Genes and Genomes (KEGG) annotates genes to pathways 
#Disease Ontology (DO) annotates genes with human disease association
#Reactome annotates gene to pathways and reactions.

#****bitr
library(clusterProfiler)
x <- c("GPX3",  "GLRX",   "LBP",   "CRYAB", "DEFB1", "HCLS1",   "SOD2",   "HSPA2",
       "ORM1",  "IGFBP1", "PTHLH", "GPC3",  "IGFBP3","TOB1",    "MITF",   "NDRG1",
       "NR1H4", "FGFR3",  "PVR",   "IL6",   "PTPRM", "ERBB2",   "NID2",   "LAMB1",
       "COMP",  "PLS3",   "MCAM",  "SPP1",  "LAMC1", "COL4A2",  "COL4A1", "MYOC",
       "ANXA4", "TFPI2",  "CST6",  "SLPI",  "TIMP2", "CPM",     "GGT1",   "NNMT",
       "MAL",   "EEF1A2", "HGD",   "TCN2",  "CDA",   "PCCA",    "CRYM",   "PDXK",
       "STC1",  "WARS",  "HMOX1", "FXYD2", "RBP4",   "SLC6A12", "KDELR3", "ITM2B")
eg = bitr(x, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
head(eg)
library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)
ids <- bitr(x, fromType="SYMBOL", toType=c("UNIPROT", "ENSEMBL"), OrgDb="org.Hs.eg.db")
head(ids)
#*********************1. GO Analysis*************************************
#****GO classification
data(geneList, package="DOSE")
gene <- names(geneList)[abs(geneList) > 2]
gene.df <- bitr(gene, fromType = "ENTREZID",
                toType = c("ENSEMBL", "SYMBOL"),
                OrgDb = org.Hs.eg.db)
head(gene.df)
ggo <- groupGO(gene     = gene,
               OrgDb    = org.Hs.eg.db,
               ont      = "CC",
               level    = 3,
               readable = TRUE)

head(ggo)
barplot(ggo, drop=TRUE, showCategory=12)

#****GO over-representation test
ego <- enrichGO(gene          = gene,
                universe      = names(geneList),
                OrgDb         = org.Hs.eg.db,
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
head(ego)

barplot(ego, showCategory=8)
dotplot(ego)
enrichMap(ego)
cnetplot(ego, categorySize="pvalue", foldChange=geneList)
plotGOgraph(ego)

ego2 <- enrichGO(gene         = gene.df$ENSEMBL,
                 OrgDb         = org.Hs.eg.db,
                 keyType       = "ENSEMBL",
                 ont           = "CC",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05)
ego2 <- setReadable(ego2, OrgDb = org.Hs.eg.db)
head(ego2)

#****GO Gene Set Enrichment Analysis
ego3 <- gseGO(geneList     = geneList,
              OrgDb        = org.Hs.eg.db,
              ont          = "CC",
              nPerm        = 1000,
              minGSSize    = 100,
              maxGSSize    = 500,
              pvalueCutoff = 0.05,
              verbose      = FALSE)
head(ego3)


#*********************2. KEGG analysis*************************************

search_kegg_organism('mmu', by='kegg_code')

#****KEGG over-representation test

kk <- enrichKEGG(gene         = gene,
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)
head(kk)
browseKEGG(kk, 'hsa04110')

#KEGG Module over-representation test
mkk <- enrichMKEGG(gene = gene,
                   organism = 'hsa')

#****KEGG Gene Set Enrichment Analysis

kk2 <- gseKEGG(geneList     = geneList,
               organism     = 'hsa',
               nPerm        = 1000,
               minGSSize    = 120,
               pvalueCutoff = 0.05,
               verbose      = FALSE)
head(kk2)
gseaplot(kk2, geneSetID = "hsa04145")
#KEGG Module Gene Set Enrichment Analysis
mkk2 <- gseMKEGG(geneList = geneList,
                 organism = 'hsa')

#3. Disease analysis
#4. Reactome pathway analysis
#********************************5. DAVID functional analysis********************
source("https://bioconductor.org/biocLite.R")
biocLite(c("RDAVIDWebService","ReactomePA"))
david <- enrichDAVID(gene = gene,
                     idType = "ENTREZ_GENE_ID",
                     annotation = "KEGG_PATHWAY",
                     david.user = "clusterProfiler@hku.hk")

#******************************6. Universal enrichment analysis*****************
gmtfile <- system.file("extdata", "c5.cc.v5.0.entrez.gmt", package="clusterProfiler")
c5 <- read.gmt(gmtfile)

egmt <- enricher(gene, TERM2GENE=c5)
head(egmt)

egmt2 <- GSEA(geneList, TERM2GENE=c5, verbose=FALSE)
head(egmt2)




library("pathview")
hsa04110 <- pathview(gene.data  = geneList,
                     pathway.id = "hsa04110",
                     species    = "hsa",
                     limit      = list(gene=max(abs(geneList)), cpd=1))

#***************************7. Biological theme comparison**************
data(gcSample)
lapply(gcSample, head)
ck <- compareCluster(geneCluster = gcSample, fun = "enrichKEGG")
head(as.data.frame(ck))
dotplot(ck)
#****Formula interface of compareCluster

mydf <- data.frame(Entrez=names(geneList), FC=geneList)
mydf <- mydf[abs(mydf$FC) > 1,]
mydf$group <- "upregulated"
mydf$group[mydf$FC < 0] <- "downregulated"
mydf$othergroup <- "A"
mydf$othergroup[abs(mydf$FC) > 2] <- "B"
library(org.Hs.eg.db)
formula_res <- compareCluster(Entrez~group+othergroup, data=mydf, fun="enrichGO", OrgDb         = org.Hs.eg.db)

head(as.data.frame(formula_res))

dotplot(formula_res)
dotplot(formula_res, x=~group) + ggplot2::facet_grid(~othergroup)

sessionInfo()
