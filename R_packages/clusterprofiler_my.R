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

setwd("C:/Users/Administrator/Desktop")
library(clusterProfiler)
library(org.Mm.eg.db)
d2 <- read.csv("BvsCD4.csv")
geneList2 <- d2[,2]
names(geneList2) = as.character(d2[,1])
geneList2 = sort(geneList2, decreasing = TRUE)
head(geneList2)
class(geneList2)
length(geneList2)
geneList2_none <- geneList2[abs(geneList2) > 1]


d <- read.csv("BvsCD4_ENTREID.csv")
#制作geneList
geneList <- d[,2]
names(geneList) = as.character(d[,1])
geneList = sort(geneList, decreasing = TRUE)
head(geneList)
class(geneList)
length(geneList)
#转化格式
gene <- names(geneList)[abs(geneList) > 1]
library(org.Mm.eg.db)
keytypes(org.Mm.eg.db)
gene.df <- bitr(gene, fromType ="ENTREZID",
                toType =c("ENSEMBL", "SYMBOL"),
                OrgDb = org.Mm.eg.db)
head(gene.df)
write.csv(gene.df,file="E:/5. Coding/R/R_Conductor/clusterprofiler/gene.df.csv")

#*********************1. GO Analysis*************************************
#####
#****GO classification (测试通过,使用原始ENSEMBL数据)
getGOLevel("CC", 3)
ggo <- groupGO(gene     = names(geneList2_none),
               keytype       = "ENSEMBL",
               OrgDb    = org.Mm.eg.db,
               ont      = "BP",
               level    = 3,
               readable = T)
write.table(ggo,file="E:/5. Coding/R/R_Conductor/clusterprofiler/ggo.txt")
head(ggo)
barplot(ggo, drop=TRUE, showCategory=12)

#****GO over-representation test (把差异比较大的基因都拿出来计算,但是差异小的就忽略了) (测试通过,使用原始ENSEMBL数据)
ego <- enrichGO(gene          = names(geneList2_none),
                universe      = names(geneList2),
                keyType       = "ENSEMBL",
                OrgDb         = org.Mm.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = T)
gofilter(ego, level = 4)
write.table(ego,file="E:/5. Coding/R/R_Conductor/clusterprofiler/ego.txt")
head(ego)
#simplify(ego, cutoff = 0.7, by = "p.adjust",
#         select_fun = min, measure = "Wang", semData = NULL)
#barplot(ego, showCategory=8)
dotplot(ego)
enrichMap(ego)
cnetplot(ego, categorySize="pvalue", foldChange=geneList)
plotGOgraph(ego)

#****GO Gene Set Enrichment Analysis (弥补上述的缺陷,每一个基因都用来分析,因为有一些情况基因变化小但是有一致性)(测试通过,使用原始ENSEMBL数据>1)
ego2 <- gseGO(geneList     = geneList2_none,
              keyType      = "ENSEMBL",
              OrgDb        = org.Mm.eg.db,
              ont          = "BP",
              nPerm        = 1000,
              minGSSize    = 100,
              maxGSSize    = 500,
              pvalueCutoff = 0.05,
              verbose      = T)
head(ego2)
write.table(ego2,file="E:/5. Coding/R/R_Conductor/clusterprofiler/ego2.txt")
dotplot(ego2)
#*********************2. KEGG analysis*************************************
#####
search_kegg_organism('mmu', by='kegg_code')

#****KEGG over-representation test (测试通过,转化后ENTREZID数据>1)

kk <- enrichKEGG(gene         = gene,
                 organism     = 'mmu',
                 qvalueCutoff = 0.2,
                 universe      = names(geneList),
                 pvalueCutoff = 0.05)
head(kk)
write.table(kk,file="E:/5. Coding/R/R_Conductor/clusterprofiler/kk.txt")
browseKEGG(kk, 'mmu05012')
dotplot(kk)

#KEGG Module over-representation test
mkk <- enrichMKEGG(gene = gene,
                   organism = 'mmu')
write.table(mkk,file="E:/5. Coding/R/R_Conductor/clusterprofiler/mkk.txt")
dotplot(mkk)
#****KEGG Gene Set Enrichment Analysis (未通过:ENSEMBL is not supported for mmu ...)

kk2 <- gseKEGG(geneList     = geneList,
               
               organism     = 'mmu',
               nPerm        = 1000,
               minGSSize    = 120,
               pvalueCutoff = 0.05,
               verbose      = T)
head(kk2)
dotplot(kk2)
write.table(kk2,file="E:/5. Coding/R/R_Conductor/clusterprofiler/kk2.txt")
gseaplot(kk2, geneSetID = "mmu05012")
#KEGG Module Gene Set Enrichment Analysis
mkk2 <- gseMKEGG(geneList = geneList,
                 organism = 'mmu')
#####
#******************************4. Universal enrichment analysis*****************
#####
library("pathview")
hsa04110 <- pathview(gene.data  = geneList,
                     pathway.id = "mmu05012",
                     species    = "mmu",
                     limit      = list(gene=max(abs(geneList)), cpd=1))
#d导出文件
#####
#***************************5. Biological theme comparison**************
#****Formula interface of compareCluster
#****enrichGO
mydf <- data.frame(Entrez=names(geneList2), FC=geneList2)
mydf <- mydf[abs(mydf$FC) > 1.065,]
mydf$group <- "upregulated"
mydf$group[mydf$FC < 0] <- "downregulated"
mydf$othergroup <- "A (FC > 0)"
mydf$othergroup[abs(mydf$FC) > 1] <- "B (FC > 1)"
mydf$othergroup[abs(mydf$FC) > 2] <- "C (FC > 2)"
mydf$othergroup[abs(mydf$FC) > 3] <- "D (FC > 3)"
mydf$othergroup[abs(mydf$FC) > 4] <- "E (FC > 4)"
library(org.Mm.eg.db)
formula_res <- compareCluster(Entrez~group, 
                              data     = mydf, 
                              fun      ='enrichGO', 
                              ont = "BP",
                              universe      = names(geneList2),
                              keyType      = "ENSEMBL",
                              OrgDb    ='org.Mm.eg.db',
                              readable = T)
head(as.data.frame(formula_res))
dotplot(formula_res,title="GO enrichment analysis (BP) of differential genes (top20%), BvsCD4")
dotplot(formula_res, x=~group,title="GO enrichment analysis (BP) of differential genes, BvsCD4") + ggplot2::facet_grid(~othergroup)


#****enrichPathway
de <- names(geneList)[abs(geneList) > 1.5]
head(de)

x <- enrichPathway(gene=de,pvalueCutoff=0.05, qvalueCutoff = 0.05,readable=T, organism = "mouse")
head(as.data.frame(x))

dotplot(x, showCategory=15)


mydf <- data.frame(entrez=names(geneList), FC=geneList)
mydf <- mydf[abs(mydf$FC) > 1,]
mydf$group <- "upregulated"
mydf$group[mydf$FC < 0] <- "downregulated"
mydf$othergroup <- "A"
mydf$othergroup[abs(mydf$FC) > 2] <- "B"

res <- compareCluster(Entrez~group+othergroup, data=mydf, fun="enrichPathway",organism = "mouse")
plot(res)


mydf <- read.xlsx("confidence_interval_filtered_list.xlsx",4)
mydf <- read.xlsx("E:/Lu_RNA-seq/microarray/0531/microarray_diff_gene.xlsx",1)

formula_res <- compareCluster(gene~group, 
                              data     = mydf, 
                              fun      ='enrichGO', 
                              ont = "BP",
                              keyType      = "ENTREZID",
                              OrgDb    ='org.Mm.eg.db',
                              readable = T)
ggplot(formula_res,aes(Cluster, Description), showCategory=8) +
  geom_point(aes(color=qvalue, size=GeneRatio)) + scale_color_continuous(low='red', high='green')



p <- dotplot(formula_res, title="GO enrichment (BP) Microarray",showCategory=10)
p2 <- p + scale_color_continuous(low='red', high='green') 
p2



mydf <- read.xlsx("confidence_interval_filtered_list.xlsx",1)
gene.df <- bitr(mydf$MA9_3, fromType ="ENSEMBL",
                toType =c("ENTREZID", "SYMBOL"),
                OrgDb = org.Mm.eg.db)

ego <- enrichGO(gene          = mydf$gene,
                OrgDb         = org.Mm.eg.db,
                keyType       = "ENTREZID",
                ont           = "BP",
                pvalueCutoff  = 0.05,
                pAdjustMethod = "BH",
                qvalueCutoff  = 0.2,
                minGSSize     = 10, 
                maxGSSize     = 500, 
                readable      = T)
dotplot(ego,showCategory=50,title="MA9_3h ENTREZID",x='Count')

















