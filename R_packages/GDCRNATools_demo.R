setwd("F:/Dropbox/1_Yenlab/1_8_Arch/Microarray/Microarry/GDCRNATools_demo")
library(GDCRNATools)
library(DT)
# ==============================================================================
# 1. Quick start
# ==============================================================================
data(rnaCounts)
data(mirCounts)
rnaExpr <- gdcVoomNormalization(counts = rnaCounts, filter = FALSE)
mirExpr <- gdcVoomNormalization(counts = mirCounts, filter = FALSE)

####### Parse and filter RNAseq metadata #######
metaMatrix.RNA <- gdcParseMetadata(project.id = 'TCGA-CHOL',
                                   data.type  = 'RNAseq', 
                                   write.meta = FALSE)

metaMatrix.RNA <- gdcFilterDuplicate(metaMatrix.RNA)
metaMatrix.RNA <- gdcFilterSampleType(metaMatrix.RNA)
datatable(as.data.frame(metaMatrix.RNA[1:5,]), extensions = 'Scroller',
          options = list(scrollX = TRUE, deferRender = TRUE, scroller = TRUE))


DEGAll <- gdcDEAnalysis(counts     = rnaCounts, 
                        group      = metaMatrix.RNA$sample_type, 
                        comparison = 'PrimaryTumor-SolidTissueNormal', 
                        method     = 'limma')
datatable(as.data.frame(DEGAll), 
          options = list(scrollX = TRUE, pageLength = 5))

### All DEGs
deALL <- gdcDEReport(deg = DEGAll, gene.type = 'all')

### DE long-noncoding
deLNC <- gdcDEReport(deg = DEGAll, gene.type = 'long_non_coding')

### DE protein coding genes
dePC <- gdcDEReport(deg = DEGAll, gene.type = 'protein_coding')

ceOutput <- gdcCEAnalysis(lnc         = rownames(deLNC), 
                          pc          = rownames(dePC), 
                          lnc.targets = 'starBase', 
                          pc.targets  = 'starBase', 
                          rna.expr    = rnaExpr, 
                          mir.expr    = mirExpr)
datatable(as.data.frame(ceOutput), 
          options = list(scrollX = TRUE, pageLength = 5))

ceOutput2 <- ceOutput[ceOutput$hyperPValue<0.01 
                      & ceOutput$corPValue<0.01 & ceOutput$regSim != 0,]
### Export edges
edges <- gdcExportNetwork(ceNetwork = ceOutput2, net = 'edges')
datatable(as.data.frame(edges), 
          options = list(scrollX = TRUE, pageLength = 5))

### Export nodes
nodes <- gdcExportNetwork(ceNetwork = ceOutput2, net = 'nodes')
datatable(as.data.frame(nodes), 
          options = list(scrollX = TRUE, pageLength = 5))

# ==============================================================================
# 2. Case study: TCGA-CHOL
# ==============================================================================
project <- 'TCGA-CHOL'
rnadir <- paste(project, 'RNAseq', sep='/')
mirdir <- paste(project, 'miRNAs', sep='/')
clinicaldir <- paste(project, 'Clinical', sep='/')

####### Download RNAseq data #######
gdcRNADownload(project.id     = 'TCGA-CHOL', 
               data.type      = 'RNAseq', 
               write.manifest = FALSE,
               method         = 'gdc-client',
               directory      = rnadir)

####### Download mature miRNA data #######
gdcRNADownload(project.id     = 'TCGA-CHOL', 
               data.type      = 'miRNAs', 
               write.manifest = FALSE,
               method         = 'gdc-client',
               directory      = mirdir)

####### Download clinical data #######
gdcClinicalDownload(project.id     = 'TCGA-CHOL', 
                    write.manifest = FALSE,
                    method         = 'gdc-client',
                    directory      = clinicaldir)

####### Parse RNAseq metadata #######
metaMatrix.RNA <- gdcParseMetadata(project.id = 'TCGA-CHOL', data.type  = 'RNAseq', write.meta = FALSE)

####### Filter duplicated samples in RNAseq metadata #######
metaMatrix.RNA <- gdcFilterDuplicate(metaMatrix.RNA)

####### Filter non-Primary Tumor and non-Solid Tissue Normal samples in RNAseq metadata #######
metaMatrix.RNA <- gdcFilterSampleType(metaMatrix.RNA)

####### Parse miRNAs metadata #######
metaMatrix.MIR <- gdcParseMetadata(project.id = 'TCGA-CHOL', data.type  = 'miRNAs', write.meta = FALSE)

####### Filter duplicated samples in miRNAs metadata #######
metaMatrix.MIR <- gdcFilterDuplicate(metaMatrix.MIR)

####### Filter non-Primary Tumor and non-Solid Tissue Normal samples in miRNAs metadata #######
metaMatrix.MIR <- gdcFilterSampleType(metaMatrix.MIR)

####### Merge RNAseq data #######
rnaCounts <- gdcRNAMerge(metadata  = metaMatrix.RNA, 
                         path      = rnadir, # the folder in which the data stored
                         organized = FALSE, # if the data are in separate folders
                         data.type = 'RNAseq')

####### Merge miRNAs data #######
mirCounts <- gdcRNAMerge(metadata  = metaMatrix.MIR,
                         path      = mirdir, # the folder in which the data stored
                         organized = FALSE, # if the data are in separate folders
                         data.type = 'miRNAs')

####### Merge clinical data #######
clinicalDa <- gdcClinicalMerge(path = clinicaldir, key.info = TRUE)
clinicalDa[1:6,5:10]

####### Normalization of RNAseq data #######
rnaExpr <- gdcVoomNormalization(counts = rnaCounts, filter = FALSE)

####### Normalization of miRNAs data #######
mirExpr <- gdcVoomNormalization(counts = mirCounts, filter = FALSE)

DEGAll <- gdcDEAnalysis(counts     = rnaCounts, 
                        group      = metaMatrix.RNA$sample_type, 
                        comparison = 'PrimaryTumor-SolidTissueNormal', 
                        method     = 'limma')

data(DEGAll)

### All DEGs
deALL <- gdcDEReport(deg = DEGAll, gene.type = 'all')

### DE long-noncoding
deLNC <- gdcDEReport(deg = DEGAll, gene.type = 'long_non_coding')

### DE protein coding genes
dePC <- gdcDEReport(deg = DEGAll, gene.type = 'protein_coding')

gdcVolcanoPlot(DEGAll)
gdcBarPlot(deg = deALL, angle = 45, data.type = 'RNAseq')
degName = rownames(deALL)
gdcHeatmap(deg.id = degName, metadata = metaMatrix.RNA, rna.expr = rnaExpr)

ceOutput <- gdcCEAnalysis(lnc         = rownames(deLNC), 
                          pc          = rownames(dePC), 
                          lnc.targets = 'starBase', 
                          pc.targets  = 'starBase', 
                          rna.expr    = rnaExpr, 
                          mir.expr    = mirExpr)

### load miRNA-lncRNA interactions
data(lncTarget)

### load miRNA-mRNA interactions
data(pcTarget)
pcTarget[1:3]

ceOutput <- gdcCEAnalysis(lnc         = rownames(deLNC), 
                          pc          = rownames(dePC), 
                          lnc.targets = lncTarget, 
                          pc.targets  = pcTarget, 
                          rna.expr    = rnaExpr, 
                          mir.expr    = mirExpr)

ceOutput2 <- ceOutput[ceOutput$hyperPValue<0.01 & 
                        ceOutput$corPValue<0.01 & ceOutput$regSim != 0,]

edges <- gdcExportNetwork(ceNetwork = ceOutput2, net = 'edges')
nodes <- gdcExportNetwork(ceNetwork = ceOutput2, net = 'nodes')

write.table(edges, file='edges.txt', sep='\t', quote=F)
write.table(nodes, file='nodes.txt', sep='\t', quote=F)

gdcCorPlot(gene1    = 'ENSG00000251165', 
           gene2    = 'ENSG00000091831', 
           rna.expr = rnaExpr, 
           metadata = metaMatrix.RNA)

shinyCorPlot(gene1    = rownames(deLNC), 
             gene2    = rownames(dePC), 
             rna.expr = rnaExpr, 
             metadata = metaMatrix.RNA)

####### CoxPH analysis #######
survOutput <- gdcSurvivalAnalysis(gene     = rownames(deALL), 
                                  method   = 'coxph', 
                                  rna.expr = rnaExpr, 
                                  metadata = metaMatrix.RNA)

####### KM analysis #######
survOutput <- gdcSurvivalAnalysis(gene     = rownames(deALL), 
                                  method   = 'KM', 
                                  rna.expr = rnaExpr, 
                                  metadata = metaMatrix.RNA, 
                                  sep      = 'median')

gdcKMPlot(gene     = 'ENSG00000136193',
          rna.expr = rnaExpr,
          metadata = metaMatrix.RNA,
          sep      = 'median')

shinyKMPlot(gene = rownames(deALL), rna.expr = rnaExpr, 
            metadata = metaMatrix.RNA)

enrichOutput <- gdcEnrichAnalysis(gene = rownames(deALL), simplify = TRUE)

data(enrichOutput)
gdcEnrichPlot(enrichOutput, type = 'bar', category = 'GO', num.terms = 10)

gdcEnrichPlot(enrichOutput, type='bubble', category='GO', num.terms = 10)
library(pathview)

deg <- deALL$logFC
names(deg) <- rownames(deALL)

pathways <- as.character(enrichOutput$Terms[enrichOutput$Category=='KEGG'])
pathways

shinyPathview(deg, pathways = pathways, directory = 'pathview')
