chrM <- read.table("hg38.50bp.chrM_cutfreqs.bed", sep = "\t", header = T)
colnames(chrM) <- gsub("Human_|_ATACseq|_chrM","",colnames(chrM))
scaled_chrM <- cbind(chrM[,1:3], scale(chrM[,4:ncol(chrM)]))

pdf("Human.pdf")
circos.clear()
circos.par(start.degree = 90)
circos.initializeWithIdeogram(species ="hg38", sort.chr = TRUE, chromosome.index="chrM",
                              plotType = c("labels", "axis"))

gene_bed = read.table("Mouse_chrM.gtf", sep = "\t", header = F)

circos.genomicTrack(gene_bed, ylim = c(0,0.2),
                    panel.fun = function(region, value, ...) {
                      circos.genomicRect(region, value, col = "red", border = "white", ...)
                    }, track.height = 0.05, bg.border = NA)

col_fun = colorRamp2(c(-3, 0, 3), rev(brewer.pal(n = 3, name = "RdBu")))
circos.genomicHeatmap(scaled_chrM, col = col_fun, side = "inside", border = NA, heatmap_height = 0.4,)

circos.clear()
dev.off()





