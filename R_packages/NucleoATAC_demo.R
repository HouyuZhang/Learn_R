library(NucleoATACR)

#Read in nucleosome and nfr positions:
nucs <- readNucs("test.nucmap_combined.bed.gz")
nfrs <- readNFRs("test.nfrpos.bed.gz")

#Read in nucleoatac signal track for particular locus:
signal <- readBedgraph("test.nucleoatac_signal.bedgraph.gz", "chr1", 70006551, 70770500)

#Read in vplot and plot:
v <- read_vplot("test.VMat")
plotV(v)
#Get +1 and -1 nucleosomes:
fake_tss <- GenomicRanges::GRanges("chr1", IRanges::IRanges(start = c(707200,707500), width = 1), 
                                     strand = c("+","-")) 
p1 <- get_p1_nuc(nucs.ranges = nucs, sites = fake_tss)
m1 <- get_m1_nuc(nucs.ranges = nucs, sites = fake_tss)
