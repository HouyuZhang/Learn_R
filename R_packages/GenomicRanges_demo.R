source("https://bioconductor.org/biocLite.R")
biocLite("GenomicRanges")
library(GenomicRanges)


#******************1. GRanges: Genomic Ranges********************
gr <- GRanges(
  seqnames = Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
  ranges = IRanges(101:110, end = 111:120, names = head(letters, 10)),
  strand = Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)),
  score = 1:10,
  GC = seq(1, 0, length=10))
gr

#提取gr非注释部分
granges(gr)

#提取gr的染色体名字，范围，链信息
seqnames(gr)
ranges(gr)
range(gr)
strand(gr)

#提取注释信息为df
mcols(gr)
mcols(gr)$score

#可以对序列长度赋值，以人为例子
seqlengths(gr) <- c(249250621, 243199373, 198022430)
seqlengths(gr)
names(gr)
length(gr)

#**********1.1 Splitting and combining GRanges objects
sp <- split(gr, rep(1:2, each=5))
sp
#合并在一起
c(sp[[1]], sp[[2]])

#**********1.2 Subsetting GRanges objects
gr[2:3]
#只提取注释的GC列
gr[2:3, "GC"]

#对gr重新赋值
singles <- split(gr, names(gr))
grMod <- gr
grMod[2] <- singles[[1]]
head(grMod, n=3)

# 重复，颠倒，选定操作
rep(singles[[2]], times = 3)
rev(gr)
head(gr,2)
tail(gr,n=2)
window(gr, start=2,end=4)
gr[IRanges(start=c(2,7), end=c(3,9))]

#**********1.3 Basic interval operations for GRanges objects
#****Intra-range methods 
g <- gr[1:3]
g <- append(g, singles[[10]])
start(g)
end(g)
width(g)

# 向上游扩增
flank(g, 10)
#向上游一定一定碱基数
shift(g, 5)
#把范围变成30 碱基
resize(g, 30)

#****Inter-range methods 
#merge overlapping的区间
reduce(g)
#取间隙
gaps(g)
#把所有的有overlap的区间重新组织成为没有overlap
disjoin(g)

coverage(g)

#**********1.4 Interval set operations for GRanges objects
g2 <- head(gr, n=2)
union(g, g2)
intersect(g, g2)
setdiff(g, g2)

g3 <- g[1:2]
ranges(g3[1]) <- IRanges(start=105, end=112)
punion(g2, g3)
pintersect(g2, g3)
psetdiff(g2, g3)

methods(class="GRanges")

#******************2. GRangesList: Groups of Genomic Ranges********************














