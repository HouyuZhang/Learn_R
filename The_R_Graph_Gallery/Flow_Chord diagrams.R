#https://www.r-graph-gallery.com/chord-diagram
library(tidyverse)
library(circlize)

# ==============================================================================
# 1. INTRODUCTION TO CIRCLIZE
# ==============================================================================
###1.1 Introduction to the circlize package ---- 
data = data.frame(
  factor = sample(letters[1:8], 1000, replace = TRUE),
  x = rnorm(1000), 
  y = runif(1000)
)

# Step1: Initialise the chart giving factor and x-axis.
circos.initialize( factors=data$factor, x=data$x )

# Step 2: Build the regions. 
circos.trackPlotRegion(factors = data$factor, y = data$y, panel.fun = function(x, y) {
  circos.axis()
})

# Step 3: Add points
circos.trackPoints(data$factor, data$x, data$y, col = "blue", pch = 16, cex = 0.5) 

###1.2 Circular chart customization for the circlize R package ----
par(
  mar = c(1, 1, 1, 1),           # Margin around chart
  bg = rgb(0.4,0.1,0.7,0.05)     # background color
) 
circos.par("track.height" = 0.6) # track hight, 0.6 = 60% of total height

circos.initialize( factors=data$factor, x=data$x )
circos.trackPlotRegion(factors = data$factor, y = data$y, panel.fun = function(x, y) {
  circos.axis(
    h="top",                   # x axis on the inner or outer part of the track?
    labels=TRUE,               # show the labels of the axis?
    major.tick=TRUE,           # show ticks?
    labels.cex=0.5,            # labels size (higher=bigger)
    labels.font=1,             # labels font (1, 2, 3 , 4)
    direction="outside",       # ticks point to the outside or inside of the circle ?
    minor.ticks=4,             # Number of minor (=small) ticks
    major.tick.percentage=0.1, # The size of the ticks in percentage of the track height
    lwd=2                      # thickness of ticks and x axis.
  )
})
circos.trackPoints(data$factor, data$x, data$y, col = "#69b3a2", pch = 16, cex = 0.5)

###1.3 Available circular chart types with circlize ----
circos.par("track.height" = 0.4)

#Circular Line chart
circos.initialize(factors=data$factor, x=data$x)
circos.trackPlotRegion(factors = data$factor, y = data$y, panel.fun = function(x, y) {circos.axis()})
circos.trackLines(data$factor, data$x[order(data$x)], data$y[order(data$x)], 
                  col = rgb(0.1,0.5,0.8,0.3), lwd=2)

#Vertical ablines
circos.initialize( factors=data$factor, x=data$x )
circos.trackPlotRegion(factors = data$factor, y = data$y, panel.fun = function(x, y) {circos.axis()})
circos.trackLines(data$factor, data$x[order(data$x)], data$y[order(data$x)], 
                  col = rgb(0.1,0.5,0.8,0.3), lwd=2, type="h")

#Circular histogram
circos.initialize( factors=data$factor, x=data$x )
circos.trackHist(data$factor, data$x, bg.col = "white", col = "#69b3a2")

###1.4 Circular chart with several tracks ----
circos.clear()

par(mar = c(1, 1, 1, 1) ) 
circos.initialize(factors = data$factor, x = data$x )

# Build the regions of track #1 
circos.trackPlotRegion(factors = data$factor, y=data$y, panel.fun = function(x, y) {
  circos.axis(labels.cex=0.5, labels.font=1, lwd=0.8)})
# --> Add a scatterplot on it:
circos.trackPoints(data$factor, data$x, data$y, col = rgb(0.1,0.5,0.8,0.3), pch=20)

# Build the regions of track #2:
circlize::circos.trackPlotRegion(factors = data$factor, y=data$y, panel.fun = function(x, y) {
  circos.axis(labels=FALSE, major.tick=FALSE)})
# --> Add a scatterplot on it
circos.trackPoints(data$factor, data$x, data$y, col = rgb(0.9,0.5,0.8,0.3), pch=20, cex=2)

# Add track #3 --> don't forget you can custom the height of tracks!
circos.par("track.height" = 0.4)
circos.trackPlotRegion(factors = data$factor, y=data$y, panel.fun = function(x, y) {
  circos.axis(labels=FALSE, major.tick=FALSE)})
circos.trackLines(data$factor, data$x, data$y, col = rgb(0.9,0.5,0.1,0.3), pch=20, cex=2, type="h")
# and continue as long as needed!

###1.5 Display several chord diagrams on same figure ----
# Arrange the layout
layout(matrix(1:9, 3, 3)) 

# A loop to create 9 circular plots
for(i in 1:9) {
  par(mar = c(0.5, 0.5, 0.5, 0.5), bg=rgb(1,1,1,0.1) )
  factors = 1:8
  circos.par(cell.padding = c(0, 0, 0, 0)) 
  circos.initialize(factors, xlim = c(0, 1)) 
  circos.trackPlotRegion(ylim = c(0, 1), track.height = 0.05, bg.col = rand_color(8), bg.border = NA ) 
  
  # add links
  for(i in 1:20) {
    se = sample(1:8, 2)
    circos.link(se[1], runif(2), se[2], runif(2), col = rand_color(1, transparency = 0.4))}
  circos.clear()
}

###1.6 Draw part of the circular chart only
# Create data
factors <- letters[1:4]
x1 <- runif(100)
y1 <- runif(100)

# general parameter of the plot. 
# With canvas.xlim and canvas.ylim we kind of "zoom on a part of the plot:
par(mar = c(1, 2, 0.1, 0.1) )
circos.par("track.height" = 0.7, "canvas.xlim" = c(0, 1), 
           "canvas.ylim" = c(0, 1), "gap.degree" = 0, "clock.wise" = FALSE)

# Make the usual plot, but with no border
circos.initialize(factors = factors, xlim = c(0, 1)) 
circos.trackPlotRegion(factors = factors, ylim = c(0, 1), bg.border = NA ) 

# Finally we plot only the firs sector, so let's change its border to "black" to see it
circos.updatePlotRegion(sector.index = "a", bg.border = "grey" , bg.lwd=0.2)

# Now we can add a plot in this section! You can repeat these steps to add several regions
circos.lines(x1, y1, pch = 16, cex = 0.5, type="h" , col="#69b3a2" , lwd=3)

# Add axis
circos.axis(h="bottom" , labels.cex=0.4, direction = "inside" )

circos.clear()
# ==============================================================================
# 2. STATIC CHORD DIAGRAM
# ==============================================================================
###2.1 Introduction to chord diagram ----
set.seed(123)
data = data.frame(
  factor = sample(letters[1:8], 1000, replace = TRUE),
  x = rnorm(1000), 
  y = runif(1000)
)

# Initialize the plot.
par(mar = c(1, 1, 1, 1) ) 
circos.initialize(factors = data$factor, x = data$x )

# Build the regions of track #1
circos.trackPlotRegion(factors = data$factor, y=data$y , bg.col = rgb(0.1,0.1,seq(0,1,0.1),0.4) , bg.border = NA)

# Add a link between a point and another
circos.link("a", 0, "b", 0, h = 0.4)

# Add a link between a point and a zone
circos.link("e", 0, "g", c(-1,1), col = "green", lwd = 2, lty = 2, border="black" )

# Add a link between a zone and another
circos.link("c", c(-0.5, 0.5), "d", c(-0.5,0.5), col = "red", border = "blue", h = 0.2)


###2.2 Advanced chord diagram with R and circlize ----
#Chord diagram from adjacency matrix

library(circlize)
library(migest)
library(dplyr)

m <- data.frame(order = 1:6,
                country = c("Ausralia", "India", "China", "Japan", "Thailand", "Malaysia"),
                V3 = c(1, 150000, 90000, 180000, 15000, 10000),
                V4 = c(35000, 1, 10000, 12000, 25000, 8000),
                V5 = c(10000, 7000, 1, 40000, 5000, 4000),
                V6 = c(7000, 8000, 175000, 1, 11000, 18000),
                V7 = c(70000, 30000, 22000, 120000, 1, 40000),
                V8 = c(60000, 90000, 110000, 14000, 30000, 1),
                r = c(255,255,255,153,51,51),
                g = c(51, 153, 255, 255, 255, 255),
                b = c(51, 51, 51, 51, 51, 153),
                stringsAsFactors = FALSE)
df1 <- m[, c(1,2, 9:11)]
m <- m[,-(1:2)]/1e04
m <- as.matrix(m[,c(1:6)])
dimnames(m) <- list(orig = df1$country, dest = df1$country)
#Sort order of data.frame and matrix for plotting in circos
df1 <- arrange(df1, order)
df1$country <- factor(df1$country, levels = df1$country)
m <- m[levels(df1$country),levels(df1$country)]


### Define ranges of circos sectors and their colors (both of the sectors and the links)
df1$xmin <- 0
df1$xmax <- rowSums(m) + colSums(m)
n <- nrow(df1)
df1$rcol<-rgb(df1$r, df1$g, df1$b, max = 255)
df1$lcol<-rgb(df1$r, df1$g, df1$b, alpha=200, max = 255)

### Plot sectors (outer part)
par(mar=rep(0,4))
circos.clear()

### Basic circos graphic parameters
circos.par(cell.padding=c(0,0,0,0), track.margin=c(0,0.15), start.degree = 90, gap.degree =4)

### Sector details
circos.initialize(factors = df1$country, xlim = cbind(df1$xmin, df1$xmax))

### Plot sectors
circos.trackPlotRegion(ylim = c(0, 1), factors = df1$country, track.height=0.1,
                       #panel.fun for each sector
                       panel.fun = function(x, y) {
                         #select details of current sector
                         name = get.cell.meta.data("sector.index")
                         i = get.cell.meta.data("sector.numeric.index")
                         xlim = get.cell.meta.data("xlim")
                         ylim = get.cell.meta.data("ylim")
                         
                         #text direction (dd) and adjusmtents (aa)
                         theta = circlize(mean(xlim), 1.3)[1, 1] %% 360
                         dd <- ifelse(theta < 90 || theta > 270, "clockwise", "reverse.clockwise")
                         aa = c(1, 0.5)
                         if(theta < 90 || theta > 270)  aa = c(0, 0.5)
                         
                         #plot country labels
                         circos.text(x=mean(xlim), y=1.7, labels=name, facing = dd, cex=0.6,  adj = aa)
                         
                         #plot main sector
                         circos.rect(xleft=xlim[1], ybottom=ylim[1], xright=xlim[2], ytop=ylim[2], 
                                     col = df1$rcol[i], border=df1$rcol[i])
                         
                         #blank in part of main sector
                         circos.rect(xleft=xlim[1], ybottom=ylim[1], xright=xlim[2]-rowSums(m)[i], ytop=ylim[1]+0.3, 
                                     col = "white", border = "white")
                         
                         #white line all the way around
                         circos.rect(xleft=xlim[1], ybottom=0.3, xright=xlim[2], ytop=0.32, col = "white", border = "white")
                         
                         #plot axis
                         circos.axis(labels.cex=0.6, direction = "outside", major.at=seq(from=0,to=floor(df1$xmax)[i],by=5), 
                                     minor.ticks=1, labels.away.percentage = 0.15)
                       })

### Plot links (inner part)
### Add sum values to df1, marking the x-position of the first links
### out (sum1) and in (sum2). Updated for further links in loop below.
df1$sum1 <- colSums(m)
df1$sum2 <- numeric(n)

### Create a data.frame of the flow matrix sorted by flow size, to allow largest flow plotted first
df2 <- cbind(as.data.frame(m),orig=rownames(m),  stringsAsFactors=FALSE)
df2 <- reshape(df2, idvar="orig", varying=list(1:n), direction="long",
               timevar="dest", time=rownames(m),  v.names = "m")
df2 <- arrange(df2,desc(m))

### Keep only the largest flows to avoid clutter
df2 <- subset(df2, m > quantile(m,0.6))

### Plot links
for(k in 1:nrow(df2)){
  #i,j reference of flow matrix
  i<-match(df2$orig[k],df1$country)
  j<-match(df2$dest[k],df1$country)
  
  #plot link
  circos.link(sector.index1=df1$country[i], point1=c(df1$sum1[i], df1$sum1[i] + abs(m[i, j])),
              sector.index2=df1$country[j], point2=c(df1$sum2[j], df1$sum2[j] + abs(m[i, j])),
              col = df1$lcol[i])
  
  #update sum1 and sum2 for use when plotting the next link
  df1$sum1[i] = df1$sum1[i] + abs(m[i, j])
  df1$sum2[j] = df1$sum2[j] + abs(m[i, j])
}

# Create an adjacency matrix: 
# a list of connections between 20 origin nodes, and 5 destination nodes:
numbers <- sample(c(1:1000), 100, replace = T)
data <- matrix( numbers, ncol=5)
rownames(data) <- paste0("orig-", seq(1,20))
colnames(data) <- paste0("dest-", seq(1,5))

# Make the circular plot
chordDiagram(data, transparency = 0.5)

# Create an edge list: a list of connections between 10 origin nodes, and 10 destination nodes:
origin <- paste0("orig ", sample(c(1:10), 20, replace = T))
destination <- paste0("dest ", sample(c(1:10), 20, replace = T))
data <- data.frame(origin, destination)

# Transform input data in a adjacency matrix
adjacencyData <- with(data, table(origin, destination))

# Make the circular plot
chordDiagram(adjacencyData, transparency = 0.5)

###2.3 
library(viridis)
library(patchwork) #devtools::install_github("thomasp85/patchwork")
library(hrbrthemes)
library(chorddiag)  #devtools::install_github("mattflor/chorddiag")

# Load dataset from github
data <- read.table("https://raw.githubusercontent.com/holtzy/data_to_viz/master/Example_dataset/13_AdjacencyDirectedWeighted.csv", header=TRUE)

# short names
colnames(data) <- c("Africa", "East Asia", "Europe", "Latin Ame.",   "North Ame.",   "Oceania", "South Asia", "South East Asia", "Soviet Union", "West.Asia")
rownames(data) <- colnames(data)

# I need a long format
data_long <- data %>%
  rownames_to_column %>%
  gather(key = 'key', value = 'value', -rowname)

# parameters
circos.clear()
circos.par(start.degree = 90, gap.degree = 4, track.margin = c(-0.1, 0.1), points.overflow.warning = FALSE)
par(mar = rep(0, 4))

# color palette
mycolor <- viridis(10, alpha = 1, begin = 0, end = 1, option = "D")
mycolor <- mycolor[sample(1:10)]

# Base plot
chordDiagram(
  x = data_long, 
  grid.col = mycolor,
  transparency = 0.25,
  directional = 1,
  direction.type = c("arrows", "diffHeight"), 
  diffHeight  = -0.04,
  annotationTrack = "grid", 
  annotationTrackHeight = c(0.05, 0.1),
  link.arr.type = "big.arrow", 
  link.sort = TRUE, 
  link.largest.ontop = TRUE)

# Add text and axis
circos.trackPlotRegion(
  track.index = 1, 
  bg.border = NA, 
  panel.fun = function(x, y) {
    
    xlim = get.cell.meta.data("xlim")
    sector.index = get.cell.meta.data("sector.index")
    
    # Add names to the sector. 
    circos.text(
      x = mean(xlim), 
      y = 3.2, 
      labels = sector.index, 
      facing = "bending", 
      cex = 0.8
    )
    
    # Add graduation on axis
    circos.axis(
      h = "top", 
      major.at = seq(from = 0, to = xlim[2], by = ifelse(test = xlim[2]>10, yes = 2, no = 1)), 
      minor.ticks = 1, 
      major.tick.percentage = 0.5,
      labels.niceFacing = FALSE)
  }
)

# ==============================================================================
# 3. INTERACTIVITY
# ==============================================================================
library(chorddiag)

# Create dummy data
m <- matrix(c(11975,  5871, 8916, 2868,
              1951, 10048, 2060, 6171,
              8010, 16145, 8090, 8045,
              1013,   990,  940, 6907),
            byrow = TRUE,
            nrow = 4, ncol = 4)

# A vector of 4 colors for 4 groups
haircolors <- c("black", "blonde", "brown", "red")
dimnames(m) <- list(have = haircolors,
                    prefer = haircolors)
groupColors <- c("#000000", "#FFDD89", "#957244", "#F26223")

# Build the chord diagram:
p <- chorddiag(m, groupColors = groupColors, groupnamePadding = 20)
p

# save the widget
# library(htmlwidgets)
# saveWidget(p, file=paste0( getwd(), "/HtmlWidget/chord_interactive.html"))
