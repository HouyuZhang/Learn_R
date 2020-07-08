#=========================================================================================
# This file contains common used functions in Tn5 project
#=========================================================================================
library(easypackages)
libraries("ggthemes", "Biostrings",
          "Logolas","tidyverse",#"DNAshapeR","fields",
          "RColorBrewer","ComplexHeatmap","circlize")

nuc_color <- c('#109648','#255C99', '#D62839','#F7B32B') #ACTG

YlOrRd_colors <- colorRampPalette(brewer.pal(n = 8, name = "YlOrRd"))(9)
Dark2_colors <- colorRampPalette(brewer.pal(n = 8, name = "Dark2"))(8)
#"#1B9E77","#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666"
RdBu_colors <- colorRampPalette(brewer.pal(n = 8, name = "RdBu"))(11)

Inter_bp_shape <- c("Rise","Shift","Slide","HelT","Roll","Tilt") #long_shape
Intra_bp_shape <- c("ProT","Stretch","Buckle","Shear","Opening","Stagger") #short_shape
groove_shape <- c("MGW","EP")
All_shapes <- c(Inter_bp_shape,Intra_bp_shape,groove_shape)
#=========================================================================================
#1. Miscellaneous functions for control R envs
#=========================================================================================
#This function list all variables in R env, ranked by Size
.ls.objects <- function (pos = 1, pattern, order.by = "Size", decreasing=TRUE, head = TRUE, n = 10) {
  napply <- function(names, fn) sapply(names, function(x) fn(get(x, pos = pos)))
  names <- ls(pos = pos, pattern = pattern)
  obj.class <- napply(names, function(x) as.character(class(x))[1])
  obj.mode <- napply(names, mode)
  obj.type <- ifelse(is.na(obj.class), obj.mode, obj.class)
  obj.size <- napply(names, object.size) / 10^6 # megabytes
  obj.dim <- t(napply(names, function(x) as.numeric(dim(x))[1:2]))
  vec <- is.na(obj.dim)[, 1] & (obj.type != "function")
  obj.dim[vec, 1] <- napply(names, length)[vec]
  out <- data.frame(obj.type, obj.size, obj.dim)
  names(out) <- c("Type", "Size", "Rows", "Columns")
  out <- out[order(out[[order.by]], decreasing=decreasing), ]
  if (head)
    out <- head(out, n)
  out
}

#=========================================================================================
#2. This function is for plotting nucleotide frequency given frequency table
#=========================================================================================
#! This function need to reform for ggplot
#plot_Nuc(file, file_ran, flank_window = 50, norm = F, plot_type = "all")
plot_Nuc_freqs_around <- function(file, file_ran, flank_window = 5, norm = F, plot_type = "all"){
  
  #get plotting region from given flank_window 
  freq <- read.table(file)
  prefix <- sub(".fa_nuc.txt","",file)
  n <- ncol(freq)
  st <- (n+1)/2 - flank_window
  ed <- (n+1)/2 + flank_window
  freq <- freq[,c(st:ed)]
  
  #get plotting region from given flank_window (random file)
  freq_ran <- read.table(file_ran)
  freq_ran <- cbind(freq_ran,freq_ran,freq_ran[,1])
  n_ran <- ncol(freq_ran)
  st_ran <- (n_ran-1)/2 - flank_window
  ed_ran <- (n_ran-1)/2 + flank_window
  freq_ran <- freq_ran[,c(st_ran:ed_ran)]
  
  #Normalize ACTG unbanlance due to N (This solved in "Fasta_GC_content.R")
  if (norm){
    for (i in seq(1:ncol(freq))){
      sf <- 1/sum(freq[,i])
      freq[,i]<- freq[,i]*sf
    }
  }
  
  n <- ncol(freq)
  rownames(freq) <- c("A","C","T","G")
  colnames(freq) <- paste0(seq(-(n-1)/2,(n-1)/2))
  
  lw = 1.5
  if (flank_window < 100){lw = 2}
  else if (flank_window > 1000){lw = 1}
  
  ###plot Nuc 
  if (plot_type =="all"){
    
    pdf(paste0(prefix,"_f",flank_window,"bp_Nuc.pdf"))
    plot(as.numeric(freq[1,]),type = "l",col = nuc_color[1],ylim = c(0,0.5),lwd=lw,
         xlab = "Position releative to Tn5 cut sites (bp)",
         ylab = "Nucleotide frequency",
         xaxt="n", cex.axis=1.2, cex.lab = 1.5)
    
    abline(v=(n-1)/2+1,lwd=2,lty=3,col="grey")
    axis(1,at=seq(1,n),labels = seq(-(n-1)/2,(n-1)/2),cex.axis=1)
    
    lines(as.numeric(freq[4,]),type = "l",col = nuc_color[4],lwd=lw)
    lines(as.numeric(freq[2,]),type = "l",col = nuc_color[2],lwd=lw)
    lines(as.numeric(freq[3,]),type = "l",col = nuc_color[3],lwd=lw)
    
    lines(as.numeric(freq_ran[4,]),type = "l",lty=2,col = alpha(nuc_color[4], 0.3),lwd=lw)
    lines(as.numeric(freq_ran[1,]),type = "l",lty=2,col = alpha(nuc_color[1], 0.3),lwd=lw)
    lines(as.numeric(freq_ran[2,]),type = "l",lty=2,col = alpha(nuc_color[2], 0.3),lwd=lw)
    lines(as.numeric(freq_ran[3,]),type = "l",lty=2,col = alpha(nuc_color[3], 0.3),lwd=lw)
    
    legend(1, 0.5, legend = c("A","C","T","G"),col = nuc_color,lty = "solid",cex = 1.1,lwd=lw)
    dev.off()
  }
  
  ###Plot GC
  pdf(paste0(prefix,"_f",flank_window,"bp_GC.pdf"))
  plot(as.numeric(freq[1,])+as.numeric(freq[3,]),type = "l",col = nuc_color[3], ylim = c(0.2,0.8),lwd=lw,
       xlab = "Position releative to Tn5 cut sites (bp)",
       ylab = "Nucleotide frequency",
       xaxt="n",cex.axis=1.2,cex.lab = 1.5)
  abline(v=(n-1)/2+1,lwd=2,lty=3,col="grey")
  axis(1,at=seq(1,n),labels = seq(-(n-1)/2,(n-1)/2),cex.axis=1)
  lines(as.numeric(freq[2,])+as.numeric(freq[4,]),type = "l",col = nuc_color[2], lwd=lw)
  
  lines(as.numeric(freq_ran[1,])+as.numeric(freq_ran[3,]),type = "l",lty=2,col = alpha(nuc_color[3], 0.3),lwd=lw)
  lines(as.numeric(freq_ran[2,])+as.numeric(freq_ran[4,]),type = "l",lty=2,col = alpha(nuc_color[2], 0.3),lwd=lw)
  legend(1, 0.8, legend = c("random AT","random CG","AT","CG"),
         col = rep(nuc_color[3:2],2) ,lty = c(2,2,1,1),cex = 1.1,lwd=lw)
  dev.off()
}

#=========================================================================================
#3. This function is for plotting MEME result - using exported 'Probability Matrix' file
#=========================================================================================
#plot_MEME_result(file,flank_window = 5)
plot_MEME_PSSM <- function(file, flank_window = 5){
  freqs <- read.table(file)
  prefix <- sub(".txt","",file)
  
  colnames(freqs) <- c("A","C","G","T")
  freqs <- t(freqs)
  
  st <- (ncol(freqs) + 1)/2 - flank_window
  ed <- (ncol(freqs) + 1)/2 + flank_window
  
  freqs <- freqs[,c(st:ed)]
  colnames(freqs) <- seq(-flank_window,flank_window)
  
  #GetConsensusSeq(freqs)
  pdf(paste0(prefix,"_MEME_motif_bits.pdf"),height = 4,width = 8)
  logomaker(freqs, type = "Logo", color_type = "per_row",
            logo_control = list(main_fontsize=20, yscale_change=F,
                                xaxis_fontsize=10, xlab_fontsize=15, y_fontsize=15,
                                xlab = "Position relative to Tn5 cut sites (bp)",ylab = "Bits"),
            colors = c(nuc_color[c(4,2,3,1)]))
  dev.off()
  
  # pdf(paste0(prefix,"_MEME_motif_prob.pdf"),height = 4,width = 8)
  # logomaker(freqs, type = "Logo", color_type = "per_row",
  #           logo_control = list(main_fontsize=20, ic.scale=F,
  #                               xaxis_fontsize=10, xlab_fontsize=15, y_fontsize=15,
  #                               xlab = "Position relative to Tn5 cut sites (bp)",ylab = "Probability"),
  #           colors = c(nuc_color[c(4,2,3,1)]))
  # dev.off()
}

#plot_MEME_result_txt(meme_file, species="Mouse")
plot_MEME_result_txt <- function(meme_file, species=NULL){
  #Parse meme.txt result
  prefix <- dirname(meme_file)
  bn <- basename(dirname(meme_file))
  
  motif <- readLines(meme_file)
  
  num <- grep("MEME-1 position-specific probability matrix",motif)
  w <- as.integer(str_extract(str_extract(motif[(num+2)],"w= [0-9]+"),"[0-9]+"))
  
  if (w < 19){
    
    num <- grep("MEME-2 position-specific probability matrix",motif)
    w <- as.integer(str_extract(str_extract(motif[(num+2)],"w= [0-9]+"),"[0-9]+"))
  }
  
  if (w < 19){
    num <- grep("MEME-3 position-specific probability matrix",motif)
    w <- as.integer(str_extract(str_extract(motif[(num+2)],"w= [0-9]+"),"[0-9]+"))
  }
  
  Tn5_motif <- motif[(num+3):(num+2+w)]
  
  freqs <- matrix(0,nrow = w,ncol = 4)
  for (line in seq(1:length(Tn5_motif))){
    x <- Tn5_motif[line]
    freqs[line,] <- matrix(scan(text = x,quiet=T),nrow = 1,byrow = TRUE)[1,]
  }
  
  write.table(freqs,paste0(bn,"_Tn5_motif_1_freqs.txt"),quote = F,sep = "\t",col.names = F, row.names = F)
  
  colnames(freqs) <- c("A","C","G","T")
  freqs <- t(freqs)
  flank <- w - (ncol(freqs) + 1)/2
  colnames(freqs) <- seq(-flank,flank)
  
  nuc_color <- c('#109648','#255C99', '#D62839','#F7B32B') #ACTG
  
  if (!missing(species)){
    if (species == "Fish"){bg <- c(0.317, 0.183, 0.183, 0.317)}
    else if (species == "Fly"){bg <- c(0.29, 0.21, 0.21, 0.29)}
    else if (species == "huffia"){bg <- c(0.403, 0.097, 0.097, 0.403)}
    else if (species == "Human"){bg <- c(0.296, 0.204, 0.204, 0.296)}
    else if (species == "maize"){bg <- c(0.266, 0.234, 0.234, 0.266)}
    else if (species == "Mouse"){bg <- c(0.292, 0.208, 0.208, 0.292)}
    else if (species == "Plant"){bg <- c(0.32, 0.18, 0.18, 0.32)}
    else if (species == "Worm"){bg <- c(0.323, 0.177, 0.177, 0.323)}
  } else {bg <- c(0.25, 0.25, 0.25, 0.25)}
  names(bg) <- c("A", "C", "G", "T")
  
  #GetConsensusSeq(freqs)
  pdf(paste0(bn,"_MEME_motif_1_bits.pdf"), height = 4,width = 8)
  logomaker(freqs, type = "Logo", color_type = "per_row", bg=bg,
            logo_control = list(main_fontsize=20, yscale_change=F,
                                xaxis_fontsize=10, xlab_fontsize=15, y_fontsize=15,
                                xlab = "Position relative to Tn5 cut sites (bp)", ylab = "Bits"),
            colors = c(nuc_color[c(4,2,3,1)]))
  dev.off()
  
  # pdf(paste0(bn,"_MEME_motif_1_prob.pdf"),height = 4,width = 8)
  # logomaker(freqs, type = "Logo", color_type = "per_row", bg=bg,
  #           logo_control = list(main_fontsize=20, ic.scale=F,
  #                               xaxis_fontsize=10, xlab_fontsize=15, y_fontsize=15,
  #                               xlab = "Position relative to Tn5 cut sites (bp)", ylab = "Probability"),
  #           colors = c(nuc_color[c(4,2,3,1)]))
  # dev.off()
}

#=========================================================================================
#4. This function is for plotting DNA shape around
#=========================================================================================
#Plot_DNA_shape(path, flank_window = 50)
Plot_DNA_shape <- function(path, flank_window = 50){
  files <- list.files(path, ".rds")
  
  for (file in files){
    #file <- files[1]
    cat("Processing",file,"...\n")
    tmp <- readRDS(paste0(path,"/",file))
    shapes <- names(tmp)
    
    NUMC <- ncol(tmp[[shapes[1]]])
    mean_shape <- matrix(0, nrow = length(shapes), ncol = NUMC) %>% as.data.frame()
    rownames(mean_shape) <- shapes
    
    for (shape in shapes){
      #shape <- shapes[1]
      if (shape %in% Inter_bp_shape){
        mean_shape[shape,1:(NUMC-1)] <- colMeans(tmp[[shape]], na.rm = T)
      } else {
        mean_shape[shape,] <- colMeans(tmp[[shape]], na.rm = T)
      }
    }
    rownames(mean_shape) <- paste0(rownames(mean_shape), "_" ,file)
    assign(paste0(file,"_mean_shape"), mean_shape)
  }
  
  v <- paste0(files,"_mean_shape")
  final_matrix <- do.call(rbind,mget(v))
  
  #Plot
  n <- (ncol(final_matrix)-1)/2 + 1
  final_matrix <- t(final_matrix[,(n-flank_window):(n+flank_window)]) %>% as.data.frame()
  final_matrix$position <- -flank_window:flank_window
  
  final_matrix_melt <- reshape2::melt(final_matrix,c("position"))
  
  final_matrix_melt$variable <- gsub(".*_mean_shape.|_ext.*","",final_matrix_melt$variable)
  final_matrix_melt$shape <- gsub("_.*","",final_matrix_melt$variable)
  final_matrix_melt$sample <- sub(".*?_","",final_matrix_melt$variable)
  
  pdf(paste0(path,"/DNA shapes around +-", flank_window,"bp.pdf"), height = 9, width = 14)
  
  p <- ggplot(final_matrix_melt, aes(position, value, color = sample)) + 
    #geom_point(size = 0.1, alpha = 0.4) +
    geom_smooth(method="loess", se=T, fullrange=F, level=0.95) + 
    facet_wrap(~shape, ncol = 4, scales = "free_y") +
    theme_bw() +
    scale_fill_brewer(palette="Dark2") +
    labs(x="Relative position to dyad (bp)", y = "") +
    theme(
      plot.title = element_text(color="black", size=20, face="bold"),
      axis.title.x = element_text(color="black", size=14, face="bold"),
      axis.text.x = element_text(color="black", size=11, face="bold"),
      axis.title.y = element_text(color="black", size=14, face="bold"),
      axis.text.y = element_text(color="black", size=11, face="bold"),
      legend.title = element_text(color="black", size=14, face="bold"), 
      legend.text = element_text(color="black", size=12, face="bold"),
      strip.text.x = element_text(size = 14, face="bold")
    ) 
  plot(p)
  dev.off()
}

#=========================================================================================
#5. This function is for calculating relative importance of given predictors
#=========================================================================================
relweights <- function(data,top=10){
  R <- cor(data)
  which(is.na(colSums(R)))
  nvar <- ncol(R)
  rxx <- R[2:nvar, 2:nvar]
  rxy <- R[2:nvar, 1]
  svd <- eigen(rxx)
  evec <- svd$vectors
  ev <- svd$values + 1e-10
  delta <- diag(sqrt(ev))
  lambda <- evec %*% delta %*% t(evec)
  lambdasq <- lambda ^ 2
  beta <- solve(lambda) %*% rxy
  rsquare <- colSums(beta ^ 2)
  rawwgt <- lambdasq %*% beta ^ 2
  import <- (rawwgt / rsquare) * 100
  import <- as.data.frame(import)
  row.names(import) <- colnames(data)[-1]
  names(import) <- "Weights"
  import <- import[order(import),1, drop=FALSE]
  
  import$predictors <- rownames(import)
  import$predictors <- factor(import$predictors, levels = import$predictors[order(import$Weights)])
  
  import <- import[(nrow(import)-top+1):nrow(import),]
  p <- ggplot(import,aes(predictors,Weights)) + 
    geom_bar(stat = "identity") +
    theme_bw() +
    coord_flip() +
    ggtitle("Relative Importance of Predictors") +
    ylab("Relative importance") +
    theme(
      plot.title = element_text(color="black", size=16, face="bold"),
      axis.title.x = element_text(size=14, face="bold"),
      axis.text.x = element_text(size=12,face="bold"),
      axis.text.y = element_text(size=12,face="bold"),
      axis.title.y = element_text(size=14, face="bold")
    )
  print(p)
  return(import)
}

varImp_my <- function(object, lambda = NULL, ...) {
  beta <- predict(object, s = lambda, type = "coef")
  if(is.list(beta)) {
    out <- do.call("cbind", lapply(beta, function(x) x[,1]))
    out <- as.data.frame(out)
  } else out <- data.frame(coefficients = beta[,1])
  out <- abs(out[rownames(out) != "(Intercept)",,drop = FALSE])
  out$predictors <- rownames(out)
  out <- out[order(out$coefficients,decreasing = T),]
  out$predictors <- factor(out$predictors, levels = out$predictors[order(out$coefficients)])
  out
}

varImp_noabs <- function(object, lambda = NULL, ...) {
  beta <- predict(object, s = lambda, type = "coef")
  if(is.list(beta)) {
    out <- do.call("cbind", lapply(beta, function(x) x[,1]))
    out <- as.data.frame(out)
  } else out <- data.frame(coefficients = beta[,1])
  out <- out[rownames(out) != "(Intercept)",,drop = FALSE]
  out$predictors <- rownames(out)
  out <- out[order(out$coefficients,decreasing = T),]
  out <- out %>% arrange(abs(coefficients))
  out$predictors <- factor(out$predictors, levels = out$predictors)
  #out$predictors <- factor(out$predictors, levels = out$predictors[order(out$coefficients)])
  out
}


