
# LOAD DATA
# interp <- read.table("simple_names.interp.txt", col.names = c("metric", "better", "x"))
# long.names <- read.table("names_underscores.txt", col.names = c("metric"))

## LIBRARIES
library("lattice")




get.split.indexes <- function(z, stringlen, i=1, maxdiff=20){
  # x is output of gregexpr() 
  # i is the ith element of z - the element you want to analyze
  # maxdiff is is the maximum difference allowed between splits
  ## ASSUMES first indexes are in general < maxdiff apart
  zlen <- length(z[[i]])
  splits <- c()
  if (zlen > 1){
    a <- 0
    for (j in 1:zlen){
      b <- z[[i]][j]
      if(b - a > maxdiff){
        splits <- c(splits, z[[i]][j-1])
        a <- z[[i]][j-1]
        if(b - a > maxdiff){
          splits <- c(splits, b)
          a <- b
        }
      }
    }
    if(stringlen - a > maxdiff){
      splits <- c(splits, b)
    }
  }
  splits
}

add.linebreaks.to.strings <- function(x, split.char=" ", max.line.len=20){
  # x is character vector
  indexes <- gregexpr(pattern = split.char, x)
  xlen <- length(x)
  for (i in 1:xlen){
    string <- x[i]
    stringlen <- nchar(string)
    if (stringlen > max.line.len){
      stringsplit <- strsplit(string, "")[[1]]
      splindex <- get.split.indexes(z = indexes, stringlen = stringlen, i = i, maxdiff = max.line.len)
      newstring <- paste(substring(text = string, first = c(1,(splindex+1)), last = c(splindex-1,stringlen)), sep="", collapse="\n")
      x[i] <- newstring
    }
  }
  x
}


trim.numbers.in.strings <- function(x, digits=6){
  ##x is character vector
  N <- length(x)
  for (i in 1:N){
    ## if string can be converted into numeric, then try trimming
    y <- suppressWarnings( as.numeric(x[i]) )
    if ( !(is.na( y ) )){
      y <- round(y, digits)
      x[i] <- as.character(y)
    }
  }
  #RETURN
  x
}

read.longtable.values <- function(x, digits=6, split.char=" ", perform.splitting=FALSE, replacement = "\n"){
  data <- read.table(x, col.names = c("score"), colClasses = c("character"))
  data$score <- trim.numbers.in.strings(data$score, digits)
  data$score <- gsub(pattern = '\\', replacement = '', x = data$score, fixed=TRUE)
  data$score <- gsub(pattern = '_', replacement = ' ', x = data$score, fixed=TRUE)
  if(perform.splitting){
    data$score <- gsub(pattern = split.char, replacement = replacement, x = data$score, fixed=TRUE)
  }
  data
}

varname.as.string <- function(v1) {
  deparse(substitute(v1))
}

argparse <- function(...){
  as.list(sys.call())
}
generate.longtable <- function(metric.names, ...){
  args <- as.list(sys.call())
  out <- cbind(metric.names, ...)
  N <- length(args)
  colnames(out) <- c("metric",args[3:N]) 
  out
}

## FUNCTIONS
addpoly <- function(x,y,step, y0=0, ...){
  polygon(x = c(x,x,x+step,x+step,x), y = c(y0,y,y,y0,y0), ...)
}

asm.compare <- function(s1, s2){
  s3 <- matrix(nrow = length(s1), ncol = 3)
  for (i in 1:length(s1)){
    s3[i,1:2] <- rank(c(s1[i], s2[i]))
    s3[i,3] <- s3[i,1]-s3[i,2]
  }
  s3
}

asm.comparison.report <- function(s3){
  ## an idea of what is better -- more neg or more pos?
  summed.scores <- sum(s3[,3])
  mean.score <- mean(s3[,3])
  
  # pct metrics that are better in new update
  pct.better <- sum(s3[,3]==1)/length(s3[,3])
  pct.worse <- sum(s3[,3]==-1)/length(s3[,3])
  pct.same <- sum(s3[,3]==0)/length(s3[,3])
  
  
  ## RETURN
  list(summed.scores=summed.scores, mean.score=mean.score, pct.worse=pct.worse, pct.same=pct.same, pct.better=pct.better)
}

asm.comparison.table <- function(interp,s1,s2,s3=c(NA),rows=c(NA), col.names = c(NA), reinterp=FALSE, perform.splitting=TRUE, split.char="_", replace.char = " "){
  ## see as table
  if(is.na(s3[1])){
    s3<-asm.compare(s1,s2)
  }
  if(reinterp){re <- interp$x}
  else{ re <- rep(1, length(interp$x))}
  if(is.na(rows[1])){
    tab <- cbind(interp[,1:2],data.frame(s1=s1*re,s2=s2*re,s3=s3))
  } else{
    tab <- cbind(interp[rows,1:2],data.frame(s1=s1[rows]*re[rows],s2=s2[rows]*re[rows],s3=s3[rows,]))
  }
  if (!(is.na(col.names[1])) & length(col.names) == dim(tab)[2]){
    colnames(tab) <- col.names
  }
  if(perform.splitting){
    tab[,1] <- gsub(pattern = split.char, replacement = replace.char, x = tab[,1], fixed = TRUE)
  }
  #RETURN
  tab
}

asm.comparison.dotplot <- function(s3,col=c(NA)){
  if(is.na(col[1])){col<-c(rep("red",5),"blue",rep("dark green",3), rep("brown",5), "purple",rep("red",2),rep("pink",11),rep("red3",7), rep("dark cyan",2),rep("cyan",11), rep("blue",4),rep("dark blue",6), rep("orange",4),rep("purple",4), rep("grey",6))}
  plot(1:length(s3[,3]), s3[,3], type="n", las=1, ylab="", xlab="",yaxt="n"); abline(h=0, lty=2)
  axis(side=2, at = c(-1,0,1), labels = c("Worse","Same", "Better"), las=1)
  abline(v=c(5.5,6.5,9.5,14.5,15.5,17.5,28.5,35.5,37.5,48.5,52.5,58.5,62.5,66.5,69.5), lty=3)
  points(1:length(s3[,3]), s3[,3],col=col)
}

asm.comparison.pctplot <- function(s3){
  worse <- 100*sum(s3[,3]==-1)/length(s3[,3])
  same <- 100*sum(s3[,3]==0)/length(s3[,3])
  better <- 100*sum(s3[,3]==1)/length(s3[,3])
  x <- seq(0,5,1)
  plot(x,x, ylim=c(0,max(worse,same,better)), type="n", ylab="Percent of Metrics",xlab="",las=1,xaxt="n")
  axis(side=1, at = c(1,2.5,4), labels = c("Worse","Same", "Better"), las=1)
  addpoly(0.5, y = worse, step = 1)
  addpoly(2, y = same, step = 1)
  addpoly(3.5, y = better, step = 1)
}

asm.comparison.pct.categs <- function(s3, rows="all", categs="all"){
  #rows is which rows from s3 to look at -- e.g. some vars I have used are: all, basic, basic.plus
  #categs is a vector of strings as long as rows that categorizes each score
  # --> choose 'all' for it to be defined inside this function
  basic <- c(3,6,7,10,15,16,17,36,37)
  basic.plus <- c(basic, 49:52, 60:61, 64:65, 67:68)
  most <- c(1:28,36:68) ## excludes reapr size stats -- needed this for times when one asm was aggressively broken and another was not
  all <- 1:68
  if(rows == "all"){rows <- all}
  else if(rows == "basic"){rows <- basic}
  else if(rows == "basic.plus"){rows <- basic.plus}
  else if(rows == "most"){rows <- most}
  if (categs == "all"){categs <- c(rep("Contiguity",5), "BUSCO", rep("Illumina",42), rep("BioNano",10), rep("ONT",4),rep("PacBio",4), "ONT", "PacBio")}
  else if(categs == "basic"){categs <- c("size", "genes", rep("Illumina",7))}
  else if(categs == "basic.plus"){categs <- c("size", "genes", rep("Illumina",7), rep("BioNano",4), rep("ONT",2),rep("PacBio",2), "ONT", "PacBio")}
  else if(categs == "most"){categs <- c(rep("Contiguity",5), "BUSCO", rep("Illumina",35), rep("BioNano",10), rep("ONT",4),rep("PacBio",4), "ONT", "PacBio")}
  uniq.categs <- unique(categs)
  nuniq <- length(uniq.categs)
  scoreslist <- list(categs=vector(length = nuniq+1), worse=vector(length = nuniq+1), same=vector(length=nuniq+1), better=vector(length=nuniq+1))
  for(i in 1:nuniq){
    categ <- uniq.categs[i]
    scoreslist$categs[i] <- categ
    scores <- s3[rows,3][categs == categ]
    scoreslist$worse[i] <- 100*sum(scores == -1)/length(scores)
    scoreslist$same[i] <- 100*sum(scores==0)/length(scores)
    scoreslist$better[i] <- 100*sum(scores==1)/length(scores)
  } 
  #MEANS
  i <- i+1
  scoreslist$categs[i] <- "Mean"
  scoreslist$worse[i] <- mean(scoreslist$worse[1:nuniq])
  scoreslist$same[i] <- mean(scoreslist$same[1:nuniq])
  scoreslist$better[i] <- mean(scoreslist$better[1:nuniq])
  #RETURN
  scoreslist
}

asm.comparison.pct.categs.legend <- function(legend = c("Worse", "Same", "Better"), cex=1){
  plot(1:10, type="n", xaxt="n", yaxt="n", ylab="", xlab="", bty="n")
  legend(x=1,y=8, legend = legend, fill = c("dark blue", "grey", "dark red"), cex=cex)  
}

asm.comparison.pct.categs.plot <- function(scoreslist, stacked=FALSE, main="All Metrics", ...){
  ## scoreslist is output of asm.comparison.pct.categs() 
  nuniq <- length(scoreslist$categs)
  if(!(stacked)){
    x <- seq(0,nuniq*6-3,1)
    plot(x,x, ylim=c(0,100), main=main, type="n", ylab="Percent of Metrics",xlab="",las=1,xaxt="n", ...)
    for(i in 1:nuniq){
      xs <- 6*(i-1)
      axis(side=1, at = c(xs,xs+1.5,xs+3), labels = c("Worse","Same", "Better"), las=2)
      addpoly(xs-0.5, y = scoreslist$worse[i], step = 1, col="dark blue")
      addpoly(xs+1, y = scoreslist$same[i], step = 1, col="grey")
      addpoly(xs+2.5, y = scoreslist$better[i], step = 1, col="dark red") 
      mtext(scoreslist$categs[i], side = 1, at = xs+1.5, line=3.5)
    }
  } else {
    ## alt polygon stacks
    x <- seq(0,nuniq*2-1,1)
    plot(x,x, ylim=c(0,100), main=main, type="n", ylab="Percent of Metrics",xlab="",las=1,xaxt="n", ...)
    for(i in 1:nuniq){
      xs <- 2*(i-1)
      axis(side=1, at = xs+0.5, labels = scoreslist$categ[i], las=1)
      worse <- scoreslist$worse[i]
      same <- scoreslist$same[i]
      better <- scoreslist$better[i]
      addpoly(xs, y = worse, step = 1, col="dark blue")
      addpoly(xs, y = same+worse, y0=worse, step = 1, col="grey")
      addpoly(xs, y = better+same+worse, y0=worse+same, step = 1, col="dark red") 
    }
    abline(h=50, lty=2)
  }
}

asm.compare.pipeline <- function(s1,s2,rows="all",main="All Metrics", ...){
  s3<-asm.compare(s1,s2)
  scoreslist <- asm.comparison.pct.categs(s3, rows = rows, categs = rows)
  asm.comparison.pct.categs.plot(scoreslist, stacked = TRUE, main=main, ...)
}

asm.compare.pipeline.3way <- function(s1,s2, main,asm1=NA, asm2=NA,most.instead.of.all=FALSE, legcex=1.5, ...){
  par(mfrow=c(2,2), mar=c(3,4,1,1))
  legend <- c("Worse", "Same", "Better")
  if (!(is.na(asm1) | is.na(asm2))){
    main <- paste0(asm1, " vs. ", asm2)
    legend <- c(asm1, "Same", asm2)
  }
  asm.compare.pipeline(s1,s2,main = paste0(main,": Basic"), rows = "basic", ...)
  asm.compare.pipeline(s1,s2,main = paste0(main,": Basic+"), rows = "basic.plus", ...)
  if (most.instead.of.all){asm.compare.pipeline(s1,s2,main = paste0(main,": Most"), rows = "most", ...)}
  else{asm.compare.pipeline(s1,s2,main = paste0(main,": All"), ...)}
  asm.comparison.pct.categs.legend(legend=legend, cex=legcex)
  par(mfrow=c(1,1))
}



### DUDS

# asm.comparison.heat <- functon(s3){
#   ### ALT VIZ
#   ##
#   lattice.options(axis.padding=list(factor=0.5)) ## Gets rid of white padding on sides of matrix. Use tck = c(1,0) in scales list to get rid of tick marks on top and right.
#   xmax <- length(s3)
#   x<-xmax:1
# #   levelplot(as.matrix(s3), scales=list(x=list(rot=0), y=list(labels=x, at=x), tck = c(1,0)), pretty=TRUE, 
# #             col.regions=colorRampPalette(c("blue", "black","red")), ylab="", xlab="", at=seq(-1,1,0.001), 
# #             colorkey=list(space="right", col=colorRampPalette(c("blue", "black","red")), at=seq(-1,1,0.001), 
# #                           raster=TRUE, tck=0))
#   
#   levelplot(as.matrix(s3[,3]), pretty=TRUE, col.regions=colorRampPalette(c("blue", "white","red")))
# #   levelplot(t(as.matrix(s3[,3])), pretty=TRUE, col.regions=colorRampPalette(c("blue", "white","red")))
# }
