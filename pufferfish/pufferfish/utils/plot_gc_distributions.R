
args <- commandArgs(trailingOnly=TRUE)

gcdist <- args[1] ## Only argument is GC distribution file -- 2 columns: %GC, comma-separated values seen for that %GC

if (length(args) == 2){
  fname <- args2
} else {
  fnamevec <- unlist(strsplit(args[1], split='/'))
  basename <- fnamevec[length(fnamevec)]
  fname <- paste0(basename,".pdf")
}
print(fname)

data <- read.table(gcdist, col.names=c('gc','dist'), colClasses=c('numeric','character'), stringsAsFactors = FALSE)

pdf(fname)

for (gc in data$gc){
  d <- as.numeric( unlist( strsplit(data$dist[data$gc == gc], split = ',') ) )
  n <- length(d)
  mu <- mean(d)
  med <- median(d)
  std <- sd(d)
  mad <- median(abs(d-med))
  
  
  rmu <- rnorm(n = n, mean = mu, sd = std)
  rmed <- rnorm(n = n, mean = med, sd = mad)
  breaks <- max(5, round(n/5))
  hd <- hist(d, breaks=breaks, plot = FALSE)
  hrmu <- hist(rmu, breaks=breaks, plot=FALSE)
  hrmed <- hist(rmed, breaks=breaks, plot=FALSE)
  
  plot(hd$mids, hd$counts, type='l', main=paste0('GC = ', gc))
  lines(hrmu$mids, hrmu$counts, col='red')
  lines(hrmed$mids, hrmed$counts, col='blue')
  abline(v=gc, lwd=2)
  abline(v=mu, lwd=2, col='blue')
  abline(v=med, lwd=2, col='red')
}

  
dev.off()




