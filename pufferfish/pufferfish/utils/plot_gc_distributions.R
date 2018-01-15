
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
  mn <- min(d)
  mx <- max(d)
  
  rmu <- rnorm(n = n, mean = mu, sd = std)
  rmed <- rnorm(n = n, mean = med, sd = mad)

  RNG <- range(d,rmu,rmed)
  l.out <- min(c(500, n/10))
  breaks <- seq(RNG[1]-mad, RNG[2]+mad, length.out=l.out)
  hd <- hist(d, breaks=breaks, plot = FALSE)
  hrmu <- hist(rmu, breaks=breaks, plot=FALSE)
  hrmed <- hist(rmed, breaks=breaks, plot=FALSE)
  
  plot(hd$mids, hd$counts, type='l', main=paste0('GC = ', gc))
  lines(hrmu$mids, hrmu$counts, col='red')
  lines(hrmed$mids, hrmed$counts, col='blue')
  points(mu, 0, col='blue')
  points(med, 0,col='red')
}

  
dev.off()




