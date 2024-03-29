
args <- commandArgs(trailingOnly=TRUE)

gcdist <- args[1] ## Only argument is GC distribution file -- 2 columns: %GC, comma-separated values seen for that %GC

PSEUDOCOUNT <- 1

if (length(args) == 2){
  fname <- args2
} else {
  fnamevec <- unlist(strsplit(args[1], split='/'))
  basename <- fnamevec[length(fnamevec)]
  fname <- paste0(basename,".gamma.pdf")
}
print(fname)

data <- read.table(gcdist, col.names=c('gc','dist'), colClasses=c('numeric','character'), stringsAsFactors = FALSE)

pdf(fname)

DIST <- c()
N <- 0
SUMTOTAL <- 0
for (gc in data$gc){
  d <- as.numeric( unlist( strsplit(data$dist[data$gc == gc], split = ',') ) )
  d <- d+PSEUDOCOUNT
  n <- length(d)
  mu <- mean(d)
  med <- median(d)
  std <- sd(d)
  mad <- median(abs(d-med))
  shape.mu <- mu^2 / std^2  ## mu^2 / var
  scale.mu <- std^2 / mu ## var / mu

  shape.med <- med^2 / mad^2
  scale.med <- mad^2 / med

  shape.robust <- med^2 / (mad*1.4826)^2
  scale.robust <- (mad*1.4826)^2 / med

  mn <- min(d)
  mx <- max(d)
  
  DIST <- c(DIST, d)
  N <- N+n

  rmu <- rgamma(n = n, shape = shape.mu, scale=scale.mu)
  rmed <- rgamma(n = n, shape = shape.med, scale=scale.med)
  rrobust <- rgamma(n = n, shape = shape.robust, scale = scale.robust)

  RNG <- range(d,rmu,rmed,rrobust)
  l.out <- round(min(c(500, n/10)))
  if (l.out < 2) {l.out <- 2}
  breaks <- seq(RNG[1]-mad, RNG[2]+mad, length.out=l.out)
  print(c(gc, RNG[1], RNG[2], l.out))

  hd <- hist(d, breaks=breaks, plot = FALSE)
  hrmu <- hist(rmu, breaks=breaks, plot=FALSE)
  hrmed <- hist(rmed, breaks=breaks, plot=FALSE)
  hrobust <- hist(rrobust, breaks=breaks, plot=FALSE)  

  plot(hd$mids, hd$counts, type='l', main=paste0('GC = ', gc), lwd=3, col="grey")
  lines(hrmu$mids, hrmu$counts, col='red', lty=3)
  lines(hrmed$mids, hrmed$counts, col='blue', lty=3)
  lines(hrobust$mids, hrobust$counts, col='orange', lty=3, lwd=2)
  points(mu, 0, col='blue')
  points(med, 0,col='red')
}

## ALL
mu <- mean(DIST)
med <- median(DIST)
std <- sd(DIST)
mad <- median(abs(DIST-med))
shape.mu <- mu^2 / std^2  ## mu^2 / var
scale.mu <- std^2 / mu ## var / mu
shape.med <- med^2 / mad^2
scale.med <- mad^2 / med
shape.robust <- med^2 / (mad*1.4826)^2
scale.robust <- (mad*1.4826)^2 / med
mn <- min(DIST)
mx <- max(DIST)
rmu <- rgamma(n = N, shape = shape.mu, scale=scale.mu)
rmed <- rgamma(n = N, shape = shape.med, scale=scale.med)
rrobust <- rgamma(n = N, shape = shape.robust, scale = scale.robust)

RNG <- range(DIST,rmu,rmed,rrobust)
l.out <- round(min(c(500, N/10)))
if (l.out < 2) {l.out <- 2}
breaks <- seq(RNG[1]-mad, RNG[2]+mad, length.out=l.out)
print(c("ALL", RNG[1], RNG[2], l.out))
hd <- hist(DIST, breaks=breaks, plot = FALSE)
hrmu <- hist(rmu, breaks=breaks, plot=FALSE)
hrmed <- hist(rmed, breaks=breaks, plot=FALSE)
hrobust <- hist(rrobust, breaks=breaks, plot=FALSE)

plot(hd$mids, hd$counts, type='l', main=paste0('GC = ALL'), lwd=3, col="grey")
lines(hrmu$mids, hrmu$counts, col='red', lty=3)
lines(hrmed$mids, hrmed$counts, col='blue', lty=3)
lines(hrobust$mids, hrobust$counts, col='orange', lty=3, lwd=2)
points(mu, 0, col='blue')
points(med, 0,col='red')

q <- quantile(x=DIST, probs = 0.99)  
plot(hd$mids[hd$mids<q], hd$counts[hd$mids<q], type='l', main=paste0('GC = ALL up to 99_pctile'), lwd=3, col="grey")
lines(hrmu$mids, hrmu$counts, col='red', lty=3)
lines(hrmed$mids, hrmed$counts, col='blue', lty=3)
lines(hrobust$mids, hrobust$counts, col='orange', lty=3, lwd=2)
points(mu, 0, col='blue')
points(med, 0,col='red')

  
dev.off()




