### qPCR Analysis





#############


reactionEfficiency <- function(slope){
  ## reaction efficiency is 10^(-1/slope) -1
  return(10^(-1/slope) - 1)
}

slopeFromRxnEff <- function(RxnEff){
  return(-1/log10(RxnEff + 1))
}

efficiencyTest <- function(slope){
  ## returns logical of pass or fail (T/F)
  ## slope of -3.58 to -3.10 == rxnEff of 90% to 110% -- the acceptable range
  return((slope >= -3.58 & slope <= -3.10))
}

dilutionSeries <- function(startingAmount=25, dilutionFactor=10, numPoints=5, highToLow=TRUE){
  ##startingAmount is ng/ul of 1:1 -- we usually start with 25 ng/ul
  ## dilutionFactor is how many fold you dilute each in series -- default is 10
  ## numPoints is number of samples/concentrations in standard curve -- we usually do 5
    ## 1, 1/10, 1/100, 1/1000, 1/10000
  dilSer <- c()
  for(i in 1:numPoints){
    nextDil <- startingAmount/(dilutionFactor^(i-1))
    dilSer <- c(dilSer, nextDil)
  }
  ## if data has Cts reported in lowToHigh conc. order
  if(!highToLow){dilSer <- dilSer[length(dilSer):1]}
  return(dilSer)
}


meanTm <- function(qPCRdata){
  ## will give mean Tm for a primer pair in qPCR data -- will use all instances (Standard, Known, etc) where it is named
  
  # Grab a list of all detectors
  primerPairs <- levels(qPCRdata$Detector)
  
  # initialize output dataframe
  prodTms <- data.frame(primerPair = primerPairs, meanTm = 0, medTm=0, sdTm = 0, distTm = 0)
 
  # For each detector/primerPair, calculate the primer stats and say whether it passes validation or not
  for(primer in primerPairs){
    isNA <- is.na(subset(qPCRdata, Detector == primer)$Tm)
    allTm <- subset(qPCRdata, Detector == primer)$Tm[!isNA]
    allTm <- sort(allTm)
    meanTm <- mean(allTm)
    medTm <- median(allTm)
    sdTm <- sd(allTm)
    TmString <- ""
    for(i in 1:(length(allTm)-1)){
      TmString <- paste0(TmString, allTm[i], ",")  
    }
    TmString <- paste0(TmString, allTm[length(allTm)])
    
    prodTms[prodTms$primerPair == primer,]$meanTm <- meanTm
    prodTms[prodTms$primerPair == primer,]$medTm <- medTm
    prodTms[prodTms$primerPair == primer,]$sdTm <- sdTm
    prodTms[prodTms$primerPair == primer,]$distTm <- TmString
  }
  return(prodTms)
}


primerValidation <- function(qPCRdata, numTechReps=2, startingAmount=25, dilutionFactor=10, numPoints=4, highToLow=TRUE){
  ## determine slope, efficiency, and R^2 for all
  ## qPCRdata is a table with columns (class): Well (character); Detector(factor); Task (factor); Ct (numeric); Tm (numeric)
  
  # Grab a list of all detectors
  primerPairs <- levels(qPCRdata$Detector)
  
  # initialize output dataframe
  primerStats <- data.frame(primerPair = primerPairs, slope = 0, efficiency = 0, Rsqr = 0, pass=0)
  
  # Create vector of concentrations for standard curve 
  logNg <- rep(log10(dilutionSeries(startingAmount=startingAmount, dilutionFactor=dilutionFactor, numPoints=numPoints, highToLow=highToLow)), each=numTechReps)
  
  # For each detector/primerPair, calculate the primer stats and say whether it passes validation or not
  for(primer in primerPairs){
    isNA <- is.na(subset(qPCRdata, Detector == primer & Task == "Standard")$Ct)
    Cts <- subset(qPCRdata, Detector == primer & Task == "Standard")$Ct[!isNA]
    dilSer <- logNg[!isNA]
    linMod <- lm(Cts ~ dilSer)
    slope <- linMod$coeff[2]
    Rsqr <- (cor(dilSer, Cts))^2
    efficiency <- reactionEfficiency(slope)
    ## Pass or fail? Slope/eff cutoffs are from ABI manual; Rsqr min is from Marcela
    pass <- slope <= -3.1 & slope >= -3.58 & efficiency >= 0.9 & efficiency <= 1.1 & Rsqr >= 0.97
    primerStats[primerStats$primerPair == primer,]$slope <- slope
    primerStats[primerStats$primerPair == primer,]$efficiency <- efficiency
    primerStats[primerStats$primerPair == primer,]$Rsqr <- Rsqr
    primerStats[primerStats$primerPair == primer,]$pass <- pass
  }
  return(primerStats)
}



## retired just use read.csv or read.table where relevant
## header=T, stringsAsFactors=F
# readqPCRdata <- function(qPCRTxtFile){
#   ## qPCRTxtFile Copy the following columns into a separate tab on the qPCR data spreadsheet:
#     ## Well, Detector, Task, Ct, Tm
#     ## Have no empty columns
#     ## then just copy/paste it into a txt file in nano
#     ## load that into here -- that is the qPCR data
#     qPCRdata <- read.table(qPCRTxtFile, header=T, stringsAsFactors=F)
#     suppressWarnings(qPCRdata$Ct <- as.numeric(qPCRdata$Ct))
#     return(qPCRdata)
# }

plotStdCurve <- function(qPCRdata, sampleName, numberTechReps=2, startingAmount=25, dilutionFactor=10, numPoints=4, ...){
  ## for plotting 1 std curve
  ## you may be interested in plotting std curves for all primer pairs on same plot (as well as have points of unknown)
  logNg <- log10(rep(dilutionSeries(startingAmount, dilutionFactor, numPoints), each=numberTechReps))
  isNA <- is.na(subset(qPCRdata, Detector == sampleName & qPCRdata$Task == "Standard")$Ct)
  Ctvals <- subset(qPCRdata, qPCRdata$Detector==sampleName & qPCRdata$Task=="Standard")$Ct[!isNA]
  logNg <- logNg[!isNA]
  plot(logNg, Ctvals, ...)
  fit <- lm(Ctvals ~ logNg)
  lines(logNg, fit$fitted, col="red")
  abline(v=0,col="black")
  grid()
  yint <- paste0("y-int = ", fit$coeff[1])
  slope <- paste0("Slope = ", fit$coeff[2])
  R <- cor(logNg, Ctvals)  
  Rsq <- paste0("R^2 = ", R^2)
  RxnEff <- paste0("Reaction Efficiency = ", reactionEfficiency(fit$coeff[2]))
  text(0, 24, slope, cex=0.8)
  text(0, 23, yint, cex=0.8)
  text(0, 22, Rsq, cex=0.8)
  text(0, 21, RxnEff, cex=0.8)
  return(fit)
}

## split() -- can split dataframe up by sampleName or w/e you want
## dont forget %in%


dCtFE <- function(){
  
}


ddCtFE <- function(sampleSOICt, sampleNormCt, calibratorSOICt, calibratorNormCt){
  ## SOI == site of interest, i.e. where NS emit from
  ## Norm is the site/gene used as interna normalization
  ## Sample is the treatment condition - e.g. ChIP or NSE
  ## Control/calibrator is the non-treatment condition -- e.g. gDNA
  ## the ddCt method is an improvement over dCt method
  ## It compares results from the experimental sample with both a calibrator and a normalizer
  ##    whereas dCt only compares it to a calibrator (does not internally control/normalize)
  ## FE = 2^-ddCt
  ## ddCt = dCt_sample - dCt_calibrator
  ## dCt_sample = sampleSOICt - sampleNormCt
  ## dCt_calibrator = calibratorSOICt - calibratorNormCt
  ## This assumes that the efficiencies for SOI and Norm primer pairs are 'identical'
  ## Std Curve still needed to check consistency of dCt
  dCt_sample = sampleSOICt - sampleNormCt
  dCt_calibrator = calibratorSOICt - calibratorNormCt
  ddCt = dCt_sample - dCt_calibrator
  FE = 2^-ddCt 
  return(FE)
  ## Note: same as (2^sampleSOICt/2^sampleNormCt)/(2^calibratorSOICt/2^calibratorNormCt)
  ## Which is the same as: (2^(sampleSOICt-sampleNormCt))/(2^(calibratorSOICt-calibratorNormCt))
  ## Which is same as: 2^((sampleSOICt-sampleNormCt) - (calibratorSOICt-calibratorNormCt))
  ## Which is same as: 2^(dCt_sample-dCt_calibrator)
  ## Which is same as: 2^ddCt ...QED?...where does the minus sign come in (2^-ddCt)?
  ## -1 comes in b/c higher enrichment is smaller Ct -- i.e. 25 is higher enriched that 27 by 2^-(25-27) = 2^-(-2) = 2^2 = 4-fold
  ## thus, the negative should have been there from the beginning of the derivation
  ## alternatively it could written as normalizer-target (or normalizer/target) -- and the -1 would not be necessary
  ##
  ## This method requires that the efficiencies for both target and normalizer are identical.
  ## What range of deviation is acceptable?
  ## see relativeEfficiencyTable below
}

ddCtFE2 <- function(sampleSOICt, sampleNormCt, calibratorSOICt, calibratorNormCt){
  ## diff way of implementing
  dCt_sample = 2^(sampleSOICt) / 2^(sampleNormCt)
  dCt_calibrator = 2^(calibratorSOICt) / 2^(calibratorNormCt)
  ddCt = dCt_sample / dCt_calibrator
  FE = 1/ddCt ## b/c if sample more enriched it will be smaller number e.g. 0.5 for 2X
  return(FE)
}



alt.ddCtFE <- function(sampleSOICt, sampleNormCt, calibratorSOICt, calibratorNormCt){
  ## This instead does:
  dCt_SOI = sampleSOICt - calibratorSOICt
  dCt_normalizer = sampleNormCt - calibratorNormCt
  ddCt = dCt_SOI - dCt_normalizer
  FE = 2^-ddCt
  return(FE)
}

allFE <- function(qPCRdata, normalizer, numTechReps=3, task="Unknown", method=ddCtFE, testSampFirst=TRUE){
  primerPairs <- levels(qPCRdata$Detector)
  if(testSampFirst){
    testSample <- 1:numTechReps
    calibratorSample <- (numTechReps+1):(2*numTechReps) 
  } else {
    calibratorSample <- 1:numTechReps
    testSample <- (numTechReps+1):(2*numTechReps)
  }
  testNormalizerCt <- mean(subset(qPCRdata, Detector == normalizer & Task == task)$Ct[testSample], na.rm=TRUE)
  calibratorNormalizerCt <- mean(subset(qPCRdata, Detector == normalizer & Task == task)$Ct[calibratorSample], na.rm=TRUE)
  feTable <- data.frame(primerPair = primerPairs, FE = 0)
  for(primer in primerPairs){
    primerdata <- subset(qPCRdata, Detector == primer & Task == task)
    testTargetCt <- mean(primerdata$Ct[testSample], na.rm=TRUE) 
    calibratorTargetCt <- mean(primerdata$Ct[calibratorSample], na.rm=TRUE)
    FE <- method(sampleSOICt=testTargetCt, sampleNormCt=testNormalizerCt, calibratorSOICt=calibratorTargetCt, calibratorNormCt=calibratorNormalizerCt)
    ##add to df
    feTable[feTable$primerPair == primer,]$FE <- FE
  }
  feTable
}

cov <- function(qPCRdata, numTechReps=3, task="Unknown", testSampFirst=TRUE, covcutoff=1, doremoval=FALSE){
  primerPairs <- levels(qPCRdata$Detector)
  if(testSampFirst){
    testSample <- 1:numTechReps
    calibratorSample <- (numTechReps+1):(2*numTechReps) 
  } else {
    calibratorSample <- 1:numTechReps
    testSample <- (numTechReps+1):(2*numTechReps)
  }
#   testNormalizerCt <- mean(subset(qPCRdata, Detector == normalizer & Task == task)$Ct[testSample], na.rm=TRUE)
#   calibratorNormalizerCt <- mean(subset(qPCRdata, Detector == normalizer & Task == task)$Ct[calibratorSample], na.rm=TRUE)
  primerpair <- c()
  sampletype <- c()
  mu <- c()
  sigma <- c()
  covs <- c()
  pass <- c()
  rmc <- c()
  for(primer in primerPairs){
    primerdata <- subset(qPCRdata, Detector == primer & Task == task)
    Ct <- primerdata$Ct[testSample]
    wells <- primerdata$Well[testSample]
    primerpair <- c(primerpair, primer)
    sampletype <- c(sampletype, "test")
    m <- mean(Ct, na.rm = TRUE)
    mu <- c(mu, m)
    s <- sd(Ct, na.rm = TRUE)
    sigma <- c(sigma, s)
    cv <- 100*s/m
    covs <- c(covs, cv)
    if (cv > covcutoff) {
      pass <- c(pass, 0) 
      removal.cand <- wells[which.max(abs(Ct-m))]
      if(doremoval && sum(!(is.na(Ct))) >= 2){qPCRdata$Ct[qPCRdata$Well == removal.cand] <- NA}
    }
    else {pass <- c(pass, 1); removal.cand <- "-"}
    rmc <- c(rmc, removal.cand)
    Ct <- primerdata$Ct[calibratorSample]
    wells <- primerdata$Well[calibratorSample]
    primerpair <- c(primerpair, primer)
    sampletype <- c(sampletype, "calibrator")
    m <- mean(Ct, na.rm = TRUE)
    mu <- c(mu, m)
    s <- sd(Ct, na.rm = TRUE)
    sigma <- c(sigma, s)
    cv <- 100*s/m
    covs <- c(covs, cv)
    if (cv > covcutoff) {
      pass <- c(pass, 0) 
      removal.cand <- wells[which.max(abs(Ct-m))]
      if(doremoval && sum(!(is.na(Ct))) >= 2){qPCRdata$Ct[qPCRdata$Well == removal.cand] <- NA}
    }
    else {pass <- c(pass, 1); removal.cand <- "-"}
    rmc <- c(rmc, removal.cand)
  }
  covtable <- data.frame(primerPair = primerpair, sampleType = sampletype, mean = mu, stdev = sigma, cov = covs, pass = pass, rm.suggest = rmc)
  if(doremoval){qPCRdata}
  else {covtable}
}

coefplot <- function(data, covcutoff=2){
  X <- round(dim(data)[1]/3)
  x <- 0:(X+1)
  y <- seq(10, 35, length.out = length(x))
  plot(x,y,type="n", las=1,ylab="Ct values", xlab="samples")
  abline(h=c(15,20,25,30), lty=3)
  for(i in seq(3,X*3,3)){
    d <- data$Ct[(i-2):i]
    wells <- data$Well[(i-2):i]
    abline(v = i/3, lty=3)
    points(rep(i/3,3), c(d))
    #   
#     a <- c(i, d, sd(d, na.rm = TRUE), mean(d, na.rm = TRUE))
#     print(a)
    m <- mean(d, na.rm = TRUE)
    s <- sd(d, na.rm = TRUE)
    coef <- 100*s/m
    removal.cand <- wells[which.max(abs(d-m))]
    if(coef > covcutoff){
      y <- 18+rnorm(1,0,2)
      col <- c("black","red","blue", "dark blue", "dark cyan")[i%%5+1]
      text(x = i/3, y, labels = data$Detector[i], cex=0.75, col=col)
      text(x = i/3, y-2, labels = round(coef,digits = 3), cex=0.75, col=col)
      text(x = i/3, y-4, labels = removal.cand, cex=0.75, col=col)
    }
  }
}
  

groupQPCRdata <- function(qPCRdata){
  primerPairs <- levels(qPCRdata$Detector)
  for(primer in primerPairs){
    print(subset(qPCRdata, Detector == primer))
  }
}

avgCTs <- function(qPCRdata, numTechReps=2, startingAmount=25, dilutionFactor=10, numPoints=4, highToLow=TRUE, task="Standard"){
  # Grab a list of all detectors
  primerPairs <- levels(qPCRdata$Detector)
  
  # Create vector of concentrations for standard curve 
  Ng <- dilutionSeries(startingAmount=startingAmount, dilutionFactor=dilutionFactor, numPoints=numPoints, highToLow=highToLow)
  logNgLevels <- log10(Ng)
  logNg  <- rep(logNgLevels, each=numTechReps)
  
  # initialize output dataframe
  # make each primer pair repeat as many times as there are numPoints
  repPrimerPairs <- rep(primerPairs, each=numPoints)
  # make each primer be associated with each ng level
  repLogNgLevels <- rep(logNgLevels, times=length(primerPairs))
  ## df
  ctAvgs <- data.frame(primerPair = repPrimerPairs, logNg = repLogNgLevels, avgCt = 0, sdCt=0, coeffVar=0)
  
  # For each detector/primerPair, calculate the primer stats and say whether it passes validation or not
  for(primer in primerPairs){
    isNA <- is.na(subset(qPCRdata, Detector == primer & Task == task)$Ct)
    Cts <- subset(qPCRdata, Detector == primer & Task == task)$Ct[!isNA]
    dilSer <- logNg[!isNA]
    for(ngAmt in logNgLevels){
      relevantCts <- dilSer == ngAmt
      avgCt <- mean(Cts[relevantCts])#, na.rm=TRUE)  <-- should not need that as NAs removed as a processing step
      sdCt <- sd(Cts[relevantCts], na.rm=TRUE)
      coeffVar <- sdCt/avgCt
      ctAvgs[ctAvgs$primerPair == primer & ctAvgs$logNg == ngAmt,]$avgCt <- avgCt
      ctAvgs[ctAvgs$primerPair == primer & ctAvgs$logNg == ngAmt,]$sdCt <- sdCt
      ctAvgs[ctAvgs$primerPair == primer & ctAvgs$logNg == ngAmt,]$coeffVar <- coeffVar
      ## Z-scores not interesting -- for 2 tech reps, it is always: 0.707106781186547,-0.707106781186547
      #Zstring <- character()
      #z.scores <- (Cts[relevantCts] - avgCt)/sdCt
      #for(i in 1:(numTechReps-1)){
      #  Zstring <- paste0(Zstring, z.scores[i], ",")
      #}
      #Zstring <- paste0(Zstring, z.scores[numTechReps])
      #ctAvgs[ctAvgs$primerPair == primer & ctAvgs$logNg == ngAmt,]$zscores = Zstring
    }
  }
  return(ctAvgs)
  
}

relativeEfficiencyTable <- function(normalizerAvgCts, targetAvgCts, dilSeries){
  ## The ddCtFE method requires that the efficiencies for both target and normalizer are identical.
  ## What range of deviation is acceptable?
  ## The way to determine this is to generate a standard curve for botht the target and normalizer using the same samples.
  ## The average dCt between the normalizer and target is obtained for each dilution.
  ## The value itself is not important; it is the consistency of that value across each dilution in series that matters.
  ## Make table of dilution series:
  ## Input amount (ng) | Normalizer Avg Ct | Target Avg Ct | dCt (Normalizer - Target)
  ##
  ## normalizerAvgCts = vector of average Ct values for each dilution series point 
  ## targetAvgCts = vector of average Ct values for each dilution series point 
  ## dilSeries = vector of ng quantities (not log10) in dilution series starting at 1:1 
  ##    -- can be output of dilutionSeries() 
  RelEffTable <- data.frame(inputNg = dilSeries, normAvgCt = normalizerAvgCts, targetAvgCt = targetAvgCts, dCt = normalizerAvgCts-targetAvgCts, relCopyNum = 2^(normalizerAvgCts-targetAvgCts))
  return(RelEffTable)
}

relativeEfficiencyPlot <- function(RelEffTable, ...){
  ## takes output of relativeEfficiencyTable
  xVals <- 1:length(RelEffTable$inputNg)
  plot(xVals, RelEffTable$dCt, ylab="dCt", xlab="Input Amount (ng)", col="red", type="b", xaxt="n", xlim=c(0,(length(RelEffTable$inputNg)+1)), ylim=c((min(RelEffTable$dCt)-2), (max(RelEffTable$dCt)+2)), ...)
  axis(side=1, at=xVals, labels=RelEffTable$inputNg)
  grid()
  fit <- lm(RelEffTable$dCt ~ xVals)
  lines(xVals, fit$fitted, col="black", type="b")
  slope <- fit$coeff[2] 
  ## Note that it is not strictly ng vs. dCt -- they are all equally spaced (mapped to integers) despite not really being equally spaced
  ##  This is exactly how done in Life Tech book and produces same results with same data
  ## ideally slope = 0 if both have identical efficiencies 
  ## but anything < 0.1 is acceptable for the ddCt method --- currently this puts out negative slopes though so > -0.1 to 0
  yint <- fit$coeff[1]
  line <- paste0("Y = ", slope, "x + ", yint)
  legend("topright", legend=c("dCt",line), fill=c("red", "black"))
  if(slope <= 0.1 & slope >= 0){
    return(paste0("Efficiencies of Target and Normalizer pass consistency test with slope = ", slope))
  }
  else{return(paste0("Efficiencies of Target and Normalizer FAILED consistency test with slope = ", slope))}
}



fullRelativeEfficiencyTable <- function(avgCtTable, normalizer, startingAmount=25, dilutionFactor=10, numPoints=4, highToLow=TRUE){
  ## takes in avgCtTable, the output of avgCTs() -- mutliple avgCTs() outputs may even be bound together with rbind() beforehand
  ## normalizer can be any primerPair name (i.e. Detector) -- typical for MYC locus is MSF1-1 or JU-9-1 (e.g. normalizer="JU-9-1")
  primerPairs <- levels(avgCtTable$primerPair)
  normCts <- subset(avgCtTable, primerPair == normalizer)$avgCt
  fullRelTable <- data.frame(primerPair = NULL, inputNg = NULL, normAvgCt = NULL, targetAvgCt = NULL, dCt = NULL)
  for(primer in primerPairs){
    dilSer <- dilutionSeries(startingAmount=startingAmount, dilutionFactor=dilutionFactor, numPoints=numPoints, highToLow=highToLow)
    targetCts <- subset(avgCtTable, primerPair == primer)$avgCt
    primerPair <- rep(primer, numPoints)
    newTable <- relativeEfficiencyTable(normalizerAvgCts=normCts, targetAvgCts=targetCts, dilSeries=dilSer)
    newTable <- cbind(primerPair, newTable)
    fullRelTable <- rbind(fullRelTable, newTable)
  }
  return(fullRelTable)
}

relativeEfficiencyTest <- function(fullRelTable){
  ## takes in output table from fullRelativeEfficiencyTable()
  primerPairs <- levels(fullRelTable$primerPair)
  relEffTestTable <- data.frame(primerPair = primerPairs, slope = 0, pass = 0)
  for(primer in primerPairs){
    yVals <- subset(fullRelTable, primerPair == primer)$dCt
    xVals <- 1:length(yVals)
    fit <- lm(yVals ~ xVals)
    slope <- fit$coeff[2] 
    ## Note that it is not strictly ng vs. dCt -- they are all equally spaced (mapped to integers) despite not really being equally spaced
    ##  This is exactly how done in Life Tech book and produces same results with same data
    ## ideally slope = 0 if both have identical efficiencies 
    ## but any |slope| < 0.1 is acceptable for the ddCt method --- currently this puts out negative slopes though so > -0.1 to 0
    pass <- as.numeric(abs(slope) >= 0 & abs(slope) <= 0.1)
    relEffTestTable[relEffTestTable$primerPair == primer, ]$slope = slope
    relEffTestTable[relEffTestTable$primerPair == primer, ]$pass = pass
  }
  return(relEffTestTable)
}


mockRelativeCopyNumber <- function(dCt){
  ##Takes in dCt, outputs 'FE'
  return(2^-dCT)
}
  
  
qPCRrecipe <- function(singleRxnVol, DNAvol, primerMol=10, numRxns=1){
  ## primerMol should be reported as uM (e.g. 10 uM or 1 uM)
  ## singleRxnVol is the total volume of a single reaction -- e.g. 15 ul, 20 ul, or 25 ul
  ## target primer conc is 100 nM
  ## e.g. 10 uMol/L * (25/100)ul / 25 ul = 10*25/100*25 uMol/L = 1/10 uM = 100 nM
  ## primerMol*X/rxnVol = 100 nM
  ## X = 0.1 uM * rxnVol / primerMol
  Fwd = 0.1*singleRxnVol/primerMol ## same as vol1() below
  Rev = Fwd
  SYBRGreen2X = singleRxnVol/2
  UPW = singleRxnVol - Fwd - Rev - SYBRGreen2X - DNAvol
  cat(DNAvol, " ul DNA per reaction (not part of master mix) \n", numRxns*UPW, " ul Ultra Pure Water \n", numRxns*SYBRGreen2X, " ul 2X SYBR Green \n", numRxns*Fwd, " ul Fwd \n", numRxns*Rev, " ul Rev \n")
  cat()
  cat("Use ", SYBRGreen2X+Fwd+Rev+UPW, " ul of the master mix per reaction. Then add the ", DNAvol, " ul of DNA to that.")
}

vol1 <- function(C2,V2,C1){
  ## C1*V1 = C2*V2
  ## what volume of a given [C1] solution is needed to make V2 vol of [C2] solution?
  ## e.g. what volume of 10 uM stock is needed to make 100 ul of 1 uM stock?
  ## 10uM * V1 = 1 uM * 100 ul
  ## V1 = 1uM * 100ul / 10uM
  ## V1 = C2*V2/C1
  V1 <- C2*V2/C1
  volWater <- V2 - V1
  return(c(V1, volWater))
}


compactAvgCts <- function(fullRelTable){
  primerPairs <- levels(fullRelTable$primerPair)
  avgCtsCompact <- data.frame(primerPair = primerPairs, avgCts = 0)
  for(primer in primerPairs){
    CTs <- subset(fullRelTable, primerPair == primer)$targetAvgCt
    CTstring <- ""
    for(i in 1:(length(CTs)-1)){
      CTstring <- paste0(CTstring, CTs[i], ",")  
    }
    CTstring <- paste0(CTstring, CTs[length(CTs)])
    avgCtsCompact[avgCtsCompact$primerPair == primer, ]$avgCts = CTstring
  }
  return(avgCtsCompact)
}







#TODO/TOFINISH
stdCurveMethFD <- function(E_target, E_Norm, calibratorSOICt, sampleSOICt, calibratorNormCt, sampleNormCt){
  ##### THIS IS POTENTIALLY NOT IMPLEMENTED CORRECTLY ######
  ## this improves on ddCt by adding an reaction efficiency (between SOI and normalizer sites) correction
  ##life tech pg 44
  ## fold difference = (E_target*dCt_target)/(E_normalizer*dCt_normalizer)
  ## E = efficiency from standard curve = 10^(-1/slope) = reactionEfficiency(slope)
  ## dCt_target = calibratorSOICt - sampleSOICt
  ## dCt_normalizer = calibratorNormCt - sampleNormCt
  ##
  ## Recall:
  ##  calibrator sample is the control sample such as gDNA control
  ##  normalizer locus is the internal control locus
  dCt_target <- calibratorSOICt - sampleSOICt
  dCt_normalizer <- calibratorNormCt - sampleNormCt
  FD <- (E_target*dCt_target)/(E_Norm*dCt_normalizer)
  
  ## above seeems wrong ... im going to try it as efficiency adjusted ddCt
  dCt_sample = E_target*sampleSOICt - E_Norm*sampleNormCt
  dCt_calibrator = E_target*calibratorSOICt - E_Norm*calibratorNormCt
  effAdjddCt = dCt_sample - dCt_calibrator
  FE = 2^-effAdjddCt
  ## this seems wrong too.....
  return(c(FD, FE))
}

