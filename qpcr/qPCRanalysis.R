### qPCR Analysis





#############


reactionEfficiency <- function(slope,logbase=10){
  ## reaction efficiency is 10^(-1/slope) -1
  return(logbase^(-1/slope) - 1)
}

slopeFromRxnEff <- function(RxnEff, logbase=10){
  return(-1/log(RxnEff + 1, base=logbase))
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
  if(length(dilutionFactor) == 1){
    dilSer <- c()
    for(i in 1:numPoints){
      nextDil <- startingAmount/(dilutionFactor^(i-1))
      dilSer <- c(dilSer, nextDil)
    }
  } else if (length(dilutionFactor)>1){
    dilSer <- vector(mode = "numeric", length=numPoints)
    dilSer[1] <- startingAmount
    for(i in 2:numPoints){
      dilSer[i] <- dilSer[i-1]/dilutionFactor[i-1]
    }
  }
  ## if data has Cts reported in lowToHigh conc. order
  if(!highToLow){dilSer <- dilSer[length(dilSer):1]}
  return(dilSer)
}


getlogNg <- function(dilutionFactor,startingAmount,numPoints,highToLow, numTechReps,logbase){
  ## THIS WAS PORTED TO dlutionSeries()
  # if(length(dilutionFactor)==1){
  #   relAmnts <- dilutionSeries(startingAmount=startingAmount, dilutionFactor=dilutionFactor, numPoints=numPoints, highToLow=highToLow)
  # } else if (length(dilutionFactor)>1){
  #   relAmnts <- vector(mode = "numeric", length=numPoints)
  #   relAmnts[1] <- startingAmount
  #   for(i in 2:numPoints){
  #     relAmnts[i] <- relAmnts[i-1]/dilutionFactor[i-1]
  #   }
  # }
  relAmnts <- dilutionSeries(startingAmount=startingAmount, dilutionFactor=dilutionFactor, numPoints=numPoints, highToLow=highToLow)
  logNg <- rep(log(relAmnts, base=logbase), each=numTechReps)
  return(logNg)
}

meanTm <- function(qPCRdata){
  ## will give mean Tm for a primer pair in qPCR data -- will use all instances (Standard, Known, etc) where it is named
  
  # Grab a list of all detectors
  #primerPairs <- levels(qPCRdata$Detector)
  primerPairs <- unique(as.character(qPCRdata$Detector))
  
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


primerValidation <- function(qPCRdata, numTechReps=2, startingAmount=2, dilutionFactor=4, numPoints=4, highToLow=TRUE, single=FALSE, qpcrtask="Standard", logbase=10){
  ## determine slope, efficiency, and R^2 for all
  ## qPCRdata is a table with columns (class): Well (character); Detector(factor); Task (factor); Ct (numeric); Tm (numeric)
  
  # Grab a list of all detectors
  #primerPairs <- levels(qPCRdata$Detector)
  primerPairs <- unique(as.character(qPCRdata$Detector))
  # print(primerPairs)
  # initialize output dataframe
  primerStats <- data.frame(primerPair = primerPairs, slope = 0, efficiency = 0, Rsqr = 0, pass=0)
  
  # Create vector of concentrations for standard curve 
  logNg <- getlogNg(dilutionFactor,startingAmount,numPoints,highToLow, numTechReps,logbase)
  

  # For each detector/primerPair, calculate the primer stats and say whether it passes validation or not
  for(primer in primerPairs){
    # print(primer)
    isNA <- is.na(subset(qPCRdata, Detector == primer & Task == qpcrtask)$Ct)
    # print(qPCRdata)
    # print(isNA)
    Cts <- subset(qPCRdata, Detector == primer & Task == qpcrtask)$Ct[!isNA]
    dilSer <- logNg[!isNA]
    # print(Cts)
    # print(dilSer)
    linMod <- lm(Cts ~ dilSer)
    # print(linMod)
    slope <- linMod$coeff[2]
    Rsqr <- (cor(dilSer, Cts))^2
    efficiency <- reactionEfficiency(slope, logbase = logbase)
    # print(slope)
    # print(Rsqr)
    # print(efficiency)
    ## Pass or fail? Slope/eff cutoffs are from ABI manual; Rsqr min is from Marcela
    pass <- slope <= -3.1 & slope >= -3.58 & efficiency >= 0.9 & efficiency <= 1.1 & Rsqr >= 0.97
    primerStats[primerStats$primerPair == primer,]$slope <- slope
    primerStats[primerStats$primerPair == primer,]$efficiency <- efficiency
    primerStats[primerStats$primerPair == primer,]$Rsqr <- Rsqr
    primerStats[primerStats$primerPair == primer,]$pass <- pass
  }
  if(single){return(list(stats=primerStats, model=linMod))}
  else {return(primerStats)}
}

autoValidatePrimers <- function(data, controls=NA, numTechReps=2, startingAmount=2,dilutionFactor=4,numPoints=4,highToLow=TRUE, qpcrtask="Standard", logbase=10){
  if(is.na(controls[1])){print("Must add vector of control detectors for relative efficiency test -- e.g., c(C1) or c(C1,C2)"); return()}
#   print(1)
  primerVal <- primerValidation(data, numTechReps, startingAmount, dilutionFactor, numPoints, highToLow, qpcrtask = qpcrtask, logbase = logbase)
  avgCtTable <- avgCTs(data, numTechReps, startingAmount, dilutionFactor, numPoints, highToLow, qpcrtask, logbase = logbase)
#   print(2)
#   print(avgCtTable)
  frt <- fullRelativeEfficiencyTable(avgCtTable=avgCtTable, normalizer=controls[1], startingAmount=startingAmount, dilutionFactor=dilutionFactor, numPoints=numPoints, highToLow=highToLow)
#   print(3)
#   print(frt)
  rett1 <- relativeEfficiencyTest(fullRelTable=frt)
#   print(4)
#   print(rett1)
  all <- rett1[,c(1,3)]
#   print(5)
#   print(controls)
  if(length(controls) > 1){
    for (control in controls[2:length(controls)]){
      #     print(control)
      frt <- fullRelativeEfficiencyTable(avgCtTable=avgCtTable, normalizer=control, startingAmount=startingAmount, dilutionFactor=dilutionFactor, numPoints=numPoints, highToLow=highToLow)
      rett <- relativeEfficiencyTest(fullRelTable=frt)
#       print(all)
      all <- cbind(all, rett[,3])
      #     print(all)
#         print(6)
#       print(all)
#       print(controls)
#       print(dim(all))
      K <- dim(all)[2]-1
      colnames(all) <- c("pair",controls[1:K])
      rownames(all) <- all[,1]
#         print(7)
#         print(all)
      rtest1_controlscores <- apply(X = all[,2:dim(all)[2]], MARGIN = 2, FUN = sum)
      #transform counts into percents
      rtest1_controlscores <- rtest1_controlscores/length(primerVal$primerPair)
      rtest1 <- apply(X = all[,2:dim(all)[2]], MARGIN = 1, FUN = sum)
      rscores <- c()
      # Be sure it is in same order
      for(primer in primerVal$primerPair){rscores <- c(rscores, rtest1[primer])}
#         print(8)
      # transform counts into percents
      rscores <- rscores/length(controls)
      final1 <- cbind(primerVal,rscores)
    } 
  } else {
    ## for now it assumes 2 or more
    rscores <- rett1[,3]
    final1 <- cbind(primerVal,rscores)
    rtest1_controlscores <- rep(1, length(controls))
  }
#   print(final1)
#   cat("\n")
#   print(rtest1)
#   print(rscores)
#   cat("\n")
#   print(rtest1_controlscores)
# print(9)
  list(scores=final1, controlscores=rtest1_controlscores)
}

autoCorrectPrimers <- function(data, scores, numTechReps=2, startingAmount=2, dilutionFactor=4, numPoints=4, highToLow=TRUE, qpcrtask="Standard", logbase=10){
  ## scores is from output of autoValidatePrimers -- ans$scores
  if(dilutionFactor <= 4){xlim <- c(-2,2)}
  else{xlim <- c(-4,1)}
  n <- length(scores$primerPair)
  notes <- ""
  for(i in 1:n){
    needsfixing <- scores$pass[i] == 0
    if(needsfixing){
      fix <- as.character(scores$primerPair[i])
      notes <- paste0(notes, "Removed from ", fix, ": ")
      cat("Fixing ", fix, "...\n")
      cat("Observe data below and plot...\n")
      print(subset(data, Detector == fix))
#       print("HEY"); print(c(numTechReps,startingAmount, dilutionFactor))
      plotStdCurve(qPCRdata = data, sampleName = fix, numTechReps, startingAmount, dilutionFactor, numPoints, qpcrtask, ylim=c(15,33), xlim=xlim, col="blue", lwd=3, main = paste0(fix, ", Initial condition"), logbase = logbase) ## std curve
      del <- 0
#       print(needsfixing)
#       print(del <= 4)
      while(needsfixing & del <= 4){
        del <- del+1
#         print(del)
        subtable <- subset(data, subset = Detector == fix)
        subtable$Detector <- factor(subtable$Detector) ## reset factor levels on new subtable to be only "fix"
        ## Get model
        fmod <- primerValidation(subtable, numTechReps, startingAmount, dilutionFactor, numPoints, highToLow, single=TRUE, qpcrtask = qpcrtask, logbase = logbase)
        ## Get residuals
        fres <- residuals(fmod$model)
        ## find index of largest residual
        fresmax <- which.max(fres)
        ## NA the Ct
        data$Ct[data$Detector == fix][fresmax] = NA
        #NOTES
        notes <- paste0(notes, data$Well[data$Detector == fix][fresmax], " ")
        subtable <- subset(data, subset = Detector == fix)
        subtable$Detector <- factor(subtable$Detector) ## reset factor levels on new subtable to be only "fix"
        ## Get pass info
        finfo <- primerValidation(subtable, numTechReps, startingAmount, dilutionFactor, numPoints, highToLow, single=TRUE, qpcrtask=qpcrtask, logbase = logbase)
        ## REPLOT
        cat(fix, "Fix", del,": Observe data below and plot...\n")
        print(subset(data, Detector == fix))
        plotStdCurve(qPCRdata = data, sampleName = fix, numTechReps, startingAmount, dilutionFactor, numPoints, qpcrtask, ylim=c(15,33), xlim=xlim, col="blue", lwd=3, main = paste0(fix, ", Fix", del), logbase = logbase) ## std curve
        ## needs fixing ?
        needsfixing <- finfo$stats$pass == 0
#         cat("END",needsfixing,"\n")
      }
      notes <- paste0(notes, "\n")
    }
  }
  cat("\n","Summary: \n", notes, "\n")
  return(list(data=data, notes=notes))
}

validateMyPrimers <- function(filename, controls=NA, numTechReps=2, startingAmount=2,dilutionFactor=4,numPoints=4,highToLow=TRUE,conditional=FALSE, qpcrtask="Standard", logbase=10){
  ## read in data - ADD /path/to/clean-data.tsv to the double quotes below
  data <- read.table(filename, header=T, colClasses=c("character", "factor", "factor", "character", "numeric"))
  ## Make Ct variable numeric instead of character
  data$Ct[data$Ct == "Undetermined"] = NA
  class(data$Ct) = "numeric"
  ## Get validation scores
  print("Begin Auto-validation")
  final1 <- autoValidatePrimers(data,controls,numTechReps, startingAmount,dilutionFactor,numPoints,highToLow, qpcrtask=qpcrtask, logbase=logbase)
  print(final1)
  print("Begin Auto-Correction")
  cordata <- autoCorrectPrimers(data, final1$scores, numTechReps, startingAmount, dilutionFactor, numPoints, highToLow, qpcrtask=qpcrtask, logbase=logbase)
#   print("A")
  final2 <- autoValidatePrimers(cordata$data,controls,numTechReps, startingAmount,dilutionFactor,numPoints,highToLow, qpcrtask=qpcrtask, logbase=logbase)
  print(final2)
  print("Finishing up....")
  ## ALL TOGETHER NOW
  prvalscores <- final1$scores$pass+final2$scores$pass
#   print("ass")
#   print(final1$scores$rscores)
#   print(final2$scores$rscores)
  final3 <- cbind(final1$scores$pass,final2$scores$pass,prvalscores,final1$scores$rscores,final2$scores$rscores)
  colnames(final3) <- c("test1","test2","testsum","rel1", "rel2")
  final4 <- cbind(final1$controlscores, final2$controlscores)
print("Finished...")
  return(list(eff=final3, rel=final4, notes=cordata$notes))
}


validateMyPrimers.BioRad <- function(filename, controls=NA, numTechReps=2, startingAmount=2,dilutionFactor=4,numPoints=4,highToLow=TRUE,conditional=FALSE, qpcrtask="Standard", logbase=10){
  ## read in data - ADD /path/to/clean-data.tsv to the double quotes below
  data <- read.table(filename, header=T, as.is = TRUE, colClasses=c("character", "character", "factor", "numeric", "numeric"))
  data <- data[data$Task != "NaN", ]
  levels(data$Detector) <- unique(data$Detector)
  ## Make Ct variable numeric instead of character
  data$Ct[data$Ct == "Undetermined"] = NA
  class(data$Ct) = "numeric"
  ## Get validation scores
  print("Begin Auto-validation")
  final1 <- autoValidatePrimers(data,controls,numTechReps, startingAmount,dilutionFactor,numPoints,highToLow, qpcrtask=qpcrtask, logbase=logbase)
  print(final1)
  print("Begin Auto-Correction")
  cordata <- autoCorrectPrimers(data, final1$scores, numTechReps, startingAmount, dilutionFactor, numPoints, highToLow, qpcrtask=qpcrtask, logbase=logbase)
  #   print("A")
  final2 <- autoValidatePrimers(cordata$data,controls,numTechReps, startingAmount,dilutionFactor,numPoints,highToLow, qpcrtask=qpcrtask, logbase=logbase)
  print(final2)
  print("Finishing up....")
  ## ALL TOGETHER NOW
  prvalscores <- final1$scores$pass+final2$scores$pass
  #   print("ass")
  #   print(final1$scores$rscores)
  #   print(final2$scores$rscores)
  final3 <- cbind(final1$scores$pass,final2$scores$pass,prvalscores,final1$scores$rscores,final2$scores$rscores)
  colnames(final3) <- c("test1","test2","testsum","rel1", "rel2")
  final4 <- cbind(final1$controlscores, final2$controlscores)
  print("Finished...")
  return(list(eff=final3, rel=final4, notes=cordata$notes))
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

plotStdCurve <- function(qPCRdata, sampleName, numberTechReps=2, startingAmount=25, dilutionFactor=10, numPoints=4, qpcrtask="Standard", logbase=10,  highToLow=TRUE, ...){
  ## for plotting 1 std curve
  ## you may be interested in plotting std curves for all primer pairs on same plot (as well as have points of unknown)
  ##logNg <- log(rep(dilutionSeries(startingAmount, dilutionFactor, numPoints), each=numberTechReps),base=logbase)
  logNg <- getlogNg(dilutionFactor,startingAmount,numPoints,highToLow, numTechReps,logbase)
  isNA <- is.na(subset(qPCRdata, Detector == sampleName & qPCRdata$Task == qpcrtask)$Ct)
  Ctvals <- subset(qPCRdata, qPCRdata$Detector==sampleName & qPCRdata$Task==qpcrtask)$Ct[!isNA]
  logNg <- logNg[!isNA]
  par(mfrow=c(2,1))
  plot(logNg, Ctvals, ...)
  fit <- lm(Ctvals ~ logNg)
  lines(logNg, fit$fitted, col="red")
  abline(v=0,col="black")
  grid()
  yint <- paste0("y-int = ", fit$coeff[1])
  slope <- paste0("Slope = ", fit$coeff[2])
  R <- cor(logNg, Ctvals)  
  Rsq <- paste0("R^2 = ", R^2)
  RxnEff <- paste0("Reaction Efficiency = ", reactionEfficiency(fit$coeff[2],logbase=logbase))
  text(median(logNg), 26, slope, cex=0.8)
  text(median(logNg), 24, yint, cex=0.8)
  text(median(logNg), 22, Rsq, cex=0.8)
  text(median(logNg), 20, RxnEff, cex=0.8)
  ## Plot Deltas
  avgCtTable <- avgCTs(qPCRdata = data, qpcrtask = qpcrtask, logbase = logbase, numTechReps = numberTechReps, startingAmount = startingAmount, dilutionFactor = dilutionFactor, numPoints = numPoints)
  avg <- avgCtTable$avgCt[avgCtTable$primerPair == sampleName]
  deltas <- avg[2:length(avg)] - avg[1:(length(avg)-1)]
  DFs <- 2^deltas
  x <- avgCtTable$logNg[avgCtTable$primerPair == sampleName]
  mDF <- max(DFs, na.rm = TRUE)
  plot(x, rep(mDF,length(x)), ylim=c(0,mDF), type="n", xlab="", ylab="DF between dilutions")
  for(i in 2:length(x)){
    segments(x0 = x[i-1], x1 = x[i], y0 = DFs[i-1], y1 = DFs[i-1])
    text(x=mean(c(x[i-1],x[i])), y=DFs[i-1]-1, labels = round(DFs[i-1], digits=3))
  }
  par(mfrow=c(1,1))
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
  #primerPairs <- levels(qPCRdata$Detector)
  primerPairs <- unique(as.character(qPCRdata$Detector))
  if(testSampFirst){
    testSample <- 1:numTechReps
    calibratorSample <- (numTechReps+1):(2*numTechReps) 
  } else {
    calibratorSample <- 1:numTechReps
    testSample <- (numTechReps+1):(2*numTechReps)
  }
  
  testNormalizerCt <- mean(subset(qPCRdata, Detector == normalizer & Task == task)$Ct[testSample], na.rm=TRUE)
  calibratorNormalizerCt <- mean(subset(qPCRdata, Detector == normalizer & Task == task)$Ct[calibratorSample], na.rm=TRUE)
  print(testNormalizerCt)
  print(calibratorNormalizerCt)
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

allFE.multiple.controls <- function(qPCRdata, normalizers, numTechReps=3, task="Unknown", method=ddCtFE, testSampFirst=TRUE){
  fet <- allFE(qPCRdata, normalizers[1], numTechReps, task, method, testSampFirst)
  fe <- data.frame(FE1=fet$FE)
  n <- length(normalizers)
  if(n>=1){
    for(i in 2:n){
    newfet <- allFE(qPCRdata, normalizers[i], numTechReps, task, method, testSampFirst)
    fe[[paste0("FE",i)]] <- newfet$FE
    }
  }
  mu <- apply(fe, 1, mean)
  stdv <- apply(fe, 1, sd)
  fet <- cbind(primerPair=fet$primerPair, fe, mu=mu, stdv=stdv)
  fet
}

allFE.multiple.controls.multiple.reps.error.prop <- function(fetlist){
  n <- length(fetlist)
  mu <- data.frame(mu1=fetlist[[1]]$mu)
  stdv <- data.frame(sd1=fetlist[[1]]$stdv)
  if (n>=2){
    for(i in 2:n){
      mu[[paste0("mu",i)]] <- fetlist[[i]]$mu
      stdv[[paste0("sd",i)]] <- fetlist[[i]]$stdv
    }
  }
  biorepsFE <- data.frame(primerPair=fetlist[[1]]$primerPair, mu=apply(mu, 1, mean), stdv=sqrt(apply(stdv^2, 1, sum)))
  return(biorepsFE)
}

plot.biorepsFE <- function(biorepsFE, ord=NA, barcol="dark cyan", addh1=TRUE, addplot=FALSE, ymin=NA, ymax=NA, bglines=c(5,10,15), plot_type="bar", sqsub=c(1), ylab="FE", LOG2=FALSE, fevar="mu", ...){
  if(LOG2){transform <- log2}else{transform <- identity}
  if(is.na(ymax)){ymax <- max(transform(biorepsFE$mu+biorepsFE$stdv))}
  if(is.na(ymin)){ymin <- min(transform(biorepsFE$mu-biorepsFE$stdv))}
  plotFE(biorepsFE, ord, barcol, addh1, addplot, ymin, ymax, bglines, plot_type, sqsub, ylab, LOG2, fevar, ...)
  n <- nrow(biorepsFE)
  y0 <- transform(biorepsFE$mu-biorepsFE$stdv)
  y1 <- transform(biorepsFE$mu+biorepsFE$stdv)
  segments(x0 = 1:n, x1 = 1:n, y0 = y0, y1 = y1)
}

allstats <- function(qPCRdata, normalizer, numTechReps=3, task="Unknown", method=ddCtFE, testSampFirst=TRUE){
  ## include mean FE, sd FE
  ## include p-values and corrected p-values for a couple tests (e.g. t-test and rank sum)
  ## include mean Ct and sd Ct values -- for each site for each sample
  ## these will be used in a future function that takes in stats tables for bioreps and computes stats from there...
  #primerPairs <- levels(qPCRdata$Detector)
  primerPairs <- unique(as.character(qPCRdata$Detector))
  if(testSampFirst){
    testSample <- 1:numTechReps
    calibratorSample <- (numTechReps+1):(2*numTechReps) 
  } else {
    calibratorSample <- 1:numTechReps
    testSample <- (numTechReps+1):(2*numTechReps)
  }
  testNormalizerCt.mu <- mean(subset(qPCRdata, Detector == normalizer & Task == task)$Ct[testSample], na.rm=TRUE)
  testNormalizerCt.sd <- sd(subset(qPCRdata, Detector == normalizer & Task == task)$Ct[testSample], na.rm=TRUE)  
  calibratorNormalizerCt.mu <- mean(subset(qPCRdata, Detector == normalizer & Task == task)$Ct[calibratorSample], na.rm=TRUE)
  calibratorNormalizerCt.sd <- sd(subset(qPCRdata, Detector == normalizer & Task == task)$Ct[calibratorSample], na.rm=TRUE)
  statTable <- data.frame(primerPair = primerPairs, FE = 0)
  for(primer in primerPairs){
    primerdata <- subset(qPCRdata, Detector == primer & Task == task)
    testTargetCt <- mean(primerdata$Ct[testSample], na.rm=TRUE) 
    calibratorTargetCt <- mean(primerdata$Ct[calibratorSample], na.rm=TRUE)
    FE <- method(sampleSOICt=testTargetCt, sampleNormCt=testNormalizerCt, calibratorSOICt=calibratorTargetCt, calibratorNormCt=calibratorNormalizerCt)
    ##add to df
    statsTable[feTable$primerPair == primer,]$FE <- FE
  }
  statsTable
}


cov <- function(qPCRdata, numTechReps=3, task="Unknown", testSampFirst=TRUE, covcutoff=1, doremoval=FALSE){
  #primerPairs <- levels(qPCRdata$Detector)
  primerPairs <- unique(as.character(qPCRdata$Detector))
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
  N <- c()
  ctvals <- c()
  for(primer in primerPairs){
    primerdata <- subset(qPCRdata, Detector == primer & Task == task)
    Ct <- primerdata$Ct[testSample]
    ctval <- paste(round(Ct, digits = 1), collapse=",")
    ctvals <- c(ctvals, ctval)
    wells <- primerdata$Well[testSample]
    primerpair <- c(primerpair, primer)
    sampletype <- c(sampletype, "test")
    n <- sum(!(is.na(Ct)))
    N <- c(N, n)
    m <- mean(Ct, na.rm = TRUE)
    mu <- c(mu, m)
    s <- sd(Ct, na.rm = TRUE)
    sigma <- c(sigma, s)
    cv <- 100*s/m
    if(n < 2){ cv <- 100} ##Cannot use only 1 data point, so assume large CV
    covs <- c(covs, cv)
    if (cv > covcutoff) {
      pass <- c(pass, 0) 
      if(sum(!(is.na(Ct))) > 2){removal.cand <- wells[which.max(abs(Ct-m))]}
      else {removal.cand <- "NA"}
      if(doremoval && sum(!(is.na(Ct))) > 2){qPCRdata$Ct[qPCRdata$Well == removal.cand] <- NA} # only remove outlier if more than 2 datapoints -- if only 2, they are equally far from mean 
    }
    else {pass <- c(pass, 1); removal.cand <- "-"}
    rmc <- c(rmc, removal.cand)
    Ct <- primerdata$Ct[calibratorSample]
    ctval <- paste(round(Ct, digits = 1), collapse=",")
    ctvals <- c(ctvals, ctval)
    wells <- primerdata$Well[calibratorSample]
    primerpair <- c(primerpair, primer)
    sampletype <- c(sampletype, "calibrator")
    n <- sum(!(is.na(Ct)))
    N <- c(N, n)
    m <- mean(Ct, na.rm = TRUE)
    mu <- c(mu, m)
    s <- sd(Ct, na.rm = TRUE)
    sigma <- c(sigma, s)
    cv <- 100*s/m
    if(n < 2){ cv <- 100} ##Cannot use only 1 data point, so assume large CV
    covs <- c(covs, cv)
    if (cv > covcutoff) {
      pass <- c(pass, 0) 
      if(sum(!(is.na(Ct))) > 2){removal.cand <- wells[which.max(abs(Ct-m))]}
      else {removal.cand <- "NA"}
      if(doremoval && sum(!(is.na(Ct))) > 2){qPCRdata$Ct[qPCRdata$Well == removal.cand] <- NA} # only remove outlier if more than 2 datapoints -- if only 2, they are equally far from mean
    }
    else {pass <- c(pass, 1); removal.cand <- "-"}
    rmc <- c(rmc, removal.cand)
  }
  covtable <- data.frame(primerPair = primerpair, sampleType = sampletype, N=N, ctvals=ctvals, mean = mu, stdev = sigma, cov = covs, pass = pass, rm.suggest = rmc)
  if(doremoval){qPCRdata}
  else {covtable}
}

coefplot <- function(data, covcutoff=2, ylim=c(10,40), ...){
  X <- round(dim(data)[1]/3)
  x <- 0:(X+1)
  y <- seq(10, 35, length.out = length(x))
  plot(x,y,type="n", las=1, ylab="Ct values", xlab="", xaxt="n", ylim=ylim, ...)
  abline(h=c(15,20,25,30), lty=3)
  for(i in seq(3,X*3,3)){
    d <- data$Ct[(i-2):i]
    wells <- data$Well[(i-2):i]
    axis(side=1, at = i/3, labels = FALSE) 
    #lablist <- paste(data$Detector[i], wells[1],wells[2],wells[3])#, las=2, cex.axis=0.5)
    #text(i/3-1.5, par("usr")[3] - 2.5, labels = lablist, srt = 45, pos = 1, xpd = TRUE, cex=0.75)
    lablist <- paste(data$Detector[i], wells[1],wells[2],wells[3], sep = "\n")
    axis(side=1, at = i/3, labels = lablist, las=1, cex.axis=0.75, pos = 6, tick = FALSE) 
    abline(v = i/3, lty=3)
    #points(rep(i/3,3), c(d))
    text(x = rep(i/3,3), y = c(d), labels = wells)
    #   
#     a <- c(i, d, sd(d, na.rm = TRUE), mean(d, na.rm = TRUE))
#     print(a)
    m <- mean(d, na.rm = TRUE)
    s <- sd(d, na.rm = TRUE)
    coef <- 100*s/m
    removal.cand <- wells[which.max(abs(d-m))]
    if(sum(!(is.na(d))) < 2){coef<-100} ## cause it to FAIL
    if(coef > covcutoff){
      y <- 18+rnorm(1,0,2)
      col <- c("black","red","blue", "dark blue", "dark cyan")[i%%5+1]
      text(x = i/3, y, labels = data$Detector[i], cex=0.75, col=col)
      if(sum(!(is.na(d))) <= 2){
        text(x = i/3, y-2, labels = "Too few data points", cex=0.75, col=col)
        text(x = i/3, y-3, labels = "for automated removal.", cex=0.75, col=col)
      }
      else {
        text(x = i/3, y-2, labels = round(coef,digits = 3), cex=0.75, col=col)
        text(x = i/3, y-4, labels = removal.cand, cex=0.75, col=col)
      }
    }
    if (sum(!(is.na(d))) < 2) {
      text(x = i/3, y-6, labels = c(paste0("There is only ", sum(!(is.na(d)))," datapoint.")), cex=0.75, col=col)
      text(x = i/3, y-7, labels = c(paste0("Consider discarding all")), cex=0.75, col=col)
      text(x = i/3, y-8, labels = c(paste0("data for this site.")), cex=0.75, col=col)
    }
  }
}
  

groupQPCRdata <- function(qPCRdata){
  #primerPairs <- levels(qPCRdata$Detector)
  primerPairs <- unique(as.character(qPCRdata$Detector))
  for(primer in primerPairs){
    print(subset(qPCRdata, Detector == primer))
  }
}

avgCTs <- function(qPCRdata, numTechReps=2, startingAmount=25, dilutionFactor=10, numPoints=4, highToLow=TRUE, qpcrtask="Standard", logbase=10){
  # Grab a list of all detectors
  #primerPairs <- levels(qPCRdata$Detector)
  primerPairs <- unique(as.character(qPCRdata$Detector))
  
  # Create vector of concentrations for standard curve 
  Ng <- dilutionSeries(startingAmount=startingAmount, dilutionFactor=dilutionFactor, numPoints=numPoints, highToLow=highToLow)
  logNgLevels <- log(Ng, base=logbase)
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
    isNA <- is.na(subset(qPCRdata, Detector == primer & Task == qpcrtask)$Ct)
    Cts <- subset(qPCRdata, Detector == primer & Task == qpcrtask)$Ct[!isNA]
    dilSer <- logNg[!isNA]
#     print(dilSer)
#     print(logNgLevels)
    for(ngAmt in logNgLevels){
      relevantCts <- dilSer == ngAmt
      avgCt <- mean(Cts[relevantCts], na.rm=TRUE)  #<-- should not need that as NAs removed as a processing step
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

relativeEfficiencyPlot <- function(RelEffTable, rounddigits=4, ...){
  ## takes output of relativeEfficiencyTable
  xVals <- 1:length(RelEffTable$inputNg)
  plot(xVals, RelEffTable$dCt, ylab="dCt", xlab="Input Amount (ng)", col="red", type="p", pch=19, cex=1.5, lwd=2, xaxt="n", xlim=c(0,(length(RelEffTable$inputNg)+1)), ylim=c((min(RelEffTable$dCt)-2), (max(RelEffTable$dCt)+2)), ...)
  axis(side=1, at=xVals, labels=RelEffTable$inputNg)
  grid()
  fit <- lm(RelEffTable$dCt ~ xVals)
  lines(xVals, fit$fitted, col="black", type="l", lwd=2) # type="b"
  slope <- round(fit$coeff[2], digits = rounddigits)
  ## Note that it is not strictly ng vs. dCt -- they are all equally spaced (mapped to integers) despite not really being equally spaced
  ##  This is exactly how done in Life Tech book and produces same results with same data
  ## ideally slope = 0 if both have identical efficiencies 
  ## but anything < 0.1 is acceptable for the ddCt method --- currently this puts out negative slopes though so > -0.1 to 0
  yint <- round(fit$coeff[1], digits = rounddigits)
  line <- paste0("Y = ", slope, "x + ", yint)
  legend("topright", legend=c("dCt",line), fill=c("red", "black"))
  if(slope <= 0.1 & slope >= 0){
    return(paste0("Efficiencies of Target and Normalizer pass consistency test with slope = ", slope))
  }
  else{return(paste0("Efficiencies of Target and Normalizer FAILED consistency test with slope = ", slope))}
}



fullRelativeEfficiencyTable <- function(avgCtTable, normalizer, startingAmount=2, dilutionFactor=4, numPoints=4, highToLow=TRUE, returnfrt=TRUE, plottables=FALSE){
  ## takes in avgCtTable, the output of avgCTs() -- mutliple avgCTs() outputs may even be bound together with rbind() beforehand
  ## normalizer can be any primerPair name (i.e. Detector) -- typical for MYC locus is MSF1-1 or JU-9-1 (e.g. normalizer="JU-9-1")
  primerPairs <- levels(avgCtTable$primerPair)
  normCts <- subset(avgCtTable, primerPair == normalizer)$avgCt
  fullRelTable <- data.frame(primerPair = NULL, inputNg = NULL, normAvgCt = NULL, targetAvgCt = NULL, dCt = NULL)
  for(primer in primerPairs){
    dilSer <- dilutionSeries(startingAmount=startingAmount, dilutionFactor=dilutionFactor, numPoints=numPoints, highToLow=highToLow)
#     print(dilSer)
    targetCts <- subset(avgCtTable, primerPair == primer)$avgCt
#     print(targetCts)
    primerPair <- rep(primer, numPoints)
#     print(primerPair)
# print(normCts)
    newTable <- relativeEfficiencyTable(normalizerAvgCts=normCts, targetAvgCts=targetCts, dilSeries=dilSer)
#     print(newTable)
    if(plottables){relativeEfficiencyPlot(newTable, main=paste0(primerPair[1]," vs. ",normalizer))}
    newTable <- cbind(primerPair, newTable)
    fullRelTable <- rbind(fullRelTable, newTable)
  }
  if(returnfrt){return(fullRelTable)}
}

allRelativeEfficiencyPlots <- function(avgCtTable, normalizer, startingAmount=2, dilutionFactor=4, numPoints=4, highToLow=TRUE, returnfrt=TRUE, plottables=FALSE){
  fullRelativeEfficiencyTable(avgCtTable, normalizer, startingAmount, dilutionFactor, numPoints, highToLow, returnfrt=FALSE, plottables=TRUE)
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

plotAvgCts <- function(avgCtTable, numPoints=4, label=c("1:1", "1:4", "1:16", "1:64"), highlight=c(NA), ...){
  idx <- c(1:(numPoints-1),0)
  #abcd <- matrix(nrow = numPoints, ncol = length(levels(avgCtTable$primerPair)))
  abcd <- matrix(nrow = numPoints, ncol = length(unique(as.character(avgCtTable$primerPair))))
#   print(abcd)
  for(i in 1:numPoints){
#     avgCtTable[(1:length(avgCtTable$avgCt))%%numPoints == i, c(1,3)]
    #a <- as.matrix(avgCtTable[(1:length(avgCtTable$avgCt))%%numPoints == idx[i], c(3)])
    a <- as.matrix(avgCtTable$avgCt[(1:length(avgCtTable$avgCt))%%numPoints == idx[i]])
    dim(a)
#     print(dim(a))
    abcd[i,] <- a
#     print(abcd)
  }
#   print(abcd)
  rand <- rnorm(length(a)*4, 0, 0.05)
  x <- rep(1:numPoints, each=length(a))+rand
  y <- as.vector(t(abcd))
  l <- avgCtTable[(1:length(avgCtTable$avgCt))%%numPoints == 1, c(1)]
  cex=0.6
  plot(x, y, pch=1, cex=cex, xlab="Dilution", ylab="Average Ct", xaxt="n", yaxt="n", las=2, type="n", xlim=c(0.5,4.5), ...)
  text(x,y,labels = l, cex=cex)
  axis(side=1, at=1:numPoints, label=label)
  axis(side=2, at=21:31, label=21:31, las=2)
  if(!(is.na(highlight[1]))){
    sub <- (1:length(avgCtTable$avgCt))%%numPoints == 1
    cols <- c("blue", "red", "green", "orange", "purple", "dark cyan", "pink", "dark blue", "dark green", "dark red", "gold", "silver", "bronze")
    colcnt <- 0
    for (high in highlight){
      colcnt <- colcnt+1
      bools <- avgCtTable$primerPair[sub] %in% c(high)
      text(x[bools], y[bools], labels=high, col=cols[colcnt], cex=cex) 
    }
  }
}

pairwisePlotAvgCts <- function(avgCtTable, numPoints=4, label=c("1:1", "1:4", "1:16", "1:64"), highlight=c(NA), control){
  ## NEED TO PROVIDE CONTROL
  pairs <- as.character(unique(avgCtTable$primerPair))
  n <- length(pairs)-1
  nr <- ceiling(sqrt(n))
  nc <- nr
  par(mfrow=c(nr,nc), mar=c(2.5,4,1,1))
  for (test in pairs){
    if (test != control){
      subTable <- subset(x = avgCtTable, subset = avgCtTable$primerPair %in% c(control, test))
      plotAvgCts(avgCtTable = subTable, numPoints = numPoints, highlight = highlight, main=paste0(test, " vs. ", control))
    }
  }
  par(mfrow=c(1,1))
}

##Currently working on....Jul 2016
plotFE <- function(fetable, ord=NA, barcol="dark cyan", addh1=TRUE, addplot=FALSE, ymin=NA, ymax=NA, bglines=c(5,10,15), plot_type="bar", sqsub=c(1), ylab="FE", LOG2=FALSE, fevar="FE", ...){
  #fetable is output of allFE
  #ord is order to show genes or genomic sites in -- given ordered primer pair names
  if(is.character(barcol)){barcol <- c(barcol)}
  if(length(barcol) == 1){barcol <- rep(barcol[1], length(fetable$primerPair))}
  if(length(sqsub) == 1){sqsub <- rep(sqsub[1], length(fetable$primerPair))}
  if (is.na(ord[1])){ord <- fetable$primerPair}
  #x <- 1:length(ord)
  x <- 0:(length(ord)+1)
  cap<-1
  if(LOG2){
    #fetable$FE <- log2(fetable$FE);
    fetable[[fevar]] <- log2(fetable[[fevar]])
    cap<-0; 
    #if(is.na(ymin)){ymin <- -1*max(abs(fetable$FE))}
    if(is.na(ymin)){ymin <- -1*max(abs(fetable[[fevar]]))}
  } else {
    if(is.na(ymin)){ymin <- 0}
  }
  #if (is.na(ymax[1])){ymax <- max(abs(fetable$FE))}
  if (is.na(ymax[1])){ymax <- max(abs(fetable[[fevar]]))}
  if (!(addplot)){
    plot(x, seq(ymin, ymax, length.out = length(x)), type="n", xlab="", ylab=ylab, las=1, xaxt="n", ...)
    axis(side=1, at = x[2:(length(ord)+1)], labels = ord, las=2, ...)
    abline(h=bglines,lty=3)
  }
  j<-0
  for(i in 1:length(ord)){
    j <- j+1
    if (ord[i] %in% fetable$primerPair){
      # l = i-0.32
      # r = i+0.52
      l = i-0.5
      r = i+0.5
      if(plot_type == "bar"){
        #t = fetable$FE[fetable$primerPair == ord[i]]
        t = fetable[[fevar]][fetable$primerPair == ord[i]]
        b = 0 
      } else if (plot_type == "square"){
        #p = fetable$FE[fetable$primerPair == ord[i]]
        p = fetable[[fevar]][fetable$primerPair == ord[i]]
        t = p+sqsub[j]
        b = p-sqsub[j]
      }
      polygon(x = c(l,l,r,r,l), y = c(b,t,t,b,b), col = barcol[j])
    } 
  }
  if (addh1) {abline(h=cap,lty=3,lwd=2)}
}



read_biorad <- function(filename, samplepresent=TRUE){
  if(samplepresent){data <- read.table(filename, header=T, as.is = TRUE, colClasses=c("character", "character", "factor", "factor", "numeric", "numeric"))}
  else{data <-read.table(filename, header=T, as.is = TRUE, colClasses=c("character", "character", "factor", "numeric", "numeric"))}
  options("scipen"=100, "digits"=10)
  data$Task <- qpcrtask
  data <- data[data$Task != "NaN", ]
  levels(data$Detector) <- unique(data$Detector)
  # Make Ct variable numeric instead of character
  data$Ct[data$Ct == 0] = NA
  data$Ct[data$Ct == 100] = NA ## as part of making clean.tsv I now make NaN Ct values into 100 -- I did not include this line in original pass through this analysis (some results may change slightly)
  class(data$Ct) = "numeric"
  return(data)
}

read_appbio_quantstudio3 <- function(filename){
  data <- read.table(filename, header=T, as.is = TRUE, colClasses=c("character", "character", "factor", "factor", "numeric", "numeric", "numeric", "numeric"))
  return(data)
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

