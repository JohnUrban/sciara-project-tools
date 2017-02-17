## Set working directory for R inside double quotes
setwd("")

## qPCR performed ... DATE
## Dissociation notes -- ADD NOTES
## Amp plot notes -- ADD NOTES
## Data generated from ABI must be processed into a cleaned tab separated text file (file-clean-data.tsv)
## Assumes columns: Well, Detector, Task, Ct, Tm
## Assumes these are in the header.
## Typically does not use "Task" column - so this can be anything.
## "Tm" column can be anything for most purposes as well - though you might want it for Tm information.
## If your columns differ, make sure the file header reflects that and make sure Well, Detector, and Ct columns are present.
## If your columns differ, you will also need to change the colClasses below to reflect it.



## FOR BOTH AUTOMATED AND MANUAL
## Add /path/to/qPCRanalysis.R in double quotes
qpcrscripts <- "qPCRanalysis.R"  ## EXAMPLE -- change as needed
source(qpcrscripts)
## ADD control Detector names into the vector
## Best to include at least two control names here (if you only have target and control loci, then just add both)
controls <- c("JU-control-1", "JU-control-2", "JU-control-3", "JU-control-5") ## EXAMPLE -- change as needed
## ADD /path/to/clean-data.tsv to the double quotes below
filename <- "example-validation-clean-data.tsv" ## EXAMPLE -- change as needed




### COMPLETELY AUTOMATED VERSION ####
validateMyPrimers(filename, controls)

## IMPORTANT NOTE:
## Autocorrecting for conditional passing is still only lightly tested.
## Going through manual version below also recommended.

# CONCLUSIONS:
# Pass:
# all except JU-target-9
# 
# Fail:
# JU-target-9 ((conditionally passes -- see below))


######################### MANUAL ###############################

### MANUALLY GO THROUGH - an example of things you can do#####
## read in data - ADD /path/to/clean-data.tsv to ABOVE filename variable
data <- read.table(filename, header=T, colClasses=c("character", "factor", "factor", "character", "numeric"))

## observe data
head(data)

## check classes
sapply(data, class)

# Make Ct variable numeric instead of character
data$Ct[data$Ct == "Undetermined"] = NA
class(data$Ct) = "numeric"

## check Ct values to show "undetermined" as NA
data$Ct

## confirm classes again
sapply(data, class)

## look at range of Ct values
range(data$Ct, na.rm = TRUE)
h1 <- hist(data$Ct, breaks=20, xlim=c(15,35))

## Look at the 4 Ct bins from histogram with highest counts -- likely reflects Ct for the 4 concentrations in dilution series
## These should generally be in Dilution Factor intervals away from each other (if DF=4, then deltaCt=2; if DF=10, then deltaCT ~= 3.3)
sort(h1$mids[order(h1$counts, decreasing = TRUE)][1:4])


## determine slope, efficiency, and R^2 for all
primerVal <- primerValidation(data, numTechReps=2, startingAmount=2, dilutionFactor=4, numPoints=4, highToLow=TRUE)
primerVal


## RELATIVE EFFICIENCIES
## pick your favorite normalizer Detector for norm variable - put inside quotes
norm = controls[1] ## EXAMPLE - change as needed
## What is the highest ng amount in your series? 
startingAmount=2 ## EXAMPLE - change as needed
## What was your dilution factor in the series?
DF = 4 ## EXAMPLE - change as needed
## How many concentrations are there in the series (typically 4)?
num = 4 ## EXAMPLE - change as needed

# Calculate average CT values, their std devs, and coeff of Variations (plot later)
avgCtTable <- avgCTs(data)

## BEFORE calculating this -- FILL IN ABOVE VALUES
frt <- fullRelativeEfficiencyTable(avgCtTable=avgCtTable, normalizer=norm, startingAmount=startingAmount, dilutionFactor=DF, numPoints=num, highToLow=TRUE)
frt

## What is the range of the copy numbers relative to the normalizer?
range(frt$relCopyNum)

## What loci have >2 fold difference in CN than normalizer?
frt[frt$relCopyNum < 0.5 | frt$relCopyNum > 2,] ## all with >2 fold diff in copy number

## Make the Relative Efficiency Table -- NEED "frt" from ABOVE
rett1 <- relativeEfficiencyTest(fullRelTable=frt)
rett1


## Look at average CTs etc 
avgCtTable
## plot avgCt table -- you can choose to highlight various detectors 
plotAvgCts(avgCtTable, highlight = c("JU-target-9")) # highlights 1
plotAvgCts(avgCtTable, highlight = c("JU-target-9","JU-target-8")) ## highlights 2
## print out avg Cts for each in compact manner
compactAvgCts(frt)


## GET VALUES FOR ALL POSSIBLE NORMALIZERS
all <- rett1[,c(1,3)]
for (control in controls[2:length(controls)]){
  frt <- fullRelativeEfficiencyTable(avgCtTable=avgCtTable, normalizer=control, startingAmount=2, dilutionFactor=4, numPoints=4, highToLow=TRUE)
  rett <- relativeEfficiencyTest(fullRelTable=frt)
  all <- cbind(all, rett[,3])
}
colnames(all) <- c("pair",controls)
rownames(all) <- all[,1]
all

## HOW MANY NORMALIZERS DOES EACH DETECTOR PASS WITH? WHAT PROPORTION? i.e. DOES EACH PASS?
relpass <- apply(X = all[,2:dim(all)[2]], MARGIN = 1, FUN = sum); relpass
## WHAT PROPORTION?
relpassp <- relpass/length(controls); relpassp

## HOW MANY DETECTORS DOES EACH NORMALIZER WORK WITH? i.e. WHAT NORMALIZERS ARE BEST?
normpass <- apply(X = all[,2:dim(all)[2]], MARGIN = 2, FUN = sum); normpass
## WHAT PROPORTION?
normpassp <- normpass/length(primerVal$primerPair); normpassp


## REVIEW
## Efficiency and Relative Efficency together
cbind(primerVal,relpassp)
## Best normalizers?
sort(normpassp, decreasing = TRUE)


##################### FINDING CONTIONAL PASSES #######################################
## CAN THE BAD ONES BE FIXED UP?
## GO THROUGH AND REMOVE BAD DATA POINTS FROM FAILED DETECTORS
## IF THEY PASS AFTER REMOVAL, GIVE THEM A CONDITIONAL PASS AND RE-TEST
## IF THEY STILL FAIL, OPTIONALLY RETEST OR DESIGN NEW PRIMERS
## In this example, only 1 primer needs changing. Copy/paste its template for more as needed.

## JU-target-9 -- EXAMPLE
fix <- "JU-target-9"
## Identify the offending well by looking at data and plotting 
subset(data, Detector == fix)
plotStdCurve(qPCRdata = data, sampleName = fix, numberTechReps = 2, startingAmount = 1, dilutionFactor = 4, ylim=c(15,33), xlim=c(-2,1), col="blue", lwd=3) ## std curve
offender <- "G2"
## NA the offending well
data[data$Detector == fix & data$Well == offender,]$Ct = NA
## Review by looking at the data and plotting again
subset(data, Detector == fix)
plotStdCurve(qPCRdata = data, sampleName = fix, numberTechReps = 2, startingAmount = 1, dilutionFactor = 4, ylim=c(15,33), xlim=c(-2,1), col="blue", lwd=3) ## std curve



## CHECK ALL AGAIN
primerVal2 <- primerValidation(data, numTechReps=2, startingAmount=2, dilutionFactor=4, numPoints=4, highToLow=TRUE)
primerVal2

## Get average CTs ... etc... 
# Calculate average CT values, their std devs, and coeff of Variations (plot later)
avgCtTable2 <- avgCTs(data)
avgCtTable2  ## MAY SEE NAs in STDDEV and COV columns -- if only 1 data point, there is no stdev. If no std dev, then no COV.

## BEFORE calculating this -- FILL IN ABOVE VALUES
frt2 <- fullRelativeEfficiencyTable(avgCtTable=avgCtTable2, normalizer=norm, startingAmount=startingAmount, dilutionFactor=DF, numPoints=num, highToLow=TRUE)
frt2

## What is the range of the copy numbers relative to the normalizer?
range(frt2$relCopyNum)

## What loci have >2 fold difference in CN than normalizer?
frt2[frt2$relCopyNum < 0.5 | frt2$relCopyNum > 2,] ## all with >2 fold diff in copy number

## Make the Relative Efficiency Table -- NEED "frt" from ABOVE
rett2 <- relativeEfficiencyTest(fullRelTable=frt2)
rett2


## Look at average CTs etc 
avgCtTable2
## plot avgCt table -- you can choose to highlight various detectors 
plotAvgCts(avgCtTable2, highlight = c("JU-target-9")) # highlights 1
plotAvgCts(avgCtTable2, highlight = c("JU-target-9","JU-target-8")) ## highlights 2
## print out avg Cts for each in compact manner
compactAvgCts(frt2)

## Check as many more detectors as possible normalizers as you'd like
all2 <- rett2[,c(1,3)]
for (control in controls[2:length(controls)]){
  frt2 <- fullRelativeEfficiencyTable(avgCtTable=avgCtTable2, normalizer=control, startingAmount=2, dilutionFactor=4, numPoints=4, highToLow=TRUE)
  rett <- relativeEfficiencyTest(fullRelTable=frt2)
  all2 <- cbind(all2, rett[,3])
}
all2
colnames(all2) <- c("pair",controls)
rownames(all2) <- all2[,1]
all2

## HOW MANY NORMALIZERS DOES EACH DETECTOR PASS WITH? WHAT PROPORTION? i.e. DOES EACH PASS?
relpass2 <- apply(X = all2[,2:dim(all2)[2]], MARGIN = 1, FUN = sum); relpass2
## WHAT PROPORTION?
relpass2p <- relpass2/length(controls); relpass2p

## HOW MANY DETECTORS DOES EACH NORMALIZER WORK WITH? i.e. WHAT NORMALIZERS ARE BEST?
normpass2 <- apply(X = all2[,2:dim(all2)[2]], MARGIN = 2, FUN = sum); normpass2
## WHAT PROPORTION?
normpass2p <- normpass2/length(primerVal$primerPair); normpass2p


## ALL TOGETHER NOW
prvalscores <- primerVal$pass+primerVal2$pass
final <- cbind(primerVal$pass,primerVal2$pass,prvalscores,relpassp,relpass2p) 
colnames(final) <- c("test1","test2","testsum","rel1", "rel2") #, "norm1", "norm2")
final
cbind(normpassp,normpass2p)

# CONCLUSIONS:
# Pass:
# all except JU-target-9
# 
# 
# Conditional Pass:
# JU-target-9
# 
# 
# Fail:
# None