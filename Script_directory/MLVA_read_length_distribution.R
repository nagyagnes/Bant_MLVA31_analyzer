#======================================================================
# *** Computing descriptive statistics of the length distribution ***
#           of the amplicons created by in silico PCR
#     Aim: Identifying and removing amplicons with outlier length value
#
# This script is part of Bant_MLVA31_analyzer pipeline
#
# Author: Peter Sály (saly.peter@ecolres.hu)
#
# version: 2022-09-04
#======================================================================

#############################
### SET WORKING DIRECTORY ###
#############################
pathToWD <- readLines("workdir2.txt")
pathToWD
setwd(pathToWD)

#############################
######### INPUT DATA ########
#############################
# Get the name of the input txt file
dir(pattern="workfile.txt")
dataFileName <- readLines("workfile.txt")
dataFileName

############################################
########### READ THE INPUT FILE ############
############################################
dir(pattern=dataFileName)
d <- read.table(dataFileName, header=F)
#fix(d)
#str(d)
# Set colname:
colnames(d)[1] <- "readLength"
head(d) # OK

##################################################################
#### ANALYSIS OF STATISTICAL DISTRIBUTION OF AMPLICON LENGTHS ####
##################################################################
## summary(d)
## stripchart(d$readLength, method="jitter", vertical=T)
## boxplot(d$readLength)
bp <- boxplot(d$readLength, plot=F)

# Number of outliers in the data:
nOutliers <- length(bp$out)
#nOutliers

##########################################
############# EXCLUDE OUTLIERS ###########
##########################################
# Removing outliers
# If there are no outliers, data not change:
if(length(bp$out)>0) { # There are outliers
    outliers <- which(d$readLength %in%  bp$out)
    outliers # Identified outlier(s)
    d2 <- d[-outliers,,drop=F] # Data without outliers
    } else {
        d2 <- d } # Original data
# Statistical distribution of output data:
##summary(d2)
##boxplot(d2$readLength)

################################################
################### OUTPUT #####################
##### DESCRIPTIVE STATISTICS OF OUTPUT DATA ####
################################################
resTab <- data.frame("Min"    = numeric(length=1),
                     "Max"    = numeric(length=1),
                     "Median" = numeric(length=1),
                     "Mean"   = numeric(length=1))
resTab
resTab[,"Min"]    <- min(d2$readLength, na.rm=T)
resTab[,"Max"]    <- max(d2$readLength, na.rm=T)
resTab[,"Median"] <- median(d2$readLength, na.rm=T)
resTab[,"Mean"]   <- mean(d2$readLength, na.rm=T)
resTab

######################
## OUTPUT FILE NAME ##
######################
#dataFileName
temp <- strsplit(dataFileName, split = "_ref_unit.sizes.txt")
#temp
strainName  <- temp[[1]]
outFileName <- paste0(strainName, "_ref_unit.stat.txt")
#outFileName

############################
##### SAVE OUTPUT FILE #####
############################
getwd()
write.table(resTab, file=outFileName, sep=";", row.names=F, quote=F)

#======================================================================
# *** EoF ***
