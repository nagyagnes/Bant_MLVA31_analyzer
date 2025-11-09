#======================================================================
# Compute Repeat Unit Number from amplicon length for Bacillus anthracis MLVA31
#
# This script is part of Bant_MLVA31_analyzer pipeline
#
# Author: Peter Sály (saly.peter@ecolres.hu)
#
# version3: 2024-05-06
#======================================================================

Sys.setenv(LANG="en_US.UTF-8")
Sys.setenv(LC_ALL="en_US.UTF-8")

###############################################################
### SET WORKING DIRECTORY for "Calculation_template.csv" ###
###############################################################
# Path to working directory are temporarily stored in a txt file named 'workdir.txt'.
# Getting and setting the path:
pathToWD01 <- getwd() #readLines("workdir.txt")
pathToWD01
setwd(pathToWD01)

#############################
### INPUT REFERENCE TABLE ###
#############################
# Input reference table is stored in script directory'.
dir()
refTab <- read.table("Calculation_template.csv",
                sep=";",
                header=T,
                as.is=1)
##fix(refTab)
##str(refTab)

###############################################
### SET WORKING DIRECTORY for "workdir.txt" ###
###############################################
# Path to working directory are temporarily stored in a txt file named 'workdir.txt'.
# Getting and setting the path:
pathToWD02 <- readLines("workdir.txt")
pathToWD02
setwd(pathToWD02)

##################
### INPUT DATA ###
##################
# Get the name of the input txt file
inputFileName <- dir(pattern="MLVA_sizes.txt")
inputFileName

# Read the input txt file
d <-  read.table(inputFileName,
                 sep=" ",
                 header=F,
                 as.is=1)
##fix(d)
##str(d)

# Get the name of the comment file
commentFileName <- dir(pattern="comment.txt")
commentFileName

# Read the comment txt file
commFile <- read.table(commentFileName,
                       sep="",
                       header=F,
                       as.is=1:2)
##fix(commFile)
##str(commFile)

##############################################################
### DATA CUT AND PASTE INTO REFERENCE TABLE (SampleLength) ###
##############################################################
# Dealing with data
#...Check set of items
refTab$RepeatName %in% d$V1 # OK

#...Write sample length into refTab
refTab$SampleLength <- d$V2[match(x=refTab$RepeatName, table=d$V1)]

# Dealing with comments
#...Check set of items
refTab$RepeatName %in% commFile$V1 # OK

#...Write comments into refTab
refTab$Comment <- commFile$V2[match(x=refTab$RepeatName, table=commFile$V1)]

#################################################
### COMPUTE REPEAT UNIT NUMBER (RepeatUnitNr) ###
#################################################
head(refTab)

# Function for computation:
RepeatUnitNumber <- function(tab){
    refTab$RepeatUnitNr <-
        tab$ReferenceUnitNr - ((tab$ReferenceLength-tab$SampleLength)/tab$RepeatLength)
    return(refTab)
    }

# Computation
resTab <- RepeatUnitNumber(refTab)
resTab

# Round RepeatUnitNr to integer
resTab$RepeatUnitNr <- round(resTab$RepeatUnitNr, digits = 1)
resTab

round2 <- function(x, n) {
  posneg = sign(x)
  z = abs(x)*10^n
  z = z + 0.5
  z = trunc(z)
  z = z/10^n
  z * posneg
}

resTab$RepeatUnitNr <- round2(resTab$RepeatUnitNr, 0)
resTab

###############################
### WRITE OUT RESULT OUTPUT ###
###############################
# Create name for the outputfile
filename <- inputFileName
filename
pat = gregexpr("[MLVA]",filename)
outputFileName <- paste0(substr(filename,1,max(pat[[1]])), "_unit.csv")
#outputFileName <- paste0(substr(inputFileName,start=1, stop=7), "_unit.csv")

# Write out
write.table(x=resTab[,c("RepeatName","SampleLength","RepeatUnitNr", "Comment")],
            file=outputFileName,quote=F,
            sep=";", row.names=F)

#======================================================================
# *** EoF ***
