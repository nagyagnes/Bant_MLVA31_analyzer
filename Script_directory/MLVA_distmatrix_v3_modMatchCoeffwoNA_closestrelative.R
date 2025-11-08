#=========================================================================================================================
# Determining closest relatives of an unknown Bacillus anthracis strain by calculating dissimilarity matrix of MLVA31 data
#
# This script is part of Bant_MLVA31_analyzer pipeline
#
# Author: Peter Sály (saly.peter@ecolres.hu)
#
# version: 2025-11-08
#=========================================================================================================================

# Loading required packages:

library(permute, lib.loc="/usr/lib/R/library/")
library(lattice, lib.loc="/usr/lib/R/library/")
library(vegan, lib.loc="/usr/lib/R/library/")
library(ape, lib.loc="/usr/lib/R/library/")
library(cluster, lib.loc="/usr/lib/R/library/")

################################################
### SET WORKING DIRECTORY for "workdir3.txt" ###
################################################
# Path to working directory are temporarily stored in a txt file named 'workdir3.txt'.
# Getting and setting the path:
pathToWD <- readLines("workdir3.txt")
pathToWD
setwd(pathToWD)

##################
### INPUT DATA ###
##################
# Get the name of the input txt file (filename saved to workdir3.txt)

inputFileName <- readLines("filename.txt")
inputFileName

MLVA31 <- read.table(inputFileName, header = TRUE, sep = ";", as.is = c(1:3))
#fix(MLVA31)
str(MLVA31)
dim(MLVA31) # OK

#################################################################
# Function for calculating pairwise similarity from MLVA31 data #
#################################################################
#-------------------------------------------------------------------------------
# Start of Function:...

# Similarity between strains is calculated excluding NAs (repeat number not determined)
#...Creating similarity coefficient: modified matching coefficient without NA
# modMatch = (A) / (A+B)
# A: the copy number of a given repeat is equal in two strains
# B: the copy number of a given repeat is different in two strains
# C: there is no data (NA) for a given repeat in one of the strains
# D: there is no data (NA) for a given repeat in both strains

# Data matrix: strains (objects) in rows, repeats in columns
# Calculation:
# If the copy number of a given repeat is equal in two strains, the two objects are identical to this repeat, the similarity value is 1
# If the copy number of a given repeat is different in two strains, the two objects are different to this repeat, the similarity value is 0
# If there is no data (NA) for a given repeat in one or both of the strains, the given repeat is excluded from the calculation
# The similarity value for two strains is the ratio of number of identical repeats to number of all repeats

modMatchNA = function(x) {
  # Number of objects (strains)
  nSTRAIN = nrow(x)
  # Repeat numbers
  nREP = ncol(x)
  # Result storage:
  resMat = matrix(nrow=nSTRAIN, ncol=nSTRAIN, dimnames=list(rownames(x), rownames(x)))
  # Loop
  for(i in 1:nSTRAIN) {
    for(j in 1:nSTRAIN) {
      temp = x[c(i,j),]
      # Repeats with NA value:
      nNA = sum( apply(temp,2,function(x) any(is.na(x))) )
      #...Loop for each repeat
      temp2 = vector(mode="logical", length=nREP) # temporary result storage
      for(r in 1:nREP){
        # Excluding NAs
        temp3 = temp[,r, drop=F]
        if((is.na(temp3[1,]) & !is.na(temp3[2,])) | 
           (!is.na(temp3[1,]) & is.na(temp3[2,])) |
           (is.na(temp3[1,]) & is.na(temp3[2,]))) {
           temp2[r] = NA
        }
        
        # If there are no NAs
        if(!is.na(temp3[1,]) & !is.na(temp3[2,])) {
          temp2[r] = ifelse(temp3[1,] == temp3[2,], TRUE, FALSE)
        }
        
      }
      # Modified matching coefficient without NAs
      mN = sum(temp2, na.rm=T) / (nREP-nNA)
      #
      # Writing results
      resMat[i,j] = mN
    }
  }
  #resMat
  return(resMat)
} # EoF
#...End of Function.

###############################################
### Calculating distance matrix of strains ###
###############################################
# Distance = 1 - similarity
dim(MLVA31)
str(MLVA31)
#...NAs
any(is.na(MLVA31))
NAk = which(is.na(MLVA31), arr.ind=T)
NAk
MLVA31[unique(NAk[,1]), unique(NAk[,2])]
head(MLVA31)
#...Names of strains to rownames:
MLVA31cimkekkel = MLVA31[,-1]
rownames(MLVA31cimkekkel) = MLVA31[,1]
head(MLVA31cimkekkel)
#...Distance matrix
dissimMat =   1 - modMatchNA(x=MLVA31cimkekkel)
head(dissimMat)
dim(dissimMat)
#fix(dissimMat)

##############################################
############## Cluster analysis ##############
##############################################

dmat = as.dist(dissimMat)
summary(dmat)
osztalyozas <- hclust(d=dmat, method="average")
str(osztalyozas)

#########################################################
## Determining closest neighbour of an object (strain) ##
#########################################################

MLVA31$Strain
length(MLVA31$Strain)
temp = as.matrix(dmat)
dimnames(temp) = list(MLVA31$Strain,MLVA31$Strain)
head(temp)
rownames(temp)
sampleName <- readLines("samplename.txt")
sampleName
sampleName %in% rownames(temp)
res = as.matrix(sort(temp[sampleName,])[1:21]) # (The first 20 closest neigbours)
res

#####################################
#####CREATE AND SAVE OUTPUT FILE ####
#####################################
#...Rownames to columns
res2 = data.frame(Strain = rownames(res), Dist = res[,1])
res2

#getwd()
write.table(res2, file=paste0(sampleName,"_close_rel.csv"), sep=";", row.names=F)
#======================================================================
# *** EoF ***
