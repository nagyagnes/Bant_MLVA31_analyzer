# Installing and loading required R packages for Bant_MLVA31_analyzer supporting R scripts

# This script is part of Bant_MLVA31_analyzer pipeline installed with apptainer

myPaths <- .libPaths()
myPaths <- c(myPaths, '/usr/lib/R/library/')
.libPaths(myPaths)  # add new path 

if(!require(permute)){
  install.packages("permute")
  library(permute)
}
if(!require(lattice)){
  install.packages("lattice")
  library(lattice)
}
if(!require(vegan)){
  install.packages("vegan")
  library(vegan)
}
if(!require(ape)){
  install.packages("ape")
  library(ape)
}
if(!require(cluster)){
  install.packages("cluster")
  library(cluster)
}