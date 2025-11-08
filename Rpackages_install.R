# Installing and loading required R packages for Bant_MLVA31_analyzer supporting R scripts

# This script is part of Bant_MLVA31_analyzer pipeline installed with apptainer

if(!require(permute)){
  install.packages("permute", lib="/usr/lib/R/library/")
  library(permute, lib.loc="/usr/lib/R/library/")
}
if(!require(lattice)){
  install.packages("lattice", lib="/usr/lib/R/library/")
  library(lattice, lib.loc="/usr/lib/R/library/")
}
if(!require(vegan)){
  install.packages("vegan", lib="/usr/lib/R/library/")
  library(vegan, lib.loc="/usr/lib/R/library/")
}
if(!require(ape)){
  install.packages("ape", lib="/usr/lib/R/library/")
  library(ape, lib.loc="/usr/lib/R/library/")
}
if(!require(cluster)){
  install.packages("cluster", lib="/usr/lib/R/library/")
  library(cluster, lib.loc="/usr/lib/R/library/")
}