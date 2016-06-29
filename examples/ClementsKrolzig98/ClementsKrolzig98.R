library(rmrs)
setwd("C:\\Users\\chiyahn\\Dropbox\\Work\\June08\\rmrs\\examples\\ClementsKrolzig98")
data.use <- read.csv("RealGDP.csv") # read the data (assuming it's in the same dir)
data.use <- data.use$GDPGROWTH
MRSMLEIndep(y = data.use, M = 3, p = 4)

