###################################################################
#  Rachel Ferrill
#  6 Jul 2015
#  Respiration Graphs
####################################################################

setwd("~/GitHub/StarvationTraits")
rm(list=ls())

# Define Functions
sem <- function(x){sd(x)/sqrt(3)}

# Import Data Files (Output from GrowthCurve Analysis)
plate1 <- read.csv("./data/GrowthCurves/output/plate1.txt", header=T)
plate2 <- read.csv("./data/GrowthCurves/output/plate2.txt", header=T)
plate3 <- read.csv("./data/GrowthCurves/output/plate3.txt", header=T)
plate4 <- read.csv("./data/GrowthCurves/output/plate4.txt", header=T)
plate5 <- read.csv("./data/GrowthCurves/output/plate5.txt", header=T)
plate6 <- read.csv("./data/GrowthCurves/output/plate6.txt", header=T)

# Import names file
st.names <- read.csv("./data/GrowthCurves/OrganismNames.txt", header=T)

# Merge into dataframe
growth.data <- as.data.frame(matrix(NA, dim(st.names)[1], 4))
colnames(growth.data) <- c("ID", "Isolate", "Rep", "G.Rate")
growth.data[,1] <- st.names$ID
for (i in 1:dim(growth.data)[1]){
  growth.data[i,2] <- strsplit(as.character(st.names$ID), "rep")[[i]][1]
  growth.data[i,3] <- strsplit(as.character(st.names$ID), "rep")[[i]][2]
}
for (i in 1:dim(growth.data)[1]){
  plate <- get(as.character(st.names$Plate[i]))
  well  <- as.character(st.names$Well[i])
  if (length(which(plate$Curve == well)) == 1){
  growth.data[i,4] <- plate$umax[which(plate$Curve == well)]}
  else {growth.data[i,4] <- "NA"}
}

