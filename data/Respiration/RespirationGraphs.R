###################################################################
#  Rachel Ferrill
#  6 Jul 2015
#  Respiration Graphs
####################################################################

setwd("~/GitHub/StarvationTraits/data/Respiration")
rm(list=ls())

# Define Functions
sem <- function(x){sd(x)/sqrt(3)}

# Import Data Files
resp1 <- read.csv("20150701_BacterialRespiration_a_RNF_Output.txt")
resp2 <- read.csv("20150701_BacterialRespiration_b_RNF_Output.txt")
resp3 <- read.csv("20150701_BacterialRespiration_c_RNF_Output.txt")


# Seperate Based on Isolate
A7022 <- resp1[1:3,]
A7022_mean <- mean(A7022$Rate..uM.O2.Hr.1.)
A7022_sem <- sem(A7022$Rate..uM.O2.Hr.1.)
B7022 <- resp1[4:6,]
B7022_mean <- mean(B7022$Rate..uM.O2.Hr.1.)
B7022_sem <- sem(B7022$Rate..uM.O2.Hr.1.)
C7022 <- resp1[7:9,]
C7022_mean <-mean(C7022$Rate..uM.O2.Hr.1.)
C7022_sem <- sem(C7022$Rate..uM.O2.Hr.1.)
D7022 <- resp1[10:12,]
D7022_mean <-mean(D7022$Rate..uM.O2.Hr.1.)
D7022_sem <-sem(D7022$Rate..uM.O2.Hr.1.)
A7025 <- resp1[13:15,]
A7025_mean <-mean(A7025$Rate..uM.O2.Hr.1.)
A7025_sem <-sem(A7025$Rate..uM.O2.Hr.1.)
B7025 <- resp1[16:18,]
B7025_mean <-mean(B7025$Rate..uM.O2.Hr.1.)
B7025_sem <-sem(B7025$Rate..uM.O2.Hr.1.)
C7025 <- resp1[19:21,]
C7025_mean <-mean(C7025$Rate..uM.O2.Hr.1.)
C7025_sem <-sem(C7025$Rate..uM.O2.Hr.1.)
mean(C7025[,4])
sem(C7025[,4])




D7025 <- resp2[1:3,]
D7025_mean <-mean(D7025$Rate..uM.O2.Hr.1.)
D7025_sem <-sem(D7025$Rate..uM.O2.Hr.1.)
A7026 <- resp2[4:6,]
A7026_mean <-mean(A7026$Rate..uM.O2.Hr.1.)
A7026_sem <-sem(A7026$Rate..uM.O2.Hr.1.)
B7026 <- resp2[7:9,]
B7026_mean <-mean(B7026$Rate..uM.O2.Hr.1.)
B7026_sem <-sem(B7026$Rate..uM.O2.Hr.1.)
C7026 <- resp2[10:12,]
C7026_mean <-mean(C7026$Rate..uM.O2.Hr.1.)
C7026_sem <-sem(C7026$Rate..uM.O2.Hr.1.)
D7026 <- resp2[13:15,]
D7026_mean <-mean(D7026$Rate..uM.O2.Hr.1.)
D7026_sem <-sem(D7026$Rate..uM.O2.Hr.1.)
A7031 <- resp2[16:18,]
A7031_mean <-mean(A7031$Rate..uM.O2.Hr.1.)
A7031_sem <-sem(A7031$Rate..uM.O2.Hr.1.)
B7031 <- resp2[19:21,]
B7031_mean <-mean(B7031$Rate..uM.O2.Hr.1.)
B7031_sem <-sem(B7031$Rate..uM.O2.Hr.1.)



A7032 <- resp3[1:3,]
A7032_mean <-mean(A7032$Rate..uM.O2.Hr.1.)
A7032_sem <-sem(A7032$Rate..uM.O2.Hr.1.)
B7032 <- resp3[4:6,]
B7032_mean <-mean(B7032$Rate..uM.O2.Hr.1.)
B7032_sem <-sem(B7032$Rate..uM.O2.Hr.1.)
C7032 <- resp3[7:9,]
C7032_mean <-mean(C7032$Rate..uM.O2.Hr.1.)
C7032_sem <-sem(C7032$Rate..uM.O2.Hr.1.)
D7032 <- resp3[10:12,]
D7032_mean <-mean(D7032$Rate..uM.O2.Hr.1.)
D7032_sem <-sem(D7032$Rate..uM.O2.Hr.1.)
A7034 <- resp3[13:15,]
A7034_mean <-mean(A7034$Rate..uM.O2.Hr.1.)
A7034_sem <-sem(A7034$Rate..uM.O2.Hr.1.)

means <- c(A7022_mean, B7022_mean, C7022_mean, D7022_mean, A7025_mean, B7025_mean, C7025_mean, 
           D7025_mean, A7026_mean, B7026_mean, C7026_mean, D7026_mean, A7031_mean, B7031_mean, 
           A7032_mean, B7032_mean, C7032_mean, D7032_mean, A7034_mean)
sems  <- c(A7022_sem, B7022_sem, C7022_sem, D7022_sem, A7025_sem, B7025_sem, C7025_sem, D7025_sem,
           A7026_sem, B7026_sem, C7026_sem, D7026_sem, A7031_sem, B7031_sem, A7032_sem, B7032_sem,
           C7032_sem, D7032_sem, A7034_sem)

bp <-barplot(means), xlab = "Isolate",ylab "Cellular Respiration (uM/Hr)", names.arg = c("A7022", "B7022", "C7022", "D7022", "A7025", "B7025", "C7025", "D7025", "A7026", "B7026", "C7026", "D7026", "A7031", "B7031", "A7032", "B7032", "C7032", "D7032", "A7034")
arrows(x0=bp, y0=means, y1=means-sems, angle=90, length=0.1, 1wd=1)
arrows(x0=bp, y0=means, y1=means+sems, angle=90, length=0.1, 1wd=1)