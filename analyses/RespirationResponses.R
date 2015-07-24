###################################################################
#  Rachel Ferrill
#  23 Jul 2015
#  Respiration Responses
####################################################################

setwd("~/GitHub/StarvationTraits")
rm(list=ls())

# Define Functions
sem <- function(x){sd(na.omit(x))/sqrt(length(na.omit(x)))}

# Import Data Files (Output from PreSens Respiration Analysis)
resp1  <- read.csv("./data/Respiration/20150701_BacterialRespiration_a_RNF_Output.txt")
resp2  <- read.csv("./data/Respiration/20150701_BacterialRespiration_b_RNF_Output.txt")
resp3  <- read.csv("./data/Respiration/20150701_BacterialRespiration_c_RNF_Output.txt")
resp4  <- read.csv("./data/Respiration/20150707_BacterialRespiration_d_RNF_Output.txt")
resp5  <- read.csv("./data/Respiration/20150707_BacterialRespiration_e_RNF_Output.txt")
resp6  <- read.csv("./data/Respiration/20150710_BacterialRespiration_f_RNF_Output.txt")
resp7  <- read.csv("./data/Respiration/20150710_BacterialRespiration_g_RNF_Output.txt")
resp8  <- read.csv("./data/Respiration/20150714_BacterialRespiration_h_RNF_Output.txt")
resp9  <- read.csv("./data/Respiration/20150715_BacterialRespiration_i_RNF_Output.txt")
resp10 <- read.csv("./data/Respiration/20150715_BacterialRespiration_j_RNF_Output.txt")
resp11 <- read.csv("./data/Respiration/20150716_BacterialRespiration_k_RNF_Output.txt")
resp12 <- read.csv("./data/Respiration/20150716_BacterialRespiration_l_RNF_Output.txt")

resp.data <- rbind(resp1, resp2, resp3, resp3, resp4, resp5, resp6, resp7, resp8,
               resp9, resp10, resp11, resp12)

colnames(resp.data) <- c("Sample", "Start", "End", "Rate", "R2", "P")

# Remove Blanks
resp.data <- resp.data[resp.data$Sample != "blank" & resp.data$Sample != "Blank", ]

for (i in 1:dim(resp.data)[1]){
  resp.data$organism[i] <- strsplit(as.character(resp.data$Sample), "rep")[[i]][1]
  resp.data$rep[i] <- strsplit(as.character(resp.data$Sample), "rep")[[i]][2]
}

# Import Cell Count Data
counts <- read.csv("./data/Respiration/RespirationCounts.txt", header=T)
counts$conc <- counts$NumColonies * 10^(-counts$plate) * 10

# Create Data Frame
resp <- as.data.frame(matrix(NA, dim(counts)[1],12))
colnames(resp) <- c("ID", "Organism", "Evol", "Conc", "Rep1_raw", "Rep2_raw",
                         "Rep3_raw", "Rep1_cor", "Rep2_cor", "Rep3_cor",
                         "Resp_avg", "Resp_sem")
resp$ID <- counts$ID
resp$Conc <- counts$conc

for (i in 1:dim(resp)[1]){
  resp$Organism[i] <- regmatches(resp$ID[i], gregexpr(".{3}", resp$ID[i]))[[1]][1]
  if (nchar(as.character(resp$ID[i])) == 3) {
    resp$Evol[i] <- "Ancestor"
  }else{
    resp$Evol[i] <- "Derived"
    }
}

for (i in 1:dim(resp)[1]){
  rep1 <- resp.data[resp.data$rep == 1, ]
  rep2 <- resp.data[resp.data$rep == 2, ]
  rep3 <- resp.data[resp.data$rep == 3, ]
  resp$Rep1_raw[i] <- rep1$Rate[match(resp$ID[i], rep1$organism)]
  resp$Rep2_raw[i] <- rep2$Rate[match(resp$ID[i], rep2$organism)]
  resp$Rep3_raw[i] <- rep3$Rate[match(resp$ID[i], rep3$organism)]
}


# change Units of Respiration to pico molar
resp$Rep1_raw <- resp$Rep1_raw * 10^6
resp$Rep2_raw <- resp$Rep2_raw * 10^6
resp$Rep3_raw <- resp$Rep3_raw * 10^6

resp$Rep1_cor <- resp$Rep1_raw / resp$Conc
resp$Rep2_cor <- resp$Rep2_raw / resp$Conc
resp$Rep3_cor <- resp$Rep3_raw / resp$Conc

# Calculate Average and SEM
resp$Resp_avg <- round(apply(resp[,8:10], 1, mean), 3)
resp$Resp_sem <- round(apply(resp[,8:10], 1, sem), 3)


kbs701 <- resp[resp$Organism == "701",]
kbs702 <- resp[resp$Organism == "702",]

# Combine means and sems based on isolate
means701 <-c(A701_mean, A7011_mean, B7011_mean, C7011_mean, D7011_mean, A7011b_mean, B7011b_mean,
           C7011b_mean, D7011b_mean, A7013_mean, B7013_mean, C7013_mean, D7013_mean)
means702 <-c(A702_mean, A7022_mean, B7022_mean, C7022_mean, D7022_mean, A7025_mean, B7025_mean,
            C7025_mean, D7025_mean, A7026_mean, B7026_mean, C7026_mean, D7026_mean)
means703 <-c(A703_mean, A7031_mean, B7031_mean, C7031_mean, D7031_mean, A7032_mean, B7032_mean,
             C7032_mean, D7032_mean, A7034_mean, B7034_mean, C7034_mean, D7034_mean)
means710 <-c(A710_mean, A7101_mean, B7101_mean, C7101_mean, D7101_mean, A7102_mean, B7102_mean,
             C7102_mean, D7102_mean, A7103_mean, B7103_mean, C7103_mean, D7103_mean)
means723 <-c(A723_mean, A7231_mean, B7231_mean, C7231_mean, D7231_mean, A7232_mean, B7232_mean,
             C7232_mean, D7232_mean, A7233_mean, B7233_mean, C7233_mean, D7233_mean)
means724 <-c(A724_mean, A7241_mean, B7241_mean, C7241_mean, D7241_mean, A7242_mean, B7242_mean,
             C7242_mean, D7242_mean, A7243_mean, B7243_mean, C7243_mean, D7243_mean)

sems701 <-c(A701_sem, A7011_sem, B7011_sem, C7011_sem, D7011_sem, A7011b_sem, B7011b_sem,
           C7011b_sem, D7011b_sem, A7013_sem, B7013_sem, C7013_sem, D7013_sem)
sems702 <-c(A702_sem, A7022_sem, B7022_sem, C7022_sem, D7022_sem, A7025_sem, B7025_sem, C7025_sem,
            D7025_sem,A7026_sem, B7026_sem, C7026_sem, D7026_sem)
sems703 <-c(A703_sem, A7031_sem, B7031_sem, C7031_sem, D7031_sem, A7032_sem, B7032_sem, C7032_sem,
            D7032_sem, A7034_sem, B7034_sem, C7034_sem, D7034_sem)
sems710 <-c(A710_sem, A7101_sem, B7101_sem, C7101_sem, D7101_sem, A7102_sem, B7102_sem,
            C7102_sem, D7102_sem, A7103_sem, B7103_sem, C7103_sem, D7103_sem)
sems723 <-c(A723_sem, A7231_sem, B7231_sem, C7231_sem, D7231_sem, A7232_sem, B7232_sem,
            C7232_sem, D7232_sem, A7233_sem, B7233_sem, C7233_sem, D7233_sem)
sems724 <-c(A724_sem, A7241_sem, B7241_sem, C7241_sem, D7241_sem, A7242_sem, B7242_sem,
            C7242_sem, D7242_sem, A7243_sem, B7243_sem, C7243_sem, D7243_sem)

# Set Default Plot Parameters
par(mar=c(6, 5, 1, 1) + 0.1)

# Plot Respiration Responses by Organism
bp701 <-barplot(kbs701$Resp_avg, las=2, ylim = c(0, 1.2*max(kbs701$Resp_avg+kbs701$Resp_sem)),
                ylab=expression(paste("Respiration (pM O"^2," Hr"^-1,")")),
                las = 2, names.arg=c("A701", "A7011", "B7011", "C7011", "D7011",
                "A7011b", "B7011b", "C7011b", "D7011b", "A7013", "B7013",
                "C7013", "D7013"),
                col = c("black", rep("gray48", 4), rep("gray73", 4), rep("gray98", 4)))
mtext("Isolate", side=1, line = 4.5)
arrows(x0=bp701, y0=kbs701$Resp_avg, y1=means701-sems701, angle=90, length=0.1, lwd=1)
arrows(x0=bp701, y0=means701, y1=means701+sems701, angle=90, length=0.1, lwd=1)
arrows(x0=bp701[1], y0=means701[1], y1=means701[1]-sems701[1], angle=90,
       length=0.1, lwd=1, col = "white")
legend("topright", legend = c("Ancestor", "Tube 1", "Tube 1b", "Tube 3"),
       fill=c("Black", "gray48", "gray73","gray98"), bty="n")

bp702 <-barplot(means702, las=2, ylim = c(0, 1.6*max(means702+sems702)),
                ylab=expression(paste("Respiration (",mu,"M O"^2," Hr"^-1,")")),
                las = 2, names.arg=c("A702", "A7022", "B7022", "C7022", "D7022",
                                     "A7025", "B7025", "C7025", "D7025", "A7026",
                                     "B7026", "C7026", "D7026"),
                col = c("black", rep("gray48", 4), rep("gray73", 4), rep("gray98", 4)))
mtext("Isolate", side=1, line = 4.5)
arrows(x0=bp702, y0=means702, y1=means702-sems702, angle=90, length=0.1, lwd=1)
arrows(x0=bp702, y0=means702, y1=means702+sems702, angle=90, length=0.1, lwd=1)
arrows(x0=bp702[1], y0=means702[1], y1=means702[1]-sems702[1], angle=90,
       length=0.1, lwd=1, col = "white")
legend("topright", legend = c("Ancestor", "Tube 1", "Tube 2", "Tube 3"),
       fill=c("Black", "gray48", "gray73","gray98"), bty="n")


bp703 <-barplot(means703, las=2, ylim = c(0, 1.7*max(means703+sems703)),
                ylab=expression(paste("Respiration (",mu,"M O"^2," Hr"^-1,")")),
                las = 2, names.arg=c("A703", "A7031", "B7031", "C7031", "D7031",
                                     "A7032", "B7032", "C7032", "D7032", "A7034",
                                     "B7034", "C7034", "D7034"),
                col = c("black", rep("gray48", 4), rep("gray73", 4), rep("gray98", 4)))
mtext("Isolate", side=1, line = 4.5)
arrows(x0=bp703, y0=means703, y1=means703-sems703, angle=90, length=0.1, lwd=1)
arrows(x0=bp703, y0=means703, y1=means703+sems703, angle=90, length=0.1, lwd=1)
arrows(x0=bp703[1], y0=means703[1], y1=means703[1]-sems703[1], angle=90,
       length=0.1, lwd=1, col = "white")
legend("topright", legend = c("Ancestor", "Tube 1", "Tube 2", "Tube 4"),
       fill=c("Black", "gray48", "gray73","gray98"), bty="n")

bp710 <-barplot(means710, las=2, ylim = c(0, 1.55*max(means710+sems710)),
                ylab=expression(paste("Respiration (",mu,"M O"^2," Hr"^-1,")")),
                las = 2, names.arg=c("A710", "A7101", "B7101", "C7101", "D7101",
                                     "A7102", "B7102", "C7102", "D7102", "A7103",
                                     "B7103", "C7103", "D7103"),
                col = c("black", rep("gray48", 4), rep("gray73", 4), rep("gray98", 4)))
mtext("Isolate", side=1, line = 4.5)
arrows(x0=bp710, y0=means710, y1=means710-sems710, angle=90, length=0.1, lwd=1)
arrows(x0=bp710, y0=means710, y1=means710+sems710, angle=90, length=0.1, lwd=1)
arrows(x0=bp710[1], y0=means710[1], y1=means710[1]-sems710[1], angle=90,
       length=0.1, lwd=1, col = "white")
legend("topright", legend = c("Ancestor", "Tube 1", "Tube 2", "Tube 3"),
       fill=c("Black", "gray48", "gray73","gray98"), bty="n")

bp723 <-barplot(means723, las=2, ylim = c(0, 1.65*max(means723+sems723)),
                ylab=expression(paste("Respiration (",mu,"M O"^2," Hr"^-1,")")),
                las = 2, names.arg=c("A723", "A7231", "B7231", "C7231", "D7231",
                                     "A7232", "B7232", "C7232", "D7232", "A7233",
                                     "B7233", "C7233", "D7233"),
                col = c("black", rep("gray48", 4), rep("gray73", 4), rep("gray98", 4)))
mtext("Isolate", side=1, line = 4.5)
arrows(x0=bp723, y0=means723, y1=means723-sems723, angle=90, length=0.1, lwd=1)
arrows(x0=bp723, y0=means723, y1=means723+sems723, angle=90, length=0.1, lwd=1)
arrows(x0=bp723[1], y0=means723[1], y1=means723[1]-sems723[1], angle=90,
       length=0.1, lwd=1, col = "white")
legend("topright", legend = c("Ancestor", "Tube 1", "Tube 2", "Tube 3"),
       fill=c("Black", "gray48", "gray73","gray98"), bty="n")

bp724 <-barplot(means724, las=2, ylim = c(0, 1.2*max(means724+sems724)),
                ylab=expression(paste("Respiration (",mu,"M O"^2," Hr"^-1,")")),
                las = 2, names.arg=c("A724", "A7241", "B7241", "C7241", "D7241",
                                     "A7242", "B7242", "C7242", "D7242", "A7243",
                                     "B7243", "C7243", "D7243"),
                col = c("black", rep("gray48", 4), rep("gray73", 4), rep("gray98", 4)))
mtext("Isolate", side=1, line = 4.5)
arrows(x0=bp724, y0=means724, y1=means724-sems724, angle=90, length=0.1, lwd=1)
arrows(x0=bp724, y0=means724, y1=means724+sems724, angle=90, length=0.1, lwd=1)
arrows(x0=bp724[1], y0=means724[1], y1=means724[1]-sems724[1], angle=90,
       length=0.1, lwd=1, col = "white")
legend("topright", legend = c("Ancestor", "Tube 1", "Tube 2", "Tube 3"),
       fill=c("Black", "gray48", "gray73","gray98"), bty="n")




# Statistics
