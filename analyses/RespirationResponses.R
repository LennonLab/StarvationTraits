###################################################################
#  Rachel Ferrill
#  24 Jul 2015
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
resp <- as.data.frame(matrix(NA, dim(counts)[1],13))
colnames(resp) <- c("ID", "Organism", "Evol", "Tube", "Conc", "Rep1_raw", "Rep2_raw",
                         "Rep3_raw", "Rep1_cor", "Rep2_cor", "Rep3_cor",
                         "Resp_avg", "Resp_sem")
resp$ID <- counts$ID
resp$Conc <- counts$conc

for (i in 1:dim(resp)[1]){
  resp$Organism[i] <- regmatches(resp$ID[i], gregexpr(".{3}", resp$ID[i]))[[1]][1]
  if (nchar(as.character(resp$ID[i])) == 3) {
    resp$Evol[i] <- "Ancestor"
    resp$Tube[i] <- "Ancestor"
  }else{
    resp$Evol[i] <- "Derived"
    resp$Tube[i] <- sub("^...(.).*", "\\1", resp$ID[i])[[1]][1]
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
resp$Resp_avg <- round(apply(resp[,9:11], 1, mean), 3)
resp$Resp_sem <- round(apply(resp[,9:11], 1, sem), 3)


kbs701 <- resp[resp$Organism == "701",]
kbs702 <- resp[resp$Organism == "702",]
kbs703 <- resp[resp$Organism == "703",]
kbs710 <- resp[resp$Organism == "710",]
kbs723 <- resp[resp$Organism == "723",]
kbs724 <- resp[resp$Organism == "724",]



# Set Default Plot Parameters
par(mar=c(6, 5, 1, 1) + 0.1)

# Plot Respiration Responses by Organism
bp701 <-barplot(kbs701$Resp_avg, las=2, ylim = c(0, 1.2*max(kbs701$Resp_avg+kbs701$Resp_sem)),
                ylab=expression(paste("Respiration (pM O"^2," Hr"^-1,")")),
                las = 2, names.arg=c("701", "7011A", "7011B", "7011C", "7011D",
                "7011bA", "7011bB", "7011bC", "7011bD", "7013A", "7013B",
                "7013C", "7013D"),
                col = c("black", rep("gray48", 4), rep("gray73", 4), rep("gray98", 4)))
mtext("Isolate", side=1, line = 4.5)
arrows(x0=bp701, y0=kbs701$Resp_avg, y1=kbs701$Resp_avg-kbs701$Resp_sem, angle=90, length=0.1, lwd=1)
arrows(x0=bp701, y0=kbs701$Resp_avg, y1=kbs701$Resp_avg+kbs701$Resp_sem, angle=90, length=0.1, lwd=1)
arrows(x0=bp701[1], y0=kbs701$Resp_avg[1], y1=kbs701$Resp_avg[1]-kbs701$Resp_sem[1], angle=90,
       length=0.1, lwd=1, col = "white")
legend("topright", legend = c("Ancestor", "Tube 1", "Tube 1b", "Tube 3"),
       fill=c("Black", "gray48", "gray73","gray98"), bty="n")

bp702 <-barplot(kbs702$Resp_avg, las=2, ylim = c(0, 1.6*max(kbs702$Resp_avg+kbs702$Resp_sem)),
                ylab=expression(paste("Respiration (pM O"^2," Hr"^-1,")")),
                las = 2, names.arg=c("702", "7022A", "7022B", "7022C", "7022D",
                                     "7025A", "7025B", "7025C", "7025D", "7026A",
                                     "7026B", "7026C", "7026D"),
                col = c("black", rep("gray48", 4), rep("gray73", 4), rep("gray98", 4)))
mtext("Isolate", side=1, line = 4.5)
arrows(x0=bp702, y0=kbs702$Resp_avg, y1=kbs702$Resp_avg-kbs702$Resp_sem, angle=90, length=0.1, lwd=1)
arrows(x0=bp702, y0=kbs702$Resp_avg, y1=kbs702$Resp_avg+kbs702$Resp_sem, angle=90, length=0.1, lwd=1)
arrows(x0=bp702[1], y0=kbs702$Resp_avg[1], y1=kbs702$Resp_avg[1]-kbs702$Resp_sem[1], angle=90,
       length=0.1, lwd=1, col = "white")
legend("topright", legend = c("Ancestor", "Tube 1", "Tube 2", "Tube 3"),
       fill=c("Black", "gray48", "gray73","gray98"), bty="n")


bp703 <-barplot(kbs703$Resp_avg, las=2, ylim = c(0, 1.7*max(kbs703$Resp_avg+kbs703$Resp_sem)),
                ylab=expression(paste("Respiration (pM O"^2," Hr"^-1,")")),
                las = 2, names.arg=c("703", "7031A", "7031B", "7031C", "7031D",
                                     "7032A", "7032B", "7032C", "7032D", "7034A",
                                     "7034B", "7034C", "7034D"),
                col = c("black", rep("gray48", 4), rep("gray73", 4), rep("gray98", 4)))
mtext("Isolate", side=1, line = 4.5)
arrows(x0=bp703, y0=kbs703$Resp_avg, y1=kbs703$Resp_avg-kbs703$Resp_sem, angle=90, length=0.1, lwd=1)
arrows(x0=bp703, y0=kbs703$Resp_avg, y1=kbs703$Resp_avg+kbs703$Resp_sem, angle=90, length=0.1, lwd=1)
arrows(x0=bp703[1], y0=kbs703$Resp_avg[1], y1=kbs703$Resp_avg[1]-kbs703$Resp_sem[1], angle=90,
       length=0.1, lwd=1, col = "white")
legend("topright", legend = c("Ancestor", "Tube 1", "Tube 2", "Tube 4"),
       fill=c("Black", "gray48", "gray73","gray98"), bty="n")

bp710 <-barplot(kbs710$Resp_avg, las=2, ylim = c(0, 1.55*max(kbs710$Resp_avg+kbs710$Resp_sem)),
                ylab=expression(paste("Respiration (pM O"^2," Hr"^-1,")")),
                las = 2, names.arg=c("710", "7101A", "7101B", "7101C", "7101D",
                                     "7102A", "7102B", "7102C", "7102D", "7103A",
                                     "7103B", "7103C", "7103D"),
                col = c("black", rep("gray48", 4), rep("gray73", 4), rep("gray98", 4)))
mtext("Isolate", side=1, line = 4.5)
arrows(x0=bp710, y0=kbs710$Resp_avg, y1=kbs710$Resp_avg-kbs710$Resp_sem, angle=90, length=0.1, lwd=1)
arrows(x0=bp710, y0=kbs710$Resp_avg, y1=kbs710$Resp_avg+kbs710$Resp_sem, angle=90, length=0.1, lwd=1)
arrows(x0=bp710[1], y0=kbs710$Resp_avg[1], y1=kbs710$Resp_avg[1]-kbs710$Resp_sem[1], angle=90,
       length=0.1, lwd=1, col = "white")
legend("topright", legend = c("Ancestor", "Tube 1", "Tube 2", "Tube 3"),
       fill=c("Black", "gray48", "gray73","gray98"), bty="n")

bp723 <-barplot(kbs723$Resp_avg, las=2, ylim = c(0, 1.65*max(kbs723$Resp_avg+kbs723$Resp_sem)),
                ylab=expression(paste("Respiration (pM O"^2," Hr"^-1,")")),
                las = 2, names.arg=c("723", "7231A", "7231B", "7231C", "7231D",
                                     "7232A", "7232B", "7232C", "7232D", "7233A",
                                     "7233B", "7233C", "7233D"),
                col = c("black", rep("gray48", 4), rep("gray73", 4), rep("gray98", 4)))
mtext("Isolate", side=1, line = 4.5)
arrows(x0=bp723, y0=kbs723$Resp_avg, y1=kbs723$Resp_avg-kbs723$Resp_sem, angle=90, length=0.1, lwd=1)
arrows(x0=bp723, y0=kbs723$Resp_avg, y1=kbs723$Resp_avg+kbs723$Resp_sem, angle=90, length=0.1, lwd=1)
arrows(x0=bp723[1], y0=kbs723$Resp_avg[1], y1=kbs723$Resp_avg[1]-kbs723$Resp_sem[1], angle=90,
       length=0.1, lwd=1, col = "white")
legend("topright", legend = c("Ancestor", "Tube 1", "Tube 2", "Tube 3"),
       fill=c("Black", "gray48", "gray73","gray98"), bty="n")

bp724 <-barplot(kbs724$Resp_avg, las=2, ylim = c(0, 1.2*max(kbs724$Resp_avg+kbs724$Resp_sem)),
                ylab=expression(paste("Respiration (pM O"^2," Hr"^-1,")")),
                las = 2, names.arg=c("724", "7241A", "7241B", "7241C", "7241D",
                                     "7242A", "7242B", "7242C", "7242D", "7243A",
                                     "7243B", "7243C", "7243D"),
                col = c("black", rep("gray48", 4), rep("gray73", 4), rep("gray98", 4)))
mtext("Isolate", side=1, line = 4.5)
arrows(x0=bp724, y0=kbs724$Resp_avg, y1=kbs724$Resp_avg-kbs724$Resp_sem, angle=90, length=0.1, lwd=1)
arrows(x0=bp724, y0=kbs724$Resp_avg, y1=kbs724$Resp_avg+kbs724$Resp_sem, angle=90, length=0.1, lwd=1)
arrows(x0=bp724[1], y0=kbs724$Resp_avg[1], y1=kbs724$Resp_avg[1]-kbs724$Resp_sem[1], angle=90,
       length=0.1, lwd=1, col = "white")
legend("topright", legend = c("Ancestor", "Tube 1", "Tube 2", "Tube 3"),
       fill=c("Black", "gray48", "gray73","gray98"), bty="n")




# Statistics
m.701 <- melt(kbs701[,c(1:4, 9:11)])
fit.701 <- aov(value ~ Evol + Tube + ID, data = m.701)
summary(fit.701)
TukeyHSD(fit.701)
