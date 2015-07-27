###################################################################
#  Rachel Ferrill
#  23 Jul 2015
#  Growth Rate Responses
####################################################################

setwd("~/GitHub/StarvationTraits")
rm(list=ls())

# Define Functions
sem <- function(x){sd(na.omit(x))/sqrt(length(na.omit(x)))}

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

growth.data$G.Rate <- as.numeric(growth.data$G.Rate)

growth <- as.data.frame(matrix(NA, length(levels(as.factor(growth.data$Isolate))), 9))
colnames(growth) <- c("ID", "Organism", "Evol", "Tube", "Rep1", "Rep2", "Rep3",
                      "G.Avg", "G.SEM")

growth$ID <- levels(as.factor(growth.data$Isolate))

for (i in 1:dim(growth)[1]){
  growth$Organism[i] <- regmatches(growth$ID[i], gregexpr(".{3}", growth$ID[i]))[[1]][1]
  if (nchar(as.character(growth$ID[i])) == 3) {
    growth$Evol[i] <- "Ancestor"
    growth$Tube[i] <- "Ancestor"
  }else {
    if (nchar(as.character(growth$ID[i])) == 5) {
      growth$Evol[i] <- "Derived"
      growth$Tube[i] <- sub("^...(.).*", "\\1", growth$ID[i])[[1]][1]
    } else {
      growth$Evol[i] <- "Derived"
      growth$Tube[i] <- sub("^...(..).*", "\\1", growth$ID[i])[[1]][1]
    }}
}

for (i in 1:dim(growth)[1]){
  rep1 <- growth.data[growth.data$Rep == 1, ]
  rep2 <- growth.data[growth.data$Rep == 2, ]
  rep3 <- growth.data[growth.data$Rep == 3, ]
  growth$Rep1[i] <- rep1$G.Rate[match(growth$ID[i], rep1$Isolate)]
  growth$Rep2[i] <- rep2$G.Rate[match(growth$ID[i], rep2$Isolate)]
  growth$Rep3[i] <- rep3$G.Rate[match(growth$ID[i], rep3$Isolate)]
}

# Calculate Average and SEM
growth$G.Avg <- round(apply(growth[,5:7], 1, mean, na.rm=T), 3)
growth$G.SEM <- round(apply(growth[,5:7], 1, sem), 3)

# Seperate By Organism
kbs701 <- growth[growth$Organism == "701",]
kbs702 <- growth[growth$Organism == "702",]
kbs703 <- growth[growth$Organism == "703",]
kbs710 <- growth[growth$Organism == "710",]
kbs723 <- growth[growth$Organism == "723",]
kbs724 <- growth[growth$Organism == "724",]


# Set Default Plot Parameters
par(mar=c(6, 5, 1, 1) + 0.1)

# Plot Respiration Responses by Organism
bp701g <-barplot(kbs701$G.Avg, las=2, ylim = c(0, 1.2*max(kbs701$G.Avg+kbs701$G.SEM)),
                ylab=expression(paste(mu, "max (day"^-1, ")")), cex.lab= 1.5,
                las = 2, names.arg=c("701", "7011A", "7011B", "7011C", "7011D",
                                     "7011bA", "7011bB", "7011bC", "7011bD", "7013A", "7013B",
                                     "7013C", "7013D"),
                col = c("black", rep("gray48", 4), rep("gray73", 4), rep("gray98", 4)))
mtext("Isolate", side=1, line = 4.5, cex = 1.2)
arrows(x0=bp701g, y0=kbs701$G.Avg, y1=kbs701$G.Avg-kbs701$G.SEM, angle=90, length=0.1, lwd=1)
arrows(x0=bp701g, y0=kbs701$G.Avg, y1=kbs701$G.Avg+kbs701$G.SEM, angle=90, length=0.1, lwd=1)
arrows(x0=bp701g[1], y0=kbs701$G.Avg[1], y1=kbs701$G.Avg[1]-kbs701$G.SEM[1], angle=90,
       length=0.1, lwd=1, col = "white")
legend("topright", legend = c("Ancestor", "Tube 1", "Tube 1b", "Tube 3"),
       fill=c("Black", "gray48", "gray73","gray98"), bty="n")

bp702g <-barplot(means702g, las=2, ylim = c(0, 1.6*max(means702g+sems702g)),
                ylab="Growth Rate",
                las = 2, names.arg=c("702", "7022A", "7022B", "7022C", "7022D",
                                     "7025A", "7025B", "7025C", "7025D", "7026A",
                                     "7026B", "7026C", "7026D"),
                col = c("black", rep("gray48", 4), rep("gray73", 4), rep("gray98", 4)))
mtext("Isolate", side=1, line = 4.5)
arrows(x0=bp702g, y0=means702g, y1=means702g-sems702g, angle=90, length=0.1, lwd=1)
arrows(x0=bp702g, y0=means702g, y1=means702g+sems702g, angle=90, length=0.1, lwd=1)
arrows(x0=bp702g[1], y0=means702g[1], y1=means702g[1]-sems702g[1], angle=90,
       length=0.1, lwd=1, col = "white")
legend("topright", legend = c("Ancestor", "Tube 1", "Tube 2", "Tube 3"),
       fill=c("Black", "gray48", "gray73","gray98"), bty="n")


bp703g <-barplot(means703g, las=2, ylim = c(0, 1.7*max(means703g+sems703g)),
                ylab="Growth Rate",
                las = 2, names.arg=c("703", "7031A", "7031B", "7031C", "7031D",
                                     "7032A", "7032B", "7032C", "7032D", "7034A",
                                     "7034B", "7034C", "7034D"),
                col = c("black", rep("gray48", 4), rep("gray73", 4), rep("gray98", 4)))
mtext("Isolate", side=1, line = 4.5)
arrows(x0=bp703g, y0=means703g, y1=means703g-sems703g, angle=90, length=0.1, lwd=1)
arrows(x0=bp703g, y0=means703g, y1=means703g+sems703g, angle=90, length=0.1, lwd=1)
arrows(x0=bp703g[1], y0=means703g[1], y1=means703g[1]-sems703g[1], angle=90,
       length=0.1, lwd=1, col = "white")
legend("topright", legend = c("Ancestor", "Tube 1", "Tube 2", "Tube 4"),
       fill=c("Black", "gray48", "gray73","gray98"), bty="n")

bp710g <-barplot(means710g, las=2, ylim = c(0, 1.55*max(means710g+sems710g)),
                ylab="Growth Rate",
                las = 2, names.arg=c("710", "7101A", "7101B", "7101C", "7101D",
                                     "7102A", "7102B", "7102C", "7102D", "7103A",
                                     "7103B", "7103C", "7103D"),
                col = c("black", rep("gray48", 4), rep("gray73", 4), rep("gray98", 4)))
mtext("Isolate", side=1, line = 4.5)
arrows(x0=bp710g, y0=means710g, y1=means710g-sems710g, angle=90, length=0.1, lwd=1)
arrows(x0=bp710g, y0=means710g, y1=means710g+sems710g, angle=90, length=0.1, lwd=1)
arrows(x0=bp710g[1], y0=means710g[1], y1=means710g[1]-sems710g[1], angle=90,
       length=0.1, lwd=1, col = "white")
legend("topright", legend = c("Ancestor", "Tube 1", "Tube 2", "Tube 3"),
       fill=c("Black", "gray48", "gray73","gray98"), bty="n")

bp723g <-barplot(means723g, las=2, ylim = c(0, 1.65*max(means723g+sems723g)),
                ylab="Growth Rate",
                las = 2, names.arg=c("723", "7231A", "7231B", "7231C", "7231D",
                                     "7232A", "7232B", "7232C", "7232D", "7233A",
                                     "7233B", "7233C", "7233D"),
                col = c("black", rep("gray48", 4), rep("gray73", 4), rep("gray98", 4)))
mtext("Isolate", side=1, line = 4.5)
arrows(x0=bp723g, y0=means723g, y1=means723g-sems723g, angle=90, length=0.1, lwd=1)
arrows(x0=bp723g, y0=means723g, y1=means723g+sems723g, angle=90, length=0.1, lwd=1)
arrows(x0=bp723g[1], y0=means723g[1], y1=means723g[1]-sems723g[1], angle=90,
       length=0.1, lwd=1, col = "white")
legend("topright", legend = c("Ancestor", "Tube 1", "Tube 2", "Tube 3"),
       fill=c("Black", "gray48", "gray73","gray98"), bty="n")

bp724g <-barplot(means724g, las=2, ylim = c(0, 1.2*max(means724g+sems724g)),
                ylab="Growth Rate",
                las = 2, names.arg=c("724", "7241A", "7241B", "7241C", "7241D",
                                     "7242A", "7242B", "7242C", "7242D", "7243A",
                                     "7243B", "7243C", "7243D"),
                col = c("black", rep("gray48", 4), rep("gray73", 4), rep("gray98", 4)))
mtext("Isolate", side=1, line = 4.5)
arrows(x0=bp724g, y0=means724g, y1=means724g-sems724g, angle=90, length=0.1, lwd=1)
arrows(x0=bp724g, y0=means724g, y1=means724g+sems724g, angle=90, length=0.1, lwd=1)
arrows(x0=bp724g[1], y0=means724g[1], y1=means724g[1]-sems724g[1], angle=90,
       length=0.1, lwd=1, col = "white")
legend("topright", legend = c("Ancestor", "Tube 1", "Tube 2", "Tube 3"),
       fill=c("Black", "gray48", "gray73","gray98"), bty="n")



# Statistics
m.701 <- melt(kbs701[,1:7])
fit.701 <- aov(value ~ Evol + Tube + ID, data = m.701)
summary(fit.701)
TukeyHSD(fit.701)

m.702 <- melt(kbs702[,1:7])
fit.702 <- aov(value ~ Evol + Tube + ID, data = m.702)
summary(fit.702)
TukeyHSD(fit.702)

m.703 <- melt(kbs703[,1:7])
fit.703 <- aov(value ~ Evol + Tube + ID, data = m.703)
summary(fit.703)
TukeyHSD(fit.703)

m.710 <- melt(kbs710[,1:7])
fit.710 <- aov(value ~ Evol + Tube + ID, data = m.710)
summary(fit.710)
TukeyHSD(fit.710)

m.723 <- melt(kbs723[,1:7])
fit.723 <- aov(value ~ Evol + Tube + ID, data = m.723)
summary(fit.723)
TukeyHSD(fit.723)

m.724 <- melt(kbs724[,1:7])
fit.724 <- aov(value ~ Evol + Tube + ID, data = m.724)
summary(fit.724)
TukeyHSD(fit.724)





