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

# Remove 7011 because of contamination
growth <- growth[-c(2,3,8,9), ]

# Seperate By Organism
kbs701 <- growth[growth$Organism == "701",]
kbs702 <- growth[growth$Organism == "702",]
kbs703 <- growth[growth$Organism == "703",]
kbs710 <- growth[growth$Organism == "710",]
kbs723 <- growth[growth$Organism == "723",]
kbs724 <- growth[growth$Organism == "724",]

# Remove Some Replicates Due to Issues
kbs723 <- kbs723[kbs723$ID != "7231C", ]
kbs710 <- kbs710[kbs710$ID != "7102B", ]
kbs701 <- kbs701[kbs701$ID != "7011bB", ]

# Set Default Plot Parameters
par(mar=c(6, 5, 1, 1) + 0.1)

# Plot Respiration Responses by Organism
bp701g <- barplot(kbs701$G.Avg, las=2, ylim = c(0, 1.2*max(kbs701$G.Avg+kbs701$G.SEM)),
                  ylab=expression(paste(mu, "max (day"^-1, ")")), cex.lab= 1.5,
                  las = 2, names.arg=kbs701$ID,
                  col = c("black", rep("gray48", 4), rep("gray98", 4)))
          mtext("Isolate", side=1, line = 4.5, cex = 1.2)
          arrows(x0=bp701g, y0=kbs701$G.Avg, y1=kbs701$G.Avg-kbs701$G.SEM,
                 angle=90, length=0.1, lwd=1)
          arrows(x0=bp701g, y0=kbs701$G.Avg, y1=kbs701$G.Avg+kbs701$G.SEM,
                 angle=90, length=0.1, lwd=1)
          arrows(x0=bp701g[1], y0=kbs701$G.Avg[1], y1=kbs701$G.Avg[1]-kbs701$G.SEM[1],
                 angle=90, length=0.1, lwd=1, col = "white")
          legend("topright", legend = c("Ancestor",  "Tube 1b", "Tube 3"),
                 fill=c("Black", "gray48", "gray98"), bty="n")

bp702g <- barplot(kbs702$G.Avg, las=2, ylim = c(0, 1.6*max(kbs702$G.Avg+kbs702$G.SEM)),
                  ylab=expression(paste(mu, "max (day"^-1, ")")), cex.lab= 1.5,
                  las = 2, names.arg=kbs702$ID,
                  col = c("black", rep("gray48", 4), rep("gray73", 4), rep("gray98", 4)))
          mtext("Isolate", side=1, line = 4.5)
          arrows(x0=bp702g, y0=kbs702$G.Avg, y1=kbs702$G.Avg-kbs702$G.SEM,
                 angle=90, length=0.1, lwd=1)
          arrows(x0=bp702g, y0=kbs702$G.Avg, y1=kbs702$G.Avg+kbs702$G.SEM,
                 angle=90, length=0.1, lwd=1)
          arrows(x0=bp702g[1], y0=kbs702$G.Avg[1], y1=kbs702$G.Avg[1]-kbs702$G.SEM[1],
                 angle=90, length=0.1, lwd=1, col = "white")
          legend("topright", legend = c("Ancestor", "Tube 1", "Tube 2", "Tube 3"),
                 fill=c("Black", "gray48", "gray73","gray98"), bty="n")


bp703g <- barplot(kbs703$G.Avg, las=2, ylim = c(0, 1.7*max(kbs703$G.Avg+kbs703$G.SEM)),
                  ylab=expression(paste(mu, "max (day"^-1, ")")), cex.lab= 1.5,
                  las = 2, names.arg=kbs703$ID,
                  col = c("black", rep("gray48", 4), rep("gray73", 4), rep("gray98", 4)))
          mtext("Isolate", side=1, line = 4.5)
          arrows(x0=bp703g, y0=kbs703$G.Avg, y1=kbs703$G.Avg-kbs703$G.SEM,
                 angle=90, length=0.1, lwd=1)
          arrows(x0=bp703g, y0=kbs703$G.Avg, y1=kbs703$G.Avg+kbs703$G.SEM,
                 angle=90, length=0.1, lwd=1)
          arrows(x0=bp703g[1], y0=kbs703$G.Avg[1], y1=kbs703$G.Avg[1]-kbs703$G.SEM[1],
                 angle=90, length=0.1, lwd=1, col = "white")
          legend("topright", legend = c("Ancestor", "Tube 1", "Tube 2", "Tube 4"),
                 fill=c("Black", "gray48", "gray73","gray98"), bty="n")

bp710g <- barplot(kbs710$G.Avg, las=2, ylim = c(0, 1.55*max(kbs710$G.Avg+kbs710$G.SEM)),
                  ylab=expression(paste(mu, "max (day"^-1, ")")), cex.lab= 1.5,
                  las = 2, names.arg=kbs710$ID,
                  col = c("black", rep("gray48", 4), rep("gray73", 4), rep("gray98", 4)))
          mtext("Isolate", side=1, line = 4.5)
          arrows(x0=bp710g, y0=kbs710$G.Avg, y1=kbs710$G.Avg-kbs710$G.SEM,
                 angle=90, length=0.1, lwd=1)
          arrows(x0=bp710g, y0=kbs710$G.Avg, y1=kbs710$G.Avg+kbs710$G.SEM,
                 angle=90, length=0.1, lwd=1)
          arrows(x0=bp710g[1], y0=kbs710$G.Avg[1], y1=kbs710$G.Avg[1]-kbs710$G.SEM[1],
                         angle=90,length=0.1, lwd=1, col = "white")
          legend("topright", legend = c("Ancestor", "Tube 1", "Tube 2", "Tube 3"),
                 fill=c("Black", "gray48", "gray73","gray98"), bty="n")

bp723g <- barplot(kbs723$G.Avg, las=2, ylim = c(0, 1.65*max(kbs723$G.Avg+kbs723$G.SEM)),
                  ylab=expression(paste(mu, "max (day"^-1, ")")), cex.lab= 1.5,
                  las = 2, names.arg=kbs723$ID,
                  col = c("black", rep("gray48", 4), rep("gray73", 4), rep("gray98", 4)))
          mtext("Isolate", side=1, line = 4.5)
          arrows(x0=bp723g, y0=kbs723$G.Avg, y1=kbs723$G.Avg-kbs723$G.SEM,
                 angle=90, length=0.1, lwd=1)
          arrows(x0=bp723g, y0=kbs723$G.Avg, y1=kbs723$G.Avg+kbs723$G.SEM,
                 angle=90, length=0.1, lwd=1)
          arrows(x0=bp723g[1], y0=kbs723$G.Avg[1], y1=kbs723$G.Avg[1]-kbs723$G.SEM[1],
                 angle=90, length=0.1, lwd=1, col = "white")
          legend("topright", legend = c("Ancestor", "Tube 1", "Tube 2", "Tube 3"),
                 fill=c("Black", "gray48", "gray73","gray98"), bty="n")

bp724g <- barplot(kbs724$G.Avg, las=2, ylim = c(0, 1.2*max(kbs724$G.Avg+kbs724$G.SEM)),
                  ylab=expression(paste(mu, "max (day"^-1, ")")), cex.lab= 1.5,
                  las = 2, names.arg=kbs724$ID,
                  col = c("black", rep("gray48", 4), rep("gray73", 4), rep("gray98", 4)))
          mtext("Isolate", side=1, line = 4.5)
          arrows(x0=bp724g, y0=kbs724$G.Avg, y1=kbs724$G.Avg-kbs724$G.SEM,
                 angle=90, length=0.1, lwd=1)
          arrows(x0=bp724g, y0=kbs724$G.Avg, y1=kbs724$G.Avg+kbs724$G.SEM,
                 angle=90, length=0.1, lwd=1)
          arrows(x0=bp724g[1], y0=kbs724$G.Avg[1], y1=kbs724$G.Avg[1] - kbs724$G.SEM[1],
                 angle=90, length=0.1, lwd=1, col = "white")
          legend("topright", legend = c("Ancestor", "Tube 1", "Tube 2", "Tube 3"),
                 fill=c("Black", "gray48", "gray73","gray98"), bty="n")


################################################################################
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
test.701 <- lm(value ~ Evol + Evol/Tube, data = m.701)
tests.701 <- lmer(value ~ Evol + Tube + (1|Tube/ID), data = m.701)
aic(fit.703)

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

require (lme4)

require(nlme)

fit.702 <- aov(value ~ Evol + Tube + ID, data = m.702)
fit.702.int <- aov(value ~ Evol * Tube * ID, data = m.702)
lm.702 <- lm(value ~ Evol + Evol/Tube, data = m.702)
lmer.702 <- lmer(value ~ Evol + Tube + (1|Tube/ID), data = m.702)
nlmer.702 <- nlmer(value ~ Evol + Tube + (1|Tube/ID), data = m.702)

AIC(lmer.702, nlmer.702)


gr.res <- data.frame(m.702$ID, m.702.r$ID, m.702$value,m.702.r$value)
plot(gr.res[,3]~gr.res[,4])

