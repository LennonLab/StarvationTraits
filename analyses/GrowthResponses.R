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


#Plate 1
A7011g <- growth.data[1:3,]
A7011g_mean <- mean(A7011g$G.Rate)
A7011g_sem <- sem(A7011g$G.Rate)
B7011g <- growth.data[4:6,]
B7011g_mean <- mean(B7011g$G.Rate)
B7011g_sem <- sem(B7011g$G.Rate)
C7011g <- growth.data[7:9,]
C7011g_mean <- mean(C7011g$G.Rate)
C7011g_sem <- sem(C7011g$G.Rate)
D7011g <- growth.data[10:12,]
D7011g_mean <- mean(D7011g$G.Rate)
D7011g_sem <- sem(D7011g$G.Rate)
A7022g <- growth.data[13:15,]
A7022g_mean <- mean(A7022g$G.Rate)
A7022g_sem <- sem(A7022g$G.Rate)
B7022g <- growth.data[16:18,]
B7022g_mean <- mean(B7022g$G.Rate)
B7022g_sem <- sem(B7022g$G.Rate)
C7022g <- growth.data[19:21,]
C7022g_mean <- mean(C7022g$G.Rate)
C7022g_sem <- sem(C7022g$G.Rate)
D7022g <- growth.data[22:24,]
D7022g_mean <- mean(D7022g$G.Rate)
D7022g_sem <- sem(D7022g$G.Rate)
A7025g <- growth.data[25:27,]
A7025g_mean <- mean(A7025g$G.Rate)
A7025g_sem <- sem(A7025g$G.Rate)
B7025g <- growth.data[28:30,]
B7025g_mean <- mean(B7025g$G.Rate)
B7025g_sem <- sem(B7025g$G.Rate)
C7025g <- growth.data[31:33,]
C7025g_mean <- mean(C7025g$G.Rate)
C7025g_sem <- sem(C7025g$G.Rate)
D7025g <- growth.data[34:36,]
D7025g_mean <- mean(D7025g$G.Rate)
D7025g_sem <- sem(D7025g$G.Rate)
A7031g <- growth.data[37:39,]
A7031g_mean <- mean(A7031g$G.Rate)
A7031g_sem <- sem(A7031g$G.Rate)
B7031g <- growth.data[40:42,]
B7031g_mean <- mean(B7031g$G.Rate)
B7031g_sem <- sem(B7031g$G.Rate)
C7031g <- growth.data[43:45,]
C7031g_mean <- mean(C7031g$G.Rate, na.rm=T)
C7031g_sem <- sem(C7031g$G.Rate)

#Plate 2
B7032g <- growth.data[46:48,]
B7032g_mean <- mean(B7032g$G.Rate)
B7032g_sem <- sem(B7032g$G.Rate)
A7032g <- growth.data[49:51,]
A7032g_mean <- mean(A7032g$G.Rate)
A7032g_sem <- sem(A7032g$G.Rate)
D7031g <- growth.data[52:54,]
D7031g_mean <- mean(D7031g$G.Rate)
D7031g_sem <- sem(D7031g$G.Rate)
D7026g <- growth.data[55:57,]
D7026g_mean <- mean(D7026g$G.Rate)
D7026g_sem <- sem(D7026g$G.Rate)
C7026g <- growth.data[58:60,]
C7026g_mean <- mean(C7026g$G.Rate)
C7026g_sem <- sem(C7026g$G.Rate)
B7026g <- growth.data[61:63,]
B7026g_mean <- mean(B7026g$G.Rate)
B7026g_sem <- sem(B7026g$G.Rate)
A7026g <- growth.data[64:66,]
A7026g_mean <- mean(A7026g$G.Rate)
A7026g_sem <- sem(A7026g$G.Rate)
D7013g <- growth.data[67:69,]
D7013g_mean <- mean(D7013g$G.Rate)
D7013g_sem <- sem(D7013g$G.Rate)
C7013g <- growth.data[70:72,]
C7013g_mean <- mean(C7013g$G.Rate)
C7013g_sem <- sem(C7013g$G.Rate)
B7013g <- growth.data[73:75,]
B7013g_mean <- mean(B7013g$G.Rate)
B7013g_sem <- sem(B7013g$G.Rate)
A7013g <- growth.data[76:78,]
A7013g_mean <- mean(A7013g$G.Rate)
A7013g_sem <- sem(A7013g$G.Rate)
C7011bg <- growth.data[79:81,]
C7011bg_mean <- mean(C7011bg$G.Rate)
C7011bg_sem <- sem(C7011bg$G.Rate)
B7011bg <- growth.data[82:84,]
B7011bg_mean <- mean(B7011bg$G.Rate, na.rm=TRUE)
B7011bg_sem <- sem(B7011bg$G.Rate)
A7011bg <- growth.data[85:87,]
A7011bg_mean <- mean(A7011bg$G.Rate, na.rm=TRUE)
A7011bg_sem <- sem(A7011bg$G.Rate)

#Plate3
C7032g <- growth.data[88:90,]
C7032g_mean <- mean(C7032g$G.Rate)
C7032g_sem <- sem(C7032g$G.Rate)
D7032g <- growth.data[91:93,]
D7032g_mean <- mean(D7032g$G.Rate)
D7032g_sem <- sem(D7032g$G.Rate)
A7034g <- growth.data[94:96,]
A7034g_mean <- mean(A7034g$G.Rate)
A7034g_sem <- sem(A7034g$G.Rate)
B7034g <- growth.data[97:99,]
B7034g_mean <- mean(B7034g$G.Rate)
B7034g_sem <- sem(B7034g$G.Rate)
C7034g <- growth.data[100:102,]
C7034g_mean <- mean(C7034g$G.Rate)
C7034g_sem <- sem(C7034g$G.Rate)
D7034g <- growth.data[103:105,]
D7034g_mean <- mean(D7034g$G.Rate)
D7034g_sem <- sem(D7034g$G.Rate)
A7101g <- growth.data[106:108,]
A7101g_mean <- mean(A7101g$G.Rate)
A7101g_sem <- sem(A7101g$G.Rate)
B7101g <- growth.data[109:111,]
B7101g_mean <- mean(B7101g$G.Rate)
B7101g_sem <- sem(B7101g$G.Rate)
C7101g <- growth.data[112:114,]
C7101g_mean <- mean(C7101g$G.Rate)
C7101g_sem <- sem(C7101g$G.Rate)
D7101g <- growth.data[115:117,]
D7101g_mean <- mean(D7101g$G.Rate)
D7101g_sem <- sem(D7101g$G.Rate)
A7102g <- growth.data[118:120,]
A7102g_mean <- mean(A7102g$G.Rate)
A7102g_sem <- sem(A7102g$G.Rate)
B7102g <- growth.data[121:123,]
B7102g_mean <- mean(B7102g$G.Rate)
B7102g_sem <- sem(B7102g$G.Rate)
C7102g <- growth.data[124:126,]
C7102g_mean <- mean(C7102g$G.Rate)
C7102g_sem <- sem(C7102g$G.Rate)
D7102g <- growth.data[127:129,]
D7102g_mean <- mean(D7102g$G.Rate)
D7102g_sem <- sem(D7102g$G.Rate)
A7103g <- growth.data[130:132,]
A7103g_mean <- mean(A7103g$G.Rate, na.rm=T)
A7103g_sem <- sem(A7103g$G.Rate)

#Plate4
B7103g <- growth.data[133:135,]
B7103g_mean <- mean(B7103g$G.Rate)
B7103g_sem <- sem(B7103g$G.Rate)
C7103g <- growth.data[136:138,]
C7103g_mean <- mean(C7103g$G.Rate)
C7103g_sem <- sem(C7103g$G.Rate)
D7103g <- growth.data[139:141,]
D7103g_mean <- mean(D7103g$G.Rate)
D7103g_sem <- sem(D7103g$G.Rate)
A7231g <- growth.data[142:144,]
A7231g_mean <- mean(A7231g$G.Rate)
A7231g_sem <- sem(A7231g$G.Rate)
B7231g <- growth.data[145:147,]
B7231g_mean <- mean(B7231g$G.Rate)
B7231g_sem <- sem(B7231g$G.Rate)
C7231g <- growth.data[148:150,]
C7231g_mean <- mean(C7231g$G.Rate)
C7231g_sem <- sem(C7231g$G.Rate)
D7231g <- growth.data[151:153,]
D7231g_mean <- mean(D7231g$G.Rate)
D7231g_sem <- sem(D7231g$G.Rate)
A7232g <- growth.data[154:156,]
A7232g_mean <- mean(A7232g$G.Rate)
A7232g_sem <- sem(A7232g$G.Rate)
B7232g <- growth.data[157:159,]
B7232g_mean <- mean(B7232g$G.Rate)
B7232g_sem <- sem(B7232g$G.Rate)
C7232g <- growth.data[160:162,]
C7232g_mean <- mean(C7232g$G.Rate)
C7232g_sem <- sem(C7232g$G.Rate)
D7232g <- growth.data[163:165,]
D7232g_mean <- mean(D7232g$G.Rate)
D7232g_sem <- sem(D7232g$G.Rate)
A7233g <- growth.data[166:168,]
A7233g_mean <- mean(A7233g$G.Rate)
A7233g_sem <- sem(A7233g$G.Rate)
B7233g <- growth.data[169:171,]
B7233g_mean <- mean(B7233g$G.Rate)
B7233g_sem <- sem(B7233g$G.Rate)
C7233g <- growth.data[172:174,]
C7233g_mean <- mean(C7233g$G.Rate)
C7233g_sem <- sem(C7233g$G.Rate)
D7233g <- growth.data[175:177,]
D7233g_mean <- mean(D7233g$G.Rate, na.rm=T)
D7233g_sem <- sem(D7233g$G.Rate)

#Plate5
A7241g <- growth.data[178:180,]
A7241g_mean <- mean(A7241g$G.Rate)
A7241g_sem <- sem(A7241g$G.Rate)
B7241g <- growth.data[181:183,]
B7241g_mean <- mean(B7241g$G.Rate)
B7241g_sem <- sem(B7241g$G.Rate)
C7241g <- growth.data[184:186,]
C7241g_mean <- mean(C7241g$G.Rate)
C7241g_sem <- sem(C7241g$G.Rate)
D7241g <- growth.data[187:189,]
D7241g_mean <- mean(D7241g$G.Rate)
D7241g_sem <- sem(D7241g$G.Rate)
A7242g <- growth.data[190:192,]
A7242g_mean <- mean(A7242g$G.Rate)
A7242g_sem <- sem(A7242g$G.Rate)
B7242g <- growth.data[193:195,]
B7242g_mean <- mean(B7242g$G.Rate)
B7242g_sem <- sem(B7242g$G.Rate)
C7242g <- growth.data[196:198,]
C7242g_mean <- mean(C7242g$G.Rate)
C7242g_sem <- sem(C7242g$G.Rate)
D7242g <- growth.data[199:201,]
D7242g_mean <- mean(D7242g$G.Rate)
D7242g_sem <- sem(D7242g$G.Rate)
A7243g <- growth.data[202:204,]
A7243g_mean <- mean(A7243g$G.Rate)
A7243g_sem <- sem(A7243g$G.Rate)
B7243g <- growth.data[205:207,]
B7243g_mean <- mean(B7243g$G.Rate)
B7243g_sem <- sem(B7243g$G.Rate)
C7243g <- growth.data[208:210,]
C7243g_mean <- mean(C7243g$G.Rate)
C7243g_sem <- sem(C7243g$G.Rate)
D7243g <- growth.data[211:213,]
D7243g_mean <- mean(D7243g$G.Rate)
D7243g_sem <- sem(D7243g$G.Rate)
A701g <- growth.data[214:216,]
A701g_mean <- mean(A701g$G.Rate)
A701g_sem <- sem(A701g$G.Rate)
A702g <- growth.data[217:219,]
A702g_mean <- mean(A702g$G.Rate)
A702g_sem <- sem(A702g$G.Rate)
A703g <- growth.data[220:222,]
A703g_mean <- mean(A703g$G.Rate, na.rm=T)
A703g_sem <- sem(A703g$G.Rate)

#Plate6
D7011bg <- growth.data[223:225,]
D7011bg_mean <- mean(D7011bg$G.Rate)
D7011bg_sem <- sem(D7011g$G.Rate)
A724g <- growth.data[226:228,]
A724g_mean <- mean(A724g$G.Rate)
A724g_sem <- sem(A724g$G.Rate)
A723g <- growth.data[229:231,]
A723g_mean <- mean(A723g$G.Rate, na.rm=T)
A723g_sem <- sem(A723g$G.Rate)
A710g <- growth.data[232:234,]
A710g_mean <- mean(A710g$G.Rate, na.rm=T)
A710g_sem <- sem(A710g$G.Rate)


# Combine means and sems based on isolate
means701g <-c(A701g_mean, A7011g_mean, B7011g_mean, C7011g_mean, D7011g_mean, A7011bg_mean, B7011bg_mean,
             C7011bg_mean, D7011bg_mean, A7013g_mean, B7013g_mean, C7013g_mean, D7013g_mean)
means702g <-c(A702g_mean, A7022g_mean, B7022g_mean, C7022g_mean, D7022g_mean, A7025g_mean, B7025g_mean,
             C7025g_mean, D7025g_mean, A7026g_mean, B7026g_mean, C7026g_mean, D7026g_mean)
means703g <-c(A703g_mean, A7031g_mean, B7031g_mean, C7031g_mean, D7031g_mean, A7032g_mean, B7032g_mean,
             C7032g_mean, D7032g_mean, A7034g_mean, B7034g_mean, C7034g_mean, D7034g_mean)
means710g <-c(A710g_mean, A7101g_mean, B7101g_mean, C7101g_mean, D7101g_mean, A7102g_mean, B7102g_mean,
             C7102g_mean, D7102g_mean, A7103g_mean, B7103g_mean, C7103g_mean, D7103g_mean)
means723g <-c(A723g_mean, A7231g_mean, B7231g_mean, C7231g_mean, D7231g_mean, A7232g_mean, B7232g_mean,
             C7232g_mean, D7232g_mean, A7233g_mean, B7233g_mean, C7233g_mean, D7233g_mean)
means724g <-c(A724g_mean, A7241g_mean, B7241g_mean, C7241g_mean, D7241g_mean, A7242g_mean, B7242g_mean,
             C7242g_mean, D7242g_mean, A7243g_mean, B7243g_mean, C7243g_mean, D7243g_mean)

sems701g <-c(A701g_sem, A7011g_sem, B7011g_sem, C7011g_sem, D7011g_sem, A7011bg_sem, B7011bg_sem,
            C7011bg_sem, D7011bg_sem, A7013g_sem, B7013g_sem, C7013g_sem, D7013g_sem)
sems702g <-c(A702g_sem, A7022g_sem, B7022g_sem, C7022g_sem, D7022g_sem, A7025g_sem, B7025g_sem, C7025g_sem,
            D7025g_sem,A7026g_sem, B7026g_sem, C7026g_sem, D7026g_sem)
sems703g <-c(A703g_sem, A7031g_sem, B7031g_sem, C7031g_sem, D7031g_sem, A7032g_sem, B7032g_sem, C7032g_sem,
            D7032g_sem, A7034g_sem, B7034g_sem, C7034g_sem, D7034g_sem)
sems710g <-c(A710g_sem, A7101g_sem, B7101g_sem, C7101g_sem, D7101g_sem, A7102g_sem, B7102g_sem,
            C7102g_sem, D7102g_sem, A7103g_sem, B7103g_sem, C7103g_sem, D7103g_sem)
sems723g <-c(A723g_sem, A7231g_sem, B7231g_sem, C7231g_sem, D7231g_sem, A7232g_sem, B7232g_sem,
            C7232g_sem, D7232g_sem, A7233g_sem, B7233g_sem, C7233g_sem, D7233g_sem)
sems724g <-c(A724g_sem, A7241g_sem, B7241g_sem, C7241g_sem, D7241g_sem, A7242g_sem, B7242g_sem,
            C7242g_sem, D7242g_sem, A7243g_sem, B7243g_sem, C7243g_sem, D7243g_sem)

# Set Default Plot Parameters
par(mar=c(6, 5, 1, 1) + 0.1)

# Plot Respiration Responses by Organism
bp701g <-barplot(means701g, las=2, ylim = c(0, 1.2*max(means701g+sems701g)),
                ylab="Growth Rate",
                las = 2, names.arg=c("701", "7011A", "7011B", "7011C", "7011D",
                                     "7011bA", "7011bB", "7011bC", "7011bD", "7013A", "7013B",
                                     "7013C", "7013D"),
                col = c("black", rep("gray48", 4), rep("gray73", 4), rep("gray98", 4)))
mtext("Isolate", side=1, line = 4.5)
arrows(x0=bp701g, y0=means701g, y1=means701g-sems701g, angle=90, length=0.1, lwd=1)
arrows(x0=bp701g, y0=means701g, y1=means701g+sems701g, angle=90, length=0.1, lwd=1)
arrows(x0=bp701g[1], y0=means701g[1], y1=means701g[1]-sems701g[1], angle=90,
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
