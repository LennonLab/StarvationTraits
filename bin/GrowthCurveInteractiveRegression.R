################################################################################
#                                                                              #
#  Analysis of SynergyMX Growth Curve Data                                     #
#   Version 2.0                                                                #
#  Written By: Mario Muscarella                                                #
#  Last Update: 16 Jul 2015                                                    #
#                                                                              #
#  This analysis uses the R package - rpanel which does interactive regression #
#  The goal is to use the interactive view to pick the area to analyze         #
#  rpanel - Interactive Regression                                             #
#  Multiple trend chart version; Annual anomaly data                           #
#  Added Support for csv and txt                                               #
#  Added Support for Row and Column format inputs                              #
#  Turned Source into Function                                                 #
#                                                                              #
################################################################################

Growth.Expo <- function(infile = " ", outfile = " ", in.format = "Rows"){
# Load Dependencies
  require("rpanel")||install.packages("rpanel");require("rpanel")
  source("./bin/ReadSynergy.R")

# Global Options
  options(digits=6)
  script <- "Interactive Regression Analysis"

# Import Data
  data.in <- read.synergy(input = infile, skip = "32")

# Exponential Transformation
  data.log <- matrix(NA, nrow = dim(data.in)[1], ncol = dim(data.in)[2])
  colnames(data.log) <- colnames(data.in)
  data.log[,1:2] <- as.matrix(data.in[,1:2])
  for (i in 3:dim(data.in)[2]){
    data.log[,i] <- round(sapply(data.in[,i], log1p), 3)
  }

  data.log <- as.data.frame(data.log)
  data.log <- data.in # no transformation - testing

# Creat Output
  titles <- c("Sample", "Start", "End", "Rate", "R2", "P-value")
  # write.table(as.matrix(t(titles)), file=outfile, append=T, row.names=F, col.names=F, sep=",", quote=FALSE)

# Select Samples
  samples <- as.factor(colnames(data.log)[3:26])

# Create Plotting Window
  windows()
  par(las=1)
  par(fig=c(0,1,0, 1), new = F)
  par(ps=9); par(cex.axis=c(0.9)); par(cex.lab=c(0.9)); par(oma=c(3,1,1,0.5)); par(mar=c(4,4,2,1))

# Attach Data
  attach(data.log)

### rpanel function ################
  draw <- function(panel) {
    if (panel$end > panel$start)
      {start <- panel$start
      end <- panel$end}
    else
      {start <- min(Time)
      end <- max(Time)}

    data.log$samp <- panel$samp
    name <- panel$sample.name

# Par Settings
    par(mar=c(5,5,1,1))

# Text placement points
  xaxis_pt <- max(Time) - 0.1*(max(Time)-min(Time))
  yaxis_pt <- 375

# Make data subset based on start & end yrs
  sub <- subset(data.log, Time >= start)
  sub <- subset(sub, Time <= end)

# Linear trend line lm & coefs / stats
  trend <- lm(sub$samp ~ sub$Time)
  a <- as.numeric(coef(trend)[1]);  b <- as.numeric(coef(trend)[2])
  r2 <- round(summary(trend)$r.squared, 3)
  p <- round(anova(trend)$'Pr(>F)'[1], 4)
  p <- ifelse (p == 0, "<0.001", p)
  start.2 <- signif(start, digits = 3)
  end.2 <- signif(end, digits = 3)
  rate <- signif(-b,3)

# Make Basic plot
  plot(log1p(panel$samp) ~ Time, type = "b", col = "darkgrey", xlab = "Time (Hrs)",
    ylab = expression(paste("Log OD (uM O "[2],")")), log="y",
    par(bty="n"),xlim=c(0,max(Time)+0.5),# ylim=c(0, 1),
    xaxs = "i", yaxs = "i", axes = FALSE, cex.main = 1.25, cex.lab = 1.5,
    main = "Interactive Regression of Exponential Growth Rate Data")
  axis(1, col = "grey"); axis(2, col = "grey")

# Calc vals for start-end regression line & add line
  x_vals = c(panel$start, panel$end)
  y_vals = c(a+b*panel$start, a+b*panel$end)
	lines(x_vals, y_vals, col= "red")
  points(sub$Time, sub$samp, pch = 19, col = "red", type = "p")
  text(xaxis_pt, yaxis_pt, paste("Sample: ", name))
  text(xaxis_pt, yaxis_pt - 15, paste("Period: ", start.2, " to " , end.2, "Hrs"))
  text(xaxis_pt, yaxis_pt - 30, bquote(Rate == .(rate) ~ Hr^-1), cex=1)
  text(xaxis_pt, yaxis_pt - 45, bquote(R^2 == .(r2)), cex=1)
  text(xaxis_pt, yaxis_pt - 60, paste("P-value = ", p), cex=1)

### Outer Margin Annotation
	my_date <- format(Sys.time(), "%m/%d/%y")
	mtext(script, side = 1, line = .75, cex=0.8, outer = T, adj = 0)
	mtext(my_date, side = 1, line =.75, cex = 0.8, outer = T, adj = 1)
	data.out <- list(Start=start, End=end, Rate=rate, R2=r2, Pvalue=p)
	panel
  }

  collect.data <- function(panel) {
    data.log$samp <- panel$samp
    name <- panel$sample.name
    start <- panel$start
    end <- panel$end
    # Text placement points
    xaxis_pt <- max(Time) - 0.1*(max(Time)-min(Time))
    yaxis_pt <- 375
    # Make data subset based on start & end yrs
    sub <- subset(data.lot, Time >= start)
    sub <- subset(sub, Time <= end)
    # Linear trend line lm & coefs / stats
    trend <- lm(sub$samp ~ sub$Time)
    a <- as.numeric(coef(trend)[1]);  b <- as.numeric(coef(trend)[2])
    r2 <- round(summary(trend)$r.squared, 3)
    p <- round(anova(trend)$'Pr(>F)'[1], 4)
    p <- ifelse (p == 0, "<0.001", p)
    start.2 <- signif(start, digits = 3)
    end.2 <- signif(end, digits = 3)
    rate <- signif(-b,3)
    data.sample <- name
    data.start <- start.2
    data.end <- end.2
    data.rate <- rate
    data.r2 <- r2
    data.p <- p
    data.out <- c(data.sample, data.start, data.end, data.rate, data.r2, data.p)
    #write.table(as.matrix(t(data.out)), file=outfile, append=T, row.names=F, col.names=F, sep=",", quote=FALSE)
    panel
    }

  end.session <- function(panel) {
    dev.off(2)
    print.noquote("Good-Bye: Computer Will Now Self-Destruct")
    panel
    }

# rpanel controls - enter start and end yrs for portion of full data set to be used
  rpplot <- rp.control(title="Interactive Regression", start=0, end = max(Time), initval = samples[1])
  rp.listbox(rpplot, variable = samp, vals = "Samples", labels = samples, action = draw)
  rp.slider(rpplot, start, action = draw, from = 0, to =  max(Time))
  rp.slider(rpplot, end, action = draw, from = 0, to =  max(Time))
  rp.doublebutton(rpplot, var = start, step = 0.05, title = "Start Fine Adjustment", action = draw)
  rp.doublebutton(rpplot, var = end, step = 0.05, title = "End Fine Adjustment", action = draw)
  rp.textentry(rpplot, var = sample.name, action = draw, labels = "Sample Name", initval = "")
  rp.button(rpplot, title = "save", action = collect.data)
  rp.button(rpplot, title = "quit", action = end.session, quitbutton=T)

  }
