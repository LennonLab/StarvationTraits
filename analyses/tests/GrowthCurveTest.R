################################################################################
#                                                                              #
#  Test Script for Analysis of Exponential Growth Curve Data                   #
#   This uses ModifiedGomp.R Version 2.0                                       #
#  Written By: Mario Muscarella                                                #
#  Last Update: 29 Jan 2015                                                    #
#                                                                              #
#  Use this file to analyze Synergy MX Growth Curve data                       #
#                                                                              #
################################################################################

setwd("~/GitHub/StarvationTraits/")
rm(list=ls())

# Inport the function from source file
source("./bin/GrowthCurveInteractiveRegression.R")
source("./bin/ModifiedGomp.R")

# Create Directory For Output
dir.create("./data/GrowthCurves/output", showWarnings = FALSE)

################################################################################
# Example ######################################################################
################################################################################

# Run Example with Test Data
growth.modGomp("./data/GrowthCurves/GrowthCurveExample.txt", "test", skip=31)
