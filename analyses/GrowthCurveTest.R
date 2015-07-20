################################################################################
#                                                                              #
#  Test Script for Analysis of Exponential Growth Curve Data                   #
#   Version 2.0                                                                #
#  Written By: Rachel Ferrill & Mario Muscarella                               #
#  Last Update: 17 Jul 2015                                                    #
#                                                                              #
#  Use this file to test the PreSens.Respiration fucntion on local machine     #
#  And as a template to create your own analysis                               #
#                                                                              #
################################################################################

setwd("~/GitHub/StarvationTraits/")
rm(list=ls())

# Inport the function from source file
source("./bin/GrowthCurveInteractiveRegression.R")
source("./bin/ModifiedGomp.R")

################################################################################
# Examples #####################################################################
################################################################################

# Create Directory For Output
dir.create("./data/GrowthCurves/output", showWarnings = FALSE)

# Run Example with Test Data
growth.modGomp("./data/GrowthCurves/GrowthCurveExample.txt", "test", skip=31)

# Run Example with RNF Data
growth.modGomp("./data/GrowthCurves/RNF_GrowthCurve_20150625_plate2.txt",
               "plate2_test", skip=32)

