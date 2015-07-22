################################################################################
#                                                                              #
#  Script for Analysis of Exponential Growth Curve Data                        #
#   This uses ModifiedGomp.R Version 2.0                                       #
#  Written By: Rachel Ferrill & Mario Muscarella                               #
#  Last Update: 17 Jul 2015                                                    #
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

################################################################################
# Starvation Traits Experiment Data ############################################
################################################################################

# The following lines of code import and analyze the growth curve data for
# the Starvation Traits Experiment by Rachel Ferrill.
# The file names indicate the date the experiments were run and the plate.
# The plate information correspondes to the IDs for each experiment.
# Because of some minor temperature changes in the lab, we changed the
# temperature parameter (delta.temp) from 3 degresse to 4.

#  RNF Data
growth.modGomp("./data/GrowthCurves/RNF_GrowthCurve_20150624_plate1.txt",
               "plate1", skip=32, delta.temp = 4)

growth.modGomp("./data/GrowthCurves/RNF_GrowthCurve_20150625_plate2.txt",
               "plate2", skip=32, delta.temp = 4)

growth.modGomp("./data/GrowthCurves/RNF_GrowthCurve_20150626_plate3.txt",
               "plate3", skip=32, delta.temp = 4)

growth.modGomp("./data/GrowthCurves/RNF_GrowthCurve_20150626_plate4.txt",
               "plate4", skip=32, delta.temp = 4)

growth.modGomp("./data/GrowthCurves/RNF_GrowthCurve_20150719_plate5.txt",
               "plate5", skip=32, delta.temp = 4)

growth.modGomp("./data/GrowthCurves/RNF_GrowthCurve_20150701_plate6.txt",
               "plate6", skip=32, delta.temp = 4)

