################################################################################
#                                                                              #
#  Test Script for Analysis of PreSens Respiration Data                        #
#   Version 2.0                                                                #
#  Written By: Mario Muscarella                                                #
#  Last Update: 28 Jan 2014                                                    #
#                                                                              #
#  Use this file to test the PreSens.Respiration fucntion on local machine     #
#  And as a template to create your own analysis                               #
#                                                                              #
################################################################################

setwd("~/GitHub/StarvationTraits/")
rm(list=ls())

# Inport the function from source file
source("./bin/PreSensInteractiveRegression.R")

################################################################################
# Examples #####################################################################
################################################################################

# Example txt analysis
PreSens.Respiration(infile = "./data/Respiration/ExampleData.txt",
                    outfile = "./data/Respiration/ExampleData_Output.txt")

# Example txt analysis
PreSens.Respiration(infile = "./data/Respiration/20150701_BacterialRespiration_a_RNF_Oxygen.txt",
                    outfile = "./data/Respiration/20150701_BacterialRespiration_a_RNF_Output.txt")

