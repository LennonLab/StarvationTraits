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
#  Template for Respiration data                                               #
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

# Example txt analysis
PreSens.Respiration(infile = "./data/Respiration/20150701_BacterialRespiration_b_RNF_Oxygen.txt",
                    outfile = "./data/Respiration/20150701_BacterialRespiration_b_RNF_Output.txt")


# Example txt analysis
PreSens.Respiration(infile = "./data/Respiration/20150701_BacterialRespiration_c_RNF_Oxygen.txt",
                    outfile = "./data/Respiration/20150701_BacterialRespiration_c_RNF_Output.txt")

PreSens.Respiration(infile = "./data/Respiration/20150707_BacterialRespiration_d_RNF_Oxygen.txt",
                    outfile = "./data/Respiration/20150707_BacterialRespiration_d_RNF_Output.txt")

PreSens.Respiration(infile = "./data/Respiration/20150707_BacterialRespiration_e_RNF_Oxygen.txt",
                    outfile = "./data/Respiration/20150707_BacterialRespiration_e_RNF_Output.txt")

PreSens.Respiration(infile = "./data/Respiration/20150710_BacterialRespiration_f_RNF_Oxygen.txt",
                    outfile = "./data/Respiration/20150710_BacterialRespiration_f_RNF_Output.txt")

PreSens.Respiration(infile = "./data/Respiration/20150710_BacterialRespiration_g_RNF_Oxygen.txt",
                    outfile = "./data/Respiration/20150710_BacterialRespiration_g_RNF_Output.txt")

PreSens.Respiration(infile = "./data/Respiration/20150714_BacterialRespiration_h_RNF_Oxygen.txt",
                    outfile = "./data/Respiration/20150714_BacterialRespiration_h_RNF_Output.txt")

PreSens.Respiration(infile = "./data/Respiration/20150715_BacterialRespiration_i_RNF_Oxygen.txt",
                    outfile = "./data/Respiration/20150715_BacterialRespiration_i_RNF_Output.txt")

PreSens.Respiration(infile = "./data/Respiration/20150715_BacterialRespiration_j_RNF_Oxygen.txt",
                    outfile = "./data/Respiration/20150715_BacterialRespiration_j_RNF_Output.txt")

PreSens.Respiration(infile = "./data/Respiration/20150716_BacterialRespiration_k_RNF_Oxygen.txt",
                    outfile = "./data/Respiration/20150716_BacterialRespiration_k_RNF_Output.txt")

PreSens.Respiration(infile = "./data/Respiration/20150716_BacterialRespiration_l_RNF_Oxygen.txt",
                    outfile = "./data/Respiration/20150716_BacterialRespiration_l_RNF_Output.txt")