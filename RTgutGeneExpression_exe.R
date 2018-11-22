#############################################
# Analyse qPCR data from lightcycler 480    #
#                                           #
# EXECUTION FILE                            #
# this is the execution file which requires #
# the following libary file for running:    #
# RTgutGeneExpression_libary.R              #
#                                           #
# Authors: Hannah Schug                     #
#                                           #
#############################################

# INSTRUCTIONS
# place these two R files and the data you want to have analysed in one folder
# .txt files saved from the LightCycler 480 software serve as input files 
# example data files are provided on GitHub as well

# IMPORTANT
# please note: The script works only for files with EXACTLY the same data format as the example files


# clear the workspace
rm(list = ls())

# set working directory to source file location: Sessions --> Set wd --> to source file location
# show the files available in the working directory
getwd()
wd <- setwd(getwd())
list.files()

# The library file qPCR_library.R contains all code for reading in the data from .txt, background substractions, 
# calculations of quantification cycles (Cq`s), normalisations and saving everything to an excel 
# it needs to be loaded before running the function
source("RTgutGeneExpression_library.R")


#############################################
# PLEASE MODIFY 
InName_data <- "InputData.txt" # name of data input file
InName_info <- "InputSampleInfo.txt" # name of data sample info 
title_name <- "OutputExampleDataRTgutGeneExpression" # name for titles of e.g. graphs
OutName <- "OutputExampleDataRTgutGeneExpression.xlsx" # name for output excel file


#############################################
# Do NOT modify
replicates <- 3 # Numer of technical replicates on the plate
n_std <- 5 # Number of concentrations used for the standard curve
n_samples <- 12 # Number of samples 
n_controls <- 2 # Number of controls, usually RT and H20

# running the wrapper function will execute the code
RTgutGeneExpression_wrapper(InName_data, InName_info, title_name,replicates, n_std, n_samples, n_controls, OutName)
