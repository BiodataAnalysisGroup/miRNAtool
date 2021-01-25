# clear
cat("\014")
rm(list = ls())
dev.off(dev.list()["RStudioGD"])
getwd()
analysis.path=getwd()  ;analysis.path
#analysis.path <- "D://exosomes//data" 

#setwd(analysis.path)
getwd()

# Source files
# qc.R : quality control
# diff_analysis.R : differential analysis
# fa.R : functional analysis
source("qc.R")
source("diff_analysis.R")
source("fa.R")

# Libraries
library(data.table)
library(tidyverse)
library(naniar)
library(prada)
library(gsubfn)
library(yaml)
library(saveImageHigh)
library(ComplexHeatmap)
library(multiMiR)
library(enrichR)
library(limma)

################# SETTING UP OUTPUT DIRECTORIES AND REPORT FILE ################
# Creating output folder
dir.create("output")
for (plate in plates){
  to_create <- paste('output/plate', plate, sep = '')
  dir.create(to_create)
}
dir.create("output/Differential Analysis")
dir.create("output/Functional Analysis")
dir.create("output/Tables")

# Creating the report file
report_file <- 'output/report.txt'
file.create(report_file)


############## LOADING INPUTS FROM YAML FILE ###################################
input_parameters <- read_yaml('input.yml')

# Num of plates
plates <- input_parameters$plates

# Threshold percentage (values: from 0 to 1)
# In case the percentage of NA's is higher than threshold, the sample is excluded
na_threshold_perc <- input_parameters$na_threshold
to_report <- paste('Selected NAs percentage threshold: ', na_threshold_perc, sep = '')
cat(to_report, file = report_file, sep = '\n', append = TRUE)

# RTC criterion threshold: THIS IS STANDARD
rtc_threshold <- 5
to_report <- paste('Standard RTC threshold: ', rtc_threshold, sep = '')
cat(to_report, file = report_file, sep = '\n', append = TRUE)

# Select endogenous or exogenous normalization
normalization_en_ex <- input_parameters$normalization_en_ex

# Select sign.f.table p-val criterion
sign_table_pval <- input_parameters$sign_table_pval

# Mirs: select validated or predicted
validated_or_predicted <- input_parameters$validated_or_predicted

# select go_enrich criterion
go_criterion <- input_parameters$go_criterion

# select the KEGG enrich criterion
kegg_enrich_criterion <- input_parameters$kegg_enrich_criterion



######################### MAIN CODE ############################################

list.files()

mirs<- fread("miRs_annotation_3plates.csv")
data<- fread("miRNome_data.csv")
meta<- fread("phenodata.csv")

head(meta)

data<- as.data.frame(data)
data$ID<- mirs$`miRNA ID`
data$plate<- mirs$Plate
head(data)

# Substitutes the ',' charachter with '.' and then converts strings to numeric data
# 2 stands for "iterate over columns"
data[,c(2:7)]<- apply(apply(data[,c(2:7)], 2, gsub, patt=",", replace="."), 2, as.numeric)

# Drop data with "blank" ID
data <- data[c(which(data$ID != "blank")), ]
write.csv(data, 'total_data.csv')

#TO DO: merge the three plates before the diff analysis
# Initial sample names
normalized_data <- NULL

# QC analysis
for (plate in plates){
  print(paste('Analysis for plate ', plate, sep = ''))
  norm_data <- QC(copy(data), plate, 'output', na_threshold_perc, rtc_threshold ,normalization_en_ex, report_file)
  norm_data$plate <- paste('Plate', plate, sep = ' ')
  
  if (is.null(normalized_data)){
    normalized_data <- rbind(normalized_data, norm_data)
    
  } else {
    normalized_data <- normalized_data[,which(colnames(normalized_data) %in% colnames(norm_data))]
    norm_data <- norm_data[,which(colnames(norm_data) %in% colnames(normalized_data))]
    normalized_data <- rbind(normalized_data, norm_data)
  }
}

if (!is.null(normalized_data)){
  
  # diff analysis and functional analysis
  normalized_data<- normalized_data[!duplicated(normalized_data$ID),]
  write.csv(normalized_data, 'output/Tables/normalized_data.csv')
  
  sign.table.f <- diff_analysis(normalized_data[,-dim(normalized_data)[2]], meta, 'output/Differential Analysis', sign_table_pval)
  functional_analysis(sign.table.f, validated_or_predicted, kegg_enrich_criterion, go_criterion, 'output/Functional Analysis')

  } else{
  # Todo: add some warnings - tasks here
  print("All plates were rejected from quality control analysis")
}