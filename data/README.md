# data
This is the folder where all input files are stored, for the project to run properly.

## Input data
- The input data are basically all the **.csv** files: `miRNome_data.csv`, `miRs_annotation_3plates.csv` and `phenodata.csv`.  
Hence, the input files should be of the same format, having the same structure and definetely the same names.
- Be careful with the `phenodata.csv`. This file is the metadata file, mapping samples into groups. The column names should remain as it is. The prefix of sample names should be the corresponding name of group. For example, here we have **normal** and **cancer** data and, hence, the sample names are **Normal.1**, **Cancer1.1** etc. It doens't matter whether letters are uppercase or lowercase.

## User's parameters
- The users need to specify their own parameters inside the `input.ylm` file.
