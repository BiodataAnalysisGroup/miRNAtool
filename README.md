# miRNAtool
mirsRNAtool was implemented in [R](https://www.r-project.org/). 

## Structure
miRNAStool repository consists of three folders:
- **input**: This is where all the input files are stored.
- **R**: This is where all the R scripts are stored.
- **output**: This is where all the output files are stored. This folder is created automatically while executing the project. 

## Installing
The project can be downloaded using git:
```
git clone https://github.com/BiodataAnalysisGroup/miRNAtool.git
```


## Execution:
- **Packages**: Make sure that every package loaded in **lines 20-31** inside the *script_mirs.R* is installed. The package **saveImageHigh** is basically **nikopech**'s tool, so you need to install it from [here](https://github.com/nikopech/saveImageHigh).
- **Input**: The input parameters are specified inside the **input.yml** file. No need to specify parameters inside the main script.
- **Execution**: Just run *script_mirs.R*.
- **Output**: Inside the *'output/'* directory:
  1. Directories *'plate1/'*, *'plate2/'*  and *'plate3/'* contain the outputs from quality control analysis.
  2. Directory *'Differential Analysis/'* contains the output from differential analysis.
  3. Directory *'Functional Analysis/'* contains the output from functional analysis.
  4. Tables normalized_data, data.all, multimir_results.validated@summary$target_symbol **or** multimir_results.predicted@summary$target_symbol (depends on the input), KEGG_enrich.f and GO_enrich.f are stored inside the *'Tables/'* directory as .csv files.
  5. *report.txt'*: Contains information concerning the analysis (e.g. samples that were droped out due to high NA's percentage, errors in specific plots, etc) and execution time measurements.

## License
This project is licensed under the MIT License - see the [LICENSE](https://github.com/BiodataAnalysisGroup/miRNAtool/blob/main/LICENSE) file for details.
