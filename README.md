# miRNAtool
mirsRNAtool was implemented in [R](https://www.r-project.org/). 

## Structure
miRNAStool repository consists of three folders:
- `data`: This is where all the input files are stored.
- `R`: This is where all the R scripts are stored.
- `output`: This is where all the output files are stored. This folder is created automatically while executing the project. 

## Getting started
### Dependencies
Execute the following lines to install the required packages:
- from CRAN:

```
install.packages(c("data.table", "tidyverse", "gsubfn", "yaml", "enrichR", "lubridate"))
```

- from Bioconductor:

```
BiocManager::install(c("naniar", "prada", "ComplexHeatmap", "multiMiR", "limma"))
```


### Installing
The project can be downloaded using git:
```
git clone https://github.com/BiodataAnalysisGroup/miRNAtool.git
```

### Running the project
The framework consists of six scripts:
- ```script_mirs.R```
- ```qc.R``` 
- ```diff_analysis.R```
- ```fa.R```
- ```save_as_pdf.R``` 
- ```save_image.R```

In order to run the project:
1. Set the [R](https://github.com/BiodataAnalysisGroup/miRNAtool/tree/main/R) folder as your working directory
2. Place the required input files inside the `data` directory and set the parameters inside the .yaml file. See instructions [here](https://github.com/BiodataAnalysisGroup/miRNAtool/tree/main/data). 
3. Use the following command:
```
source("script_mirs.R")
```

## Additional scripts
- `initial_script.R`: Contains the full code of the initial version of the project. There is no need to execute this file.

## Input
- For information concerning the input please move to the `data` [folder](https://github.com/BiodataAnalysisGroup/miRNAtool/tree/main/data) and checkout the correpsonding README.md file.

## Output
- For information concerning the output please move to the `output` [folder](https://github.com/BiodataAnalysisGroup/miRNAtool/tree/main/output) and checkout the correpsonding README.md file


## Acknowledgements
- [nikopech](https://github.com/nikopech/)  
Scripts `save_as_pdf.R` and `save_image.R` were taken from nikopech's project [saveImageHigh](https://github.com/nikopech/saveImageHigh)

## License
This project is licensed under the MIT License - see the [LICENSE](https://github.com/BiodataAnalysisGroup/miRNAtool/blob/main/LICENSE) file for details.
