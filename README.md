# pascPhen Package
pascPhen is a package for the analysis of Post-Acute Sequelae of SARS-CoV-2 infection (PASC) rule-based phenotypes.
It provides functions for the description of the phenotype distribution and the application of MLHO framework with the aim of extracting features that can contribute to the re-definition of the phenotype.

## Installation
You can install the package from [Github](https://github.com/rebeccamesa/pascPhen) with:

`devtools::install_github("rebeccamesa/pascPhen")`

## Analysis
2 main analysis are implemented:
* distribution of the phenotypes
* application of MLHO framework to a specific phenotype

Please run the first analysis with `pascPhen::runAnalysis_distribution(data_dir, output_dir, siteid)`, where *data_dir* is the 4CE data directory, *output_dir* is the customized directory for the results and *siteid* is the label specifying the site. If *output_dir* is null, results are saved in the `getProjectOutputDirectory()`directory. 
