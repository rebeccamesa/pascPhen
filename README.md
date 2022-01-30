# pascPhen Package
pascPhen is a package for the analysis of Post-Acute Sequelae of SARS-CoV-2 infection (PASC) rule-based phenotypes.
It provides functions for the description of the phenotype distribution and the application of MLHO framework with the aim of extracting features that can contribute to the re-definition of the phenotype.

## Installation
You can install the released version of pascPhen from [Github](https://github.com/rebeccamesa/pascPhen) with:

`devtools::install_github("rebeccamesa/pascPhen")`

## Data
You need three datasets:
* 2.1 *PatientObservations* table
* 2.1 *PatientSummary* table
* *rules* table with 4 columns:

  PASC_Phenotype | Coding | Code | Description
  :------------: | :----: | :--: | :---------:
  character | character | character | character
  
  You can load the *rules* table from the package with:
  
  `utils::data(rules)`
