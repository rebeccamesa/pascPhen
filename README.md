# pascPhen Package
pascPhen is a package for the analysis of Post-Acute Sequelae of SARS-CoV-2 infection (PASC) rule-based phenotypes.  
The main steps of the analysis are:
* Plot of the phenotypes distribution over different stages 
* Selection of the phenotypes with a high prevalence in the PASC window 
* Application of the MLHO for PheWAS framework to the selected phenotypes to find new features that can be included in the definition of the phenotypes  

The package will create automatically a html report with the results of all these steps.

## Run the analysis
To run the package:

```
if(!require(devtools)) devtools::install_github("hadley/devtools")
if(!require(mlho)) devtools::install_github("hestiri/mlho")
if(!require(pascPhen)) devtools::install_github("rebeccamesa/pascPhen")
if(!require(pacman)) install.packages("pacman")

pacman::p_load(data.table, devtools, backports, Hmisc, tidyr,dplyr,ggplot2,plyr,scales,readr,Rmisc,
               httr, DT, lubridate, tidyverse,reshape2,foreach,doParallel,caret,gbm,lubridate,praznik,
               ggridges, forcats, stats, FourCePhase2.1Data, ppcor,pascPhen, reshape2)

pascPhen::RUN_lite(data_dir = "...", # set data directory
    output_dir= ".../Output", # create an output directory where you want the html files saved
    siteid = "",
    long.thres = 60, #60 and 90
    long.perc = 0.03, #*100 (default)
    MSMR.sparsity = 0.005, #*100
    data_type = "")# 2.1, 2.2
```

### RUN_lite parameters

|**Parameter**|**Description**|
|:------------|:--------------|
|data_dir     |path to your 2.1 or 2.2 data|
|output_dir   |path to the output directory where you want the html reports saved|
|siteid       |your site id. **Important for 4CE obfuscation**
|long.thres   |PASC window. We will do 60 and 90 days|
|long.perc    |phenotype selection percentage. Default is 3%|
|MSMR.sparsity|parameter of MSMR. Set to 0.005|
|data_type    |2.1 or 2.2|


Please send your output reports via slack.


