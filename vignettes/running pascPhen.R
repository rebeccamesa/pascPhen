##implementing Rebecca's package

##some library installations
# If you don't have devtools, see here: https://www.r-project.org/nosvn/pandoc/devtools.html


if(!require(devtools)) devtools::install_github("hadley/devtools")
if(!require(mlho)) devtools::install_github("hestiri/mlho")
if(!require(pascPhen)) devtools::install_github("rebeccamesa/pascPhen")
if(!require(pacman)) install.packages("pacman")

pacman::p_load(data.table, devtools, backports, Hmisc, tidyr,dplyr,ggplot2,plyr,scales,readr,Rmisc,
               httr, DT, lubridate, tidyverse,reshape2,foreach,doParallel,caret,gbm,lubridate,praznik,
               ggridges, forcats, stats, FourCePhase2.1Data, ppcor,pascPhen, reshape2)


#
pascPhen::RUN_lite(data_dir = "...", ##set data directory
    output_dir= ".../Output", #create an output directory where you want the html files saved
    siteid = "",
    long.thres = 60, #60 to begin
    long.perc = 0.03, #*100 (default)
    MSMR.sparsity = 0.005, #*100
    data_type = "")# 2.1, 2.2, 2.2.all



