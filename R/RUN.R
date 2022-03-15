#' RUN PASC analysis
#'
#' @param data_dir directory of the 4CE tables
#' @param output_dir directory of the output report. If NULL, the result is saved in 'getProjectOutput' directory
#' @param siteid label specifying the site
#' @param long.thres "long" stage threshold
#' @param MSMR.sparsity MSMR lite parameter
#' @param long.perc percentage criterion to choose which phenotypes are interesting
#'
#' @return NULL
#' @export

# # Set input data directory
# data_dir <- "C:/Users/rebim/Desktop/Tesi/PASC Analysis/Dati ICSM 4CE"
# # Set output results directory. If NULL results are saved in 'getProjectOutput' directory
# output_dir <- "C:/Users/rebim/Desktop/Tesi/PASC Analysis/Rule based Phenotyping/Output"
#
# siteid = "icsm"
#
# long.thres = 15
# MSMR.sparsity = 0.005
# long.perc <- 0.04

RUN <- function(data_dir, output_dir, siteid, long.thres, long.perc = 0.03, MSMR.sparsity, data_type){
  data_dir <- data_dir
  output_dir <- output_dir
  if (is.null(output_dir)) {
    output_dir <- getProjectOutputDirectory()
  }
  long.thres <- long.thres
  long.perc <- long.perc
  MSMR.sparsity <- MSMR.sparsity
  data_type <- data_type

  #rmarkdown::render(paste0(output_dir,"/Report.Rmd"),output_file = paste0("report.",siteid,".", Sys.Date(),".html"))
  rmarkdown::render(system.file("rmd", "Report.Rmd", package = "pascPhen"),
                    output_file = paste0("report_",siteid,"_",long.thres,"_",substr(long.perc,3,nchar(long.perc)),"_",gsub("\\.", "",data_type),"_", Sys.Date(),".html"),
                    output_dir = output_dir)
}

