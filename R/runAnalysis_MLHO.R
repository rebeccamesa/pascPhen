#' run MLHO analysis
#'
#' @param data_dir directory of the 4CE tables
#' @param output_dir directory of the analysis results. If NULL results are saved in the `getProjectOutputDirectory()` directory
#' @param siteid label specifying the site
#' @long.thres "long" stage threshold
#' @phenotype phenotype of interest
#'
#' @return NULL.
#' @export

runAnalysis_MLHO <- function(data_dir, output_dir = NULL, siteid, long.thres, phenotype){

  PatientObservations <- utils::read.csv(file.path(data_dir, "LocalPatientObservations.csv"))
  PatientSummary <- utils::read.csv(file.path(data_dir, "LocalPatientSummary.csv"))
  utils::data("rules")

  print("Phenotyping:")
  phenotyping.output <- mlho_features(PatientObservations, PatientSummary, rules, days.long=long.thres, phenotype=phenotype, MSMR.topn=50, classifier="GLM")
  phenotyping.output$type <- rep("phenotyping", nrow(phenotyping.output))
  phenotyping.output$site <- rep(siteid, nrow(phenotyping.output))
  print(phenotyping.output)

  print("Prediction:")
  prediction.output <- mlho_features(PatientObservations, PatientSummary, rules, days.long=long.thres, phenotype=phenotype, MSMR.topn=50, classifier="GLM", prediction = TRUE)
  phenotyping.output$type <- rep("prediction", nrow(prediction.output))
  phenotyping.output$site <- rep(siteid, nrow(prediction.output))
  print(prediction.output)

  if (is.null(output_dir)) {
    output_dir <- getProjectOutputDirectory()
  }

  save(phenotyping.output,prediction.output,file = paste(paste(output_dir,"/",sep=""),paste(siteid,"results_MLHOphen.RData",sep="_"),sep=""))
}
