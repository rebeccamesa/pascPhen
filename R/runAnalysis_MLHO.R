#' run MLHO analysis
#'
#' @param data_dir directory of the 4CE tables
#' @param siteid label specifying the site
#' @long.thresh "long" stage threshold
#' @phenotype phenotype of interest
#'
#' @return NULL. Results are saved in the `getProjectOutputDirectory()` directory
#' @export

runAnalysis_MLHO <- function(data_dir, siteid, long.thresh, phenotype){

  PatientObservations <- utils::read.csv(file.path(data_dir, "LocalPatientObservations.csv"))
  PatientSummary <- utils::read.csv(file.path(data_dir, "LocalPatientSummary.csv"))
  utils::data("rules")

  print("Phenotyping:")
  phenotyping.output <- mlho_features(PatientObservations, PatientSummary, rules, days.long=long.thresh, phenotype=phenotype, MSMR.topn=50, classifier="GLM")
  phenotyping.output$type <- rep("phenotyping", nrow(phenotyping.output))
  phenotyping.output$site <- rep(siteid, nrow(phenotyping.output))
  print(phenotyping.output)

  print("Prediction:")
  prediction.output <- mlho_features(PatientObservations, PatientSummary, rules, days.long=long.thresh, phenotype=phenotype, MSMR.topn=50, classifier="GLM", prediction = TRUE)
  phenotyping.output$type <- rep("prediction", nrow(prediction.output))
  phenotyping.output$site <- rep(siteid, nrow(prediction.output))
  print(prediction.output)

  save(phenotyping.output,prediction.output,file = paste(paste(getProjectOutputDirectory(),"/",sep=""),paste(siteid,"results_MLHO.RData",sep="_"),sep=""))
}
