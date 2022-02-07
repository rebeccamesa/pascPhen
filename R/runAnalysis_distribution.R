#' run distribution analysis
#'
#' @param data_dir directory of the 4CE tables
#' @param siteid label specifying the site
#'
#' @return NULL. Results are saved in the `getProjectOutputDirectory()` directory
#' @export

runAnalysis_distribution <- function(data_dir, siteid){

  PatientObservations <- utils::read.csv(file.path(data_dir, "LocalPatientObservations.csv"))
  PatientSummary <- utils::read.csv(file.path(data_dir, "LocalPatientSummary.csv"))
  utils::data("rules")

  distribution.output <- phen_distribution(PatientObservations, PatientSummary, rules,siteid)
  tot.count <- distribution.output$tot.count
  print(tot.count)
  tot.count1stOcc <- distribution.output$first.count
  print(tot.count1stOcc)
  phen.distributionPlot <- distribution.output$phen.distributionPlot
  print(phen.distributionPlot)
  phen1stOcc.distributionPlot <- distribution.output$phen.1stDistributionPlot
  print(phen1stOcc.distributionPlot)
  phen.barPlot <- distribution.output$phen.barPlot
  print(phen.barPlot)
  phen1stOcc.barPlot <- distribution.output$phen.1stBarPlot
  print(phen1stOcc.barPlot)

  distribution.output.sex <- phen_distribution(PatientObservations, PatientSummary, rules, siteid, stratified_by = "sex")
  phen.barPlot.sex <- distribution.output.sex$phen.barPlot
  print(phen.barPlot.sex)
  phen1stOcc.barPlot.sex <- distribution.output.sex$phen.1stBarPlot
  print(phen1stOcc.barPlot.sex)

  distribution.output.age<-phen_distribution(PatientObservations, PatientSummary, rules, siteid, stratified_by = "age_group")
  phen.barPlot.age <- distribution.output.age$phen.barPlot
  print(phen.barPlot.age)
  phen1stOcc.barPlot.age <- distribution.output.age$phen.1stBarPlot
  print(phen1stOcc.barPlot.age)

  distribution.output.severity<-phen_distribution(PatientObservations, PatientSummary, rules, siteid, stratified_by = "severe")
  phen.barPlot.severity <- distribution.output.severity$phen.barPlot
  print(phen.barPlot.severity)
  phen1stOcc.barPlot.severity <- distribution.output.severity$phen.1stBarPlot
  print(phen1stOcc.barPlot.severity)

  distribution.output.race<-phen_distribution(PatientObservations, PatientSummary, rules, siteid, stratified_by = "race")
  phen.barPlot.race <- distribution.output.race$phen.barPlot
  print(phen.barPlot.race)
  phen1stOcc.barPlot.race <- distribution.output.race$phen.1stBarPlot
  print(phen1stOcc.barPlot.race)

  save(tot.count,tot.count1stOcc,phen.distributionPlot,phen1stOcc.distributionPlot,phen.barPlot,phen1stOcc.barPlot,
       phen.barPlot.sex,phen1stOcc.barPlot.sex,phen.barPlot.age,phen1stOcc.barPlot.age,
       phen.barPlot.severity,phen1stOcc.barPlot.severity,phen.barPlot.race,phen1stOcc.barPlot.race,
       file = paste(paste(getProjectOutputDirectory(),"/",sep=""),paste(siteid,"results_distribution.RData",sep="_"),sep=""))
}
