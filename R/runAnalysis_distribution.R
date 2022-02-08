#' run distribution analysis
#'
#' @param data_dir directory of the 4CE tables
#' @param output_dir directory of the analysis results. If NULL results are saved in the `getProjectOutputDirectory()` directory
#' @param siteid label specifying the site
#'
#' @return NULL.
#' @export

runAnalysis_distribution <- function(data_dir, output_dir=NULL, siteid){

  PatientObservations <- utils::read.csv(file.path(data_dir, "LocalPatientObservations.csv"))
  PatientSummary <- utils::read.csv(file.path(data_dir, "LocalPatientSummary.csv"))
  utils::data("rules")

  blur_it <- function(df, vars, blur_abs, mask_thres){
    # Obfuscate count values.
    # If blurring range is +/-3, or blur_abs = 3,
    # the count receive a small addition of a random number from -3 to 3.
    # If a count is less than mask_thres, set that count to 0.

    for (var in vars){
      var <- sym(var)
      blur_vec <- sample(seq(- blur_abs, blur_abs), nrow(df), replace = TRUE)
      df <- df %>%
        dplyr::mutate(!!var := !!var + blur_vec,
               !!var := ifelse(!!var < mask_thres, 0, !!var))
    }
    df
  }

  blur_it_matrix <- function(matrix, blur_abs, mask_thres){
    for (i in 1:dim(matrix)[3]) {
      for (j in 1:dim(matrix)[2]) {
        for (k in 1:dim(matrix)[1]) {
          matrix[k,j,i,] <- matrix[k,j,i,]+sample(seq(-blur_abs,blur_abs), 1, replace = TRUE)
          matrix[k,j,i,] <- ifelse(matrix[k,j,i,] < mask_thres, 0, matrix[k,j,i,])
        }
      }
    }
    matrix
  }

  blur_abs <- 2
  mask_thres <- FourCePhase2.1Data::getObfuscation(toupper(siteid))
  if (is_empty(mask_thres)) {
    mask_thres <- 0
  }

  distribution.output <- phen_distribution(PatientObservations, PatientSummary, rules,siteid)
  tot.count <- distribution.output$tot.count
  tot.count <- blur_it(tot.count,c("before_adm","dayN14toN1","day0to29","day30to89","day90plus"), blur_abs, mask_thres)
  print(tot.count)
  tot.count1stOcc <- distribution.output$first.count
  tot.count1stOcc <- blur_it(tot.count1stOcc,c("before_adm","dayN14toN1","day0to29","day30to89","day90plus"), blur_abs, mask_thres)
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
  tot.count.sex <- distribution.output.sex$tot.count.strat
  tot.count.sex <- blur_it_matrix(tot.count.sex, blur_abs,mask_thres)
  tot.count1stOcc.sex <- distribution.output.sex$first.count.strat
  tot.count1stOcc.sex <- blur_it_matrix(tot.count1stOcc.sex, blur_abs,mask_thres)
  phen.barPlot.sex <- distribution.output.sex$phen.barPlot
  print(phen.barPlot.sex)
  phen1stOcc.barPlot.sex <- distribution.output.sex$phen.1stBarPlot
  print(phen1stOcc.barPlot.sex)

  distribution.output.age<-phen_distribution(PatientObservations, PatientSummary, rules, siteid, stratified_by = "age_group")
  tot.count.age <- distribution.output.age$tot.count.strat
  tot.count.age <- blur_it_matrix(tot.count.age, blur_abs,mask_thres)
  tot.count1stOcc.age <- distribution.output.age$first.count.strat
  tot.count1stOcc.age <- blur_it_matrix(tot.count1stOcc.age, blur_abs,mask_thres)
  phen.barPlot.age <- distribution.output.age$phen.barPlot
  print(phen.barPlot.age)
  phen1stOcc.barPlot.age <- distribution.output.age$phen.1stBarPlot
  print(phen1stOcc.barPlot.age)

  distribution.output.severity<-phen_distribution(PatientObservations, PatientSummary, rules, siteid, stratified_by = "severe")
  tot.count.severity <- distribution.output.severity$tot.count.strat
  tot.count.severity <- blur_it_matrix(tot.count.severity, blur_abs,mask_thres)
  tot.count1stOcc.severity <- distribution.output.severity$first.count.strat
  tot.count1stOcc.severity <- blur_it_matrix(tot.count1stOcc.severity, blur_abs,mask_thres)
  phen.barPlot.severity <- distribution.output.severity$phen.barPlot
  print(phen.barPlot.severity)
  phen1stOcc.barPlot.severity <- distribution.output.severity$phen.1stBarPlot
  print(phen1stOcc.barPlot.severity)

  distribution.output.race<-phen_distribution(PatientObservations, PatientSummary, rules, siteid, stratified_by = "race")
  tot.count.race <- distribution.output.race$tot.count.strat
  tot.count.race <- blur_it_matrix(tot.count.race, blur_abs,mask_thres)
  tot.count1stOcc.race <- distribution.output.race$first.count.strat
  tot.count1stOcc.race <- blur_it_matrix(tot.count1stOcc.race, blur_abs,mask_thres)
  phen.barPlot.race <- distribution.output.race$phen.barPlot
  print(phen.barPlot.race)
  phen1stOcc.barPlot.race <- distribution.output.race$phen.1stBarPlot
  print(phen1stOcc.barPlot.race)

  if (is.null(output_dir)) {
    output_dir <- getProjectOutputDirectory()
  }

  save(tot.count,tot.count1stOcc,tot.count.sex, tot.count1stOcc.sex, tot.count.age, tot.count1stOcc.age,
       tot.count.race,tot.count1stOcc.race,tot.count.severity, tot.count1stOcc.severity,
       phen.distributionPlot,phen1stOcc.distributionPlot,phen.barPlot,phen1stOcc.barPlot,
       phen.barPlot.sex,phen1stOcc.barPlot.sex,phen.barPlot.age,phen1stOcc.barPlot.age,
       phen.barPlot.severity,phen1stOcc.barPlot.severity,phen.barPlot.race,phen1stOcc.barPlot.race,
       file = paste(paste(output_dir,"/",sep=""),paste(siteid,"results_PhenDistribution.RData",sep="_"),sep=""))
}
