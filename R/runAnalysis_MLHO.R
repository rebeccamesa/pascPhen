#' run MLHO analysis
#'
#' @param data_dir directory of the 4CE tables
#' @param siteid label specifying the site
#' @param long.thres "long" stage threshold
#' @param phenotype phenotype of interest
#' @param MSMR.sparsity MSMR lite parameter
#' @param data_type string with the type of 4CE data
#'
#' @return phenotyping output and features prevalance
#' @export

runAnalysis_MLHO <- function(data_dir, siteid, long.thres, phenotype, MSMR.sparsity, data_type){

  PatientObservations <- utils::read.csv(file.path(data_dir, "LocalPatientObservations.csv"))
  PatientSummary <- utils::read.csv(file.path(data_dir, "LocalPatientSummary.csv"))
  if (data_type != "2.1") {
    data <- filter_data(PatientObservations, PatientSummary, data_type)
    PatientObservations <- data$PatientObservations
    PatientSummary <- data$PatientSummary
  }
  utils::data("rules")

  PatientObservations$concept_code <- gsub('[.]', '', PatientObservations$concept_code)

  print(paste0("Phenotyping .....",phenotype, "  (─‿─)  "))
  phenotyping.output <- mlho_features(PatientObservations, PatientSummary, rules, days.long=long.thres, phenotype=phenotype, MSMR.topn=100, MSMR.sparsity = MSMR.sparsity)
  phenotyping.features <- phenotyping.output$model.output
  phenotyping.features$type <- rep("phenotyping", nrow(phenotyping.features))
  phenotyping.features$site <- rep(siteid, nrow(phenotyping.features))
  phenotyping.features <- as.data.frame(phenotyping.features)

  corr.matrix <- phenotyping.output$mlho.corr
  # print(phenotyping.output)

  # blur_it <- function(df, vars, blur_abs, mask_thres){
  #   # Obfuscate count values.
  #   # If blurring range is +/-3, or blur_abs = 3,
  #   # the count receive a small addition of a random number from -3 to 3.
  #   # If a count is less than mask_thres, set that count to 0.
  #
  #   for (var in vars){
  #     var <- sym(var)
  #     blur_vec <- sample(seq(- blur_abs, blur_abs), nrow(df), replace = TRUE)
  #     df <- df %>%
  #       dplyr::mutate(!!var := !!var + blur_vec,
  #                     !!var := ifelse(!!var < mask_thres, 0, !!var))
  #   }
  #   df
  # }
  #
  # blur_abs <- 2
  # mask_thres <- FourCePhase2.1Data::getObfuscation(toupper(siteid))
  # if (is_empty(mask_thres)) {
  #   mask_thres <- 0
  # }
  #
  # perc.output <- phenotyping.output$perc.output
  # perc.output.blur <- blur_it(perc.output, c("n.tot", "n"), blur_abs, mask_thres)
  # perc.output.blur <- dplyr::mutate(perc.output.blur, perc = n/n.tot)
  # perc.output$type <- rep("phenotyping", nrow(perc.output))
  # perc.output$site <- rep(siteid, nrow(perc.output))


  # print("Prediction:")
  # prediction.output <- mlho_features(PatientObservations, PatientSummary, rules, days.long=long.thres, phenotype=phenotype, MSMR.topn=50, classifier="GLM", prediction = TRUE)
  # phenotyping.output$type <- rep("prediction", nrow(prediction.output))
  # phenotyping.output$site <- rep(siteid, nrow(prediction.output))
  # print(prediction.output)

  # if (is.null(output_dir)) {
  #   output_dir <- getProjectOutputDirectory()
  # }

  return(list(
    output = phenotyping.features,
    corr = corr.matrix))

  # save(phenotyping.output,file = paste(paste(output_dir,"/",sep=""),paste(siteid,"results_MLHOphen.RData",sep="_"),sep=""))
}
