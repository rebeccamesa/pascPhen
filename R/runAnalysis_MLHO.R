#' run MLHO analysis
#'
#' @param data_dir directory of the 4CE tables
#' @param output_dir directory of the analysis results. If NULL results are saved in the `getProjectOutputDirectory()` directory
#' @param siteid label specifying the site
#' @param long.thres "long" stage threshold
#' @param phenotype phenotype of interest
#'
#' @return phenotyping output and features prevalance
#' @export

runAnalysis_MLHO <- function(data_dir, output_dir = NULL, siteid, long.thres, phenotype){

  PatientObservations <- utils::read.csv(file.path(data_dir, "LocalPatientObservations.csv"))
  PatientSummary <- utils::read.csv(file.path(data_dir, "LocalPatientSummary.csv"))
  utils::data("rules")

  phenotyping.output <- mlho_features(PatientObservations, PatientSummary, rules, days.long=long.thres, phenotype=phenotype, MSMR.topn=50, classifier="GLM")
  phenotyping.features <- phenotyping.output$model.output
  phenotyping.features$type <- rep("phenotyping", nrow(phenotyping.features))
  phenotyping.features$site <- rep(siteid, nrow(phenotyping.features))

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

  blur_abs <- 2
  mask_thres <- FourCePhase2.1Data::getObfuscation(toupper(siteid))
  if (is_empty(mask_thres)) {
    mask_thres <- 0
  }

  perc.output <- phenotyping.output$perc.output
  perc.output.blur <- blur_it(perc.output, "n.tot", "n", blur_abs, mask_thres)
  perc.output.blur <- dplyr::mutate(perc.output.blur, perc = n/n.tot)
  perc.output$type <- rep("phenotyping", nrow(perc.output))
  perc.output$site <- rep(siteid, nrow(perc.output))

  if (is.null(output_dir)) {
    output_dir <- getProjectOutputDirectory()
  }

  save(phenotyping.features, perc.output.blur, file = paste(paste(output_dir,"/",sep=""),paste(siteid,"results_MLHOphen.RData",sep="_"),sep=""))

  return(list(
    phenotyping.features = phenotyping.features,
    perc.output = perc.output))
}
