#' run distribution analysis
#'
#' @param data_dir directory of the 4CE tables
#' @param siteid label specifying the site
#' @param long.thres "long" stage threshold
#'
#' @return phenotype data and count tables
#' @export

runAnalysis_distribution <- function(data_dir, siteid, long.thres){

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

  # blur_it_matrix <- function(matrix, blur_abs, mask_thres){
  #   for (i in 1:dim(matrix)[3]) {
  #     for (j in 1:dim(matrix)[2]) {
  #       for (k in 1:dim(matrix)[1]) {
  #         matrix[k,j,i,] <- matrix[k,j,i,]+sample(seq(-blur_abs,blur_abs), 1, replace = TRUE)
  #         matrix[k,j,i,] <- ifelse(matrix[k,j,i,] < mask_thres, 0, matrix[k,j,i,])
  #       }
  #     }
  #   }
  #   matrix
  # }

  blur_abs <- 2
  mask_thres <- FourCePhase2.1Data::getObfuscation(toupper(siteid))
  if (is_empty(mask_thres)) {
    mask_thres <- 0
  }

  distribution.output <- phen_distribution(PatientObservations, PatientSummary, rules, siteid, long.thres)
  Phen.data.tot <- distribution.output$Phen.data.tot
  tot.count <- distribution.output$tot.count
  tot.count.blur <- blur_it(tot.count,c("before_adm","dayN14toN1","day0to29","day30to59", "day60to89", "day90plus"), blur_abs, mask_thres)
  post.count <- distribution.output$post.count
  post.count.blur <- blur_it(post.count, "n", blur_abs, mask_thres)
  post.prev.blur <- post.count.blur%>%dplyr::mutate(perc = 100*n/n.tot)
  # post.prev <- post.prev[order(post.prev$perc, decreasing = TRUE),]

  # Phen.data.tot$siteid <- rep(siteid, nrow(Phen.data.tot))
  # tot.count.sex<-table(Phen.data.tot$phenotype, Phen.data.tot$stage, Phen.data.tot$sex, Phen.data.tot$siteid)
  # tot.count.sex.blur <- blur_it_matrix(tot.count.sex, blur_abs,mask_thres)
  # tot.count.age<-table(Phen.data.tot$phenotype, Phen.data.tot$stage, Phen.data.tot$age_group, Phen.data.tot$siteid)
  # tot.count.age.blur <- blur_it_matrix(tot.count.age, blur_abs,mask_thres)
  # tot.count.severity<-table(Phen.data.tot$phenotype, Phen.data.tot$stage, Phen.data.tot$severe, Phen.data.tot$siteid)
  # tot.count.severity.blur <- blur_it_matrix(tot.count.severity, blur_abs,mask_thres)
  # # tot.count.race<-table(Phen.data.tot$phenotype, Phen.data.tot$stage, Phen.data.tot$race, Phen.data.tot$siteid)
  # tot.count.race <- blur_it_matrix(tot.count.race, blur_abs,mask_thres)
#
#   if (is.null(output_dir)) {
#     output_dir <- getProjectOutputDirectory()
#   }
#
#   save(tot.count.blur,post.prev.blur,tot.count.sex.blur, tot.count.age.blur,tot.count.severity.blur,
#        file = paste(paste(output_dir,"/",sep=""),paste(siteid,"results_PhenDistribution.RData",sep="_"),sep=""))

  return(list(
    phen.data = Phen.data.tot,
    tot.count = tot.count.blur,
    post.prev = post.prev.blur
  ))
}
