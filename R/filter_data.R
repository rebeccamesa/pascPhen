#' Filter data according to type
#'
#' @param PatientObservations 4CE table
#' @param PatientSummary 4CE table
#' @param data_type string with the type of 4CE data
#'
#' @return filtered 4CE tables
#' @export
#'
filter_data <- function(PatientObservations, PatientSummary, data_type){

  if (data_type == "2.2") {
    PatientSummary <- PatientSummary%>%dplyr::filter(cohort %like% "PosAdm")
  } else{
    if (data_type == "2.2.all") {
      PatientSummary<-PatientSummary%>%dplyr::filter(cohort %like% "PosAdm" | cohort %like% "PosNotAdm")
    }
  }
  pts <- unique(PatientSummary$patient_num)
  PatientObservations <- PatientObservations%>%dplyr::filter(patient_num %in% pts)

  return(list(PatientObservations = PatientObservations,
              PatientSummary = PatientSummary))
}
