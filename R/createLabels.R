#' Build labeldt for mlho_features
#' 
#' Build a table containing for each patient a label that specify whether the patient
#' matches the phenotype for the first time \code{days.long} after positive test
#' 
#' @param PatientObservations 2.1 PatientObservations table
#' @param days.long "long" stage threshold
#' @param rulesPhen table containing the rules for a specific phenotype
#' 
#' @return labeldt
#' @export

createLabels <- function(days.long,
                         rulesPhen,
                         PatientObservations){
  
  dOrder <- PatientObservations[order(PatientObservations$patient_num, PatientObservations$days_since_admission),]
  PatientObservations <- dplyr::distinct(dOrder, patient_num, concept_code, value, .keep_all = TRUE)
  PatientObservations$label <- 0
  
  PatientObservations$label[PatientObservations$days_since_admission>days.long & paste(gsub(".*-", "", PatientObservations$concept_type), PatientObservations$concept_code) %in% paste(rulesPhen$Coding, rulesPhen$Code)] <- 1
  labeldt <- PatientObservations%>%
    dplyr::select(patient_num,label)%>%
    dplyr::group_by(patient_num)%>%
    dplyr::summarize_all(max)%>%
    dplyr::mutate(label = as.factor(label))
  
  return(labeldt)
  
}
