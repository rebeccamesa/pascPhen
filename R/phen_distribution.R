#' Phenotype distribution
#'
#' Build a ridgeline plot and a box plot of the distribution of the phenotypes
#' over the dataset
#'
#' @param PatientObservations 2.1 PatientObservations table
#' @param PatientSummary 2.1 PatientSummary table
#' @param rules table containing the phenotype definitions
#' @param stratified_by label for the stratification of the bar plot
#'
#' @return two plots
#' @export

phen_distribution <- function(PatientObservations,
                              PatientSummary,
                              rules,
                              stratified_by = NULL)
{
  # Prepare data

  PatientObservations <- dplyr::select(PatientObservations, -value,-siteid)

  data <- PatientSummary %>%
    dplyr::select(patient_num, severe, deceased, sex, age_group, race, admission_date) %>%
    merge(PatientObservations, by = "patient_num")
  data$calendar_date <- as.Date(data$admission_date)+data$days_since_admission
  data <- data %>%
    dplyr::mutate(stage = dplyr::case_when(days_since_admission >= 90 ~ "day90plus",
                             days_since_admission >= 30 ~ "day30to89",
                             days_since_admission >= 0 ~ "day0to29",
                             days_since_admission >=-14 ~ "dayN14toN1",
                             TRUE ~ "before_adm"))%>%
    dplyr::mutate(stage = factor(stage, levels = c("before_adm", "dayN14toN1", "day0to29", "day30to89", "day90plus")))

  # Compute occurences

  phenotypes <- unique(rules$PASC_Phenotype)
  tot.count <- data.frame()
  Phen.data.tot <- data.frame()
  for (i in 1:length(phenotypes)) {
    Phen.rules <- dplyr::filter(rules, PASC_Phenotype == phenotypes[i])
    Phen.data <- data[paste(gsub(".*-", "", data$concept_type), data$concept_code) %in% paste(Phen.rules$Coding, Phen.rules$Code),]
    counts <- table(dplyr::select(Phen.data,stage))
    if (sum(counts)>0) {
      tot.count<-rbind(tot.count,c(phenotypes[i], counts))
      Phen.data$phenotype <- paste(phenotypes[i])
      Phen.data.tot<-rbind(Phen.data.tot, Phen.data)
    }
  }
  names(tot.count)<- c("Phenotype", "before_adm", "dayN14toN1", "day0to29", "day30to89", "day90plus")
  tot.count<-tot.count[order(tot.count$Phenotype),]

  # Ridgeline plot

  pl1<-ggplot2::ggplot(Phen.data.tot, aes(y = forcats::fct_rev(phenotype), x = as.numeric(days_since_admission), fill = stat(x))) +
    ggridges::geom_density_ridges_gradient(scale = 2, rel_min_height = 0.01, alpha = 0.8) +
    ggplot2::theme_minimal()+
    ggplot2::scale_fill_viridis_c(name = "Days since admission:",option = "G")+
    ggplot2::labs(x = "Days since admission", y = "", title = "Phenotype distribution")

  # Bar plot

  Phen.data.tot<-dplyr::mutate(Phen.data.tot, stage = factor(stage, levels = c("day90plus", "day30to89", "day0to29", "dayN14toN1", "before_adm")))
  if (is.null(stratified_by)) {
    pl2<-ggplot2::ggplot(Phen.data.tot, aes(forcats::fct_rev(phenotype), fill = stage))+
      ggplot2::geom_bar()+
      ggplot2::theme_minimal()+
      ggplot2::coord_flip()+
      ggplot2::scale_fill_brewer(name = "Stage:", labels = c("day90plus", "day30to89", "day0to29", "dayN14toN1", "before_adm"))+
      ggplot2::labs(y = "Counts", x = "", title = "Phenotype occurences", )
  } else {
    pl2<-ggplot2::ggplot(Phen.data.tot, aes(y=Phen.data.tot[,stratified_by], fill = stage, col=sex))+
        ggplot2::geom_bar()+
        ggplot2::facet_grid(phenotype~ ., switch = "y")+
        ggplot2::theme_minimal()+
        ggplot2::scale_y_discrete(position = "right")+
        ggplot2::scale_fill_brewer(name = "Stage:", labels = c("day90plus", "day30to89", "day0to29", "dayN14toN1", "before_adm"))+
        ggplot2::theme(strip.text.y.left = ggplot2::element_text(angle = 0))+
        ggplot2::labs(y = "", x = "Counts", title = paste("Phenotype occurences stratified by",stratified_by))

  }

  return(list(
    tot.count = tot.count,
    phen.distributionPlot = pl1,
    phen.barPlot = pl2)
  )


}
