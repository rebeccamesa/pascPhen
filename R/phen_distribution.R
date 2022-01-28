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
    dplyr::select(patient_num, sex, admission_date) %>%
    merge(PatientObservations, by = "patient_num")
  data$calendar_date <- as.Date(data$admission_date)+data$days_since_admission
  data <- data %>% 
    dplyr::mutate(stage = dplyr::case_when(days_since_admission > 90 ~ "long",
                             days_since_admission > 30 ~ "post",
                             days_since_admission > -14 ~ "acute",
                             TRUE ~ 'pre'))%>%
    dplyr::select(c(1,2,3,4,7,8,5,6)) %>%
    dplyr::mutate(stage = factor(stage, levels = c("pre", "acute", "post", "long")))
  
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
  names(tot.count)<- c("Phenotype", "pre", "acute", "post", "long")
  
  # Ridgeline plot
  
  pl1<-ggplot2::ggplot(Phen.data.tot, aes(y = forcats::fct_rev(forcats::fct_infreq(phenotype)), x = as.numeric(days_since_admission), fill = phenotype)) +
    ggridges::geom_density_ridges(scale = 3, rel_min_height = 0.01, alpha = 0.9) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none") +
    ggplot2::labs(x = "Days since admission", y = "", title = "Phenotype distribution")
  
  # Bar plot
  
  Phen.data.tot<-dplyr::mutate(Phen.data.tot, stage = factor(stage, levels = c("long", "post", "acute", "pre")))
  if (is.null(stratified_by)) {
    pl2<-ggplot2::ggplot(Phen.data.tot, aes(forcats::fct_rev(forcats::fct_infreq(phenotype)), fill = stage))+
      ggplot2::geom_bar()+
      ggplot2::theme_minimal()+
      ggplot2::coord_flip()+
      ggplot2::scale_fill_brewer(name = "Stage:", labels = c("Long (90, )", "Post (30, 90]", "Acute (-14, 30]", "Pre ( ,-14]"))+
      ggplot2::labs(y = "Counts", x = "", title = "Phenotype occurences")
  } else {
    pl2<-ggplot2::ggplot(Phen.data.tot, aes(y=Phen.data.tot[,stratified_by], fill = stage))+
        ggplot2::geom_bar()+
        ggplot2::facet_grid(forcats::fct_infreq(phenotype)~ ., switch = "y")+
        ggplot2::theme_minimal()+
        ggplot2::scale_y_discrete(position = "right")+
        ggplot2::scale_fill_brewer(name = "Stage:", labels = c("Long (90, )", "Post (30, 90]", "Acute (-14, 30]", "Pre ( ,-14]"))+
        ggplot2::theme(strip.text.y.left = ggplot2::element_text(angle = 0))+
        ggplot2::labs(y = "", x = "Counts", title = c("Phenotype occurences stratified by",stratified_by))
      
  }
  
  return(list(
    phen.distributionPlot = pl1,
    phen.barPlot = pl2)
  )
  
  
}