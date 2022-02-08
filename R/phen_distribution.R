#' Phenotype distribution
#'
#' Build ridgeline plots and box plots of the distribution of the phenotypes
#' over the dataset
#'
#' @param PatientObservations 2.1 PatientObservations table
#' @param PatientSummary 2.1 PatientSummary table
#' @param rules table containing the phenotype definitions
#' @param siteid label specifying the site
#' @param stratified_by label for the stratification of the bar plot
#'
#' @return four plots and count table
#' @export

phen_distribution <- function(PatientObservations,
                              PatientSummary,
                              rules,
                              siteid,
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
  tot.count<-tot.count%>%dplyr::mutate(before_adm = as.numeric(before_adm), dayN14toN1 = as.numeric(dayN14toN1),
                                       day0to29 = as.numeric(day0to29), day30to89 = as.numeric(day30to89), day90plus = as.numeric(day90plus))
  tot.count<-tot.count[order(tot.count$Phenotype),]
  dist<-data%>%dplyr::distinct(patient_num,stage)
  tot<-table(dist$stage)
  tot$Phenotype <- "TOT"
  tot.count<-rbind(tot.count, tot)
  site <- rep(siteid, nrow(tot.count))
  tot.count <- cbind(tot.count,site)

  ridgeline.plot <- function(data,siteid,first){
    # Create ridgeline plot with the distribution of the phenotypes

    pl<-ggplot2::ggplot(data, aes(y = forcats::fct_rev(phenotype), x = as.numeric(days_since_admission), fill = stat(x))) +
      ggridges::geom_density_ridges_gradient(scale = 2, rel_min_height = 0.01, alpha = 0.8) +
      ggplot2::theme_minimal()+
      ggplot2::scale_fill_viridis_c(name = "Days since admission:",option = "G")+
      ggplot2::labs(x = "Days since admission", y = "",
                    title = ifelse(first,"Phenotype distribution - first occurrence","Phenotype distribution"),
                    subtitle = siteid)
    pl
  }

  bar.plot <- function(data, stratified_by,siteid,first){
    # Create bar plot
    data<-dplyr::mutate(data, stage = factor(stage, levels = c("day90plus", "day30to89", "day0to29", "dayN14toN1", "before_adm")))
    if (is.null(stratified_by)) {
      pl<-ggplot2::ggplot(data, aes(forcats::fct_rev(phenotype), fill = stage))+
        ggplot2::geom_bar()+
        ggplot2::theme_minimal()+
        ggplot2::coord_flip()+
        ggplot2::scale_fill_brewer(name = "Stage:", labels = c("day90plus", "day30to89", "day0to29", "dayN14toN1", "before_adm"))+
        ggplot2::labs(y = "Counts", x = "",
                      title = ifelse(first,"Phenotype counts - first occurrence","Phenotype counts"),
                      subtitle = siteid)
    } else {
      pl<-ggplot(data) +
        geom_bar(aes(y = as.factor(data[,stratified_by]), fill = as.factor(data[,stratified_by]), alpha = stage)) +
        facet_grid(phenotype~., switch = "y") +
        ggplot2::theme_minimal()+
        ggplot2::scale_y_discrete(position = "right")+
        #scale_fill_manual(values = c("#ff7514","#008f39"), guide = "none")+
        scale_fill_manual(values = RColorBrewer::brewer.pal(name="Set1", n=length(unique(data[,stratified_by]))), guide="none")+
        scale_alpha_manual(values=c(0.2, 0.4, 0.6, 0.8, 1), name = "Stage", labels= c("day90plus", "day30to89", "day0to29", "dayN14toN1", "before_adm"))+
        ggplot2::theme(strip.text.y.left = ggplot2::element_text(angle = 0))+
        ggplot2::labs(y = "", x = "Counts",
                      title = ifelse(first, paste("Phenotype counts stratified by", "\'",stratified_by, "\' - first occurrence"),
                                     paste("Phenotype counts stratified by", "\'",stratified_by, "\'")),
                      subtitle = siteid)

    }
    pl
  }

  # Use all the occurences -> Phen.data.tot
  pl1<-ridgeline.plot(Phen.data.tot,siteid,first = FALSE)
  pl2<-bar.plot(Phen.data.tot,stratified_by,siteid, first = FALSE)

  Phen.data.tot$siteid <- rep(siteid,nrow(Phen.data.tot))
  tot.count.strat <- NULL
  if (!is.null(stratified_by)) {
    tot.count.strat <- table(Phen.data.tot$phenotype, Phen.data.tot$stage, Phen.data.tot[,stratified_by], Phen.data.tot$siteid)
  }

  # Use only the first time the phenotype appears -> Phen.data.first
  dOrder <- Phen.data.tot[order(Phen.data.tot$patient_num, Phen.data.tot$days_since_admission),]
  Phen.data.first <- dplyr::distinct(dOrder, patient_num, phenotype, .keep_all = TRUE)
  first.count <- as.data.frame.matrix(table(Phen.data.first$phenotype, Phen.data.first$stage))
  first.count$Phenotype<-rownames(first.count)
  first.count<-first.count[,c(6,1,2,3,4,5)]
  rownames(first.count) <- 1:dim(first.count)[1]
  first.count <- rbind(first.count,tot)
  first.count <- cbind(first.count,site)

  pl3<-ridgeline.plot(Phen.data.first,siteid,first = TRUE)
  pl4<-bar.plot(Phen.data.first,stratified_by,siteid, first = TRUE)

  Phen.data.first$siteid <- rep(siteid, nrow(Phen.data.first))
  first.count.strat<-NULL
  if (!is.null(stratified_by)) {
    first.count.strat<-table(Phen.data.first$phenotype, Phen.data.first$stage, Phen.data.first[,stratified_by], Phen.data.first$siteid)
  }

  return(list(
    tot.count = tot.count,
    tot.count.strat = tot.count.strat,
    first.count = first.count,
    first.count.strat = first.count.strat,
    phen.distributionPlot = pl1,
    phen.barPlot = pl2,
    phen.1stDistributionPlot = pl3,
    phen.1stBarPlot = pl4)
  )


}


