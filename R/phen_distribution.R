#' Phenotype distribution
#'
#'
#' @param PatientObservations 2.1 PatientObservations table
#' @param PatientSummary 2.1 PatientSummary table
#' @param rules table containing the phenotype definitions
#' @param siteid label specifying the site
#' @param long.thres "long" stage threshold
#'
#' @return count tables
#' @export

phen_distribution <- function(PatientObservations,
                              PatientSummary,
                              rules,
                              siteid,
                              long.thres)
{
  # Prepare data
  PatientObservations <- dplyr::select(PatientObservations, -value,-siteid)
  PatientObservations$concept_code <- gsub('[.]', '', PatientObservations$concept_code)

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
  #tot.count <- data.frame()
  Phen.data.tot <- data.frame()
  for (i in 1:length(phenotypes)) {
    Phen.rules <- dplyr::filter(rules, PASC_Phenotype == phenotypes[i])
    Phen.data <- data[paste(gsub(".*-", "", data$concept_type), data$concept_code) %in% paste(Phen.rules$Coding, Phen.rules$Code),]
    counts <- table(dplyr::select(Phen.data,stage))
    if (sum(counts)>0) {
      #tot.count<-rbind(tot.count,c(phenotypes[i], counts))
      Phen.data$phenotype <- paste(phenotypes[i])
      Phen.data.tot<-rbind(Phen.data.tot, Phen.data)
    }
  }
  # names(tot.count)<- c("Phenotype", "before_adm", "dayN14toN1", "day0to29", "day30to89", "day90plus")
  # tot.count<-tot.count%>%dplyr::mutate(before_adm = as.numeric(before_adm), dayN14toN1 = as.numeric(dayN14toN1),
  #                                      day0to29 = as.numeric(day0to29), day30to89 = as.numeric(day30to89), day90plus = as.numeric(day90plus))
  # tot.count<-tot.count[order(tot.count$Phenotype),]
  # dist<-data%>%dplyr::distinct(patient_num,stage)
  # tot<-table(dist$stage)
  # tot$Phenotype <- "TOT"
  # tot.count<-rbind(tot.count, tot)
  # site <- rep(siteid, nrow(tot.count))
  # tot.count <- cbind(tot.count,site)
  #Phen.data.tot$siteid <- rep(siteid,nrow(Phen.data.tot))

  # tot.count.strat <- NULL
  # if (!is.null(stratified_by)) {
  #   tot.count.strat <- table(Phen.data.tot$phenotype, Phen.data.tot$stage, Phen.data.tot[,stratified_by], Phen.data.tot$siteid)
  # }

  # Use only the first time the phenotype appears
  dOrder <- Phen.data.tot[order(Phen.data.tot$patient_num, Phen.data.tot$days_since_admission),]
  Phen.data.tot <- dplyr::distinct(dOrder, patient_num, phenotype, .keep_all = TRUE)
  tot.count <- as.data.frame.matrix(table(Phen.data.tot$phenotype, Phen.data.tot$stage))
  tot.count$Phenotype<-rownames(tot.count)
  tot.count<-tot.count[,c(6,1,2,3,4,5)]
  rownames(tot.count) <- 1:dim(tot.count)[1]
  dist<-data%>%dplyr::distinct(patient_num,stage)
  tot<-table(dist$stage)
  tot$Phenotype <- "TOT"
  tot.count <- rbind(tot.count,tot)
  site <- rep(siteid, nrow(tot.count))
  tot.count <- cbind(tot.count, site)

  post.covid <- data%>%dplyr::filter(days_since_admission >= long.thres)
  n.post <- length(unique(post.covid$patient_num))
  post.count <-dplyr::count(dplyr::distinct(Phen.data.tot%>%dplyr::filter(days_since_admission >= long.thres),patient_num,phenotype),phenotype)
  # post.prev <- post.count%>%dplyr::mutate(perc = n/n.post)
  # post.prev <- post.prev[order(post.prev$perc, decreasing = TRUE),]
  n.tot <- rep(n.post, nrow(post.count))
  site <- rep(siteid, nrow(post.count))
  post.count <- cbind(post.count, n.tot, site)

  return(list(
    Phen.data.tot = Phen.data.tot,
    tot.count = tot.count,
    post.count = post.count)
  )


}


