#' analyzing buckets of MLHO and core features
#'
#' @param data_dir directory of the 4CE tables
#' @param output_dir directory for saving the results
#' @param siteid label specifying the site
#' @param long.thres long" stage threshold
#' @param data_type string with the type of 4CE data
#' @param pascphens list of pasc phenotype to study
#' @param chart_review set TRUE if you want data for chart review
#'
#' @return NULL
#' @export

pascbuckets <- function(data_dir,
                        output_dir,
                        siteid,
                        long.thres,
                        data_type,
                        pascphens,
                        chart_review=FALSE)
{

  if (is.null(output_dir)) {
    output_dir <- getProjectOutputDirectory()
  }

  tmp.dir<-paste0(output_dir,"/",siteid,"_iterationPASC2_",stringi::stri_rand_strings(1, 4))
  if (dir.exists(tmp.dir)) {
    unlink(tmp.dir, recursive = TRUE)
  }
  dir.create(tmp.dir)


  PatientObservations <- utils::read.csv(file.path(data_dir, "LocalPatientObservations.csv"))
  PatientSummary <- utils::read.csv(file.path(data_dir, "LocalPatientSummary.csv"))
  if (data_type != "2.1") {
    data <- filter_data(PatientObservations, PatientSummary, data_type)
    PatientObservations <- data$PatientObservations
    PatientSummary <- data$PatientSummary
  }
  utils::data("rules_it2")
  rules_it2 <- subset(rules_it2,rules_it2$PASC_Phenotype %in% pascphens)

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

  PatientObservations <- dplyr::select(PatientObservations, -value,-siteid)
  PatientObservations$concept_code <- gsub('[.]', '', PatientObservations$concept_code)
  PatientObservations$concept_code <- gsub('[-]', '', PatientObservations$concept_code)

  #keep a copy for a summary
  # PatientObservations_copy <- PatientObservations
  ##lightening up
  PatientObservations <- subset(PatientObservations,PatientObservations$concept_code %in% rules_it2$Code)

  data <- PatientSummary %>%
    dplyr::select(patient_num, severe, sex, age_group, admission_date)%>%
    merge(PatientObservations, by = "patient_num")
  data$calendar_date <- as.Date(data$admission_date)+data$days_since_admission

  ###going phenotype-by-phenotype
  pasc_data <- list()
  bucketlist <- list()
  phen.feats <- list()
  bucketrows <- list()

  for (h in seq_along(pascphens)){
    tryCatch({

      feats.h <- subset(rules_it2,rules_it2$PASC_Phenotype == pascphens[h])
      dat.h.core <- subset(data,data$concept_code %in% c(as.character(subset(feats.h$Code,feats.h$feature_type == "core"))))
      sel1 <- dat.h.core %>%
        dplyr::group_by(patient_num,concept_code) %>%
        dplyr::summarise(min_days_since_admission=min(days_since_admission)) %>%
        filter(min_days_since_admission > long.thres) %>%
        mutate(key = paste0(patient_num,concept_code,min_days_since_admission) )
      dat.h.core <- dat.h.core %>%
        mutate(key = paste0(patient_num,concept_code,days_since_admission))
      dat.h.core <- subset(dat.h.core,dat.h.core$key %in% sel1$key)
      dat.h.core$feature_type <- "core"
      dat.h.core$key <- NULL
      ##data from mlho
      dat.h.mlho <- subset(data,data$concept_code %in% c(as.character(subset(feats.h$Code,feats.h$feature_type == "MLHO"))))
      dat.h.mlho$feature_type <- "MLHO"
      feats.h.mlho <- c(as.character(unique(subset(feats.h$Code,feats.h$feature_type == "MLHO"))))

      dat.h.mlho_update <- list()
      for (f in seq_along(feats.h.mlho)){
        tryCatch({
          feats.h.mlho.f <- subset(feats.h,feats.h$Code == feats.h.mlho[f])
          dat.h.mlho.f <- dat.h.mlho %>%
            filter(concept_code == feats.h.mlho[f] &
                     days_since_admission > as.integer(feats.h.mlho.f$time))

          if(feats.h.mlho.f$end_time != Inf){
            dat.h.mlho.f <- dat.h.mlho.f %>%
              filter(days_since_admission <= as.integer(feats.h.mlho.f$time))

          } else if (feats.h.mlho.f$end_time == Inf){
            dat.h.mlho.f <- dat.h.mlho.f %>%
              filter(days_since_admission <= feats.h.mlho.f$end_time)

          }

          dat.h.mlho_update[[f]] <- dat.h.mlho.f
          rm(dat.h.mlho.f,feats.h.mlho.f)
        },
        error = function(fr) {cat("ERROR :",conditionMessage(fr), "\n")})

      }
      dat.h.mlho <- do.call(rbind, lapply(dat.h.mlho_update, data.frame, stringsAsFactors=FALSE))

      temp.dat <- rbind(dat.h.mlho,dat.h.core)
      temp.dat$PASC_Phenotype <- pascphens[h]
      temp.dat$approximate.time <- floor(temp.dat$days_since_admission/10)*10
      temp.dat$site <- siteid
      ##get an obfuscated version of the features and when they were observed
      features.phen.time <- dplyr::select(temp.dat,concept_type, concept_code,feature_type,
                                          approximate.time,PASC_Phenotype,site)
      features.phen.time <- features.phen.time %>%
        dplyr::group_by(concept_type,concept_code,feature_type,approximate.time,PASC_Phenotype,site) %>%
        dplyr::summarise(count=length(concept_code))

      ##buckets
      el3 <- temp.dat %>%
        dplyr::group_by(feature_type,patient_num) %>%
        dplyr::summarise(count=length(feature_type)) %>%
        pivot_wider(names_from = feature_type, values_from = count)

      bucket4 <- subset(el3$patient_num,!is.na(el3$core) & !is.na(el3$MLHO))
      bucket3 <- subset(el3$patient_num,!is.na(el3$core) & is.na(el3$MLHO))
      bucket2 <- subset(el3$patient_num,is.na(el3$core) & !is.na(el3$MLHO))

      el3$bucket <- ifelse(el3$patient_num %in% bucket4, 4,0)
      el3$bucket <- ifelse(el3$patient_num %in% bucket3, 3,el3$bucket )
      el3$bucket <- ifelse(el3$patient_num %in% bucket2, 2,el3$bucket )
      table(el3$bucket)

 if (chart_review==TRUE) {

      bucket1.sample <- subset(data$patient_num,!(data$patient_num %in% el3$patient_num))

      bucket4 <- data.frame(sample (bucket4, size=35, replace =F))
      bucket4$bucket <- 4
      colnames(bucket4)[1] <- "patient_num"
      bucket3 <- data.frame(sample (bucket3, size=35, replace =F))
      bucket3$bucket <- 3
      colnames(bucket3)[1] <- "patient_num"
      bucket2 <- data.frame(sample (bucket2, size=20, replace =F))
      bucket2$bucket <- 2
      colnames(bucket2)[1] <- "patient_num"
      bucket1 <- data.frame(sample (bucket1.sample, size=10, replace =F))
      bucket1$bucket <- 1
      colnames(bucket1)[1] <- "patient_num"

      buckets <- rbind(bucket4,bucket3, bucket2, bucket1)
      buckets$PASC_Phenotype <- pascphens[h]
      buckets$site <- siteid

      bucketlist[[h]] <- buckets #PHI needed for chart review
      # feat.viz[[h]] <- dat.viz
}

      el3$PASC_Phenotype <- pascphens[h]
      el3$site <- siteid


      phen.feats[[h]] <- features.phen.time ##aggregate

      bucketrows[[h]] <- el3 #patient level but no PHI
      pasc_data[[h]] <- temp.dat ##PHI data
      rm(temp.dat,dat.h.mlho,dat.h.core)
    },
    error = function(fra) {cat("ERROR :",conditionMessage(fra), "\n")})
  }
  ###also look at shared patients...
  ##patient level reports for the sites
  pasc_data <- do.call(rbind, lapply(pasc_data, data.frame, stringsAsFactors=FALSE))
  pasc_data$long.thres <- long.thres
  bucketlist <- do.call(rbind, lapply(bucketlist, data.frame, stringsAsFactors=FALSE))
  bucketlist$long.thres <- long.thres
  pasc_data$patient_num <- as.character(pasc_data$patient_num)
  bucketlist$patient_num <- as.character(bucketlist$patient_num)
  ##patient level but no PHI
  bucketrows <- do.call(rbind, lapply(bucketrows, data.frame, stringsAsFactors=FALSE))
  bucketrows$long.thres <- long.thres
  bucketrows$patient_num <- as.character(bucketrows$patient_num)

  # bucketrows$patient_num <- NULL
  ##aggregated data for meta analysis
  phen.feats <- do.call(rbind, lapply(phen.feats, data.frame, stringsAsFactors=FALSE))
  phen.feats$long.thres <- long.thres

  ###summary information
  summary <- data.frame(siteid)
  summary$total_patients <- length(unique(PatientObservations$patient_num))
  summary$atleast1_PASC_34 <-  nrow(bucketrows %>%
                                      filter(bucket %in% c(3,4)) %>%
                                      dplyr::group_by(patient_num) %>%
                                      dplyr::summarise(count=length(patient_num)) )
  summary$atleast1_PASC_perc_34 <- round(summary$atleast1_PASC_34/summary$total_patients,2) * 100

  summary$multi_PASC_34  <-  nrow(bucketrows %>%
                                    filter(bucket %in% c(3,4)) %>%
                                    dplyr::group_by(patient_num) %>%
                                    dplyr::summarise(count=length(patient_num)) %>%
                                    filter(count > 1))
  summary$multi_PASC_perc_34 <- round(summary$multi_PASC_34/summary$total_patients,2) * 100

  summary$atleast1_PASC_4 <-  nrow(bucketrows %>%
                                     filter(bucket %in% c(4)) %>%
                                     dplyr::group_by(patient_num) %>%
                                     dplyr::summarise(count=length(patient_num)) )
  summary$atleast1_PASC_perc_4 <- round(summary$atleast1_PASC_4/summary$total_patients,2) * 100

  summary$multi_PASC_4  <-  nrow(bucketrows %>%
                                   filter(bucket %in% c(4)) %>%
                                   dplyr::group_by(patient_num) %>%
                                   dplyr::summarise(count=length(patient_num)) %>%
                                   filter(count > 1))
  summary$multi_PASC_perc_4 <- round(summary$multi_PASC_4/summary$total_patients,2) * 100

  ##
  pasc_countss <- bucketrows %>%
    # filter(bucket %in% c(3,4)) %>%
    dplyr::group_by(bucket,PASC_Phenotype) %>%
    dplyr::summarise(count_distinct=n_distinct(patient_num))


  if(chart_review==TRUE){
    save(pasc_data,bucketlist,file = paste0(tmp.dir,"/PHI_for chart review_",siteid,long.thres,".RData"))
    save(bucketrows,file = paste0(tmp.dir,"/PHI_for internal analysis_",siteid,long.thres,".RData"))
  }

  save(phen.feats,summary,pasc_countss,file = paste0(tmp.dir,"/Share_for meta analysis_",siteid,long.thres,".RData"))


  ####obervations from buckets


  return(print(paste0("[][][][][][][][][][][][][][][][][][][][][][][][][][][][]",
                      "Your job is completed here --congratulations! results to share are stored in :",
                      tmp.dir))

  )
}
