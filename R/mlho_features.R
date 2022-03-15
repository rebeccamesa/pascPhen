#' Apply MLHO to extract features
#'
#' @importFrom data.table ":="
#'
#' @param PatientObservations 2.1 PatientObservations table
#' @param PatientSummary 2.1 PatientSummary table
#' @param rules table containing the phenotype definitions
#' @param days.long "long" stage threshold
#' @param days.pre "pre" stage threshold
#' @param phenotype phenotype label
#' @param prediction do you want to make prediction or phenotyping?
#' @param test.sample mlho.it parameter. Percentage for testing
#' @param MSMR.sparsity MSMR.lite parameter
#' @param MSMR.topn MSMR.lite parameter
#' @param multicore do you want to parallelize the process?
#' @param iteration mlho.it parameter. Number of iterations
#' @param score.threshold MLHO score threshold
#' @param classifier mlearn parameter. Classification algorithm
#'
#' @return output features and their prevalence + MLHO features' correlation
#' @export

mlho_features <- function(PatientObservations,
                          PatientSummary,
                          rules,
                          days.long,
                          days.pre=-14,
                          phenotype,
                          prediction = FALSE,
                          test.sample=25,
                          MSMR.sparsity,
                          MSMR.topn=200,
                          multicore=FALSE,
                          iteration=5,
                          score.threshold=80,
                          classifier="GLM")
{

  # Prepare labeldt, dbmart and dems

  rules.phen <- dplyr::filter(rules, PASC_Phenotype == phenotype)

  PatientObservations <- dplyr::mutate(PatientObservations, patient_num=as.character(patient_num))

  # Let's binarize age. "other" has to be discussed
  PatientSummary <- dplyr::mutate(PatientSummary, age_group = dplyr::case_when((age_group == "80plus" | age_group == "70to79") ~ "70plus",
                                                                               TRUE ~ "00to69"))
  labeldt <- createLabels(days.long, rules.phen, PatientObservations)

  dems <- PatientSummary %>%
    dplyr::select(patient_num, sex, age_group)%>%
    dplyr::mutate(age_group = as.factor(age_group), sex = as.factor(sex))

  # PatientObservations.excl <- PatientObservations %>%
  #   dplyr::filter(!(paste(gsub(".*-", "", PatientObservations$concept_type), PatientObservations$concept_code) %in% paste(rules.phen$Coding, rules.phen$Code)))%>% # let's exclude codes used to define the phenotype
  #   dplyr::mutate(phenx = paste(concept_type,concept_code, sep = "_"))
  PatientObservations.excl <- PatientObservations %>%
    dplyr::filter(!(paste(gsub(".*-", "", PatientObservations$concept_type), PatientObservations$concept_code) %in% paste(rules.phen$Coding, rules.phen$Code)))%>% # let's exclude codes used to define the phenotype
    dplyr::mutate(stage = dplyr::case_when(days_since_admission >= 90 ~ "day90plus",
                                           days_since_admission >= 60 ~ "day60to89",
                                           days_since_admission >= 30 ~ "day30to59",
                                           days_since_admission >= 0 ~ "day0to29",
                                           days_since_admission >=-14 ~ "dayN14toN1",
                                           TRUE ~ "before_adm"))%>%
    dplyr::mutate(phenx = paste0(concept_type,"_",concept_code,"_",stage))


  if (prediction) {
    dbmart <- dplyr::filter(PatientObservations.excl, days_since_admission < days.pre) #prediction
  } else {
    #dbmart <- dplyr::filter(PatientObservations.excl, days_since_admission > days.long) #phenotyping
    dbmart <- dplyr::filter(PatientObservations.excl, days_since_admission >=days.pre)
  }
  dbmart <- dplyr::mutate(dbmart[,c('patient_num', 'phenx')], DESCRIPTION = phenx)

  pts.phen <- labeldt%>%dplyr::filter(label == 1)
  pts.phen <- pts.phen$patient_num
  pts.nophen <- unique(dbmart%>%dplyr::filter(!(patient_num %in% pts.phen))%>%dplyr::select(patient_num))
  obs.phen <- dbmart%>%dplyr::filter(patient_num %in% pts.phen)%>%dplyr::mutate(aoi = 1, n.tot = length(pts.phen))
  obs.nophen <- dbmart%>%dplyr::filter(!(patient_num %in% pts.phen))%>%dplyr::mutate(aoi = 0, n.tot = nrow(pts.nophen))
  obs <- rbind(obs.phen,obs.nophen)

  print("Implement iterative MLHO for PheWAS")

  data.table::setDT(dbmart)
  dbmart[,row :=.I]
  dbmart$value.var <- 1

  mlho.features <- mlho::mlho.it(dbmart,
                                 labels = labeldt,
                                 dems,
                                 test.sample = test.sample,
                                 MSMR.sparsity=MSMR.sparsity,
                                 MSMR.topn=MSMR.topn,
                                 mlearn.note="mlho phewas",
                                 mlearn.aoi=phenotype,
                                 multicore=multicore,
                                 iterations=iteration)



  #
  require(Rmisc)
  mlho.auroc.mean <- exp(CI(log(mlho.features$model.i.roc), ci=0.95))[2]
  mlho.auroc.up95 <- exp(CI(log(mlho.features$model.i.roc), ci=0.95))[1]
  mlho.auroc.low95 <- exp(CI(log(mlho.features$model.i.roc), ci=0.95))[3]

    # Compute MLHO confidence scores (CS)

  mlho.scores <- mlho::mlho.cs(mlho.features)
  mlho.features <- c(as.character(subset(mlho.scores$features,mlho.scores$cs >= score.threshold)))


  # Implement mlearn with no test set with a training set made of the entire population of MSMR.lite

  dbmart <- subset(dbmart,dbmart$phenx %in% mlho.features)
  setDT(dbmart)
  dbmart[,row := .I]
  dbmart$value.var <- 1
  uniqpats <- c(as.character(unique(dbmart$patient_num)))

  model.data <- mlho::MSMSR.lite(MLHO.dat=dbmart,patients = uniqpats,sparsity=NA,jmi = FALSE,labels = labeldt)
  corr <- pcor(model.data[,2:(ncol(model.data)-1)], method = "pearson")
  corr.matrix <- corr$estimate
  print("computing the multivariate model")
  model.output <- mlho::mlearn(dat.train=model.data,
                               dat.test=NULL,
                               dems=dems,
                               save.model=FALSE,
                               classifier=classifier,
                               note="extracting_coefficients",
                               aoi=phenotype,
                               multicore=multicore)

  model.output <- stats::na.omit(model.output)
  model.output$auroc.mean <- round(as.numeric(mlho.auroc.mean)[1],3)
  model.output$auroc.CI <- paste0("[ ",round(as.numeric(mlho.auroc.low95)[1],3)," - ",round(as.numeric(mlho.auroc.up95)[1],3)," ]")


  obs.output <- obs %>%dplyr::filter(phenx %in% model.output$features)
  count.output <- dplyr::count(distinct(obs.output, patient_num, phenx, aoi, .keep_all = TRUE), phenx,aoi,n.tot)
  perc.output <- count.output%>%dplyr::mutate(perc = n/n.tot)

  return(list(
    mlho.corr = corr.matrix,
    model.output = model.output,
    perc.output = perc.output))
}
