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
#' @return table containing the extracted features with their ORs
#' @export

mlho_features <- function(PatientObservations,
                          PatientSummary,
                          rules,
                          days.long,
                          days.pre=-14,
                          phenotype,
                          prediction = FALSE,
                          test.sample=25,
                          MSMR.sparsity=0.005,
                          MSMR.topn=200,
                          multicore=FALSE,
                          iteration=5,
                          score.threshold=80,
                          classifier)
{
  
  # Prepare labeldt, dbmart and dems
  
  rules.phen <- dplyr::filter(rules, PASC_Phenotype == phenotype)
  
  PatientObservations <- dplyr::mutate(PatientObservations, patient_num=as.character(patient_num))
  PatientSummary <- dplyr::mutate(PatientSummary, age_group = dplyr::case_when((age_group == "80plus" | age_group == "70to79") ~ "70plus", 
                                                                                     TRUE ~ "12to69")) #let's binarize age
  labeldt <- createLabels(days.long, rules.phen, PatientObservations)
  
  dems <- PatientSummary %>%
    dplyr::select(patient_num, sex, age_group)%>%
    dplyr::mutate(age_group = as.factor(age_group), sex = as.factor(sex))
  
  PatientObservations.excl <- PatientObservations %>%
    dplyr::filter(!(paste(gsub(".*-", "", PatientObservations$concept_type), PatientObservations$concept_code) %in% paste(rules.phen$Coding, rules.phen$Code)))%>% # let's exclude codes used to define the phenotype
    dplyr::mutate(phenx = paste(concept_type,concept_code, sep = "_"))
  
  if (prediction) {
    dbmart <- dplyr::filter(PatientObservations.excl, days_since_admission < days.pre) #prediction
  } else {
    dbmart <- dplyr::filter(PatientObservations.excl, days_since_admission > days.long) #phenotyping
  }
  dbmart <- dplyr::mutate(dbmart[,c('patient_num', 'phenx')], DESCRIPTION = phenx)
  
  # Implement iterative MLHO
  
  data.table::setDT(dbmart)
  dbmart[,row :=.I]
  dbmart$value.var <- 1
  
  mlho.features <- mlho::mlho.it(dbmart,
                           labels = labeldt,
                           dems,
                           test.sample = test.sample,
                           MSMR.sparsity=MSMR.sparsity,
                           MSMR.topn=MSMR.topn,
                           mlearn.note="mlho phewas test run",
                           mlearn.aoi=phenotype,
                           multicore=multicore,
                           iterations=iteration)
  
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
  model.output <- mlho::mlearn(dat.train=model.data,
                       dat.test=NULL,
                       dems=dems,
                       save.model=FALSE,
                       classifier=classifier,
                       note="extracting_coefficients",
                       aoi=phenotype,
                       multicore=multicore)
  
  model.output <- stats::na.omit(model.output)
  
  return(model.output)
}

