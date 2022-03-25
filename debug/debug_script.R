# # Set parameters
# data_dir <- ""
# output_dir <- ""
# siteid <- ""
# long.thres <- 60 #60,90
# long.perc <- 0.03
# MSMR.sparsity <- 0.005
# data_type <- "" #2.1, 2.2, 2.2.all

# Load packages
if(!require(devtools)) devtools::install_github("hadley/devtools")
if(!require(mlho)) devtools::install_github("hestiri/mlho")
if(!require(pascPhen)) devtools::install_github("rebeccamesa/pascPhen")
if(!require(pacman)) install.packages("pacman")

pacman::p_load(data.table, devtools, backports, Hmisc, tidyr,dplyr,ggplot2,plyr,scales,readr,Rmisc,
               httr, DT, lubridate, tidyverse,reshape2,foreach,doParallel,caret,gbm,lubridate,praznik,
               ggridges, forcats, stats, FourCePhase2.1Data, ppcor,pascPhen, reshape2)

if (is.null(output_dir)) {
  dataRepositoryUrl = "https://github.com/covidclinical/pascPhenSummariesPublic.git"
  repositoryName = gsub(x = gsub(x = dataRepositoryUrl, pattern = "https://github.com/covidclinical/",
                                 fixed = TRUE, replace = ""), pattern = ".git", fixed = TRUE,
                        replace = "")
  projectName = gsub(x = repositoryName, pattern = "SummariesPublic",
                     replace = "", fixed = TRUE)
  dirName = file.path(FourCePhase2.1Data::getContainerScratchDirectory(),
                      projectName)
  if (!dir.exists(dirName)) {
    dir.create(dirName)
  }
  output_dir <- dirName
}

tmp.dir<-paste0(output_dir,"/",stringi::stri_rand_strings(1, 8))
if (dir.exists(tmp.dir)) {
  unlink(tmp.dir, recursive = TRUE)
}
dir.create(tmp.dir)
# tmp.file <- paste0(tmp.dir,"/MLHOAnalysis.RData")


# Run analysis distribution
PatientObservations <- utils::read.csv(file.path(data_dir, "LocalPatientObservations.csv"))
PatientSummary <- utils::read.csv(file.path(data_dir, "LocalPatientSummary.csv"))
if (data_type != "2.1") {
  if (data_type == "2.2") {
    PatientSummary <- PatientSummary%>%dplyr::filter(cohort %like% "PosAdm")
  } else{
    if (data_type == "2.2.all") {
      PatientSummary<-PatientSummary%>%dplyr::filter(cohort %like% "PosAdm" | cohort %like% "PosNotAdm")
    }
  }
  pts <- unique(PatientSummary$patient_num)
  PatientObservations <- PatientObservations%>%dplyr::filter(patient_num %in% pts)
}
utils::data("rules", package = "pascPhen")


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

# Prepare data
PatientObservations <- dplyr::select(PatientObservations, -value,-siteid)
PatientObservations$concept_code <- gsub('[.]', '', PatientObservations$concept_code)

data <- PatientSummary %>%
  #dplyr::select(patient_num, severe, deceased, sex, age_group, race, admission_date) %>%
  dplyr::select(patient_num, severe, sex, age_group, admission_date)%>%
  merge(PatientObservations, by = "patient_num")
data$calendar_date <- as.Date(data$admission_date)+data$days_since_admission
data <- data %>%
  dplyr::mutate(stage = dplyr::case_when(days_since_admission >= 90 ~ "day90plus",
                                         days_since_admission >= 60 ~ "day60to89",
                                         days_since_admission >= 30 ~ "day30to59",
                                         days_since_admission >= 0 ~ "day0to29",
                                         days_since_admission >=-14 ~ "dayN14toN1",
                                         TRUE ~ "before_adm"))%>%
  dplyr::mutate(stage = factor(stage, levels = c("before_adm", "dayN14toN1", "day0to29", "day30to59", "day60to89", "day90plus")))

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

# Use only the first time the phenotype appears
dOrder <- Phen.data.tot[order(Phen.data.tot$patient_num, Phen.data.tot$days_since_admission),]
Phen.data.tot <- dplyr::distinct(dOrder, patient_num, phenotype, .keep_all = TRUE)
tot.count <- as.data.frame.matrix(table(Phen.data.tot$phenotype, Phen.data.tot$stage))
tot.count$Phenotype<-rownames(tot.count)
tot.count<-tot.count[,c(7,1,2,3,4,5,6)]
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
n.tot <- rep(n.post, nrow(post.count))
site <- rep(siteid, nrow(post.count))
post.count <- cbind(post.count, n.tot, site)


tot.count.blur <- blur_it(tot.count,c("before_adm","dayN14toN1","day0to29","day30to59", "day60to89", "day90plus"), blur_abs, mask_thres)
post.count.blur <- blur_it(post.count, "n", blur_abs, mask_thres)
post.prev.blur <- post.count.blur%>%dplyr::mutate(perc = 100*n/n.tot)

phen.data <- Phen.data.tot
tot.count <- tot.count.blur
post.prev <- post.prev.blur

tot.count.long <- pivot_longer(tot.count[-nrow(tot.count),], -c(site, Phenotype), values_to = "Count", names_to = "Stage")
tot.count.plot <- tot.count.long %>%mutate(Stage = factor(Stage, levels = c("day90plus","day60to89","day30to59","day0to29","dayN14toN1", "before_adm")))

pl1<-ggplot2::ggplot(tot.count.plot,
                     aes(x = forcats::fct_rev(Phenotype), y = Count, fill = Stage))+
  ggplot2::geom_bar(stat = "identity")+
  ggplot2::theme_minimal()+
  ggplot2::coord_flip()+
  ggplot2::scale_fill_brewer(name = "Stage:", labels = c("day90plus","day60to89","day30to59", "day0to29", "dayN14toN1", "before_adm"))+
  ggplot2::labs(y = "Counts", x = "",
                title = "Phenotype counts",
                subtitle = siteid)



phen.prev <- subset(post.prev, perc >=long.perc*100)
phen.prev <- phen.prev[order(phen.prev$perc, decreasing = TRUE),]

PASClist <- phen.prev$phenotype
phen.data.plot <- phen.data%>%
  dplyr::filter(phenotype %in% PASClist)

pl2<-ggplot2::ggplot(phen.data.plot, aes(y = forcats::fct_rev(phenotype), x = as.numeric(days_since_admission), fill = stat(x)))+
  ggridges::geom_density_ridges_gradient(scale = 2, rel_min_height = 0.01, alpha = 0.8) +
  ggplot2::theme_minimal()+
  ggplot2::scale_fill_viridis_c(name = "Days since admission:",option = "G")+
  ggplot2::labs(x = "Days since admission", y = "",
                title = "Phenotype distribution",
                subtitle = siteid)
save(pl1,phen.prev,pl2,PASClist,file = paste0(tmp.dir,"/distributions.RData"))
rm(pl1,phen.prev,pl2,phen.data.plot,
   tot.count,tot.count.long,tot.count.plot,site.country.obfuscation);gc()

# create for loop to produce the tables of new features and report rocs
for (i in seq_along(PASClist)) {
  tryCatch({
    print(paste0("running MLHO for ",PASClist[i]))
    # run MLHO on each phenotype
    phenotype <- PASClist[i]
    PatientObservations <- utils::read.csv(file.path(data_dir, "LocalPatientObservations.csv"))
    PatientSummary <- utils::read.csv(file.path(data_dir, "LocalPatientSummary.csv"))
    if (data_type != "2.1") {
      if (data_type == "2.2") {
        PatientSummary <- PatientSummary%>%dplyr::filter(cohort %like% "PosAdm")
      } else{
        if (data_type == "2.2.all") {
          PatientSummary<-PatientSummary%>%dplyr::filter(cohort %like% "PosAdm" | cohort %like% "PosNotAdm")
        }
      }
      pts <- unique(PatientSummary$patient_num)
      PatientObservations <- PatientObservations%>%dplyr::filter(patient_num %in% pts)
    }
    utils::data("rules", package = "pascPhen")

    PatientObservations$concept_code <- gsub('[.]', '', PatientObservations$concept_code)
    days.long <- long.thres
    days.pre <- -14
    phenotype <- phenotype
    prediction <-FALSE
    test.sample<-25
    MSMR.topn<-100
    multicore<-FALSE
    iteration<-5
    score.threshold<-80
    classifier<-"GLM"
    print(paste0("Phenotyping .....",phenotype, "  (─‿─)  "))

    rules.phen <- dplyr::filter(rules, PASC_Phenotype == phenotype)

    PatientObservations <- dplyr::mutate(PatientObservations, patient_num=as.character(patient_num))

    # Let's binarize age. "other" has to be discussed
    PatientSummary <- dplyr::mutate(PatientSummary, age_group = dplyr::case_when((age_group == "80plus" | age_group == "70to79") ~ "70plus",
                                                                                 TRUE ~ "00to69"))
    # Create labels
    rulesPhen <- rules.phen
    dOrder <- PatientObservations[order(PatientObservations$patient_num, PatientObservations$days_since_admission),]
    PatientObservations <- dplyr::distinct(dOrder, patient_num, concept_code, value, .keep_all = TRUE)
    PatientObservations$label <- 0

    PatientObservations$label[PatientObservations$days_since_admission>days.long & paste(gsub(".*-", "", PatientObservations$concept_type), PatientObservations$concept_code) %in% paste(rulesPhen$Coding, rulesPhen$Code)] <- 1
    labeldt <- PatientObservations%>%
      dplyr::select(patient_num,label)%>%
      dplyr::group_by(patient_num)%>%
      dplyr::summarize_all(max)%>%
      dplyr::mutate(label = as.factor(label))


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

    phenotyping.features <- model.output
    phenotyping.features$type <- rep("phenotyping", nrow(phenotyping.features))
    phenotyping.features$site <- rep(siteid, nrow(phenotyping.features))
    phenotyping.features <- as.data.frame(phenotyping.features)

    features <- phenotyping.features
    corrs <- corr.matrix
    save(features,corrs,file = paste0(tmp.dir,"/phen",i,".RData"))
    rm(features,corrs);gc()

  },
  error = function(fr) {cat("ERROR :",conditionMessage(fr), "\n")})
}

rmarkdown::render(system.file("rmd", "Report_lite.Rmd", package = "pascPhen"),
                  output_file = paste0("report_",siteid,"_",long.thres,"_",substr(long.perc,3,nchar(long.perc)),"_",gsub("\\.", "",data_type),"_", Sys.Date(),".html"),
                  output_dir = output_dir)
