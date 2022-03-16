#' RUN PASC analysis lite
#'
#' @param data_dir directory of the 4CE tables
#' @param output_dir directory of the output report. If NULL, the result is saved in 'getProjectOutput' directory
#' @param siteid label specifying the site
#' @param long.thres "long" stage threshold
#' @param MSMR.sparsity MSMR lite parameter
#' @param long.perc percentage criterion to choose which phenotypes are interesting
#'
#' @return NULL
#' @export

# # Set input data directory
# data_dir <- "C:/Users/rebim/Desktop/Tesi/PASC Analysis/Dati ICSM 4CE"
# # Set output results directory. If NULL results are saved in 'getProjectOutput' directory
# output_dir <- "C:/Users/rebim/Desktop/Tesi/PASC Analysis/Rule based Phenotyping/Output"
#
# siteid = "icsm"
#
# long.thres = 15
# MSMR.sparsity = 0.005
# long.perc <- 0.04

RUN_lite <- function(data_dir, output_dir, siteid, long.thres, long.perc = 0.03, MSMR.sparsity, data_type){

  if (is.null(output_dir)) {
    output_dir <- getProjectOutputDirectory()
  }

  distribution.output <- runAnalysis_distribution (data_dir, siteid, long.thres, data_type)
  tot.count<-distribution.output$tot.count
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
  phen.prev <- subset(distribution.output$post.prev, perc >=long.perc*100)
  phen.prev <- phen.prev[order(phen.prev$perc, decreasing = TRUE),]

  PASClist <- phen.prev$phenotype
  phen.data.plot <- distribution.output$phen.data%>%
    dplyr::filter(phenotype %in% PASClist)

  pl2<-ggplot2::ggplot(phen.data.plot, aes(y = forcats::fct_rev(phenotype), x = as.numeric(days_since_admission), fill = stat(x)))+
    ggridges::geom_density_ridges_gradient(scale = 2, rel_min_height = 0.01, alpha = 0.8) +
    ggplot2::theme_minimal()+
    ggplot2::scale_fill_viridis_c(name = "Days since admission:",option = "G")+
    ggplot2::labs(x = "Days since admission", y = "",
                  title = "Phenotype distribution",
                  subtitle = siteid)

  mlho_list<- list()
  corr_list <- list()

  # create for loop to produce the tables of new features and report rocs
  for (i in seq_along(PASClist)) {
    tryCatch({
      print(paste0("running MLHO for ",PASClist[i]))
      # run MLHO on each phenotype
      mlho_output <- runAnalysis_MLHO(data_dir, siteid, long.thres, phenotype = PASClist[i],  MSMR.sparsity, data_type)
      mlho_list[[i]] <- mlho_output$output
      corr_list[[i]] <- mlho_output$corr
    },
    error = function(fr) {cat("ERROR :",conditionMessage(fr), "\n")})
  }

  features.table <- do.call(rbind, lapply(mlho_list, data.frame, stringsAsFactors=FALSE))
  write.csv(features.table,
            file = file.path(output_dir, paste0("table_",siteid,"_",long.thres,"_",substr(long.perc,3,nchar(long.perc)),"_",gsub("\\.", "",data_type),"_", Sys.Date(),".csv")), row.names = FALSE)
  # features.table$OR <- round(features.table$OR,3)
  # features.table$low <- round(features.table$low,3)
  # features.table$high <- round(features.table$high,3)
  cor.tab <- list()
  for (j in seq_along(corr_list)) {
    dat.i<- data.frame(corr_list[j])
    dat.i$features <- rownames(dat.i)
    rownames(dat.i) <- NULL

    dat.i.wide <- dat.i %>%
      reshape2::melt(id.var="features") %>%
      dplyr::arrange(features, variable)
    colnames(dat.i.wide) <- c("var1","var2","correlation")
    dat.i.wide$correlation <- round(dat.i.wide$correlation,3)
    cor.tab[[j]] <- dat.i.wide
  }
  cor.tab <- do.call(rbind, lapply(cor.tab, data.frame, stringsAsFactors=FALSE))

  tmp.dir<-paste0(output_dir,"/tmp")
  if (dir.exists(tmp.dir)) {
    unlink(tmp.dir, recursive = TRUE)
  }
  dir.create(tmp.dir)
  tmp.file <- paste0(tmp.dir,"/MLHOAnalysis.RData")
  save(pl1,phen.prev,pl2,features.table,cor.tab,file = tmp.file)

  #rmarkdown::render(paste0(output_dir,"/Report.Rmd"),output_file = paste0("report.",siteid,".", Sys.Date(),".html"))
  rmarkdown::render(system.file("rmd", "Report_lite.Rmd", package = "pascPhen"),
                    output_file = paste0("report_",siteid,"_",long.thres,"_",substr(long.perc,3,nchar(long.perc)),"_",gsub("\\.", "",data_type),"_", Sys.Date(),".html"),
                    output_dir = output_dir)
}
