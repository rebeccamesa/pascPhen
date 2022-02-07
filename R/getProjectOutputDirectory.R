
#' Returns the name of the project's output directory. runAnalysis() should save its output here, and submitAnalysis() will expect to read from this location. This function checks to make sure the location exists, and if not, creates it.
#'
#' @keywords 4CE
#' @export

getProjectOutputDirectory <- function ()
{
  dataRepositoryUrl = getPublicSummaryRepositoryUrl()
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
  return(dirName)
}
