HemaScribeData <- new.env(parent=emptyenv())

.onAttach <- function(libname, pkgname) {
  if (!nzchar(Sys.getenv("R_CMD_CHECK"))) {
    dataDir <- tools::R_user_dir("HemaScribe_v2", which = "data")

    files <- c("HemaScribe-data.rdx", "HemaScribe-data.rdb", "StromaScribe-data.rdx", "StromaScribe-data.rdb", "HemaScape-data.rda")

    for (file in files) {
      path <- file.path(dataDir, file)
      if (!file.exists(path)) {
        repo <- if (startsWith(file, "HemaScribe") | startsWith(file, "StromaScribe")) {
          "jhfung/HemaScribe_v2"
        } else {
          "RabadanLab/HemaScribe"
        }
        download_data(file = file, repo = repo, dataDir = dataDir)
      }
    }

    if (all(file.exists(file.path(dataDir, files)))) {
      lazyLoad(file.path(dataDir, "HemaScribe-data"), envir=HemaScribeData)
      lazyLoad(file.path(dataDir, "StromaScribe-data"), envir=HemaScribeData)
      load(file.path(dataDir, "HemaScape-data.rda"), envir=HemaScribeData)
    } else {
      rlang::abort("HemaScribe cannot be loaded because reference data is not found.")
    }

    packageStartupMessage("HemaScribe reference data loaded into HemaScribeData.")
  }
}

#' Download reference data
#'
#' Downloads the reference data from the Github repo.  The user is prompted for
#' download when the package is first loaded, but this function can also be
#' called in a later session to update the reference data file.
#' @param file Name of file to be downloaded
#' @param repo GitHub repository from which file is downloaded
#' @param dataDir Path where the data should be downloaded (default: package-managed user directory)
#' @export
download_data <- function(file, repo, dataDir = NA) {
  if (is.na(dataDir)) {
    dataDir <- tools::R_user_dir("HemaScribe_v2", which = "data")
  }

  if (!dir.exists(dataDir)) {
    dir.create(dataDir, recursive = TRUE)
  }

  piggyback::pb_download(
    file=file,
    repo=repo,
    dest=dataDir,
    overwrite=TRUE
  )
}
