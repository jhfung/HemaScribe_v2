HemaScribeData <- new.env(parent=emptyenv())

.onAttach <- function(libname, pkgname) {
  if (!nzchar(Sys.getenv("R_CMD_CHECK"))) {
    dataDir <- tools::R_user_dir("HemaScribe", which = "data")
    dataIndex <- file.path(dataDir, "HemaScribe-data.rdx")
    dataBase <- file.path(dataDir, "HemaScribe-data.rdb")

    if (!all(file.exists(dataIndex, dataBase))) {
        download_data(dataDir)
    }

    if (all(file.exists(dataIndex, dataBase))) {
      lazyLoad(file.path(dataDir, "HemaScribe-data"), envir=HemaScribeData)
    } else {
      rlang::abort("HemaScribe cannot be loaded because reference data is not found.")
    }

    # Handle mapping internal data
    mappingDataFile <- file.path(dataDir, "HemaScape-data.rda")
    if (!file.exists(mappingDataFile)) {
      download_mapping_data(dataDir)
    }

    if (file.exists(mappingDataFile)) {
      load(mappingDataFile, envir = HemaScribeData)
    } else {
      rlang::warn("HemaScape mapping reference data is not found.")
    }
  }
}

#' Download reference data
#'
#' Downloads the reference data from the Github repo.  The user is prompted for
#' download when the package is first loaded, but this function can also be
#' called in a later session to update the reference data file.
#' @param dataDir Specifies where the reference data should be downloaded.  Should essentially never be altered by the user directly since the package looks for the data in a specific location.
#' @export
download_data <- function(dataDir = NA) {
  if (is.na(dataDir)) {
    dataDir <- tools::R_user_dir("HemaScribe", which = "data")
  }

  if (!dir.exists(dataDir)) {
    dir.create(dataDir, recursive = TRUE)
  }

  piggyback::pb_download(
    file="HemaScribe-data.rdx",
    repo="RabadanLab/HemaScribe",
    dest=dataDir,
    overwrite=TRUE
  )

  piggyback::pb_download(
    file="HemaScribe-data.rdb",
    repo="RabadanLab/HemaScribe",
    dest=dataDir,
    overwrite=TRUE
  )
}

#' Download mapping reference data
#'
#' Downloads the mapping reference data from the Github repo. The user is prompted for
#' download when the package is first loaded, but this function can also be
#' called in a later session to update the mapping reference data file.
#' @param dataDir Specifies where the mapping reference data should be downloaded. Should essentially never be altered by the user directly since the package looks for the data in a specific location.
#' @export
download_mapping_data <- function(dataDir = NA) {
  if (is.na(dataDir)) {
    dataDir <- tools::R_user_dir("HemaScribe", which = "data")
  }

  if (!dir.exists(dataDir)) {
    dir.create(dataDir, recursive = TRUE)
  }

  # Download mapping internal data from GitHub release using piggyback
  piggyback::pb_download(
    file = "HemaScape-data.rda",  # file name in the release
    repo = "RabadanLab/HemaScribe",   # repository
    dest = dataDir,               # destination directory
    overwrite = TRUE              # overwrite if the file already exists
  )
}












