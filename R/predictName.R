#' create the proper name for sclust files
#' @description
#' given any filename from files created by Sclust
#' it produces a string that can be used as argument
#' for the parseSCLUST function
#' @param filename (character) the filename of on of the sclust files
#' @param user.extension (dataframe) if set then the specified extensions are added to the standard list of extensions
#' @return a character used as sclust.name
#' @details
#' SClust produces
#' \itemize{
#' \item<sclust.name>_cn_summary.txt # contains general information about the purity, ploidy, contamination etc.
#' \item<sclust.name>_subclonal_cn.txt# contains the segments along with the copy number information.
#' \item<sclust.name>_allelic_states.txt # contains the information of segments with allele copies.
#' }
#'
#' ControlFreec produces
#' \itemize{
#' \item<cfreec>_info.txt
#' \item<cfreec>_CNVs
#' \item<cfreec>_ratio.txt
#' \item<cfreec>_BAF.txt
#' }
#'
#' Sequenza produces
#' \itemize{
#' \item<sequenza>_confints_CP.txt
#' \item<sequenza>_segments.txt
#' }
#'
#'
#'
#'
predictName <- function(filename, user.extension = NULL) {
  fn <- predictMethod(filename, method = FALSE, user.extension = NULL)
  if (!is.na(fn)) {
    filename <- fn
  }
  return(filename)
}


#' create the proper tool for parsing files
#'
#' given any filename from files created by tools
#' it produces a string that indicates which method has created this file
#' @param filename (character) the filename of on of the sclust files
#' @param method (logical) if set to TRUE it returns the method that corresponds to this extension. If set to FALSE it removes the extension
#' @param user.extension (dataframe) if set then the specified extensions are added to the standard list of extensions
#' @return a character
#' @description
#' Sclust produces
#' <sclust.name>_cn_summary.txt # contains general information about the purity, ploidy, contamination etc.
#' <sclust.name>_subclonal_cn.txt# contains the segments along with the copy number information.
#' <sclust.name>_allelic_states.txt # contains the information of segments with allele copies.
#' <sclust.name>_uncorr_cn.seg,
#' <sclust.name>_muts_expAF.txt
#' Standard file extensions for IGV (in SEG format)
#' <name>.seg,
#' <name>.cbs,
#' The user.extensions dataframe should have two columns
#'  extension: contains the extension of the files
#'  method: contains the name of hte method that corresponds to this extension
#'
#'
predictMethod <- function(filename, method = TRUE, user.extension = NULL) {
  extensions <- data.frame(
    "extension" = c(
      "_cn_summary.txt", "_subclonal_cn.txt", "_allelic_states.txt", "_muts_expAF.txt", "_iCN.seg",
      ".seg", ".cbs",
      "_loh.csv", "_variants.csv",
      "_info.txt", "_CNVs", "_ratio.txt", "_BAF.txt", "_CNVs.p.value.txt",
      "_subclones.txt", "_rho_and_psi.txt", "_cellularity_ploidy.txt",
      "_confints_CP.txt", "_segments.txt", "_seqzinfo.txt", "-seqzinfo.txt",
      ".results"
    ),
    "method" = c(
      rep("sclust", times = 5),
      rep("seg", times = 2),
      rep("purecn", times = 2),
      rep("freec", times = 5),
      rep("battenberg", times = 3),
      rep("sequenza", times = 4),
      rep("clonecna", times = 1)
    )
  )
  if (!is.null(user.extension)) {
    extensions <- rbind(user.extension, extensions)
  }
  # get the method based on the extension of a given file
  if (method == TRUE) {
    for (ext in 1:length(extensions)) {
      if (endsWith(filename, extensions$extension[ext])) {
        return(extensions$method[ext])
      }
    }
  }
  # remove the extension from the file.
  if (method == FALSE) {
    for (ext in extensions$extension) {
      if (endsWith(filename, ext)) {
        return(gsub(ext, "", filename))
      }
    }
  }

  return(NA)
}


#' find a file that matches the pattern provided
#'
#' @details
#' This function finds files that match a given pattern and
#' returns the first one that matches.
#'
#' @param filename (character) the filename, either a character or a regex
#' @return filename of existing file or NA
#'
findFile <- function(filename) {
  if (file.exists(filename)) {
    return(filename)
  }
  # parent directory
  parent_directory <- dirname(filename)
  # List all files in the parent directory
  matched_files <- fs::dir_info(
    path = parent_directory,
    recurse = FALSE,
    type = "any",
    glob = paste0("*", basename(filename))
  )

  # If no files match the regular expression, return NA
  if (length(matched_files) == 0) {
    return(NA)
  }

  # Get the creation times of all matched files
  creation_times <- matched_files$modification_time

  # Find the index of the most recently created file
  max_index <- which.max(creation_times)

  # Return the name of the most recently created file
  filename <- matched_files[max_index, "path", drop = TRUE]

  return(filename)
}
