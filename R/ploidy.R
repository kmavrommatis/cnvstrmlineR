#' Get the  ploidy for the sample
#'
#' Two parameters can be retrieved:
#'
#' **rounded ploidy** is the rounded ploidy (e.g. 2.1 will become 2). This is the default.
#'
#' **ploidy** is the actual value provided by the user
#'
#'
#' @param .Object Object of class cnvstrml.
#' @param rounded Logical. If TRUE, return the rounded ploidy. Otherwise return the raw ploidy.
#' @docType methods
#' @return The sample ploidy.
#' @aliases getPloidy,cnvstrml,logical-method getPloidy,cnvstrml,ANY-method
#'
setGeneric("getPloidy", function(.Object, rounded) {
  standardGeneric("getPloidy")
})

#' Get the  ploidy for the sample
#'
#' @rdname getPloidy
#' @name getPloidy
#' @export getPloidy
#'
setMethod(
  "getPloidy",
  signature(.Object = "cnvstrml", rounded = "logical"),
  function(.Object, rounded) {
    if (rounded == FALSE) {
      return(.Object@ploidy)
    }
    return(round(.Object@ploidy))
  }
)

#' Get the  ploidy for the sample
#'
#' @rdname getPloidy
#' @name getPloidy
#' @export getPloidy
#'
setMethod(
  "getPloidy",
  signature(.Object = "cnvstrml"),
  function(.Object, rounded) {
    return(getPloidy(.Object, TRUE))
  }
)


#' Calculate the ploidy
#'
#' @param segs A GRanges object that contains the segments.
#' @return Numeric effective ploidy inferred from the dominant copy-number state.
#' @details
#' The ploidy that is provided by the tools is typically the weighted average of the ploidy of all the segments
#' This function is returning the integer ploidy of the segments that account for the majority of bases
#'
calculateEffectivePloidy <- function(segs) {
  stats <- segs %>%
    as.data.frame() %>%
    dplyr::group_by(major.allele.copies, minor.allele.copies) %>%
    dplyr::summarise(n_len = sum(width, na.rm = TRUE), n_count = dplyr::n()) %>%
    dplyr::arrange(dplyr::desc(n_len)) %>%
    dplyr::mutate(ploidy = major.allele.copies + minor.allele.copies) %>%
    tidyr::drop_na(ploidy)
  ploidy <- stats$ploidy[1]

  ploidy
}
