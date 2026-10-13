#' Get the  purity for the sample
#'
#' Two parameters can be retrieved:
#'
#' **purity** is the actual value provided by the user
#'
#' @docType methods
#' @param .Object object of type cnvstrml,
#' @return the purity of the sample (0 - 1)
#'
setGeneric("getPurity", function(.Object){
  standardGeneric("getPurity")
})

#' Get the  purity for the sample
#'
#' @rdname getPurity
#' @export
#'
setMethod("getPurity",
          signature(.Object="cnvstrml"),
          function(.Object){
            return(.Object@purity)
          }
)



