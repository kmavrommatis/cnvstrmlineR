
#' Infer the sex for a sample
#'
#' If events occur on the Y chromosome the sample
#' is assigned to XY (male)
#' otherwise XX (female)
#' @param .Object object of type GRanges or GRangeslist
#' @return the sex of the sample as XX for female and XY for male
#' @docType methods
#' @section  Notes:
#' The inference of the sex depends on the identification of events on the sex chromosome
#' If any event is observed on chrY then it is assumed that the sample is male
#' The return value is M for male or F for female
#'
setGeneric("inferSex", function(.Object){
  standardGeneric("inferSex")
})


#' Infer the sex for a sample
#'
#' @rdname inferSex
#'
setMethod("inferSex",
          signature(.Object="cnvstrml"),
          function(.Object){
            seqnames=GenomeInfoDb::seqnames(.Object@segments) %>% as.vector()
            male=which( seqnames %in% c("Y",'chrY') )
            if(length(male)>0){ sex='M'}else{sex='F'}

            return(sex)
          }
)



#' Get the  sex for the sample
#'
#' @param .Object object of type cnvstrml,
#' @return the sex (M for male, F for female)
#' @export getSex
#' @docType methods
#' @name getSex
setGeneric("getSex", function(.Object){
  standardGeneric("getSex")
})

#' Get the  sex for the sample
#'
#' @rdname getSex
#'
setMethod("getSex",
          signature(.Object="cnvstrml"),
          function(.Object){
            return(.Object@sex)
          }
)

