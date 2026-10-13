#' Annotate a CNV
#'
#' annotate a CNV segment with a flag describing a Copy number aberration according to the table
#'
#' |total cn | minor cn |	type         | meaning of type                                                        |
#' |---------|----------|--------------|------------------------------------------------------------------------|
#' |0        |	0 	    | HOMDEL       | Homozygous deletion                                                    |
#' |0        |	NA 	    | HOMDEL       | Homozygous deletion                                                    |
#' |1        |	0 	    | HETDEL       | Heterozygous deletion (1 out of 2 alleles)                             |
#' |1        |	NA 	    | HETDEL       | Heterozygous deletion (this is an edge case that should not occur).    |
#' |2        |	0 	    | NEUT-LOH     | Loss of heterozygosity without change in the number of alleles         |
#' |2        |	1 	    | NEUT         | Normal                                                                 |
#' |2        |	NA 	    | NEUT/Unknown | Normal or undetermined by the algorithm                                |
#' |2+       |	0 	    | AMP-LOH      | Amplification with loss of heterozygosity                              |
#' |2+       |	1+ 	    | AMP          | Amplification                                                          |
#' |2+       |	NA 	    | AMP/LOH.     | Amplification with (possible) loss of heterozygosity. Cannot determine.|
#'
#'
#' @param major.cn Major Allele copy number (integer)
#' @param minor.cn Minor Allele copy number (integer)
#' @param total.cn Total copy number
#' @param na.minor  The default value that the minor allele will take if it is NA. Defaults to 0
#' @param ploidy Integer. The assumed ploidy for the genome. Used to set the amplification and deletion flags.
#' @section Notes:
#' The major.cn and minor.cn should add up to total.cn.
#' Either total.cn or major.cn are necessary.
#' If one of these is not provided the function
#' will calculate it assuming total.cn=major.cn + minor.cn
#' Minor.cn is not necessary, it defaults to 0 meaning that there is no minor allele
#' @section Ploidy consideration:
#' This function assumes that normal ploidy is 2.
#' @return Returns a string indicating the type of the copy number aberration.
#'
#' @export
annotateCNV <- function(major.cn = NA, minor.cn = 0, total.cn = NA, na.minor = 0, ploidy = 2) {
  # message(major.cn, " ", class(major.cn)," " ,minor.cn," ",class(minor.cn), " " , total.cn)
  if (is.nan(total.cn) && is.nan(major.cn)) {
    if (is.na(total.cn) && is.na(major.cn)) {
      stop("Please provide either the total number of alleles or the major allele copies. Received total.cn :", total.cn, " major.cn :", major.cn, " minor.cn :", minor.cn, " for ploidy: ", ploidy)
    }
  }
  major.cn <- as.numeric(major.cn)
  minor.cn <- as.numeric(minor.cn)
  total.cn <- as.numeric(total.cn)
  if (is.na(total.cn)) {
    total.cn <- ifelse(is.na(minor.cn), major.cn, major.cn + minor.cn)
  }
  if (total.cn < 0) {
    stop("total.cn (major.cn + minor.cn) cannot be a negative number")
  }
  if (is.na(minor.cn)) {
    # message("The minor.cn is not defined. This is a tricky situation in most cases, we default to 0 minor allele copies, but can be changed with the na.minor argument")
    minor.cn <- na.minor
  }
  if (minor.cn < 0) {
    stop("minor.cn (major.cn + minor.cn) cannot be a negative number")
  }
  if (ploidy != as.integer(ploidy)) {
    stop("Ploidy should be set to an integer value, but ", ploidy, " was provided.")
  }
  flag <- "UNKNOWN"
  # set the flag
  if (total.cn == 0) {
    flag <- "HOMDEL"
  }
  if ((total.cn < ploidy & total.cn > 0) & (minor.cn < total.cn | is.na(minor.cn))) {
    flag <- "HETDEL"
  } # Loss of copies with regards to Neutral
  if ((total.cn < ploidy & total.cn > 1) & (minor.cn == 0)) {
    flag <- "HETDEL-LOH"
  } # Loss of copies, and minor allele is 0. Only if the major allele >1 (otherise HETDEL)
  if (total.cn == ploidy & (!is.na(minor.cn) & minor.cn == 0)) {
    flag <- "NEUT-LOH"
  } # Neutral but minor allele is 0
  if (total.cn == ploidy & (!is.na(minor.cn) & (minor.cn > 0 & minor.cn < ploidy))) {
    flag <- "NEUT"
  } # Neutral with major and minor alleles > 0
  if (total.cn == ploidy & is.na(minor.cn)) {
    flag <- "NEUT/UNKNOWN"
  }
  if (total.cn > ploidy & (!is.na(minor.cn) & minor.cn == 0)) {
    flag <- "AMP-LOH"
  } # Amplification but minor allele is 0
  if (total.cn > ploidy & (!is.na(minor.cn) & minor.cn > 0)) {
    flag <- "AMP"
  } # Amplification with respect to Neutral
  if (total.cn > ploidy & is.na(minor.cn)) {
    flag <- "AMP/LOH"
  } # We don't have information about the minor allele, but we have an amplification

  return(flag)
}


#' Adjust the logR
#'
#' @param cnvObj A cnvstrml object whose LogR values will be adjusted.
#' @param usePurity Logical. If TRUE it uses the overall purity of the sample. Otherwise it uses the CCF per segment.
#' @return An updated cnvstrml object with the column cnv.lrrAdj.
#'
#' @description
#' Given the ploidy of genome and its purity The observed log2 ratios are treated as a mix of some fraction of tumor cells (specified by purity),
#' possibly with altered copy number, and a remainder of normal cells with neutral copy number (specified by ploidy for autosomes) by default, diploid autosomes,
#' Assuming an initial ploidy the various copy numbers (copies) are transformed to logR.
#'
#' logR= log2(copies/ploidy)
#' However taking the purity (or CCF if we want to be specific for each segment) this transformation becomes
#' logR=log2( (1 - ccf) + ccf * copies / ploidy )
#' We can reverse that equation and if we know the CCF we can extract the copy numbers and the set CCF to 1 and recalculate the logRadjusted.
#' \deqn{logR_{adjusted}=log2( (2^{logR} -(1-ccf))/ccf )}
#'
#' For the sex chromosomes. If the sex is "F" we assume ploidy = ploidy of the genome for X
#' If the sex is "M" we set ploidy =1 for the Y and ploidy -1 for X
#'
#'
adjustLogR <- function(cnvObj, usePurity = FALSE) {
  purity <- S4Vectors::mcols(cnvObj@segments)[, "cnv.cf"]
  if (usePurity == TRUE) {
    purity[] <- cnvObj@purity
  }
  S4Vectors::mcols(cnvObj@segments)["cnv.lrrAdj"] <- suppressWarnings(
    log2((2^S4Vectors::mcols(cnvObj@segments)[, "cnv.lrr"] - (1 - purity)) / purity)
  )
  S4Vectors::mcols(cnvObj@segments)[, "cnv.lrrAdj"] <- ifelse(is.nan(S4Vectors::mcols(cnvObj@segments)[, "cnv.lrrAdj"]),
    min(S4Vectors::mcols(cnvObj@segments)[, "cnv.lrrAdj"], na.rm = TRUE) - 1,
    S4Vectors::mcols(cnvObj@segments)[, "cnv.lrrAdj"]
  )
  cnvObj
}


#' Adjust the CCF
#' @param x An object of class cnvstrml or a numeric vector of cancer fractions.
#' @param ... Additional arguments passed to methods.
#' @return An updated cnvstrml object or the cancer cell fraction.
#' @description
#' Add the ccf of the segment
#' roughly speaking ccf=cf/purity
#' @export
#'
addCCF <- function(x, ...) {
  UseMethod("addCCF")
}

#' Adjust the CCF
#' @param x An object of class cnvstrml.
#' @param ... Additional arguments passed through the generic.
#' @return An updated cnvstrml object.
#' @description
#' Add the ccf of the segment
#' roughly speaking ccf=cf/purity
#' @export
#'
addCCF.cnvstrml <- function(x, ...) {
  # if(!'cnv.cf' %in% colnames(S4Vectors::mcols(x@segments))){
  #   S4Vectors::mcols(x@segments)[,'cnv.cf']=x@purity
  # }
  # ccf=S4Vectors::mcols(x@segments)[,'cnv.cf']
  # ccf[is.na(ccf)]=x@purity
  # ccf=ccf/x@purity
  ccf <- addCCF(S4Vectors::mcols(x@segments)[, "cnv.cf"], purity = x@purity)
  S4Vectors::mcols(x@segments)[, "cnv.ccf"] <- ccf
  x
}
#' Adjust the CCF
#' @param x Numeric vector of cancer fractions.
#' @param purity Numeric purity of the sample.
#' @param ... Additional arguments passed through the generic.
#' @return The CCF of the segment.
#' @description
#' Add the ccf of the segment
#' roughly speaking ccf=cf/purity
#' @export
#'
addCCF.numeric <- function(x, purity, ...) {
  if (is.na(purity) | purity == 0) {
    purity <- 1
  }
  ccf <- ifelse(is.na(x), 1, x / purity)
  ccf <- sapply(ccf, function(x) {
    min(x, 1, na.rm = TRUE)
  })
  ccf
}

#
# addCCF = function( purity, cf ){
#   if(is.na(purity) | purity ==0){
#     purity=1
#   }
#   ccf=ifelse(is.na(cf), 1, cf/purity )
#   ccf=min(1, ccf, na.rm=TRUE)
#   ccf
# }


#' Annotate an LOH
#'
#' annotate a CNV segment with a flag describing if it is an LOH or not
#'
#' @param major.cn Major Allele copy number (integer)
#' @param minor.cn Minor Allele copy number (integer)
#' @param chromosome (character) used to check if this is a segment on a sex chromosome. If not provided it is assumed an autosome.
#' @param sex (character) The sex of this sample (M or F). If not provided it is assumed F
#' @section Notes:
#' An LOH is called when the minor.cn is set to 0.
#' If this is on a sex chromosome (X, Y, chrX, chrY) and the sex is set to M(ale) then it is not considered LOH
#' @return Returns a string indicating the type of the copy number aberration.
#'
#' @export
addLOH <- function(major.cn = NA, minor.cn = 0, chromosome = NA, sex = "F") {
  # message(major.cn, " ", class(major.cn)," " ,minor.cn," ",class(minor.cn), " " , total.cn)
  major.cn <- as.numeric(major.cn)
  minor.cn <- as.numeric(minor.cn)

  if (is.na(minor.cn)) {
    # message("The minor.cn is not defined. This is a tricky situation in most cases, we default to 0 minor allele copies, but can be changed with the na.minor argument")
    minor.cn <- 0
  }
  LOH <- FALSE
  if (minor.cn == 0) {
    LOH <- TRUE
  }
  if (sex == "M" & chromosome %in% c("chrX", "chrY", "X", "Y")) {
    LOH <- FALSE
  }


  return(LOH)
}


#' Annotate segments with information from the chromosomes, telomers, centromeres and arms
#'
#' @return An annotated object with chromosome arm information.
#' @param x A cnvstrml or GRanges object to annotate.
#' @param ... Additional arguments passed to methods.
#' @details
#' Adds a column in the GRanges of the object chr.annotation which contains
#' information about the arm of the chromosome, the telomere and centromere
#' if a segment overlaps more than one of these annotations it will be divided accordingly
#'
#' @export

addChromosomeAnnotation <- function(x, ...) {
  UseMethod("addChromosomeAnnotation")
}
#' Annotate segments with information from the chromosomes, telomers, centromeres and arms
#'
#' @return An updated cnvstrml object.
#' @param x A cnvstrml object to annotate.
#' @param ... Additional arguments passed through the generic.
#' @details
#' Adds a column in the GRanges of the object chr.annotation which contains
#' information about the arm of the chromosome, the telomere and centromere
#' if a segment overlaps more than one of these annotations it will be divided accordingly
#'
#' @export
addChromosomeAnnotation.cnvstrml <- function(x, ...) {
  segs <- x@segments
  genome <- x@genome
  segs <- cnvstrmliner::addChromosomeAnnotation(segs, genome)
  x@segments <- segs
  x
}


#' Annotate segments with information from the chromosomes, telomers, centromeres and arms
#'
#' @return An updated GRanges object.
#' @param x A GRanges object to annotate.
#' @param genome Character genome build used to select the chromosome annotations.
#' @param ... Additional arguments passed through the generic.
#' @details
#' Adds a column in the GRanges of the object `arm` which contains
#' information about the arm of the chromosome, the telomere and centromere
#' if a segment overlaps more than one of these annotations it will be divided accordingly
#'
#' @export
addChromosomeAnnotation.GRanges <- function(x, genome, ...) {
  meres <- cnvstrmliner::chrom_annot[[genome]]
  ov <- GenomicRanges::findOverlaps(query = x, subject = meres)
  x <- GenomicRanges::pintersect(x[S4Vectors::queryHits(ov)], meres[S4Vectors::subjectHits(ov)])
  ov <- GenomicRanges::findOverlaps(query = x, subject = meres)
  S4Vectors::mcols(x)["arm"] <- NA
  S4Vectors::mcols(x)[S4Vectors::queryHits(ov), "arm"] <- S4Vectors::mcols(meres)[S4Vectors::subjectHits(ov), "source"]

  return(x)
}
