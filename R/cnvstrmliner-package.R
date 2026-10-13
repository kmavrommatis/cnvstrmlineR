#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom methods callNextMethod
#' @importFrom methods new
#' @importFrom data.table := .N
#' @importClassesFrom GenomicRanges GRanges
#' @name cnvstrmliner-overview
#' @title cnvstrmliner: overview
#' @details
#' CNV streamliner.
#' The aim of this package is to streamline the analysis of results produced by CNV tools.
#' It introduces the cnvstrml object and parsers for commonly used methods.
#'
#'
#'
#'
#'


#' It has some methods to compare different calls.
#' Available support for
#' \itemize{
#' \item FACETS,
#' \item PureCN,
#' \item SCLUST,
#' \item GATK CopyNumber,
#' \item ControlFreec (without subclones),
#' \item Battenberg,
#' \item SEG (Generic input file for IGV).
#' }
#' The Outputs can be
#' \itemize{
#' \item VCF (v4.2 compliant),
#' \item SEG
#' }
#'
#' This package contains methods to visualize the segments
#' of a genome or of a chromosome
#' or a cohort of samples\cr
#' It also provides metrics for the segmentation result and allows to compare
#' CNV calls among samples using standard CNV metrics.
#' Using the CNAQC package it provides a PASS/FAIL result for the sample and for
#' each segment.
#'
#' For every object we introduce a number of columns that enable to normalize the information provided by different methods:
#' \itemize{
#' \item  total.allele.copies: the total copies found in a segment (major.allele.copies + minor.allele.copies)
#' \item  major.allele.copies: the copies of the major allele
#' \item  minor.allele.copies: the copies of the minor allele
#' \item  cnv.totalmarkers: total number of markers SNPs used to test this segment. If this information is not provided by the tool it is set to NA
#' \item  cnv.markers: number of heterozygous SNPs used to test this segment
#' \item  cnv.annotation: annotation of the segment (NEUT, AMP, DEL etc). See `annotationCNV` for more information
#' \item  cnv.lrr: the log2 ratio of that segment.
#' \item  cnv.cf: This is the cell fraction that has the event. Typically the purity of the tumor sample equals the highest cnv.cf.
#' \item  cnv.cff: This is the Cancer cell fraction. Events appearing in all cancer cells (CCF = 1) are considered clonal. To find the prevalence (cn.cf) of these cells in the sample divide this by the purity.
#' \item  cnv.loh: This is set to TRUE if the segment has lost its heterozygosity. Depending on the total copy number it can be a deletion, neutral LOH or even LOH in an amplified region.
#' \item  purity: is the predicted overall purity of the sample.
#' \item  ploidy: is the ploidy of the sample
#' \item  effective_ploidy: is the ploidy of the sample rounded to the closest integer.
#' \item  sex: the sex of the sample.
#' Depending on the method either ccf or cf is typically provided. If none is provided the default value is set to purity i.e. assumed that all the cells have this segment
#'
#' }
#'
## usethis namespace: end
NULL
