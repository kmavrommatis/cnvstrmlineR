#' Chromosome Annotation Tracks
#'
#' Annotation ranges (for example telomeres and centromeres) keyed by genome build.
#'
#' @format A named list of genomic annotation objects.
#' @source Internal package data.
"chrom_annot"

#' Chromosome Arm Coordinates
#'
#' Chromosome p/q arm ranges keyed by genome build.
#'
#' @format A named list of \code{GRanges} objects.
#' @source Internal package data.
"chrom_arms"

#' Chromosome Gap Coordinates
#'
#' Gap annotations (centromere, telomere, short arms) keyed by genome build.
#'
#' @format A named list of genomic gap annotations.
#' @source Internal package data.
"chrom_gaps"

#' Chromosome Sequence Information
#'
#' Sequence-level metadata used to construct and validate genomic ranges.
#'
#' @format A named list of \code{Seqinfo} objects.
#' @source Internal package data.
"seqinfo"

#' Chromosome Sizes
#'
#' Chromosome lengths keyed by genome build.
#'
#' @format A named list of chromosome length tables.
#' @source Internal package data.
"chrom_size"
