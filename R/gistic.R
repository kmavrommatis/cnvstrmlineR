#' Export a CNV object to a GISTIC-compatible segment file
#'
#' @description
#' Exports copy number variation (CNV) segments from one or more cnvstrmliner objects
#' to a tab-delimited segment file suitable for cohort-level analysis with
#' GISTIC2.0 (Genomic Identification of Significant Targets in Cancer). Segment means
#' are reported as log2 copy ratios relative to the ploidy of each sample, which
#' removes the whole-genome-doubling offset and makes segment values comparable
#' across samples of differing ploidy.
#'
#' @details
#' The function produces a tab-delimited file with the six columns expected by
#' GISTIC2.0, in this order:
#' \describe{
#'   \item{\strong{Sample}}{The sample name, taken from the \code{sample} slot of the
#'     cnvstrmliner object}
#'   \item{\strong{Chromosome}}{Chromosome identifier, written exactly as stored in the
#'     \code{segments} slot (see the note on chromosome naming below)}
#'   \item{\strong{Start}}{Start position of the segment (1-based coordinates)}
#'   \item{\strong{End}}{End position of the segment (inclusive)}
#'   \item{\strong{Num_Markers}}{Number of heterozygous markers supporting the segment,
#'     taken from the \code{cnv.markers} field}
#'   \item{\strong{Seg.CN}}{Ploidy-relative log2 copy ratio (see below)}
#' }
#'
#' \subsection{Segment mean calculation}{
#' The segment mean is the log2 ratio of the total copy number of the segment to the
#' ploidy of the sample:
#' \deqn{Seg.CN = \log_2\left(\frac{total.allele.copies}{ploidy}\right)}
#'
#' The raw, fractional ploidy reported by the CNV caller is used
#' (\code{getPloidy(cnvObj, rounded = FALSE)}), not the integer \code{effective_ploidy}.
#' A consequence is that a segment sitting exactly at the sample ploidy yields
#' \code{Seg.CN = 0}, while a segment at the nearest integer copy number of a
#' non-integer-ploidy sample yields a small non-zero value. For a sample with
#' \code{ploidy = 2}, \code{Seg.CN} reduces to the familiar \code{log2(CN/2)}:
#' a single-copy gain gives 0.585 and a heterozygous deletion gives -1.
#'
#' If the ploidy of a sample is missing or non-positive it is replaced by 2 and a
#' warning is logged.
#' }
#'
#' \subsection{Homozygous deletions}{
#' Segments with \code{total.allele.copies == 0} would give \code{log2(0) = -Inf},
#' which GISTIC2.0 cannot parse. Such segments are floored at \code{minSegCN}
#' (default \code{-5}), a magnitude comparable to the values emitted by callers such
#' as FACETS and CNVkit for deep deletions. The floor is applied with \code{pmax}, so
#' any finite value falling below \code{minSegCN} is clamped as well.
#' }
#'
#' \subsection{Multi-sample support}{
#' When a list of cnvstrmliner objects is provided, all samples are concatenated into a
#' single output file distinguished by the \code{Sample} column. This is the expected
#' input format for GISTIC2.0, which derives its background model from the whole
#' cohort. Sample names must be unique across the list; duplicate sample names cause
#' an error, since GISTIC would otherwise treat segments from two different objects
#' as belonging to the same sample.
#' }
#'
#' \subsection{Chromosome naming}{
#' Chromosome names are written exactly as they appear in the \code{segments} slot,
#' without adding or stripping a \code{chr} prefix. Note that GISTIC2.0 matches the
#' \code{Chromosome} column against its reference genome marker file, and most
#' distributed reference files use unprefixed names (\code{1}, \code{2}, ..., \code{X},
#' \code{Y}). If the segments carry UCSC-style names (\code{chr1}, \code{chr2}) the
#' prefix may need to be removed before running GISTIC, either upstream with
#' \code{GenomeInfoDb::seqlevelsStyle} or on the resulting file.
#' }
#'
#' \subsection{Missing values}{
#' Segments whose \code{Seg.CN} cannot be computed (missing \code{total.allele.copies})
#' are dropped from the output and the number removed is logged, because GISTIC cannot
#' interpret an \code{NA} segment mean. Missing \code{Num_Markers} values are retained
#' as-is and reported with a warning; some callers (e.g. Manta) do not provide marker
#' counts, and these must be resolved before the file can be used with GISTIC.
#' }
#'
#' @param filename Character string specifying the output file path. The file will be
#'   created or overwritten if it already exists.
#' @param cnvObj A cnvstrmliner object containing CNV segments, or a list of cnvstrmliner
#'   objects for cohort-level output. When providing a list, all sample names must be
#'   unique.
#' @param minSegCN Numeric lower bound applied to the segment means. Defaults to
#'   \code{-5}. Its main purpose is to convert the \code{-Inf} produced by homozygous
#'   deletions into a finite value that GISTIC can parse.
#' @param digits Integer number of decimal places \code{Seg.CN} is rounded to,
#'   defaulting to 4. Set to \code{NULL} or \code{NA} to write full precision.
#'
#' @return Character string containing the path to the output file (\code{filename}),
#'   invisibly returned unchanged whether or not the write succeeded; check the file
#'   externally if confirmation is required.
#'
#' @note
#' \itemize{
#'   \item Column headers are written only for the first sample; subsequent samples are
#'     appended without headers.
#'   \item This function is the preset
#'     \code{saveSEG(segmean = "ploidy", minSegCN = -5, dropNA = TRUE, style = "gistic")}
#'     and shares its write path; anything documented for \code{\link{saveSEG}} applies
#'     here too.
#'   \item No purity correction is applied. Ploidy normalization removes the ploidy
#'     offset but not the dilution caused by normal-cell contamination, which are
#'     different confounders. For cohorts with widely varying purity, consider the
#'     purity-corrected log-ratio mode instead:
#'     \code{saveSEG(logR = TRUE, correct = TRUE, style = "gistic")}.
#'   \item \strong{Subclonal events are written at full strength.} Because the segment
#'     mean is derived from the fitted integer copy number, a three-copy segment
#'     present in only 30\% of tumor cells produces the same \code{Seg.CN = 0.585} as a
#'     fully clonal gain, and GISTIC cannot tell the two apart. This is the price of
#'     being caller-agnostic. If subclonal attenuation should be preserved, use
#'     \code{saveSEG(logR = TRUE, correct = TRUE, usePurity = TRUE, style = "gistic")},
#'     and see the mode comparison table in \code{\link{saveSEG}}.
#' }
#'
#' @seealso
#' \code{\link{saveSEG}}, which this function wraps, for the other segment mean modes
#' and the optional purity correction,
#' \code{\link{getPloidy}} for retrieving the ploidy used as the denominator,
#' \code{\link{writeVCF}} for VCF format output
#'
#' @examples
#' \dontrun{
#' # Export a single sample
#' exportGISTIC("sample1.gistic.seg", cnv_obj)
#'
#' # Export a cohort into a single file for GISTIC2.0
#' cnv_list <- list(sample1_cnv, sample2_cnv, sample3_cnv)
#' exportGISTIC("cohort.gistic.seg", cnv_list)
#'
#' # Use a shallower floor for homozygous deletions
#' exportGISTIC("cohort.gistic.seg", cnv_list, minSegCN = -3)
#'
#' # Drop the chr prefix if the GISTIC reference file uses unprefixed names
#' cnv_list <- lapply(cnv_list, function(X) {
#'   GenomeInfoDb::seqlevelsStyle(X@segments) <- "NCBI"
#'   X
#' })
#' exportGISTIC("cohort.gistic.seg", cnv_list)
#' }
#'
#' @export
exportGISTIC <- function(filename, cnvObj, minSegCN = -5, digits = 4) {
  saveSEG(filename, cnvObj,
    segmean = "ploidy",
    minSegCN = minSegCN,
    dropNA = TRUE,
    style = "gistic",
    digits = digits
  )
}
