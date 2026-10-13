#' Write a SEG file with the contents of a CNV object
#'
#' @description
#' Exports copy number variation (CNV) segments from one or more cnvstrmliner objects
#' to a SEG file format. The SEG format is a tab-delimited text file commonly used
#' by genomic visualization and analysis tools such as IGV (Integrative Genomics Viewer)
#' and GISTIC2.0 (Genomic Identification of Significant Targets in Cancer).
#'
#' @details
#' The function produces a tab-delimited file with the following columns:
#' \describe{
#'   \item{\strong{sample}}{The name of the sample (extracted from the cnvstrmliner object)}
#'   \item{\strong{chrom}}{Chromosome identifier}
#'   \item{\strong{start}}{Start position of the segment (1-based coordinates)}
#'   \item{\strong{end}}{End position of the segment}
#'   \item{\strong{probes}}{Number of markers/probes used by the CNV calling method for this segment}
#'   \item{\strong{segmean}}{Segment mean value, selected by \code{segmean} (see below)}
#' }
#'
#' When \code{style = "gistic"} the same six columns are written under the names
#' GISTIC2.0 uses: \code{Sample}, \code{Chromosome}, \code{Start}, \code{End},
#' \code{Num_Markers} and \code{Seg.CN}.
#'
#' \subsection{Multi-sample Support}{
#' When a list of cnvstrmliner objects is provided, the function concatenates all

#' samples into a single output file. This is the expected input format for tools
#' like GISTIC2.0 that perform cohort-level analysis. Sample names must be unique
#' across the list; duplicate sample names will cause an error.
#' }
#'
#' \subsection{Segment Mean Modes}{
#' The \code{segmean} argument selects what is written in the sixth column:
#' \describe{
#'   \item{\code{"auto"}}{The default. Derives the mode from the \code{logR} and
#'     \code{correct} flags, reproducing the behaviour of earlier versions of this
#'     function: \code{logR = FALSE} yields \code{"copies"}, \code{logR = TRUE} yields
#'     \code{"logR"}.}
#'   \item{\code{"copies"}}{Absolute total copy number (\code{total.allele.copies}).}
#'   \item{\code{"logR"}}{The log2 ratio as reported by the calling tool
#'     (\code{cnv.lrr}), or the purity-corrected \code{cnv.lrrAdj} when
#'     \code{correct = TRUE}.}
#'   \item{\code{"ploidy"}}{The log2 ratio of the total copy number to the ploidy of
#'     the sample, recomputed from the fitted allele copies:
#'     \deqn{\log_2\left(\frac{total.allele.copies}{ploidy}\right)}
#'     The raw fractional ploidy is used (\code{getPloidy(X, rounded = FALSE)}), not
#'     the integer \code{effective_ploidy}. If the ploidy of a sample is missing or
#'     non-positive it is replaced by 2 and a warning is logged.}
#' }
#'
#' The \code{"logR"} and \code{"ploidy"} modes differ in more than scale.
#' \code{"logR"} passes through the caller's own measurement, centred on whatever
#' baseline that tool chose, so a whole-genome-doubled sample usually appears flat and
#' a cohort mixing several callers mixes several baselines. \code{"ploidy"} is derived
#' from the discrete fitted copy number against an explicit ploidy, which is
#' caller-agnostic and consistent across a cohort, at the cost of discarding the
#' sub-integer signal that \code{cnv.lrr} carries.
#'
#' Because the purity correction rescales \code{cnv.lrr}, which the ploidy-relative
#' mode never reads, combining \code{correct = TRUE} with \code{segmean = "ploidy"}
#' is an error rather than a no-op. With \code{segmean = "copies"} the
#' \code{correct} flag is ignored, as in previous versions.
#' }
#'
#' \subsection{Log-Ratio Correction}{
#' When both \code{logR = TRUE} and \code{correct = TRUE}, the function applies a
#' purity/ploidy correction to rescale log2 ratios to values expected in a completely
#' pure (100\% tumor) sample. This correction removes the dilution effect caused by
#' normal cell contamination.
#'
#' The correction formula is:
#' \deqn{logR_{corrected} = \log_2\left(\frac{2^{logR_{observed}} - (1 - purity)}{purity}\right)}
#'
#' By default, the cellular fraction (CCF) of each segment is used as the purity value
#' for correction. When \code{usePurity = TRUE}, the global sample purity is used instead.
#'
#' The observed log2 ratios represent a mixture of:
#' \itemize{
#'   \item Tumor cells (fraction = purity) with potentially altered copy number
#'   \item Normal cells (fraction = 1 - purity) with neutral copy number
#' }
#' Normal cells are assumed to have diploid autosomes (CN = 2) and haploid sex
#' chromosomes depending on sample sex.
#'
#' Note that purity correction and ploidy normalization address different
#' confounders: the former removes normal-cell dilution, the latter removes the
#' whole-genome-doubling offset.
#' }
#'
#' \subsection{GISTIC Compatibility}{
#' For GISTIC2.0 analysis, prefer \code{\link{exportGISTIC}}, which is the preset
#' \code{saveSEG(segmean = "ploidy", minSegCN = -5, dropNA = TRUE, style = "gistic")}
#' and produces GISTIC's column names along with the value hygiene that tool needs.
#'
#' If you call \code{saveSEG} directly for GISTIC input, note that the default
#' \code{segmean = "auto"} with \code{logR = FALSE} writes absolute copy numbers;
#' GISTIC would read \code{CN = 2} as a large amplification. Set \code{logR = TRUE}
#' or \code{segmean = "ploidy"}.
#' }
#'
#' \subsection{Value Hygiene}{
#' Coordinates are always written as integers, so \code{write.table} cannot emit them
#' in scientific notation (e.g. \code{5e+05}), which downstream parsers reject.
#'
#' \code{minSegCN} floors the segment means. Its main purpose is the
#' \code{segmean = "ploidy"} mode, where a homozygous deletion gives
#' \code{log2(0) = -Inf}; a finite floor keeps the segment in the file as a strong
#' deletion signal. At the default of \code{-Inf} no flooring occurs.
#'
#' \code{dropNA} removes rows whose segment mean could not be computed. It defaults to
#' \code{FALSE}, so such rows are written as \code{NA} unless the caller opts in.
#'
#' Segment means are rounded to \code{digits} decimal places, 4 by default. The purity
#' correction inverts a log/exponential round trip, which leaves residue on the order of
#' 1e-16: without rounding an exactly neutral segment is written as
#' \code{-4.80513975572238e-16} rather than \code{0}. Four decimals give a resolution of
#' 1e-4 in log2 space, far finer than GISTIC's usual amplification and deletion
#' thresholds of 0.1 and 0.3, and finer than anything IGV displays. Pass
#' \code{digits = NULL} to keep full precision.
#' }
#'
#' \subsection{Choosing a Segment Mean Mode}{
#' The modes respond very differently to tumor purity and to subclonal events. The
#' table below is for a diploid sample (ploidy 2) carrying a segment of the given total
#' copy number, where \code{CCF} is the fraction of \emph{tumor} cells carrying the
#' event, so the fraction of all cells is \code{purity * CCF}:
#'
#' | CN | purity | CCF | `copies` | `ploidy` | `logR` | `correct` | `correct` + `usePurity` |
#' | --: | --: | --: | --: | --: | --: | --: | --: |
#' | 2 | any | any | 2 | 0 | 0 | 0 | 0 |
#' | 3 | 0.30 | 1.0 | 3 | 0.585 | 0.2016 | 0.585 | 0.585 |
#' | 3 | 0.80 | 1.0 | 3 | 0.585 | 0.4854 | 0.585 | 0.585 |
#' | 3 | 1.00 | 1.0 | 3 | 0.585 | 0.585 | 0.585 | 0.585 |
#' | 3 | 0.80 | 0.6 | 3 | 0.585 | 0.3103 | 0.585 | 0.3785 |
#' | 3 | 0.80 | 0.3 | 3 | 0.585 | 0.1635 | 0.585 | 0.2016 |
#' | 3 | 0.30 | 0.3 | 3 | 0.585 | 0.0635 | 0.585 | 0.2016 |
#'
#' These are the values as written at the default \code{digits = 4}.
#'
#' Reading across the table:
#' \itemize{
#'   \item \strong{Two copies against ploidy 2 is invisible in every log-ratio mode.}
#'     Two copies on a diploid background is no change, so the log ratio is 0 at every
#'     purity and every CCF, and a "subclonal" non-event has no meaning. Only
#'     \code{"copies"} distinguishes such a segment, by writing \code{2}.
#'   \item \strong{Uncorrected \code{"logR"} is compressed by purity}, following
#'     \code{log2(1 + 0.5 * purity * CCF)} for a single-copy gain. A clonal gain reads
#'     0.585 at full purity but only 0.2016 at 30\%. Comparing raw \code{"logR"}
#'     across samples of differing purity is therefore not a like-for-like comparison.
#'   \item \strong{\code{correct = TRUE} discards subclonality.} It divides by the
#'     per-segment \code{cnv.cf}, which bundles purity and CCF together, so all six
#'     three-copy rows collapse to 0.585: a gain present in 30\% of tumor cells is
#'     reported identically to a fully clonal one.
#'   \item \strong{\code{correct = TRUE} with \code{usePurity = TRUE} separates the
#'     two.} Dividing by the global purity alone leaves \code{log2(1 + 0.5 * CCF)},
#'     which depends on CCF but not on purity. Note that purity 0.80/CCF 0.3 and
#'     purity 0.30/CCF 0.3 both give 0.2016 even though their raw log ratios differ
#'     (0.1635 against 0.0635). Use this mode when subclonality matters.
#'   \item \strong{\code{"copies"} and \code{"ploidy"} ignore both effects}, being
#'     derived from the fitted integer allele copies alone. That makes them
#'     caller-agnostic and directly comparable across a cohort, but it also means a
#'     subclonal gain is written at full strength.
#' }
#' }
#'
#' @param filename Character string specifying the output file path. The file will be
#'   created or overwritten if it already exists.
#' @param cnvObj A cnvstrmliner object containing CNV segments, or a list of cnvstrmliner
#'   objects for multi-sample output. When providing a list, all sample names must be
#'   unique.
#' @param logR Logical. If \code{TRUE}, the segment mean values written to the file
#'   will be log2 ratios (\code{cnv.lrr} or \code{cnv.lrrAdj} if corrected). If \code{FALSE}
#'   (default), absolute total copy numbers (\code{total.allele.copies}) are written.
#'   Only consulted when \code{segmean = "auto"}.
#' @param correct Logical. If \code{TRUE} and the resolved mode is \code{"logR"}, applies
#'   purity/ploidy correction to the log2 ratios using \code{\link{adjustLogR}}. Ignored
#'   for \code{segmean = "copies"}, and an error for \code{segmean = "ploidy"}.
#'   Default is \code{FALSE}.
#' @param usePurity Logical. Only relevant when \code{correct = TRUE}. If \code{TRUE},
#'   uses the global sample purity for log2 ratio correction. If \code{FALSE} (default),
#'   uses the segment-specific cellular fraction (CCF) for correction.
#' @param segmean Character string selecting the value written in the sixth column:
#'   one of \code{"auto"} (default), \code{"copies"}, \code{"logR"} or \code{"ploidy"}.
#'   See the Segment Mean Modes section.
#' @param minSegCN Numeric lower bound applied to the segment means. Defaults to
#'   \code{-Inf}, which applies no floor.
#' @param dropNA Logical. If \code{TRUE}, segments whose segment mean is \code{NA} are
#'   removed from the output and the number dropped is logged. Default is \code{FALSE}.
#' @param style Character string selecting the column names: \code{"seg"} (default) for
#'   the lowercase SEG/IGV names, or \code{"gistic"} for GISTIC2.0's names.
#' @param digits Integer number of decimal places the segment means are rounded to,
#'   defaulting to 4. Set to \code{NULL} or \code{NA} to write full precision. See the
#'   Value Hygiene section for why the default rounds.
#'
#' @return Character string containing the path to the output SEG file (\code{filename}),
#'   or \code{NULL} if an error occurred during file writing.
#'
#' @note
#' \itemize{
#'   \item The function validates that all samples in a list have unique names before
#'     processing.
#'   \item For single-sample analysis, a single cnvstrmliner object can be passed directly
#'     (not wrapped in a list).
#'   \item Column headers are only written for the first sample; subsequent samples are
#'     appended without headers.
#'   \item The \code{probes} column corresponds to the \code{cnv.markers} field in the
#'     cnvstrmliner object (number of heterozygous markers).
#'   \item Missing marker counts are reported with a warning only when
#'     \code{style = "gistic"}, where GISTIC requires a value; they are unremarkable
#'     for IGV output.
#' }
#'
#' @seealso
#' \code{\link{exportGISTIC}} for the GISTIC2.0 preset of this function,
#' \code{\link{adjustLogR}} for details on the log-ratio correction algorithm,
#' \code{\link{writeVCF}} for VCF format output,
#' \code{\link{parseFACETS}} and other parse functions for creating cnvstrmliner objects
#'
#' @examples
#' \dontrun{
#' # Export a single sample with absolute copy numbers
#' saveSEG("sample1.seg", cnv_obj)
#'
#' # Export with log2 ratios for IGV visualization
#' saveSEG("sample1.seg", cnv_obj, logR = TRUE)
#'
#' # Export with corrected log2 ratios
#' saveSEG("sample1.seg", cnv_obj, logR = TRUE, correct = TRUE)
#'
#' # Use global purity instead of segment CCF for correction
#' saveSEG("sample1.seg", cnv_obj, logR = TRUE, correct = TRUE, usePurity = TRUE)
#'
#' # Multi-sample export
#' cnv_list <- list(sample1_cnv, sample2_cnv, sample3_cnv)
#' saveSEG("cohort.seg", cnv_list, logR = TRUE, correct = TRUE)
#'
#' # Ploidy-relative log2 ratios under GISTIC's column names.
#' # exportGISTIC() is the shorthand for this call.
#' saveSEG("cohort.gistic.seg", cnv_list,
#'   segmean = "ploidy", minSegCN = -5,
#'   dropNA = TRUE, style = "gistic"
#' )
#' }
#'
#' @export
#' @importFrom dplyr mutate select

saveSEG <- function(filename, cnvObj, logR = FALSE, correct = FALSE, usePurity = FALSE,
                    segmean = c("auto", "copies", "logR", "ploidy"),
                    minSegCN = -Inf,
                    dropNA = FALSE,
                    style = c("seg", "gistic"),
                    digits = 4) {
  segmean <- match.arg(segmean)
  style <- match.arg(style)

  # the legacy logR flag selects the mode unless the caller named one explicitly
  if (segmean == "auto") {
    segmean <- ifelse(logR == TRUE, "logR", "copies")
  }
  if (segmean == "ploidy" & correct == TRUE) {
    stop(
      "correct = TRUE cannot be combined with segmean = \"ploidy\": the purity ",
      "correction rescales cnv.lrr, which the ploidy-relative mode never reads."
    )
  }
  # the purity correction is only meaningful for the caller-reported log ratios
  applyCorrect <- (correct == TRUE & segmean == "logR")

  if (!is.list(cnvObj)) {
    cnvObj <- list(cnvObj)
  }
  names <- sapply(cnvObj, function(X) {
    X@sample
  })
  if (any(duplicated(names)) == TRUE) {
    stop("There are multiple instances of the same sample in this list of CNV objects")
  }
  for (i in 1:length(cnvObj)) {
    X <- cnvObj[[i]]
    stopifnot(class(X) == "cnvstrml")
    if (applyCorrect == TRUE) {
      X <- adjustLogR(X, usePurity = usePurity)
    }
    gr <- X@segments

    seg <- gr %>%
      as.data.frame(row.names = NULL) %>%
      dplyr::mutate(sample = X@sample, .before = 1)

    if (segmean == "ploidy") {
      ploidy <- getPloidy(X, rounded = FALSE)
      if (is.na(ploidy) || ploidy <= 0) {
        logger::log_warn("Sample {X@sample} has an invalid ploidy ({ploidy}). Using 2 instead.")
        ploidy <- 2
      }
    }

    segvals <- switch(segmean,
      copies = seg$total.allele.copies,
      logR = if (applyCorrect == TRUE) seg$cnv.lrrAdj else seg$cnv.lrr,
      ploidy = log2(seg$total.allele.copies / ploidy)
    )

    # floors the -Inf that a homozygous deletion produces in the ploidy mode.
    # guarded so that the default of -Inf leaves the values, and their type,
    # exactly as they were
    if (is.finite(minSegCN)) {
      segvals <- pmax(segvals, minSegCN)
    }

    # the purity correction inverts a log/exp round trip, which leaves residue of the
    # order of 1e-16 on segments that are exactly neutral. rounding last, after the
    # floor, keeps minSegCN exact and writes a plain 0 instead of -4.8e-16.
    # signif() would not help here: signif(3.2e-16, 6) is still 3.2e-16
    if (!is.null(digits) && !is.na(digits)) {
      segvals <- round(segvals, digits)
    }

    seg <- seg %>%
      dplyr::mutate(segvalue = segvals) %>%
      dplyr::select(sample, chrom = seqnames, start, end, probes = cnv.markers, segvalue) %>%
      # keep the coordinates as integers so that write.table cannot emit them in
      # scientific notation (e.g. 5e+05), which downstream parsers reject
      dplyr::mutate(
        start = as.integer(start),
        end = as.integer(end)
      )

    if (dropNA == TRUE) {
      dropped <- sum(is.na(seg$segvalue))
      if (dropped > 0) {
        logger::log_warn("Dropping {dropped} segment(s) of sample {X@sample} with a missing segment mean.")
        seg <- seg %>% dplyr::filter(!is.na(segvalue))
      }
    }
    if (style == "gistic" & any(is.na(seg$probes))) {
      logger::log_warn("Sample {X@sample} has segments without marker counts. GISTIC requires a value in Num_Markers.")
    }

    colnames(seg) <- if (style == "gistic") {
      c("Sample", "Chromosome", "Start", "End", "Num_Markers", "Seg.CN")
    } else {
      c("sample", "chrom", "start", "end", "probes", "segmean")
    }

    tryCatch(
      {
        write.table(seg,
          file = filename,
          quote = FALSE,
          sep = "\t",
          col.names = ifelse(i == 1, TRUE, FALSE),
          row.names = FALSE,
          append = ifelse(i == 1, FALSE, TRUE)
        )
      },
      warning = function(w) {
        message("Warning: ", w)
      },
      error = function(e) {
        message("Error:", e)
      }
    )
  }
  return(filename)
}


#' @rdname saveSEG
#' @export
writeSEG <- saveSEG
