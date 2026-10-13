#' Detect overfragmentation in CNV segments at the chromosomal arm level
#'
#' @description
#' Identifies chromosomal arms with statistically significant overfragmentation based on
#' segment size distribution. Uses a frequentist binomial test to assess whether the number
#' of short segments in an arm is significantly higher than expected under a null fragmentation model.
#' Originally ported from CNAqc; does not require mutation data.
#'
#' @param cnvObj A \code{cnvstrml} object containing CNV segments to test for fragmentation.
#' @param ccf_threshold Numeric clonal cancer cell fraction threshold (0–1). Only segments with
#'   \code{cnv.ccf >= ccf_threshold} are included (default: 0, all segments; 1, clonal only).
#' @param arm_percentage_cutoff Numeric segment size threshold as a fraction of arm length
#'   (0–1). Segments smaller than this proportion are considered "short." If > 1,
#'   interpreted as an absolute size cutoff in base pairs (default: 0.2, i.e., 20%).
#' @param minimum_segments_for_testing Integer minimum number of segments required per
#'   arm for statistical testing (default: 10). Arms with fewer segments are excluded.
#' @param alpha Numeric significance level for multiple-testing adjustment (default: 0.01).
#'   Bonferroni correction applied: significance threshold = \code{alpha / N_tests}.
#'
#' @return The input \code{cnvObj} is updated with a new element \code{@quality$fragmentation},
#'   a list containing:
#' \describe{
#'   \item{table}{A \code{GRanges} of chromosome arms with columns:
#'     \describe{
#'       \item{seqnames, arm}{Chromosome and arm (p or q).}
#'       \item{n_short, n_long}{Count of short and long segments per arm.}
#'       \item{jumps}{Sum of absolute copy-number changes across segments (fragmentation proxy).}
#'       \item{p_value}{Raw binomial test p-value.}
#'       \item{Bonferroni_cutoff}{Adjusted significance threshold (alpha / N_tests).}
#'       \item{significant}{Logical; \code{TRUE} if \code{p_value < Bonferroni_cutoff}.}
#'       \item{testable}{Logical; \code{TRUE} if arm meets minimum segment threshold.}
#'     }}
#'   \item{alpha}{Significance level used.}
#'   \item{N_tests}{Number of arms tested (meeting segment count threshold).}
#'   \item{arm_size_percentage_cutoff}{Arm-relative size threshold, or \code{NA}.}
#'   \item{small_segment_size}{Absolute size cutoff in bp, or \code{NA}.}
#'   \item{ccf_threshold}{CCF threshold applied.}
#'   \item{minimum_segments_for_testing}{Minimum segments per arm required.}
#' }
#'
#' @details
#' \strong{Overfragmentation concept}
#' In CNV analysis, "fragmentation" refers to excessive segmentation—many small, discrete CNV calls
#' over a short genomic region instead of one large call. Biological and technical factors can cause this:
#' \itemize{
#'   \item Over-calling by the detection algorithm.
#'   \item Subclonal CNV events creating complex patterns.
#'   \item Mapping or copy-number estimation artifacts.
#' }
#' This test flags arms with statistically excessive short segments.
#'
#' \strong{Statistical test}
#' Uses a binomial test: under the null hypothesis, short segments appear with probability
#' \code{p_cutoff_short} (derived from \code{arm_percentage_cutoff} or \code{small_segment_size}).
#' The test computes \code{P(X >= n_short | N = n_short + n_long, p = p_cutoff_short)}.
#' Bonferroni correction adjusts for multiple arms tested.
#'
#' \strong{Preprocessing}
#' \itemize{
#'   \item Requires arm annotation; \code{addChromosomeAnnotation()} is called internally if missing.
#'   \item Segments are merged across large distances (1 Mbp) regardless of logR/CCF agreement to avoid
#'     inflating short-segment counts from spurious boundaries.
#'   \item Sex chromosomes (chrX, chrY) may be included; filter externally if needed.
#' }
#'
#' \strong{Arm eligibility for testing}
#' An arm is testable if \code{n_short + n_long >= minimum_segments_for_testing}.
#' Non-testable arms still appear in the output table but with \code{testable = FALSE}
#' and \code{p_value = 1.0}.
#'
#' @examples
#' \dontrun{
#' # Detect overfragmentation with default thresholds
#' cnv <- detectOverfragmentation(
#'   cnvObj,
#'   ccf_threshold = 0,
#'   arm_percentage_cutoff = 0.2,
#'   minimum_segments_for_testing = 10,
#'   alpha = 0.01
#' )
#'
#' # Inspect results
#' frag_results <- cnv@quality$fragmentation$table
#' significant_arms <- frag_results[frag_results$significant]
#' cat("Significant fragmentation:", length(significant_arms), "arms\n")
#'
#' # Tighter test: only clonal segments, 10% size cutoff
#' cnv2 <- detectOverfragmentation(
#'   cnvObj,
#'   ccf_threshold = 0.8,
#'   arm_percentage_cutoff = 0.1,
#'   minimum_segments_for_testing = 10
#' )
#' }
#'
#' @seealso \code{\link{addChromosomeAnnotation}}, \code{\link{frequentist_test_fragmentation}},
#'   \code{\link{cnvEntropy}} for complementary fragmentation assessment
#'
#' @export
detectOverfragmentation <- function(cnvObj,
                                    ccf_threshold = 0,
                                    arm_percentage_cutoff = 0.2,
                                    minimum_segments_for_testing = 10,
                                    alpha = 0.01) {
  if (arm_percentage_cutoff < 0) {
    logger::log_warn("Cannot use negative'arm_percentage_cutoff' {arm_percentage_cutoff}", namespace = "cnvstrmliner")
    stop()
  }

  if ((arm_percentage_cutoff) > 1) {
    small_segment_size <- arm_percentage_cutoff
    arm_percentage_cutoff <- NA
  } else {
    small_segment_size <- NA
  }
  # add the chromosome arm annotation
  cnvObj <- cnvstrmliner::addChromosomeAnnotation(cnvObj)
  # merge the segments
  segs <- cnvstrmliner::mergeSegments(cnvObj@segments, minDistance = 1e6, minccfDiffPct = 1000, minlrrDiffPct = 1000)

  # we will use all the segments that have CCF above the threshold
  clonal_cna <- subset(segs, "cnv.ccf" >= ccf_threshold)
  # the size of each arm
  S4Vectors::mcols(chrom_arms[[cnvObj@genome]])["size"] <- BiocGenerics::width(chrom_arms[[cnvObj@genome]])
  # calculate the length, percentage of arm of each segment
  ov <- GenomicRanges::findOverlaps(query = clonal_cna, subject = chrom_arms[[cnvObj@genome]])
  S4Vectors::mcols(clonal_cna)["L"] <- 0
  S4Vectors::mcols(clonal_cna)["perc_length"] <- 0
  S4Vectors::mcols(clonal_cna)["smaller"] <- FALSE
  S4Vectors::mcols(clonal_cna)[S4Vectors::queryHits(ov), "L"] <- S4Vectors::mcols(chrom_arms[[cnvObj@genome]])[S4Vectors::subjectHits(ov), "size"]
  S4Vectors::mcols(clonal_cna)["perc_length"] <- BiocGenerics::width(clonal_cna) / S4Vectors::mcols(clonal_cna)[, "L"]
  if (is.na(small_segment_size)) {
    S4Vectors::mcols(clonal_cna)["smaller"] <- S4Vectors::mcols(clonal_cna)[, "perc_length"] <= arm_percentage_cutoff
  }
  if (is.na(arm_percentage_cutoff)) {
    S4Vectors::mcols(clonal_cna)["smaller"] <- S4Vectors::width(clonal_cna) <= small_segment_size
  }

  # count the number of smaller/larger segments
  counts <- clonal_cna %>%
    as.data.frame() %>%
    dplyr::filter(arm %in% c("p", "q")) %>%
    dplyr::group_by(seqnames, arm, smaller, L) %>%
    dplyr::summarise(n_short = dplyr::n(), .groups = "keep") %>%
    dplyr::ungroup()

  if (all(counts$smaller == FALSE)) {
    logger::log_warn("No short segments with these parameters.")
    return(cnvObj)
  }
  if (all(counts$smaller == TRUE)) {
    counts <- counts %>%
      tidyr::spread(smaller, n_short) %>%
      dplyr::rename(n_short = `TRUE`)
    counts$n_long <- 0
  } else {
    counts <- counts %>%
      tidyr::spread(smaller, n_short) %>%
      dplyr::rename(n_short = `TRUE`, n_long = `FALSE`)
  }
  counts$n_long[is.na(counts$n_long)] <- 0
  counts$n_short[is.na(counts$n_short)] <- 0

  jumps <- clonal_cna %>%
    as.data.frame() %>%
    dplyr::filter(arm %in% c("p", "q")) %>%
    split(., f = paste(.$seqnames, .$arm)) %>%
    lapply(., function(x) {
      if (nrow(x) <= 1) {
        return(list(seqnames = as.character(x$seqnames[1]), arm = as.character(x$arm[1]), jumps = 0))
      }
      y <- x %>%
        dplyr::pull(total.allele.copies) %>%
        dplyr::lag(., default = .[1])
      return(list(seqnames = as.character(x$seqnames[1]), arm = as.character(x$arm[1]), jumps = sum(abs(x$total.allele.copies - y))))
    }) %>%
    dplyr::bind_rows()
  counts <- dplyr::full_join(counts, jumps, by = dplyr::join_by(seqnames, arm))
  # find which arms to be tested
  testable <- counts %>%
    dplyr::mutate(total_segments = n_short + n_long) %>%
    dplyr::mutate(testable = total_segments >= minimum_segments_for_testing)
  N_tests <- sum(testable$testable)

  testable <- testable %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      size_cutoff = ifelse(is.na(small_segment_size), arm_percentage_cutoff, min(1, small_segment_size / L)),
      p_value = frequentist_test_fragmentation(n_short,
        n_long,
        testable,
        p_cutoff_short = size_cutoff,
        N_tests = N_tests,
        alpha = alpha
      ),
      Bonferroni_cutoff = alpha / N_tests,
      significant = ifelse(N_tests == 0, FALSE, p_value < (alpha / N_tests))
    ) %>%
    dplyr::arrange(p_value) %>%
    dplyr::ungroup()
  # N_sign = sum(testable$significant)


  chrom_fragmentation <- chrom_arms[[cnvObj@genome]] %>%
    as.data.frame() %>%
    dplyr::left_join(., testable, by = dplyr::join_by("seqnames", "arm")) %>%
    GenomicRanges::makeGRangesFromDataFrame(, keep.extra.columns = TRUE)


  cnvObj@quality[["fragmentation"]] <- list(
    table = chrom_fragmentation, alpha = alpha,
    N_tests = N_tests, arm_size_percentage_cutoff = arm_percentage_cutoff, small_segment_size = small_segment_size,
    ccf_threshold = ccf_threshold,
    minimum_segments_for_testing = minimum_segments_for_testing
  )
  cnvObj
}


#' Calculate the p-value for events on an arm chromosome
#' @details
#' This is a copy of the function with the same name found in CNAqc (Thank you !!!)
#'
#' @param n_short (integer) number of short segments
#' @param n_long (integer) number of long segments
#' @param testable (logical) a value indicating if the arm is testable (passing the thresholds for number of segments etc)
#' @param p_cutoff_short (float) the pvalue cutoff used in the dbinom
#' @param N_tests (integer) number of tests performed in the sample
#' @param alpha (float) threshold for adjusted pvalue
frequentist_test_fragmentation <- function(n_short, n_long, testable, p_cutoff_short = 0.2,
                                           N_tests = 1, alpha = 0.01) {
  if (!testable) {
    return(1)
  }
  N <- n_long + n_short
  p_value <- sum(sapply(n_short:N, function(k) {
    dbinom(k, N,
      prob = p_cutoff_short
    )
  }))
  significant <- p_value < (alpha / N_tests)

  return(p_value)
}
