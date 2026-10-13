# Metrics comparing two segmentation solutions genome-wide, plus their rollups.

#' Calculate pairwise similarity metrics between CNV profiles
#'
#' @description
#' This function computes comprehensive similarity metrics that quantify the concordance
#' between copy number profiles across multiple samples. It performs pairwise comparisons
#' of CNV patterns to assess how similar different samples are in terms of their genomic
#' aberrations. The function leverages the \code{CNVMetrics} package to calculate
#' overlap-based similarity measures for each CNV state (amplifications, deletions, neutral, etc.).
#'
#' **What are CNV Similarity Metrics?**
#'
#' CNV similarity metrics quantify the degree of concordance between two copy number profiles
#' by measuring how much their aberrant genomic regions overlap. These metrics are essential for:
#' \itemize{
#'   \item Comparing CNV calls from different algorithms on the same sample
#'   \item Assessing reproducibility across technical or biological replicates
#'   \item Identifying samples with similar aberration patterns in a cohort
#'   \item Validating CNV calling methods against gold standards
#'   \item Clustering samples based on CNV profiles
#' }
#'
#' **The Sørensen-Dice Coefficient:**
#'
#' This function uses the Sørensen-Dice coefficient (also known as Dice similarity or
#' F1 score in machine learning) to measure overlap between CNV segments. For two samples
#' A and B, the coefficient is calculated as:
#'
#' \deqn{S = \frac{2 \times |A \cap B|}{|A| + |B|}}
#'
#' where:
#' \itemize{
#'   \item \eqn{|A \cap B|} is the size (in base pairs) of genomic regions where both
#'         samples share the same CNV state
#'   \item \eqn{|A|} and \eqn{|B|} are the total sizes of aberrant regions in each sample
#' }
#'
#' **Key Properties of the Sørensen-Dice Coefficient:**
#' \itemize{
#'   \item **Range:** Values range from 0 to 1
#'   \item **Symmetric:** S(A,B) = S(B,A)
#'   \item **Harmonic mean:** It's the harmonic mean of precision and recall, giving
#'         balanced weight to both shared and unique segments
#'   \item **State-specific:** Calculated separately for each CNV state (AMP, DEL, NEUT, etc.)
#' }
#'
#' **Interpreting Similarity Values:**
#'
#' \itemize{
#'   \item **Similarity = 1.0 (Perfect match):** The two samples have identical CNV patterns
#'         for the given state. This indicates:
#'         \itemize{
#'           \item Perfect concordance between methods (ideal for validation)
#'           \item True biological replicates or clonal samples
#'           \item Potentially identical input data
#'         }
#'   \item **Similarity = 0.8-1.0 (High similarity):** Strong concordance with minor differences.
#'         This typically represents:
#'         \itemize{
#'           \item Good agreement between different calling methods
#'           \item Technical replicates with expected variation
#'           \item Samples from the same tumor with minor heterogeneity
#'         }
#'   \item **Similarity = 0.5-0.8 (Moderate similarity):** Partial overlap with substantial
#'         differences. This can indicate:
#'         \itemize{
#'           \item Different methods with systematic differences in sensitivity/specificity
#'           \item Related samples with some biological divergence
#'           \item Focal vs. broad CNV patterns in similar regions
#'         }
#'   \item **Similarity = 0.2-0.5 (Low similarity):** Limited overlap, suggesting:
#'         \itemize{
#'           \item Fundamentally different CNV calling approaches
#'           \item Biologically distinct samples
#'           \item Different tumor regions with high heterogeneity
#'           \item Potential quality issues in one or both samples
#'         }
#'   \item **Similarity = 0.0-0.2 (Very low/no similarity):** Minimal or no overlap,
#'         indicating:
#'         \itemize{
#'           \item Completely different biological samples
#'           \item Major quality issues or failed CNV calling
#'           \item Incompatible CNV calling parameters or methods
#'         }
#' }
#'
#' **Using Metrics in CNV Analysis:**
#'
#' \enumerate{
#'   \item **Method Comparison:** When comparing multiple CNV callers on the same sample,
#'         high similarity (>0.8) suggests robust, reliable calls. Low similarity may
#'         indicate method-specific biases or parameter tuning needs.
#'
#'   \item **Quality Control:** Comparing tumor-normal pairs or replicates can identify
#'         problematic samples. Expected high similarity that's absent suggests technical
#'         issues.
#'
#'   \item **Cohort Analysis:** Similarity matrices can reveal subgroups with shared CNV
#'         patterns, useful for cancer subtype classification or identifying convergent
#'         evolutionary paths.
#'
#'   \item **State-Specific Insights:** Analyzing metrics separately for AMP, DEL, and other
#'         states can reveal whether methods agree better on certain aberration types
#'         (e.g., good concordance on large deletions but poor on focal amplifications).
#'
#'   \item **Batch Effect Detection:** Unexpected clustering by batch in similarity
#'         matrices may indicate technical artifacts rather than biological signal.
#' }
#'
#' Running this function on a set of CNV profiles returns a table like the following:
#'   state          sample Sample_0110 Sample_0120
#'   AMP Sample_0110      0.00000000      0.86432425
#'   AMP Sample_0120      0.86432425      0.00000000
#'  NEUT Sample_0110      0.00000000      0.98330594
#'  NEUT Sample_0120      0.98330594      0.00000000
#'   DEL Sample_0110      0.00000000      0.00264313
#'   DEL Sample_0120      0.00264313      0.00000000
#' **The Simplify Parameter:**
#'
#' When \code{simplify = TRUE}, the function collapses related CNV states to reduce
#' complexity and focus on major aberration types:
#' \itemize{
#'   \item AMP-LOH → AMP (amplification with LOH treated as amplification)
#'   \item HETDEL + HOMDEL → DEL (heterozygous and homozygous deletions merged)
#'   \item NEUT-LOH → NEUT (copy-neutral LOH treated as neutral)
#'   \item INFER → NEUT (inferred states treated as neutral)
#' }
#'
#' This simplification is useful when the distinction between subtypes (e.g., HETDEL vs.
#' HOMDEL) is less important than the overall pattern of gains vs. losses. Use
#' \code{simplify = FALSE} to preserve all state distinctions for fine-grained comparisons.
#'
#' Running this function with the flag simplify set to FALSE
#' on a set of CNV profiles returns a table like the following:
#'
#'    state          sample SAMP_XX_01 SAMP_XX_01 - 90 SAMP_XX_01 - 30
#'      AMP      SAMP_XX_01  0.0000000       0.7295007       0.8443465
#'      AMP SAMP_XX_01 - 90  0.7295007       0.0000000       0.8105727
#'      AMP SAMP_XX_01 - 30  0.8443465       0.8105727       0.0000000
#'     NEUT      SAMP_XX_01  0.0000000       0.9509234       0.9681023
#'     NEUT SAMP_XX_01 - 90  0.9509234       0.0000000       0.9699873
#'     NEUT SAMP_XX_01 - 30  0.9681023       0.9699873       0.0000000
#'   HETDEL      SAMP_XX_01  0.0000000       0.9149003              NA
#'   HETDEL SAMP_XX_01 - 90  0.9149003       0.0000000              NA
#'   HETDEL SAMP_XX_01 - 30         NA              NA       0.0000000
#'  AMP-LOH      SAMP_XX_01  0.0000000       0.9845454       0.9829729
#'  AMP-LOH SAMP_XX_01 - 90  0.9845454       0.0000000       0.9898567
#'  AMP-LOH SAMP_XX_01 - 30  0.9829729       0.9898567       0.0000000
#' NEUT-LOH      SAMP_XX_01  0.0000000       0.1447130       0.9381414
#' NEUT-LOH SAMP_XX_01 - 90  0.1447130       0.0000000       0.1345986
#' NEUT-LOH SAMP_XX_01 - 30  0.9381414       0.1345986       0.0000000
#'   HOMDEL      SAMP_XX_01  0.0000000              NA              NA
#'   HOMDEL SAMP_XX_01 - 90         NA       0.0000000       0.1400067
#'   HOMDEL SAMP_XX_01 - 30         NA       0.1400067       0.0000000
#' @param cnvList A list of \code{cnvstrml} objects to compare. Each object should contain
#'   CNV segmentation data with \code{cnv.annotation} and \code{cnv.lrr} columns. The
#'   objects will be compared in a pairwise all-vs-all manner.
#'
#' @param simplify Logical. If \code{TRUE} (default), collapses CNV states into major
#'   categories: all amplifications with LOH are labeled as AMP, all deletions (HETDEL,
#'   HOMDEL) are labeled as DEL, and neutral regions with LOH (NEUT-LOH) are labeled as
#'   NEUT. If \code{FALSE}, preserves all original state distinctions. Set to \code{FALSE}
#'   when LOH information is biologically important.
#'
#' @return A data frame with the following structure:
#'   \describe{
#'     \item{state}{The CNV state for which the metric was calculated (e.g., "AMP", "DEL",
#'       "NEUT")}
#'     \item{sample}{The sample name (row identifier in the similarity matrix)}
#'     \item{sample columns}{One column per sample containing the Sørensen-Dice
#'       similarity coefficient between the row sample and that column sample}
#'   }
#'   The result is a long-format data frame where each row represents one sample, and
#'   columns represent pairwise similarities for a specific CNV state. Values range from
#'   0 (no overlap) to 1 (perfect overlap).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Compare CNV profiles from three different callers on the same sample
#' cnv_list <- list(
#'   facets = parseFACETS(facets_file),
#'   battenberg = parseBattenberg(battenberg_file),
#'   cnvkit = parseCNVkit(cnvkit_file)
#' )
#'
#' # Calculate similarity metrics with simplified states
#' metrics <- cnvMetrics(cnv_list, simplify = TRUE)
#'
#' # View amplification similarities
#' amp_metrics <- subset(metrics, state == "AMP")
#'
#' # Identify pairs with high concordance (>0.8)
#' # This helps validate that methods agree on major aberrations
#'
#' # Compare with detailed states preserved
#' detailed_metrics <- cnvMetrics(cnv_list, simplify = FALSE)
#'
#' # Check if methods agree on LOH calls
#' loh_metrics <- subset(detailed_metrics, state == "NEUT-LOH")
#' }
#'
#' @seealso
#' \code{\link{cnvEntropy}} for complexity/fragmentation metrics
#' \code{\link[CNVMetrics]{calculateOverlapMetric}} for the underlying calculation
#' \code{\link{getJaccardIndex}} for calculating the Jaccard index between CNV profiles
#'
cnvMetrics <- function(cnvList, simplify = TRUE) {
  # Validate cnvList
  if (!is.list(cnvList)) {
    stop("cnvList must be a list of cnvstrml objects")
  }

  # Check if all elements in the list are cnvstrml objects
  if (!all(sapply(cnvList, function(x) inherits(x, "cnvstrml")))) {
    stop("All elements in cnvList must be cnvstrml objects")
  }

  # CNV metrics expect a column with the name "state"
  tt <- lapply(cnvList, function(X) {
    a <- X@segments[, c("cnv.annotation", "cnv.lrr")]
    S4Vectors::mcols(a) <- S4Vectors::mcols(a) %>%
      as.data.frame() %>%
      dplyr::mutate(cnv.annotation = ifelse(cnv.annotation == "INFER", "NEUT", cnv.annotation))
    if (simplify == TRUE) {
      S4Vectors::mcols(a) <- S4Vectors::mcols(a) %>%
        as.data.frame() %>%
        dplyr::mutate(
          cnv.annotation = ifelse(cnv.annotation == "AMP-LOH", "AMP", cnv.annotation),
          cnv.annotation = ifelse(cnv.annotation == "NEUT-LOH", "NEUT", cnv.annotation),
          cnv.annotation = ifelse(cnv.annotation == "HETDEL", "DEL", cnv.annotation),
          cnv.annotation = ifelse(cnv.annotation == "HOMDEL", "DEL", cnv.annotation)
        )
    }
    colnames(S4Vectors::mcols(a)) <- c("state", "log2ratio")
    a
  })

  tt <- GenomicRanges::GRangesList(tt)
  if (!is.null(names(cnvList))) {
    names(tt) <- names(cnvList)
  } else {
    names(tt) <- sapply(cnvList, function(X) {
      X@sample
    })
  }
  states <- lapply(tt, function(X) {
    unique(X$state)
  }) %>%
    unlist() %>%
    unique()

  metrics <- CNVMetrics::calculateOverlapMetric(segmentData = tt, method = "sorensen", nJobs = 1, states = states)

  lapply(metrics, function(X) {
    X %>%
      as.dist() %>%
      as.matrix() %>%
      as_tibble(rownames = "sample")
  }) %>%
    dplyr::bind_rows(., .id = "state") %>%
    as.data.frame()
}


#' Calculate Jaccard Index (and Sorensen - Dice ) for Altered Genome Regions
#'
#' @description
#' Computes the Jaccard index (also known as Intersection over Union) to measure
#' the similarity between truth and predicted copy number altered regions at the
#' base pair level. The Jaccard index quantifies the overlap between two sets as
#' the ratio of their intersection to their union.
#' For the purposes of this function one of the cnv objects is considered the "truth"
#' (subject) and the other the "prediction" (query).
#'
#' @param comparison A comparison object (typically from \code{cnvstrmliner::compareSegmentStatus})
#'   containing a \code{table} component. The table must include the following columns:
#'   \itemize{
#'     \item \code{subject}: Copy number state in the truth/reference sample
#'     (e.g., "NEUT" for neutral, or altered states)
#'     \item \code{query}: Copy number state in the predicted sample
#'     \item \code{start}, \code{end}: Genomic coordinates of each segment
#'     \item \code{status}: Classification of the segment (e.g., "TP" for true positive)
#'   }
#'
#' @return A named list with the following numeric elements:
#'   \item{jaccard}{The Jaccard index, calculated as
#'     TP / (Truth + Predicted - TP), where TP is true positive base pairs,
#'     Truth is total altered base pairs in the truth set, and Predicted is
#'     total altered base pairs in the predicted set. Values range from 0
#'     (no overlap) to 1 (perfect agreement).}
#'  \item{sorensen}{The Sorensen-Dice coefficient, calculated as
#'   2 * TP / (Truth + Predicted). Values range from 0 (no overlap) to 1
#'    (perfect agreement).}
#'   \item{subject_altered_genome_bp}{Total number of base pairs altered in the
#'     subject/reference (subject != "NEUT").}
#'   \item{query_altered_genome_bp}{Total number of base pairs altered in
#'     the query (query != "NEUT").}
#'   \item{tp_altered_genome_bp}{Total number of true positive base pairs
#'     (both query and subject are non-neutral AND status == "TP").}
#'   \item{matched_query_altered_genome}{The fraction of predicted altered base pairs that are true positives, calculated as
#'     TP / Predicted.}
#'   \item{matched_subject_altered_genome}{The fraction of truth altered base pairs that are true positives, calculated as
#'     TP / Truth.}
#'
#' @details
#' The Jaccard index is a standard metric for comparing set similarity, defined as:
#' \deqn{J(A, B) = \frac{|A \cap B|}{|A \cup B|} = \frac{|A \cap B|}{|A| + |B| - |A \cap B|}}
#' The Sorensen-Dice coefficient is a related metric, defined as:
#' \deqn{S(A, B) = \frac{2|A \cap B|}{|A| + |B|}}
#'
#' In the context of CNV analysis:
#' \itemize{
#'   \item Neutral regions (\code{NEUT}) are excluded from both truth and prediction sets
#'   \item Only altered genome regions contribute to the calculation
#'   \item The metric is computed at base pair resolution, not segment count
#'   \item Higher values indicate better agreement between truth and prediction
#' }
#'
#' @examples
#' \dontrun{
#' # Compare truth and predicted CNV segmentations
#' comparison <- cnvstrmliner::compareSegmentStatus(truth_cnv, predicted_cnv)
#' jaccard_results <- getJaccardIndex(comparison)
#'
#' # Extract the Jaccard index
#' print(jaccard_results$jaccard)
#' }
#'
#' @seealso \code{\link{compareSegmentStatus}} for generating comparison objects
#' @seealso \code{\link{cnvMetrics}} for other CNV comparison metrics
#' @export
getJaccardIndex <- function(comparison) {
  # altered genome bp (subject=truth)
  subject_altered_genome_bp <- comparison$table %>%
    dplyr::filter(subject != "NEUT") %>%
    dplyr::mutate(width = end - start + 1) %>%
    dplyr::summarise(total_width = sum(width)) %>%
    dplyr::pull(total_width)

  query_altered_genome_bp <- comparison$table %>%
    dplyr::filter(query != "NEUT") %>%
    dplyr::mutate(width = end - start + 1) %>%
    dplyr::summarise(total_width = sum(width)) %>%
    dplyr::pull(total_width)

  tp_altered_genome_bp <- comparison$table %>%
    dplyr::filter(query != "NEUT" & subject != "NEUT" & status == "TP") %>%
    dplyr::mutate(width = end - start + 1) %>%
    dplyr::summarise(total_width = sum(width)) %>%
    dplyr::pull(total_width)

  jaccard_altered_genome_bp <- tp_altered_genome_bp / (subject_altered_genome_bp + query_altered_genome_bp - tp_altered_genome_bp)
  sorensen_altered_genome_bp <- 2 * tp_altered_genome_bp / (subject_altered_genome_bp + query_altered_genome_bp)
  return(list(
    jaccard = jaccard_altered_genome_bp,
    sorensen = sorensen_altered_genome_bp,
    subject_altered_genome_bp = subject_altered_genome_bp,
    query_altered_genome_bp = query_altered_genome_bp,
    tp_altered_genome_bp = tp_altered_genome_bp,
    matched_query_altered_genome = tp_altered_genome_bp / query_altered_genome_bp,
    matched_subject_altered_genome = tp_altered_genome_bp / subject_altered_genome_bp
  ))
}


# ===== 1. GENOME-WIDE TRUTH-VS-PREDICTION ====================================
#
# Scalar scores for one segmentation against a reference, intended for the
# simulation cohort where the reference is the known truth. They are equally
# well defined against any other reference, and nothing in them assumes the
# subject is correct.
#
# The organising idea is that raw agreement is not a quality metric. In this
# cohort the truly altered genome fraction ranges from about 0.07 to 0.83, so a
# caller that emits nothing scores a raw agreement of 0.63 on average and up to
# 0.93 -- a number that looks like skill and is not. Kappa removes that
# confound by construction; scoreNullCaller() is the regression test that keeps
# it honest; computeEventRecovery() separates the two error directions kappa
# merges into one number.
# =============================================================================


#' Chance-corrected agreement between two segmentations at bin level
#'
#' @description
#' Computes Cohen's kappa between the copy-number state a query segmentation
#' assigns to each bin and the state a reference (subject) segmentation assigns
#' to the same bin.
#'
#' Kappa is the achieved-quality metric of choice for this comparison because it
#' is *null referenced*: a caller that labels the entire genome `NEUT` scores
#' exactly 0 no matter how much of the reference is genuinely neutral. Raw bin
#' agreement does not have that property, and the size of the resulting
#' distortion is not small -- see [scoreNullCaller()], which makes it concrete
#' for a given sample.
#'
#' @section How it is computed:
#' Each bin is assigned to the single segment it overlaps most, on both sides
#' independently. The two resulting state vectors are cross-tabulated over the
#' union of their levels, and
#'
#' \deqn{\kappa = \frac{p_o - p_e}{1 - p_e}}
#'
#' where \eqn{p_o} is the weighted fraction of bins whose states match and
#' \eqn{p_e = \sum_s p_{query}(s)\,p_{subject}(s)} is the agreement expected if
#' the two label sets were assigned independently with their observed marginals.
#'
#' With `weight = "bins"` every bin counts once, which is the right choice for a
#' capture panel where each target is one unit of evidence. With
#' `weight = "width"` bins contribute in proportion to their span, which
#' reproduces a base-pair metric and is appropriate only for a uniform tiling.
#'
#' @section Choosing the bins:
#' `bins` is required and has no default on purpose. On whole-exome data a
#' uniform base-pair tiling is not a discretisation but a fabrication: most of
#' each tile is unsequenced, no caller had evidence there, and the metric
#' degenerates into a base-pair agreement measured largely over intergenic
#' space. Pass the capture target regions. On whole-genome data a
#' [GenomicRanges::tileGenome()] grid is appropriate.
#'
#' The bins also act as the region filter. To exclude the sex chromosomes --
#' which otherwise dominate the score on samples with no autosomal events --
#' pass autosome-only bins. There is deliberately no separate `drop_seqnames`
#' argument: one mechanism, not two.
#'
#' @section NA handling:
#' Bins with no assigned segment on one or both sides are dropped from the
#' denominator and reported as `n_bins_unscored`. They are never imputed to
#' `NEUT`: a caller that emitted no segment over a region has not asserted that
#' the region is neutral, and treating absence as a neutral call is precisely
#' how an empty segmentation earns a good score. Callers that genuinely want
#' gap filling should apply [cnvstrmliner::fillGaps()] first, which makes the
#' assumption explicit and visible in the object.
#'
#' A large `n_bins_unscored` relative to `n_bins_scored` means the metric is
#' answering a narrower question than it appears to, and is worth checking
#' before the score is believed.
#'
#' Empty input or zero scorable bins yields `NA_real_` with a warning. A
#' degenerate expected agreement (\eqn{p_e = 1}, reached when both sides are
#' constant and equal) yields `kappa = 0` by convention. Wrong argument types, a
#' missing state column, or an unreadable bin file are contract violations and
#' raise an error.
#'
#' @section Interpretation:
#' Landis and Koch bands, which are conventional rather than principled, and
#' should be read against [scoreNullCaller()] on the same sample:
#'
#' - **> 0.80** -- almost perfect; the caller reproduces the reference states
#' - **0.60 - 0.80** -- substantial
#' - **0.40 - 0.60** -- moderate
#' - **0.20 - 0.40** -- fair
#' - **0.00 - 0.20** -- slight; hard to distinguish from an all-neutral caller
#' - **< 0** -- worse than chance, which in practice almost always means a
#'   purity or ploidy solution has flipped rather than that segmentation failed
#'
#' @param query A `cnvstrml` object or a `GRanges` of segments -- the
#'   segmentation being scored.
#' @param subject A `cnvstrml` object or a `GRanges` of segments -- the
#'   reference, e.g. the simulation truth.
#' @param bins A `GRanges` of bin-level regions, or a character path to a BED
#'   file, which is imported with [rtracklayer::import.bed()] if that package is
#'   available. Required; see the *Choosing the bins* section.
#' @param statecol Character. Name of the `mcols` column holding the state in
#'   both `query` and `subject`. Default `"cnv.annotation"`. Set to
#'   `"total.allele.copies"` to score integer copy number instead of the
#'   annotation label.
#' @param simplify Logical. Collapse `AMP-LOH` to `AMP`, `NEUT-LOH` to `NEUT`
#'   and both `HETDEL` and `HOMDEL` to `DEL`, matching
#'   [cnvstrmliner::cnvMetrics()]. Default `FALSE`, which is deliberately the
#'   opposite of that function's default: kappa is already chance corrected and
#'   does not need the collapse to stabilise, and merging `HOMDEL` into `DEL`
#'   discards the distinction that matters most. Set `TRUE` to compare like for
#'   like with `cnvMetrics` output.
#' @param infer_as_neutral Logical. Treat the `INFER` annotation as `NEUT`.
#'   Default `TRUE`, matching `cnvMetrics`, which does it unconditionally.
#'   `INFER` marks regions with no evidence, typically produced by
#'   [cnvstrmliner::fillGaps()]; keeping it as a distinct state penalises honest
#'   gap filling.
#' @param weight One of `"bins"` (default, every bin counts once) or `"width"`
#'   (bins weighted by base pairs).
#'
#' @return A named `list`:
#' \describe{
#'   \item{kappa}{Cohen's kappa. `NA_real_` if no bin was scorable.}
#'   \item{p_o}{Weighted observed agreement.}
#'   \item{p_e}{Weighted agreement expected by chance from the marginals.}
#'   \item{n_bins_scored}{Bins with an assigned segment on both sides.}
#'   \item{n_bins_unscored}{Bins dropped for lack of an assignment on at least
#'     one side.}
#'   \item{weight}{The weighting actually used, echoed back.}
#'   \item{marginals}{A `data.frame` with columns `state`, `p_query` and
#'     `p_subject` giving the weighted marginal distribution on each side. These
#'     drive \eqn{p_e}, and inspecting them is the fastest way to see that a low
#'     kappa is a composition artefact rather than a segmentation failure.}
#' }
#'
#' @seealso
#' [scoreNullCaller()] for the all-neutral reference score kappa must beat;
#' [computeStateConfusion()] for the table kappa is read off;
#' [computeEventRecovery()] to separate the two error directions kappa merges;
#' [computeSegmentationARI()] for agreement on segment *identity* rather than
#' state; [cnvstrmliner::cnvMetrics()] for the per-state Sorensen-Dice
#' alternative.
#'
#' @examples
#' \dontrun{
#' truth <- readRDS("data/cache/truth_profiles/BRCA_XX_01.rds")
#' targets <- rtracklayer::import.bed("data/S07604514_Regions_v6_hg38_Targets.bed")
#'
#' # A segmentation against itself is the ceiling.
#' computeStateKappa(truth, truth, bins = targets)$kappa # 1
#'
#' # And the floor, for context.
#' scoreNullCaller(truth, bins = targets)$null_kappa # 0
#' }
#'
#' @export
computeStateKappa <- function(query,
                              subject,
                              bins,
                              statecol = "cnv.annotation",
                              simplify = FALSE,
                              infer_as_neutral = TRUE,
                              weight = c("bins", "width")) {
  weight <- match.arg(weight)
  query_gr <- .asSegments(query, "query", statecol)
  subject_gr <- .asSegments(subject, "subject", statecol)
  bins <- .resolveBins(bins, allow_auto = FALSE)

  empty_result <- list(
    kappa = NA_real_, p_o = NA_real_, p_e = NA_real_,
    n_bins_scored = 0L, n_bins_unscored = length(bins), weight = weight,
    marginals = data.frame(
      state = character(), p_query = numeric(),
      p_subject = numeric(), stringsAsFactors = FALSE
    )
  )

  if (length(bins) == 0 || length(query_gr) == 0 || length(subject_gr) == 0) {
    logger::log_warn("computeStateKappa: empty input, returning NA")
    return(empty_result)
  }

  query_state <- .binStateVector(bins, query_gr, statecol, simplify, infer_as_neutral)
  subject_state <- .binStateVector(bins, subject_gr, statecol, simplify, infer_as_neutral)

  bin_weight <- if (weight == "width") as.numeric(BiocGenerics::width(bins)) else NULL
  ct <- .contingency(query_state, subject_state, w = bin_weight, common_levels = TRUE)

  if (ct$n_scored == 0) {
    logger::log_warn("computeStateKappa: no bin is covered by both segmentations, returning NA")
    return(empty_result)
  }

  k <- .kappaFromTable(ct$tab)
  marginals <- k$marginals
  names(marginals) <- c("state", "p_query", "p_subject")

  list(
    kappa = k$kappa,
    p_o = k$p_o,
    p_e = k$p_e,
    n_bins_scored = ct$n_scored,
    n_bins_unscored = ct$n_dropped,
    weight = weight,
    marginals = marginals
  )
}


#' Bin-level confusion matrix between two segmentations
#'
#' @description
#' The contingency table [computeStateKappa()] is read off, plus the per-state
#' recall, precision and F1 that fall out of it.
#'
#' Kappa answers "how much better than chance", in one number. That number
#' cannot say *which* states were confused for which, and a single scalar hides
#' the difference between a caller that misses amplifications and one that
#' invents deletions. This function is exposed so that downstream code inspects
#' the table rather than rebuilding it by hand from segment coordinates, which
#' is where the bin-assignment rules quietly diverge.
#'
#' @inheritSection computeStateKappa Choosing the bins
#' @inheritSection computeStateKappa NA handling
#'
#' @inheritParams computeStateKappa
#'
#' @return A named `list`:
#' \describe{
#'   \item{counts}{A `data.frame` with one row per (`subject_state`,
#'     `query_state`) pair, and columns `n` (bins) and `w` (weight, equal to `n`
#'     when `weight = "bins"`). The full cross-product is present, including
#'     zero cells, so the shape does not depend on which confusions happened to
#'     occur.}
#'   \item{per_state}{A `data.frame` with one row per state and columns
#'     `recall` (of the bins the subject calls this state, the fraction the
#'     query agrees on), `precision` (the converse), `f1`, `support_subject` and
#'     `support_query`.}
#' }
#'
#' @seealso [computeStateKappa()], [computeEventRecovery()]
#'
#' @examples
#' \dontrun{
#' cm <- computeStateConfusion(facets_obj, truth, bins = targets)
#' subset(cm$per_state, support_subject > 0)
#' }
#'
#' @export
computeStateConfusion <- function(query,
                                  subject,
                                  bins,
                                  statecol = "cnv.annotation",
                                  simplify = FALSE,
                                  infer_as_neutral = TRUE,
                                  weight = c("bins", "width")) {
  weight <- match.arg(weight)
  query_gr <- .asSegments(query, "query", statecol)
  subject_gr <- .asSegments(subject, "subject", statecol)
  bins <- .resolveBins(bins, allow_auto = FALSE)

  empty_result <- list(
    counts = data.frame(
      subject_state = character(), query_state = character(),
      n = numeric(), w = numeric(), stringsAsFactors = FALSE
    ),
    per_state = data.frame(
      state = character(), recall = numeric(),
      precision = numeric(), f1 = numeric(),
      support_subject = numeric(), support_query = numeric(),
      stringsAsFactors = FALSE
    )
  )

  if (length(bins) == 0 || length(query_gr) == 0 || length(subject_gr) == 0) {
    logger::log_warn("computeStateConfusion: empty input, returning empty tables")
    return(empty_result)
  }

  query_state <- .binStateVector(bins, query_gr, statecol, simplify, infer_as_neutral)
  subject_state <- .binStateVector(bins, subject_gr, statecol, simplify, infer_as_neutral)

  bin_weight <- if (weight == "width") as.numeric(BiocGenerics::width(bins)) else NULL
  ct_n <- .contingency(query_state, subject_state, w = NULL, common_levels = TRUE)
  if (ct_n$n_scored == 0) {
    logger::log_warn("computeStateConfusion: no bin is covered by both segmentations")
    return(empty_result)
  }
  ct_w <- if (is.null(bin_weight)) {
    ct_n
  } else {
    .contingency(query_state, subject_state, w = bin_weight, common_levels = TRUE)
  }

  # Rows are the query, columns the subject: transpose into long form with the
  # subject first, because the reference is the thing being recovered.
  states <- rownames(ct_n$tab)
  counts <- expand.grid(
    query_state = states, subject_state = states,
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
  counts$n <- as.numeric(ct_n$tab[cbind(counts$query_state, counts$subject_state)])
  counts$w <- as.numeric(ct_w$tab[cbind(counts$query_state, counts$subject_state)])
  counts <- counts[, c("subject_state", "query_state", "n", "w")]
  counts <- counts[order(counts$subject_state, counts$query_state), , drop = FALSE]
  rownames(counts) <- NULL

  tab_w <- ct_w$tab
  agreed <- diag(tab_w)
  support_query <- rowSums(tab_w)
  support_subject <- colSums(tab_w)

  recall <- .safeDiv(agreed, support_subject)
  precision <- .safeDiv(agreed, support_query)
  # Both zero is a real, informative answer -- the state was called and missed
  # entirely -- so F1 is 0 there rather than 0/0. It stays NA only when the
  # state is absent on one side and the question does not arise.
  f1 <- .safeDiv(2 * recall * precision, recall + precision, zero = 0)
  f1[is.na(recall) | is.na(precision)] <- NA_real_

  per_state <- data.frame(
    state = states,
    recall = as.numeric(recall),
    precision = as.numeric(precision),
    f1 = as.numeric(f1),
    support_subject = as.numeric(support_subject),
    support_query = as.numeric(support_query),
    stringsAsFactors = FALSE
  )
  rownames(per_state) <- NULL

  list(counts = counts, per_state = per_state)
}


#' Score a synthetic all-neutral caller against a reference
#'
#' @description
#' Builds the caller that emits nothing -- a single neutral call over every bin
#' the reference covers -- and scores it. This is the regression test for the
#' whole labelling scheme rather than a metric of any real segmentation.
#'
#' The use is diagnostic. `null_p_o` is the raw agreement that doing no work
#' achieves on this sample, and it is the number every raw-agreement score has
#' to be read against. `null_kappa` must come out at 0; it is computed rather
#' than asserted so that a bug in the marginals shows up here instead of being
#' silently absorbed into every real score. If a candidate quality metric gives
#' this caller a good score, that metric is confounded by the composition of the
#' reference and must not be used to define whether a sample was solved.
#'
#' @inheritSection computeStateKappa Choosing the bins
#' @inheritSection computeStateKappa NA handling
#'
#' @param subject A `cnvstrml` object or `GRanges` -- the reference.
#' @param bins A `GRanges` of bin-level regions, or a BED file path. Required.
#' @param statecol Character. Name of the state column in `subject`.
#' @param simplify,infer_as_neutral Passed through to the state
#'   canonicalisation; see [computeStateKappa()].
#' @param neutral_state Character. The state the null caller emits everywhere.
#'
#' @return A named `list`:
#' \describe{
#'   \item{null_kappa}{Cohen's kappa of the all-neutral caller. Zero, up to
#'     floating point.}
#'   \item{null_p_o}{Raw bin agreement of the all-neutral caller, which is the
#'     neutral bin fraction of the reference.}
#'   \item{null_tp_bin_fraction}{Same quantity expressed as a fraction of
#'     scorable bins, kept separate so it reads as a score rather than a
#'     composition statistic.}
#'   \item{null_tp_genome_fraction}{The base-pair-weighted equivalent, for
#'     comparison with genome-fraction metrics such as those from
#'     [cnvstrmliner::compareSegmentStatus()].}
#'   \item{n_bins_scored}{Bins covered by the reference.}
#' }
#'
#' @seealso [computeStateKappa()], [computeEventRecovery()]
#'
#' @examples
#' \dontrun{
#' null_scores <- scoreNullCaller(truth, bins = targets)
#' # Any raw-agreement score below this is worse than doing nothing.
#' null_scores$null_p_o
#' }
#'
#' @export
scoreNullCaller <- function(subject,
                            bins,
                            statecol = "cnv.annotation",
                            simplify = FALSE,
                            infer_as_neutral = TRUE,
                            neutral_state = "NEUT") {
  subject_gr <- .asSegments(subject, "subject", statecol)
  bins <- .resolveBins(bins, allow_auto = FALSE)

  empty_result <- list(
    null_kappa = NA_real_, null_p_o = NA_real_,
    null_tp_bin_fraction = NA_real_, null_tp_genome_fraction = NA_real_,
    n_bins_scored = 0L
  )

  if (length(bins) == 0 || length(subject_gr) == 0) {
    logger::log_warn("scoreNullCaller: empty input, returning NA")
    return(empty_result)
  }

  subject_state <- .binStateVector(bins, subject_gr, statecol, simplify, infer_as_neutral)
  scorable <- !is.na(subject_state)
  if (!any(scorable)) {
    logger::log_warn("scoreNullCaller: the reference covers no bin, returning NA")
    return(empty_result)
  }

  null_state <- rep(neutral_state, length(bins))
  null_state[!scorable] <- NA_character_

  ct <- .contingency(null_state, subject_state, common_levels = TRUE)
  k <- .kappaFromTable(ct$tab)

  bin_width <- as.numeric(BiocGenerics::width(bins))
  is_neutral <- scorable & subject_state == neutral_state

  list(
    null_kappa = k$kappa,
    null_p_o = k$p_o,
    null_tp_bin_fraction = sum(is_neutral) / sum(scorable),
    null_tp_genome_fraction = .safeDiv(
      sum(bin_width[is_neutral]), sum(bin_width[scorable])
    ),
    n_bins_scored = ct$n_scored
  )
}


#' Recall and precision on altered regions
#'
#' @description
#' Splits the two error directions that kappa merges: the fraction of genuinely
#' altered reference territory the query recovers, and the fraction of the
#' query's altered calls that are correct.
#'
#' Both are needed because the failure modes are not interchangeable. Missing
#' real events and inventing spurious ones are different problems with different
#' remedies -- the first usually wants a more sensitive setting, the second a
#' less sensitive one -- and a single summary score can hide either behind the
#' other.
#'
#' A segment counts as recovered when at least `concordance_threshold` of its
#' bins carry the matching state, so a segment that overlaps its counterpart
#' positionally but calls the wrong copy number does not count. Segments are
#' weighted by width by default, so a chromosome-arm event is not traded off
#' against a focal one.
#'
#' @section Segments with no bins:
#' An altered segment carrying no bin at all is dropped from both the numerator
#' and the denominator, and counted in `n_truth_unscorable` or
#' `n_query_unscorable`.
#'
#' The alternative -- scoring it as a miss -- charges the caller for events
#' lying entirely between capture targets, which no caller could have seen and
#' no setting could recover. That is a property of the assay, not of the
#' segmentation, and folding it into the score makes samples with more
#' off-target truth look like worse callers.
#'
#' Note what this is *not*. It excludes segments with no evidence whatsoever,
#' which is a geometric fact about the bins. It does not exclude segments whose
#' events are too faint to detect at the sample's purity and noise level -- that
#' is a detectability question, it needs a model of the assay rather than a
#' comparison of two segmentations, and it is deliberately out of scope here.
#' The recall reported is therefore still a raw recall, and reading it as an
#' attainable one will flatter every caller on a low-purity sample.
#'
#' @inheritSection computeStateKappa Choosing the bins
#' @inheritSection computeStateKappa NA handling
#'
#' @param query A `cnvstrml` object or `GRanges` -- the segmentation scored.
#' @param subject A `cnvstrml` object or `GRanges` -- the reference.
#' @param bins A `GRanges` of bin-level regions, or a BED file path. Required.
#' @param concordance_threshold Numeric in `[0, 1]`. The per-segment state
#'   concordance at which a segment counts as recovered. Default `0.5`.
#' @param statecol Character. Name of the state column.
#' @param simplify Logical. Passed to the state canonicalisation used when
#'   scoring concordance; see [computeStateKappa()].
#' @param neutral_states Character vector of states that are not events. Note
#'   this is a coarser question than `simplify` answers, which is why the
#'   altered/not split is controlled by its own argument here.
#' @param weight One of `"width"` (default) or `"bins"`.
#'
#' @return A named `list` with `recall_altered`, `precision_altered`,
#'   `sorensen`, `n_truth_altered` (altered segments in the reference),
#'   `n_query_altered` (altered segments in the query), the two
#'   `n_*_unscorable` counts described above, and `weight`. Rates are
#'   `NA_real_` where the corresponding side has no scorable altered segment --
#'   the question does not apply rather than having been answered badly.
#'
#' @seealso [computeStateKappa()], [annotateBinStateConcordance()],
#'   [computeStateConfusion()]
#'
#' @examples
#' \dontrun{
#' rec <- computeEventRecovery(facets_obj, truth, bins = targets)
#' # A high recall with a low precision is over-calling, not under-calling.
#' c(rec$recall_altered, rec$precision_altered)
#' }
#'
#' @export
computeEventRecovery <- function(query,
                                 subject,
                                 bins,
                                 concordance_threshold = 0.5,
                                 statecol = "cnv.annotation",
                                 simplify = FALSE,
                                 neutral_states = c("NEUT", "INFER"),
                                 weight = c("width", "bins")) {
  weight <- match.arg(weight)
  query_gr <- .asSegments(query, "query", statecol)
  subject_gr <- .asSegments(subject, "subject", statecol)
  bins <- .resolveBins(bins, allow_auto = FALSE)

  empty_result <- list(
    recall_altered = NA_real_, precision_altered = NA_real_,
    sorensen = NA_real_, n_truth_altered = 0L, n_query_altered = 0L,
    n_truth_unscorable = 0L, n_query_unscorable = 0L,
    weight = weight
  )

  if (length(bins) == 0 || length(query_gr) == 0 || length(subject_gr) == 0) {
    logger::log_warn("computeEventRecovery: empty input, returning NA")
    return(empty_result)
  }

  # Each side is scored against the other, so both annotations are needed.
  query_scored <- annotateBinStateConcordance(
    query_gr, subject_gr, bins,
    statecol = statecol, simplify = simplify, infer_as_neutral = FALSE
  )
  subject_scored <- annotateBinStateConcordance(
    subject_gr, query_gr, bins,
    statecol = statecol, simplify = simplify, infer_as_neutral = FALSE
  )

  rate <- function(gr) {
    state <- .canonicalStates(
      S4Vectors::mcols(gr)[[statecol]],
      simplify = simplify, infer_as_neutral = FALSE
    )
    altered <- !is.na(state) & !(state %in% neutral_states)
    # A segment carrying no bin has no evidence either way: scoring it as a
    # failure would charge the caller for events lying entirely between capture
    # targets, which no caller could have seen. It is dropped from both the
    # numerator and the denominator and counted separately.
    has_bins <- !is.na(S4Vectors::mcols(gr)$state.n_total) &
      S4Vectors::mcols(gr)$state.n_total > 0
    scorable <- altered & has_bins

    if (!any(scorable)) {
      return(list(rate = NA_real_, n = sum(altered), n_unscorable = sum(altered & !has_bins)))
    }

    conc <- S4Vectors::mcols(gr)$state.concordance[scorable]
    conc[is.na(conc)] <- 0
    recovered <- conc >= concordance_threshold
    seg_weight <- if (weight == "width") {
      as.numeric(BiocGenerics::width(gr))[scorable]
    } else {
      as.numeric(S4Vectors::mcols(gr)$state.n_total[scorable])
    }
    seg_weight[is.na(seg_weight)] <- 0
    list(
      rate = .safeDiv(sum(recovered * seg_weight), sum(seg_weight)),
      n = sum(altered),
      n_unscorable = sum(altered & !has_bins)
    )
  }

  recall <- rate(subject_scored)
  precision <- rate(query_scored)

  # The harmonic mean of precision and recall is algebraically identical to the
  # Sorensen-Dice coefficient 2*TP / (truth + query), which is why it is named
  # `sorensen` here and matches the field of the same name returned by
  # [cnvstrmliner::getJaccardIndex()].
  sorensen <- if (is.na(recall$rate) || is.na(precision$rate)) {
    NA_real_
  } else {
    .safeDiv(2 * recall$rate * precision$rate, recall$rate + precision$rate, zero = 0)
  }

  list(
    recall_altered = recall$rate,
    precision_altered = precision$rate,
    sorensen = sorensen,
    n_truth_altered = as.integer(recall$n),
    n_query_altered = as.integer(precision$n),
    n_truth_unscorable = as.integer(recall$n_unscorable),
    n_query_unscorable = as.integer(precision$n_unscorable),
    weight = weight
  )
}


#' Jaccard index of the altered genome, counted in bins
#'
#' @description
#' The bin-level counterpart of [cnvstrmliner::getJaccardIndex()]. Same formula,
#' same field names, same definition of a true positive -- the only difference is
#' that the sets are counted in bins of the evidence grid rather than in base
#' pairs.
#'
#' Each bin is assigned the state of the segment covering most of it, on both
#' sides, by [cnvstrmliner::assignBinState()]. A bin is then:
#'
#' * **altered in the query** if its query state is not neutral,
#' * **altered in the subject** if its subject state is not neutral,
#' * a **true positive** if it is altered on both sides *and* the two states are
#'   the same label.
#'
#' Requiring the same label is what makes this the analogue of the base-pair
#' function, where a region contributes to `tp` only when it is non-neutral on
#' both sides and carries `status == "TP"`. A bin called `AMP` by one side and
#' `HETDEL` by the other is in the union and not in the intersection.
#'
#' \deqn{J = \frac{TP}{S + Q - TP}, \qquad D = \frac{2\,TP}{S + Q}}
#'
#' where \eqn{S} and \eqn{Q} are the numbers of altered bins in the subject and
#' the query.
#'
#' **Why count bins.** On capture or amplicon data most of the genome carries no
#' evidence, so a base-pair Jaccard credits a segment for the unsequenced gaps it
#' spans and behaves partly as a segment-length metric. Counting bins restricts
#' both the intersection and the union to territory the assay could see. Bins that
#' no segment covers on one or both sides are excluded from the intersection and
#' from the union alike, and are never imputed to neutral.
#'
#' @param query A `cnvstrml` object or `GRanges` -- the segmentation scored.
#' @param subject A `cnvstrml` object or `GRanges` -- the reference, treated as
#'   truth by the two directional fractions.
#' @param bins A `GRanges` of bin-level regions, or a BED file path. On capture
#'   data pass the target regions.
#' @param statecol Character. Name of the state column. Default
#'   `"cnv.annotation"`.
#' @param simplify Logical. If `TRUE`, `-LOH` suffixes are stripped and deletion
#'   states are pooled before comparison, matching
#'   [cnvstrmliner::cnvMetrics()]. Default `FALSE`, which keeps
#'   `HETDEL` and `HOMDEL` distinct.
#' @param neutral_states Character vector of states that do not count as
#'   alterations. Default `c("NEUT", "INFER")`.
#'
#' @return A named list, mirroring [cnvstrmliner::getJaccardIndex()] with bin
#'   counts in place of base pairs:
#'   \describe{
#'     \item{jaccard_bin}{`tp / (subject + query - tp)`, in `[0, 1]`.}
#'     \item{sorensen_bin}{`2 * tp / (subject + query)`, the Sorensen-Dice
#'       coefficient.}
#'     \item{subject_altered_genome_n}{Number of bins the subject calls altered.}
#'     \item{query_altered_genome_n}{Number of bins the query calls altered.}
#'     \item{tp_altered_genome_n}{Number of bins altered on both sides and given
#'       the same state.}
#'     \item{matched_query_altered_genome}{`tp / query`, a precision. Named as in
#'       [cnvstrmliner::getJaccardIndex()], since the fraction is unitless.}
#'     \item{matched_subject_altered_genome}{`tp / subject`, a recall.}
#'   }
#'   The `_n` suffix marks the three counts as numbers of bins, where the
#'   base-pair function uses `_bp`. When neither side calls a single bin altered
#'   the union is empty; `jaccard_bin` and `sorensen_bin` are then `NA_real_`
#'   rather than `NaN`, and the two directional fractions are `NA_real_` whenever
#'   their denominator is zero.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' j <- getBinJaccardIndex(query_obj, reference_obj, bins = targets)
#'
#' unlist(j[c("jaccard_bin", "sorensen_bin")])
#'
#' # The two directions the Jaccard combines. Perfect precision with poor recall
#' # is under-calling, not accuracy.
#' unlist(j[c("matched_query_altered_genome", "matched_subject_altered_genome")])
#'
#' # Same quantity in base pairs, for comparison.
#' getJaccardIndex(compareSegmentStatus(query_obj, reference_obj))$jaccard
#' }
#'
#' @seealso
#' [cnvstrmliner::getJaccardIndex()] for the base-pair version,
#' [cnvstrmliner::annotateBinJaccard()] for a per-segment rather than
#' genome-wide bin Jaccard,
#' [cnvstrmliner::computeEventRecovery()] for the event-level view, which asks
#' whether each altered segment was recovered rather than how many bins agree
#'
getBinJaccardIndex <- function(query,
                               subject,
                               bins,
                               statecol = "cnv.annotation",
                               simplify = FALSE,
                               neutral_states = c("NEUT", "INFER")) {
  bins <- .resolveBins(bins, allow_auto = FALSE)

  query_state <- assignBinState(bins, query,
    statecol = statecol, simplify = simplify
  )
  subject_state <- assignBinState(bins, subject,
    statecol = statecol, simplify = simplify
  )

  # A bin with no covering segment on either side is unscorable, and must stay
  # out of both the intersection and the union rather than counting as neutral.
  scorable <- !is.na(query_state) & !is.na(subject_state)

  query_altered <- scorable & !query_state %in% neutral_states
  subject_altered <- scorable & !subject_state %in% neutral_states
  tp_altered <- query_altered & subject_altered &
    query_state == subject_state

  n_query <- sum(query_altered)
  n_subject <- sum(subject_altered)
  n_tp <- sum(tp_altered)
  n_union <- n_subject + n_query - n_tp

  list(
    jaccard_bin = .safeDiv(n_tp, n_union),
    sorensen_bin = .safeDiv(2 * n_tp, n_subject + n_query),
    subject_altered_genome_n = as.integer(n_subject),
    query_altered_genome_n = as.integer(n_query),
    tp_altered_genome_n = as.integer(n_tp),
    matched_query_altered_genome = .safeDiv(n_tp, n_query),
    matched_subject_altered_genome = .safeDiv(n_tp, n_subject)
  )
}

#' Greedy one-to-one breakpoint matching
#'
#' @description
#' Maximum-cardinality matching between two sets of breakpoints under a distance
#' tolerance, where each breakpoint may be used at most once. Replaces the
#' many-to-many \code{.matchBreakpoints()} used by
#' \code{\link[cnvstrmliner]{computeBoundaryConcordance}}.
#'
#' @details
#' Breakpoints live on a line, so the optimal matching does not need a general
#' bipartite solver. Sorting both sets per chromosome and sweeping them with two
#' pointers -- matching the earliest pair that is within tolerance, otherwise
#' advancing whichever pointer is behind -- yields a maximum matching. The
#' exchange argument is the usual one: if an optimal solution does not match the
#' earliest matchable pair, it can be rewritten to do so without losing a pair.
#'
#' The sweep is O(n + m) after the sort, against the O(n * m) of an overlap join,
#' so this is also the cheaper of the two on dense inputs.
#'
#' @param query_bp,subject_bp \code{GRanges} of breakpoint positions, as returned
#'   by \code{\link[cnvstrmliner]{getSegmentBreakpoints}}.
#' @param tolerance_bp Numeric. Maximum distance for a pair to count as matched.
#'
#' @return A list with \code{n_matched}, the number of matched pairs. Because the
#'   assignment is one-to-one this single number serves both directions: it is the
#'   numerator of precision and of recall alike.
#'
#' @keywords internal
matchBreakpointsOneToOne <- function(query_bp, subject_bp, tolerance_bp = 3e5) {
  if (length(query_bp) == 0 || length(subject_bp) == 0) {
    return(list(n_matched = 0L))
  }

  # Only chromosomes present on both sides can contribute, and a match may never
  # cross a chromosome boundary -- hence the per-seqname sweep rather than one
  # global sort of the positions.
  q_chr <- as.character(GenomicRanges::seqnames(query_bp))
  s_chr <- as.character(GenomicRanges::seqnames(subject_bp))
  shared <- intersect(unique(q_chr), unique(s_chr))
  if (length(shared) == 0) {
    return(list(n_matched = 0L))
  }

  q_pos <- GenomicRanges::start(query_bp)
  s_pos <- GenomicRanges::start(subject_bp)

  n_matched <- sum(vapply(shared, function(ch) {
    q <- sort(q_pos[q_chr == ch])
    s <- sort(s_pos[s_chr == ch])
    i <- 1L
    j <- 1L
    n <- 0L
    nq <- length(q)
    ns <- length(s)
    while (i <= nq && j <= ns) {
      d <- q[i] - s[j]
      if (abs(d) <= tolerance_bp) {
        # Earliest matchable pair: consume both.
        n <- n + 1L
        i <- i + 1L
        j <- j + 1L
      } else if (d < 0) {
        # q[i] is too far to the left to ever match s[j] or anything after it.
        i <- i + 1L
      } else {
        j <- j + 1L
      }
    }
    n
  }, integer(1)))

  list(n_matched = as.integer(n_matched))
}


#' Breakpoint concordance with a one-to-one assignment
#'
#' @description
#' Precision, recall and F1 for breakpoint detection, counting a query
#' breakpoint as correct when a reference breakpoint lies within
#' `tolerance_bp`.
#'
#' This is the structural half of segmentation quality, deliberately blind to
#' copy-number state. A caller that finds every real transition but assigns the
#' wrong state everywhere scores 1 here and near 0 on [computeStateKappa()];
#' that contrast is the diagnostic. Overlap-based metrics conflate the two and
#' cannot make it.
#'
#' @section Choosing the tolerance:
#' The default of 100 kb is roughly one bin at the 100 kb scale these
#' analyses use, and it encodes a real claim: a breakpoint called within one
#' bin of the truth is the same breakpoint, because no exome assay localises a
#' transition more precisely than the spacing of its targets. Tighten it only if
#' the grid genuinely supports it. Note the number is not symmetric in effect --
#' a large tolerance flatters an over-segmented caller, because more of its
#' spurious breakpoints fall inside some window.
#'
#' @section Why the two match counts can differ:
#' Matching is done independently in each direction rather than as a one-to-one
#' assignment. Three query breakpoints clustered inside one tolerance window all
#' count as matched, but they match the same single reference breakpoint, so
#' `n_matched_query` is 3 while `n_matched_subject` is 1. That asymmetry is the
#' intended behaviour: precision should reward each of those three only in
#' proportion to how many distinct real transitions exist, which is what the
#' separate denominators achieve.
#'
#' @section Edge cases:
#' Two empty breakpoint sets agree perfectly and score 1 -- both segmentations
#' assert that the genome has no internal transitions, and they are right about
#' each other. An empty query against a non-empty reference scores 0. A
#' non-empty query against an empty reference scores precision 0, recall 1
#' (there are no false negatives available), F1 0.
#'
#' @section What changes:
#' We could use computeBoundaryConcordance() to match breakpoints many-to-many. It
#' asks, for each breakpoint, "is there ANY partner within tolerance", via
#'      findOverlaps(query_bp + tolerance_bp, subject_bp)
#'      n_matched_query <- length(unique(queryHits(hits)))
#' and the same in the other direction. Nothing stops one breakpoint from claiming
#' several partners, so both precision and recall are inflated whenever breakpoints
#' are dense relative to the tolerance -- which on this cohort is the normal case,
#' not an edge case:
#'    * median gap between consecutive TRUE breakpoints, in samples with more than
#'     300 of them, is 81 kb -- against a tolerance of 300 kb. 55% of neighbouring
#'     truth breakpoints are closer to each other than the tolerance, so a +/- 300 kb
#'     window routinely spans several of them.
#'   * measured on 60 FACETS/TCGA pairs, boundary_f1 as computed by the package is
#'     1.29-1.42x the value obtained from a one-to-one matching at the same
#'     tolerance (0.507 -> 0.357, 0.416 -> 0.295, 0.395 -> 0.307 by truth-breakpoint
#'     band).
#'
#' The fix is to require a one-to-one assignment: every predicted breakpoint may
#' account for at most one true breakpoint and vice versa. That is what a reader
#' assumes "precision" and "recall" mean on a set of positions.
#' With a one-to-one assignment \code{n_matched_query == n_matched_subject}, so
#' precision and recall share a numerator and differ only in their denominator --
#' the number of predicted and of true breakpoints respectively. Under the
#' many-to-many version the two numerators could differ, which is how a single
#' well-placed call in a dense cluster could lift recall without touching
#' precision.
#'
#' Note that the tolerance itself is a separate question. 300 kb is roughly four
#' times the median true segment length in fragmented samples, so a call can be
#' three or four segments away and still count as correct. This file does not
#' change the default; it only stops the double counting. Sweeping the tolerance is
#' left to the caller.
#'
#' Degenerate cases are preserved exactly as the package defines them: two empty
#' sets score 1, an empty query scores 0, and an empty subject scores
#' \code{f1 = 0} with \code{recall = 1}.
#'
#' @param query,subject \code{cnvstrml} objects. \code{query} is the prediction and
#'   \code{subject} the truth, matching the convention in \code{scripts/R/README.md}.
#' @param tolerance_bp Numeric. Matching window in base pairs. Default 10,000bp, kept
#'   for continuity; see the header of this file on why a smaller value is worth
#'   considering on fragmented truths.
#' @param exclude_chr_ends Logical. Drop the first start and last end of each
#'   chromosome, which are reference artefacts rather than caller decisions.
#'
#' @return A list with \code{boundary_f1}, \code{boundary_precision},
#'   \code{boundary_recall}, \code{n_query_breakpoints},
#'   \code{n_subject_breakpoints}, \code{n_matched_query}, \code{n_matched_subject}
#'   and \code{tolerance_bp}.
#'
#' @seealso \code{matchBreakpointsOneToOne}
#'
#' @export
computeBoundaryConcordance <- function(query, subject, tolerance_bp = 1e5,
                                       exclude_chr_ends = TRUE) {
  query_bp <- cnvstrmliner::getSegmentBreakpoints(query, exclude_chr_ends = exclude_chr_ends)
  subject_bp <- cnvstrmliner::getSegmentBreakpoints(subject, exclude_chr_ends = exclude_chr_ends)
  n_query <- length(query_bp)
  n_subject <- length(subject_bp)

  result <- function(f1, precision, recall, matched = 0L) {
    list(
      boundary_f1 = f1, boundary_precision = precision, boundary_recall = recall,
      n_query_breakpoints = as.integer(n_query),
      n_subject_breakpoints = as.integer(n_subject),
      n_matched_query = as.integer(matched),
      n_matched_subject = as.integer(matched),
      tolerance_bp = tolerance_bp
    )
  }

  if (n_query == 0 && n_subject == 0) {
    return(result(1, 1, 1))
  }
  if (n_query == 0) {
    return(result(0, 0, 0))
  }
  if (n_subject == 0) {
    return(result(0, 0, 1))
  }

  m <- matchBreakpointsOneToOne(query_bp, subject_bp, tolerance_bp)
  precision <- m$n_matched / n_query
  recall <- m$n_matched / n_subject
  f1 <- if (precision + recall == 0) 0 else 2 * precision * recall / (precision + recall)
  result(f1, precision, recall, m$n_matched)
}


#' Adjusted Rand index between two segmentations
#'
#' @description
#' Treats each segmentation as a partition of genomic bins and measures how
#' often the two agree about whether a pair of bins belongs together, corrected
#' for the agreement expected by chance.
#'
#' @section Segment identity versus copy-number state:
#' These are different questions and the answer differs sharply, so the mode is
#' explicit rather than implied.
#'
#' With `partition_by = "segment"` (the default) each bin is labelled by *which
#' segment contains it*. A caller that splits every reference segment cleanly in
#' two, assigning the correct state to both halves, scores about 0.5: it has got
#' the biology right and the granularity wrong, and that is what this mode
#' reports.
#'
#' With `partition_by = "state"` each bin is labelled by its *copy-number
#' state*. The same caller scores 1.0.
#'
#' The state mode is provided for completeness, but [computeStateKappa()] is the
#' better instrument for that question: kappa is the standard chance-corrected
#' agreement statistic and its value has a direct reading, whereas ARI over a
#' state partition is a pair-counting quantity with no comparable
#' interpretation. The identity mode is what this function is for.
#'
#' @section Choosing the bins:
#' Unlike the metrics in section 1, `bins` may be `NULL` here. A bin in this
#' context is a discretisation of the coordinate axis rather than a unit of
#' evidence, so a uniform tiling is a defensible quadrature rule even on exome
#' data. When `bins` is `NULL` the grid is tiled from the genome's sequence
#' lengths at `bin_width`, restricted to the sequences present in the inputs.
#'
#' @section Differences from the project-local computeARI:
#' Three, all intentional. Expect agreement to about three decimal places on
#' real data, not to floating point.
#'
#' The auto-tiled grid is derived from the genome, not from the coordinate range
#' of the two inputs. The project-local implementation tiles from the minimum
#' start to the maximum end of the segmentations being compared, which puts
#' every caller on its own grid and makes the resulting numbers incomparable
#' across a sweep -- the exact comparison a sweep exists to make.
#'
#' Bins are assigned to the segment they overlap most. The project-local
#' implementation labels each bin with whichever overlapping segment
#' `GenomicRanges::findOverlaps()` reports last, so a bin straddling a boundary
#' is credited to its right-hand neighbour regardless of how little of it lies
#' there, and the labelling depends on segment order in the object. This is the
#' residual source of disagreement when `bins` is passed explicitly to both.
#'
#' `ri` is the full Rand index, counting pairs the two partitions agree to
#' separate as well as pairs they agree to join. The project-local version
#' reports only the joined term, which understates agreement severely whenever
#' there are many segments -- as there always are. `ari` is unaffected, since
#' the adjustment cancels the difference.
#'
#' @inheritSection computeStateKappa NA handling
#'
#' @param query A `cnvstrml` object or `GRanges`.
#' @param subject A `cnvstrml` object or `GRanges`.
#' @param bins A `GRanges`, a BED file path, or `NULL` to auto-tile.
#' @param bin_width Numeric width of auto-tiled bins. Default `1e6`. Ignored
#'   when `bins` is supplied.
#' @param partition_by One of `"segment"` (default) or `"state"`; see the
#'   section above.
#' @param statecol Character. Name of the state column. Used only when
#'   `partition_by = "state"`.
#' @param simplify,infer_as_neutral Logical. State canonicalisation, used only
#'   when `partition_by = "state"`; see [computeStateKappa()].
#' @param genome Character genome identifier for auto-tiling. Defaults to the
#'   `genome` slot of `query`, then of `subject`.
#'
#' @return A named `list` with `ari`, `ri`, `n_bins` (bins scored), the group
#'   counts `n_query_groups` and `n_subject_groups`, and `partition_by` echoed
#'   back so a stored result cannot be misread later.
#'
#' @seealso [computeSegmentationVI()] for the information-theoretic companion,
#'   [computeFragmentation()] for the direction of the disagreement,
#'   [computeStateKappa()] for state agreement.
#'
#' @examples
#' \dontrun{
#' # Structural agreement, on a grid that does not depend on the inputs.
#' computeSegmentationARI(facets_obj, truth, bin_width = 1e6)$ari
#' }
#'
#' @export
computeSegmentationARI <- function(query,
                                   subject,
                                   bins = NULL,
                                   bin_width = 1e6,
                                   partition_by = c("segment", "state"),
                                   statecol = "cnv.annotation",
                                   simplify = FALSE,
                                   infer_as_neutral = TRUE,
                                   genome = NULL) {
  partition_by <- match.arg(partition_by)
  parts <- .partitionBins(
    query, subject, bins, bin_width, partition_by,
    statecol, simplify, infer_as_neutral, genome
  )

  if (parts$n_scored < 2) {
    logger::log_warn("computeSegmentationARI: fewer than two scorable bins, returning NA")
    return(list(
      ari = NA_real_, ri = NA_real_, n_bins = parts$n_scored,
      n_query_groups = 0L, n_subject_groups = 0L,
      partition_by = partition_by
    ))
  }

  a <- .ariFromTable(parts$tab)
  list(
    ari = a$ari,
    ri = a$ri,
    n_bins = parts$n_scored,
    n_query_groups = nrow(parts$tab),
    n_subject_groups = ncol(parts$tab),
    partition_by = partition_by
  )
}


#' Variation of information between two segmentations
#'
#' @description
#' An information-theoretic distance between two partitions of the genome:
#' how much is still unknown about one segmentation once the other is known.
#'
#' \deqn{VI(A, B) = H(A) + H(B) - 2\,I(A; B)}
#'
#' It is a true metric on partitions -- symmetric and obeying the triangle
#' inequality -- which the adjusted Rand index is not, and it decomposes: the
#' `h_query - mi` part is information the query has that the reference does not
#' (over-segmentation), the `h_subject - mi` part is the reverse
#' (under-segmentation).
#'
#' Lower is better, and 0 means the two partitions are identical.
#'
#' @section Why nvi is reported:
#' Raw VI is unbounded above and grows with the number of bins, so a VI of 1.2
#' on a 3,000-bin grid and a VI of 1.2 on a 200,000-bin grid do not mean the
#' same thing. Any cohort-level aggregation of raw VI across samples with
#' different grids is therefore meaningless. `nvi` divides by `log(n_bins)`, the
#' maximum attainable value, giving a `[0, 1]` quantity that can be averaged.
#' Prefer it for anything comparative.
#'
#' @inheritSection computeSegmentationARI Segment identity versus copy-number state
#' @inheritSection computeSegmentationARI Choosing the bins
#' @inheritSection computeStateKappa NA handling
#'
#' @inheritParams computeSegmentationARI
#'
#' @return A named `list`:
#' \describe{
#'   \item{vi}{Variation of information, in nats. Zero for identical
#'     partitions.}
#'   \item{nvi}{`vi / log(n_bins)`, in `[0, 1]`.}
#'   \item{h_query, h_subject}{Entropy of each partition.}
#'   \item{mi}{Mutual information.}
#'   \item{n_bins}{Bins scored.}
#'   \item{partition_by}{Echoed back.}
#' }
#'
#' @seealso [computeSegmentationARI()], [computeFragmentation()]
#'
#' @examples
#' \dontrun{
#' v <- computeSegmentationVI(facets_obj, truth, bin_width = 1e6)
#' # Which side carries the excess information?
#' c(over = v$h_query - v$mi, under = v$h_subject - v$mi)
#' }
#'
#' @export
computeSegmentationVI <- function(query,
                                  subject,
                                  bins = NULL,
                                  bin_width = 1e6,
                                  partition_by = c("segment", "state"),
                                  statecol = "cnv.annotation",
                                  simplify = FALSE,
                                  infer_as_neutral = TRUE,
                                  genome = NULL) {
  partition_by <- match.arg(partition_by)
  parts <- .partitionBins(
    query, subject, bins, bin_width, partition_by,
    statecol, simplify, infer_as_neutral, genome
  )

  if (parts$n_scored == 0) {
    logger::log_warn("computeSegmentationVI: no scorable bin, returning NA")
    return(list(
      vi = NA_real_, nvi = NA_real_, h_query = NA_real_,
      h_subject = NA_real_, mi = NA_real_, n_bins = 0L,
      partition_by = partition_by
    ))
  }

  v <- .viFromTable(parts$tab)
  list(
    vi = v$vi,
    nvi = v$nvi,
    h_query = v$h_a,
    h_subject = v$h_b,
    mi = v$mi,
    n_bins = parts$n_scored,
    partition_by = partition_by
  )
}


#' Direction of the segmentation disagreement
#'
#' @description
#' Whether the query splits the reference or merges it, and by how much.
#'
#' The adjusted Rand index and variation of information both say that two
#' segmentations differ without saying which way. That is the missing half of
#' the diagnosis in a sensitivity sweep, where the two directions call for
#' opposite adjustments: a caller that shatters real events wants a less
#' sensitive setting, one that merges them wants a more sensitive one. Both look
#' identical in a symmetric distance.
#'
#' The ratios are computed on the bin contingency table, so they count only
#' segments that actually carry bins -- a segment over an unsequenced region is
#' not evidence of fragmentation.
#'
#' @section Reading the numbers:
#' `mean_query_per_subject` is the average number of distinct query segments
#' overlapping one reference segment, weighted by bins. A value of 2 means the
#' query cuts each reference segment in half on average; 1 means one-to-one;
#' below 1 is not attainable, since a reference segment always meets at least
#' one query segment among its scorable bins.
#'
#' `split_ratio` and `merge_ratio` are the segment-count ratios in each
#' direction, which respond faster than the bin-weighted means but ignore how
#' much genome each segment covers. Read them together: a `split_ratio` of 5
#' with a `mean_query_per_subject` of 1.1 is a caller adding many tiny segments
#' in one place, not fragmenting throughout.
#'
#' @inheritSection computeSegmentationARI Choosing the bins
#' @inheritSection computeStateKappa NA handling
#'
#' @param query A `cnvstrml` object or `GRanges`.
#' @param subject A `cnvstrml` object or `GRanges`.
#' @param bins A `GRanges`, a BED file path, or `NULL` to auto-tile.
#' @param bin_width Numeric width of auto-tiled bins. Default `1e6`.
#' @param genome Character genome identifier for auto-tiling.
#'
#' @return A named `list` with `split_ratio`, `merge_ratio`,
#'   `n_query_segments`, `n_subject_segments` (both counted among scorable bins),
#'   `mean_query_per_subject` and `mean_subject_per_query`.
#'
#' @seealso [computeSegmentationARI()], [computeBoundaryConcordance()]
#'
#' @examples
#' \dontrun{
#' f <- computeFragmentation(facets_obj, truth)
#' f$mean_query_per_subject # > 1 means the caller is splitting
#' }
#'
#' @export
computeFragmentation <- function(query,
                                 subject,
                                 bins = NULL,
                                 bin_width = 1e6,
                                 genome = NULL) {
  parts <- .partitionBins(
    query, subject, bins, bin_width, "segment",
    "cnv.annotation", FALSE, TRUE, genome
  )

  empty_result <- list(
    split_ratio = NA_real_, merge_ratio = NA_real_,
    n_query_segments = 0L, n_subject_segments = 0L,
    mean_query_per_subject = NA_real_, mean_subject_per_query = NA_real_
  )

  if (parts$n_scored == 0) {
    logger::log_warn("computeFragmentation: no scorable bin, returning NA")
    return(empty_result)
  }

  tab <- parts$tab
  n_query <- nrow(tab)
  n_subject <- ncol(tab)

  # Number of distinct partners each segment meets, and the bin mass it holds.
  partners_per_subject <- colSums(tab > 0)
  partners_per_query <- rowSums(tab > 0)
  mass_per_subject <- colSums(tab)
  mass_per_query <- rowSums(tab)

  list(
    split_ratio = .safeDiv(n_query, n_subject),
    merge_ratio = .safeDiv(n_subject, n_query),
    n_query_segments = as.integer(n_query),
    n_subject_segments = as.integer(n_subject),
    mean_query_per_subject = .safeDiv(
      sum(partners_per_subject * mass_per_subject), sum(mass_per_subject)
    ),
    mean_subject_per_query = .safeDiv(
      sum(partners_per_query * mass_per_query), sum(mass_per_query)
    )
  )
}


# ===== 4. TRUTH-FREE CROSS-SEGMENTATION ======================================
#
# Everything above needs a reference that is known to be right, which exists
# only for simulated data. These functions compare a family of segmentations to
# each other, so they run on real samples where the answer is not available.
#
# The standing caveat applies to all three and is repeated in each: agreement
# measures typicality, not correctness. Where callers share a failure mode the
# consensus is wrong and the metric is actively misleading, and at very low
# purity every caller agrees that nothing happened -- perfect agreement, no
# information. A high score is good news only conditional on the sample being
# solvable at all, which is a separate question these functions do not answer.
# =============================================================================


#' Recurrence of breakpoints across a family of segmentations
#'
#' @description
#' Clusters the boundaries of several segmentations of the same sample and
#' reports, for each cluster, the fraction of segmentations that place a
#' boundary there.
#'
#' The output is a per-position table, which makes it usable two ways: as a map
#' of where the data forces a transition regardless of setting, and -- through
#' [annotateBreakpointConsistency()] -- as a per-segment score for how much of a
#' given segmentation rests on well-supported boundaries.
#'
#' @section How clustering works:
#' Positions on each sequence are sorted and walked in order; a new cluster
#' opens when the next position is further than `tolerance_bp` from the running
#' *median* of the open cluster. Using the median rather than the last position
#' matters: with the last position, a chain of boundaries each just inside the
#' tolerance of its predecessor would merge into one cluster of unbounded width.
#' The cluster is reported at its median position.
#'
#' Each segmentation contributes at most one vote per cluster, so a caller that
#' shatters a region into many segments does not manufacture support for it.
#'
#' @param cnv_list A named list of `cnvstrml` objects or `GRanges`. Names become
#'   the voter identities; unnamed lists are numbered.
#' @param tolerance_bp Numeric. Clustering tolerance in base pairs. Default
#'   `3e5`.
#' @param exclude_chr_ends Logical. Passed to [getSegmentBreakpoints()].
#' @param min_consistency Numeric. Drop clusters below this consistency from the
#'   output. Default `0`, keeping everything.
#'
#' @return A `data.frame` with one row per cluster and columns `seqnames`, `pos`
#'   (the cluster median), `n_methods` (segmentations voting for it),
#'   `total_methods` and `consistency` (the ratio). Sorted by position.
#'
#' @seealso [annotateBreakpointConsistency()], [computeCallerAgreement()],
#'   [computeBoundaryConcordance()] for the same idea against a known reference.
#'
#' @examples
#' \dontrun{
#' consistency <- computeBreakpointConsistency(facets_settings, tolerance_bp = 3e5)
#' # Boundaries every setting agrees on.
#' subset(consistency, consistency == 1)
#' }
#'
#' @export
computeBreakpointConsistency <- function(cnv_list,
                                         tolerance_bp = 3e5,
                                         exclude_chr_ends = TRUE,
                                         min_consistency = 0) {
  if (!is.list(cnv_list)) {
    stop("cnv_list must be a list of cnvstrml or GRanges objects")
  }

  empty <- data.frame(
    seqnames = character(), pos = numeric(), n_methods = integer(),
    total_methods = integer(), consistency = numeric(),
    stringsAsFactors = FALSE
  )

  total_methods <- length(cnv_list)
  if (total_methods == 0) {
    return(empty)
  }

  method_names <- names(cnv_list)
  if (is.null(method_names)) {
    method_names <- paste0("method_", seq_along(cnv_list))
  }

  per_method <- lapply(seq_along(cnv_list), function(i) {
    bp <- .breakpointTable(
      .asSegments(cnv_list[[i]], method_names[i]),
      exclude_chr_ends = exclude_chr_ends
    )
    if (nrow(bp) == 0) {
      return(NULL)
    }
    bp$method <- method_names[i]
    bp
  })
  per_method <- per_method[!vapply(per_method, is.null, logical(1))]

  if (length(per_method) == 0) {
    return(empty)
  }

  bp_df <- do.call(rbind, per_method)
  bp_df <- bp_df[order(bp_df$seqnames, bp_df$pos), , drop = FALSE]

  clusters <- lapply(unique(bp_df$seqnames), function(chr) {
    chr_bp <- bp_df[bp_df$seqnames == chr, , drop = FALSE]
    cluster_id <- .clusterPositions(chr_bp$pos, tolerance_bp)
    data.frame(
      seqnames = chr,
      pos = as.numeric(tapply(chr_bp$pos, cluster_id, stats::median)),
      n_methods = as.integer(tapply(
        chr_bp$method, cluster_id,
        function(m) length(unique(m))
      )),
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, clusters)
  out$total_methods <- as.integer(total_methods)
  out$consistency <- out$n_methods / total_methods
  out <- out[out$consistency >= min_consistency, , drop = FALSE]
  out <- out[order(out$seqnames, out$pos), , drop = FALSE]
  rownames(out) <- NULL
  out
}


#' Agreement among a family of callers about which genome is altered
#'
#' @description
#' Tiles the genome, asks every segmentation whether each bin is altered, and
#' summarises how much the answers coincide.
#'
#' This is the most direct truth-free proxy for "did the callers see the same
#' thing", and it is the natural companion to
#' [computeBreakpointConsistency()]: that one asks whether they cut in the same
#' places, this one whether they agree about what is there.
#'
#' @section What high agreement does and does not mean:
#' Agreement is not monotone with quality, which is why this belongs in a model
#' rather than a threshold rule. At very low purity every caller agrees that
#' nothing happened and every caller is wrong; the same maximal agreement is
#' produced by an easy sample and by an impossible one. High agreement is
#' evidence of quality only conditional on adequate purity and noise, and those
#' have to come from somewhere else.
#'
#' @section Vocabulary:
#' The pairwise quantity is called `overlap` rather than `jaccard`, and nothing
#' here is called recall or precision, because these are caller-versus-caller
#' quantities with no reference involved. Downstream truth-free guards key on
#' the name; renaming these to match the truth-based metrics would let a
#' consensus feature pass for an achieved-quality one.
#'
#' @inheritSection computeSegmentationARI Choosing the bins
#'
#' @param cnv_list A named list of `cnvstrml` objects or `GRanges`.
#' @param bins A `GRanges`, a BED file path, or `NULL` to auto-tile.
#' @param bin_width Numeric width of auto-tiled bins. Default `1e6`.
#' @param statecol Character. Name of the state column.
#' @param neutral_states Character vector of states counted as not altered.
#' @param genome Character genome identifier for auto-tiling. Defaults to the
#'   first element of `cnv_list` that carries one.
#'
#' @return A named `list`:
#' \describe{
#'   \item{pairwise}{A `data.frame` with one row per unordered caller pair and
#'     columns `method_a`, `method_b`, `overlap` (Jaccard of their altered bin
#'     sets), `n_a`, `n_b`, `n_shared`.}
#'   \item{summary}{A one-row `data.frame` with `agree_overlap_mean`,
#'     `agree_overlap_max`, `agree_overlap_sd`, `agree_altered_rate_mean`,
#'     `agree_frac_majority` (fraction of bins more than half the callers call
#'     altered), `agree_vote_entropy_mean` (per-bin vote entropy, maximal when
#'     the callers split evenly), `agree_n_methods` and `n_bins`.}
#' }
#'
#' @seealso [computeConsensusAgreement()] for the per-caller version,
#'   [computeBreakpointConsistency()], [cnvstrmliner::findConsensus()]
#'
#' @examples
#' \dontrun{
#' ag <- computeCallerAgreement(caller_list, bins = targets)
#' ag$summary$agree_overlap_mean
#' }
#'
#' @export
computeCallerAgreement <- function(cnv_list,
                                   bins = NULL,
                                   bin_width = 1e6,
                                   statecol = "cnv.annotation",
                                   neutral_states = c("NEUT", "INFER"),
                                   genome = NULL) {
  prepared <- .prepareAlteredMatrix(
    cnv_list, bins, bin_width, statecol,
    neutral_states, genome
  )
  altered <- prepared$altered
  method_names <- colnames(altered)
  n_methods <- ncol(altered)

  empty_pairwise <- data.frame(
    method_a = character(), method_b = character(), overlap = numeric(),
    n_a = integer(), n_b = integer(), n_shared = integer(),
    stringsAsFactors = FALSE
  )
  summary_cols <- c(
    "agree_overlap_mean", "agree_overlap_max", "agree_overlap_sd",
    "agree_altered_rate_mean", "agree_frac_majority",
    "agree_vote_entropy_mean", "agree_n_methods", "n_bins"
  )

  if (n_methods < 2 || nrow(altered) == 0) {
    logger::log_warn("computeCallerAgreement: fewer than two callers or no bins, returning NA")
    summary_df <- .emptyRow(summary_cols)
    summary_df$agree_n_methods <- n_methods
    summary_df$n_bins <- nrow(altered)
    return(list(pairwise = empty_pairwise, summary = summary_df))
  }

  pairs <- utils::combn(seq_len(n_methods), 2)
  pairwise <- data.frame(
    method_a = method_names[pairs[1, ]],
    method_b = method_names[pairs[2, ]],
    overlap = apply(pairs, 2, function(p) .jaccardLogical(altered[, p[1]], altered[, p[2]])),
    n_a = as.integer(colSums(altered)[pairs[1, ]]),
    n_b = as.integer(colSums(altered)[pairs[2, ]]),
    n_shared = apply(pairs, 2, function(p) sum(altered[, p[1]] & altered[, p[2]])),
    stringsAsFactors = FALSE
  )

  votes <- rowSums(altered)
  vote_fraction <- votes / n_methods
  # Binary entropy of the per-bin vote split, in nats: 0 when the callers are
  # unanimous either way, maximal when they split down the middle.
  entropy_terms <- function(p) ifelse(p > 0, -p * log(p), 0)
  vote_entropy <- entropy_terms(vote_fraction) + entropy_terms(1 - vote_fraction)

  summary_df <- data.frame(
    agree_overlap_mean = mean(pairwise$overlap),
    agree_overlap_max = max(pairwise$overlap),
    agree_overlap_sd = if (nrow(pairwise) > 1) stats::sd(pairwise$overlap) else 0,
    agree_altered_rate_mean = mean(colMeans(altered)),
    agree_frac_majority = mean(votes > n_methods / 2),
    agree_vote_entropy_mean = mean(vote_entropy),
    agree_n_methods = as.integer(n_methods),
    n_bins = as.integer(nrow(altered)),
    stringsAsFactors = FALSE
  )

  list(pairwise = pairwise, summary = summary_df)
}


#' Leave-one-out agreement of each caller with the consensus of the others
#'
#' @description
#' For every caller, how well its altered regions match the majority verdict of
#' the *other* callers.
#'
#' This is the truth-free selection signal: the output that best matches what
#' everyone else independently found is, absent a shared failure mode, more
#' likely to be right. The consensus is computed leave-one-out, so a caller is
#' never compared against a consensus it helped define -- that would reward
#' typical callers and penalise outliers inconsistently, in proportion to how
#' much each contributed to the thing they are scored against.
#'
#' @section Over- and under-calling:
#' Two of the returned quantities are directional and should be read as a pair.
#' `loo_consensus_captured` is the share of the peers' consensus this caller
#' recovers; `loo_consensus_supported` is the share of its own calls the
#' consensus backs. A caller that is high on the first and low on the second is
#' over-calling -- it finds everything the peers found and much besides. The
#' reverse is under-calling. `loo_consensus_overlap` merges both into one
#' symmetric number and cannot distinguish them.
#'
#' @section Why this does not use findConsensus:
#' [cnvstrmliner::findConsensus()] and [cnvstrmliner::ensemblCall()] merge
#' segments and fill gaps, which changes the geometry: the resulting consensus
#' is not the set of regions the peers agreed on but a post-processed
#' segmentation of it, and comparing against it would silently mix "do the
#' callers agree" with "does merging help". Leave-one-out also means building it
#' once per caller, which at capture resolution over a large sweep is minutes
#' per sample against milliseconds for a bin-vote matrix. The bin vote is used
#' instead, and deliberately no alternative code path is offered.
#'
#' @inheritSection computeCallerAgreement What high agreement does and does not mean
#' @inheritSection computeCallerAgreement Vocabulary
#' @inheritSection computeSegmentationARI Choosing the bins
#'
#' @inheritParams computeCallerAgreement
#'
#' @return A `data.frame` with one row per caller:
#' \describe{
#'   \item{method}{Caller name.}
#'   \item{loo_consensus_overlap}{Jaccard of this caller's altered bin set
#'     against the leave-one-out majority vote.}
#'   \item{loo_mean_pairwise}{Mean overlap against each other caller taken
#'     individually. Less brittle than the consensus when the peers disagree
#'     among themselves.}
#'   \item{loo_consensus_captured}{Share of the consensus altered genome this
#'     caller recovers.}
#'   \item{loo_consensus_supported}{Share of this caller's altered calls the
#'     consensus backs.}
#'   \item{altered_rate}{Fraction of bins this caller calls altered, for
#'     context: the two directional scores are uninterpretable without it.}
#'   \item{n_peer_methods}{Callers in the leave-one-out consensus.}
#' }
#'
#' @seealso [computeCallerAgreement()], [computeBreakpointConsistency()]
#'
#' @examples
#' \dontrun{
#' loo <- computeConsensusAgreement(caller_list, bins = targets)
#' loo[order(-loo$loo_consensus_overlap), c("method", "loo_consensus_overlap")]
#' }
#'
#' @export
computeConsensusAgreement <- function(cnv_list,
                                      bins = NULL,
                                      bin_width = 1e6,
                                      statecol = "cnv.annotation",
                                      neutral_states = c("NEUT", "INFER"),
                                      genome = NULL) {
  prepared <- .prepareAlteredMatrix(
    cnv_list, bins, bin_width, statecol,
    neutral_states, genome
  )
  altered <- prepared$altered
  method_names <- colnames(altered)
  n_methods <- ncol(altered)

  out_cols <- c(
    "loo_consensus_overlap", "loo_mean_pairwise",
    "loo_consensus_captured", "loo_consensus_supported",
    "altered_rate", "n_peer_methods"
  )

  if (n_methods < 2 || nrow(altered) == 0) {
    logger::log_warn("computeConsensusAgreement: fewer than two callers or no bins, returning NA")
    out <- .emptyRow(out_cols)
    out <- out[rep(1, max(n_methods, 0)), , drop = FALSE]
    if (nrow(out) == 0) {
      out <- out[0, , drop = FALSE]
    }
    out <- cbind(method = method_names, out, stringsAsFactors = FALSE)
    rownames(out) <- NULL
    return(out)
  }

  votes <- rowSums(altered)

  rows <- lapply(seq_len(n_methods), function(i) {
    self <- altered[, i]
    peer_votes <- votes - self
    n_peers <- n_methods - 1L
    # Strict majority of the peers, so a two-peer tie is not counted as
    # consensus: it is precisely the case where there is none.
    consensus <- peer_votes > n_peers / 2

    n_self <- sum(self)
    n_cons <- sum(consensus)
    n_shared <- sum(self & consensus)

    data.frame(
      method = method_names[i],
      loo_consensus_overlap = .jaccardLogical(self, consensus),
      loo_mean_pairwise = mean(vapply(
        setdiff(seq_len(n_methods), i),
        function(j) .jaccardLogical(self, altered[, j]),
        numeric(1)
      )),
      loo_consensus_captured = .safeDiv(n_shared, n_cons),
      loo_consensus_supported = .safeDiv(n_shared, n_self),
      altered_rate = mean(self),
      n_peer_methods = n_peers,
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}


# ===== 5. ROLLUPS AND THE ONE-CALL DRIVER ====================================
#
# The rollups take the per-segment tables of section 3, already coerced to
# data.frames, and reduce them to one row per group. They deliberately sit on a
# data.frame boundary with no dependency on the GRanges layer: that is what
# makes them testable on three lines of literal input and reusable against
# per-segment tables written by an earlier version of the pipeline, through
# their column-name arguments.
#
# segmentationMetrics() is the function analysis code should normally call. Its
# one hard contract is that the returned column set never depends on the input,
# so binding a sweep's worth of results together cannot produce a ragged frame.
# =============================================================================


#' Roll per-segment Jaccard columns up to one row per group
#'
#' @description
#' Summarises the per-segment overlap scores produced by
#' [annotateBinJaccard()] or [cnvstrmliner::compareSegmentOverlap()].
#'
#' @section Why five numbers per column:
#' A mean alone is a poor summary of this distribution, which is typically
#' bimodal -- segments that matched well and segments that did not match at all.
#' The median resists the tail; the weighted mean says whether the failures are
#' concentrated in small segments or large ones; `frac_ge_threshold` counts how
#' many segments are usable rather than how good the average is; and
#' `frac_no_overlap` isolates the segments with no counterpart, which are a
#' different failure from a poor match and are invisible in every other
#' statistic here.
#'
#' @section NA handling:
#' A segment with no counterpart carries `NA`. For the averages that is a
#' genuine zero -- nothing was matched -- so it is treated as one, and separately
#' reported as `frac_no_overlap`. Dropping those rows instead would score a
#' caller that invents fifty spurious segments identically to one that invents
#' none.
#'
#' @param seg_df A `data.frame` of per-segment results.
#' @param group_cols Character vector of grouping column names, e.g.
#'   `c("sample", "purity", "method")`.
#' @param value_cols Character vector of per-segment Jaccard columns to
#'   summarise. Defaults to the two this file and the package produce. Names
#'   that are absent from `seg_df` are skipped with a warning, so a table
#'   carrying only one of them still works.
#' @param weight_col Character. Column used to weight the weighted mean, usually
#'   segment width.
#' @param threshold Numeric. Cut for `frac_ge_threshold`. Default `0.5`.
#' @param prefix Character prepended to every output column, for joining several
#'   rollups of different provenance side by side.
#'
#' @return A `data.frame` with one row per group: the grouping columns, then for
#'   each value column a `median_`, `mean_`, `wmean_`, `frac_ge_` and
#'   `frac_no_overlap_` variant named after it.
#'
#' @seealso [annotateBinJaccard()], [summariseStateConcordance()]
#'
#' @examples
#' \dontrun{
#' seg_df <- as.data.frame(cmp$query)
#' summariseSegmentJaccard(seg_df, group_cols = c("sample", "method"))
#' }
#'
#' @export
summariseSegmentJaccard <- function(seg_df,
                                    group_cols,
                                    value_cols = c("jaccard.index", "bin.jaccard.index"),
                                    weight_col = "width",
                                    threshold = 0.5,
                                    prefix = "") {
  if (!is.data.frame(seg_df)) {
    stop("seg_df must be a data.frame")
  }
  missing_group <- setdiff(group_cols, colnames(seg_df))
  if (length(missing_group) > 0) {
    stop(sprintf(
      "seg_df is missing grouping columns: %s",
      paste(missing_group, collapse = ", ")
    ))
  }

  present <- intersect(value_cols, colnames(seg_df))
  absent <- setdiff(value_cols, present)
  if (length(absent) > 0) {
    logger::log_warn("summariseSegmentJaccard: skipping absent columns {paste(absent, collapse=', ')}")
  }
  if (length(present) == 0) {
    stop("none of value_cols is present in seg_df")
  }
  if (!(weight_col %in% colnames(seg_df))) {
    stop(sprintf("seg_df is missing the weight column '%s'", weight_col))
  }

  # A short suffix per value column, so jaccard.index and bin.jaccard.index do
  # not collide once prefixed.
  suffix_of <- function(x) gsub("\\.index$", "", gsub("\\.", "_", x))

  # The weight is copied to a fixed name so the across() lambdas can refer to it
  # as a bare symbol; the .data pronoun would need rlang attached, which a
  # source()-able file cannot assume.
  work <- seg_df
  work[["tmp_weight"]] <- tidyr::replace_na(as.numeric(work[[weight_col]]), 0)

  grouped <- dplyr::group_by(work, dplyr::across(dplyr::all_of(group_cols)))

  out <- dplyr::summarise(
    grouped,
    dplyr::across(
      dplyr::all_of(present),
      list(
        median = ~ stats::median(tidyr::replace_na(.x, 0), na.rm = TRUE),
        mean = ~ mean(tidyr::replace_na(.x, 0), na.rm = TRUE),
        wmean = ~ stats::weighted.mean(
          tidyr::replace_na(.x, 0),
          w = tmp_weight, na.rm = TRUE
        ),
        frac_ge = ~ mean(tidyr::replace_na(.x, 0) >= threshold, na.rm = TRUE),
        frac_no_overlap = ~ mean(is.na(.x))
      ),
      .names = "{.fn}_{.col}"
    ),
    n_segments = dplyr::n(),
    .groups = "drop"
  )
  out <- as.data.frame(out)

  metric_cols <- setdiff(colnames(out), group_cols)
  renamed <- metric_cols
  for (col in present) {
    renamed <- gsub(col, suffix_of(col), renamed, fixed = TRUE)
  }
  colnames(out)[match(metric_cols, colnames(out))] <-
    if (nzchar(prefix)) paste0(prefix, renamed) else renamed

  out
}


#' Roll per-segment state concordance up to one row per group
#'
#' @description
#' Summarises the columns [annotateBinStateConcordance()] attaches.
#'
#' @section Two accuracies, not one:
#' `pooled` sums the matched and total bin counters across segments before
#' dividing, giving the accuracy of the segmentation as a whole -- every bin
#' counts once, which is what a downstream user of the calls experiences.
#' `mean` averages the per-segment fractions, giving every segment equal say
#' regardless of size, which is what tells you whether the *segments* are
#' individually trustworthy.
#'
#' They diverge sharply and neither is a substitute for the other. A caller with
#' three enormous correct segments and forty tiny wrong ones scores well on the
#' first and badly on the second, and both readings are true.
#'
#' @param seg_df A `data.frame` of per-segment results.
#' @param group_cols Character vector of grouping column names.
#' @param matched_col,total_col,value_col Column names holding the matched-bin
#'   count, the total-bin count and the per-segment fraction. Default to the
#'   names [annotateBinStateConcordance()] writes; point them at
#'   `exon_state_n_matched` and friends to use a table from an earlier pipeline.
#' @param weight_col Character. Column used for the weighted mean, usually
#'   segment width.
#' @param prefix Character prepended to every output column.
#'
#' @return A `data.frame` with one row per group and columns
#'   `state_concordance_pooled`, `state_concordance_mean`,
#'   `state_concordance_wmean`, `state_concordance_binw`, `n_bins_scored` and
#'   `n_segments`.
#'
#' @seealso [annotateBinStateConcordance()], [computeStateKappa()]
#'
#' @examples
#' \dontrun{
#' summariseStateConcordance(seg_df, group_cols = c("sample", "purity", "method"))
#' }
#'
#' @export
summariseStateConcordance <- function(seg_df,
                                      group_cols,
                                      matched_col = "state.n_matched",
                                      total_col = "state.n_total",
                                      value_col = "state.concordance",
                                      weight_col = "width",
                                      prefix = "") {
  if (!is.data.frame(seg_df)) {
    stop("seg_df must be a data.frame")
  }
  required <- c(group_cols, matched_col, total_col, value_col, weight_col)
  missing_cols <- setdiff(required, colnames(seg_df))
  if (length(missing_cols) > 0) {
    stop(sprintf(
      "seg_df is missing required columns: %s",
      paste(missing_cols, collapse = ", ")
    ))
  }

  # Fixed internal names, so the summarise() below can use bare symbols rather
  # than the .data pronoun, which would need rlang attached.
  work <- seg_df
  work[["tmp_matched"]] <- as.numeric(work[[matched_col]])
  work[["tmp_total"]] <- as.numeric(work[[total_col]])
  work[["tmp_value"]] <- as.numeric(work[[value_col]])
  work[["tmp_weight"]] <- tidyr::replace_na(as.numeric(work[[weight_col]]), 0)

  out <- dplyr::summarise(
    dplyr::group_by(work, dplyr::across(dplyr::all_of(group_cols))),
    state_concordance_pooled = .safeDiv(
      sum(tmp_matched, na.rm = TRUE), sum(tmp_total, na.rm = TRUE)
    ),
    state_concordance_mean = mean(tmp_value, na.rm = TRUE),
    state_concordance_wmean = stats::weighted.mean(
      tmp_value,
      w = tmp_weight, na.rm = TRUE
    ),
    state_concordance_binw = stats::weighted.mean(
      tmp_value,
      w = tidyr::replace_na(tmp_total, 0), na.rm = TRUE
    ),
    n_bins_scored = sum(tmp_total, na.rm = TRUE),
    n_segments = dplyr::n(),
    .groups = "drop"
  )
  out <- as.data.frame(out)

  if (nzchar(prefix)) {
    metric_cols <- setdiff(colnames(out), group_cols)
    colnames(out)[match(metric_cols, colnames(out))] <- paste0(prefix, metric_cols)
  }
  out
}


#' Roll a breakpoint consistency table up to scalars
#'
#' @description
#' Reduces the per-cluster output of [computeBreakpointConsistency()] to a
#' single row, optionally within groups.
#'
#' @param consistency_df A `data.frame` from [computeBreakpointConsistency()].
#' @param group_cols Character vector of grouping columns. Empty by default,
#'   giving one row overall.
#' @param high_threshold Numeric. Consistency above which a cluster counts as
#'   well supported. Default `0.3`.
#' @param prefix Character prepended to every output column.
#'
#' @return A `data.frame` with `mean_consistency`, `median_consistency`,
#'   `frac_high`, `frac_unanimous` (clusters every segmentation votes for),
#'   `n_breakpoints` and `high_threshold`.
#'
#' @seealso [computeBreakpointConsistency()], [annotateBreakpointConsistency()]
#'
#' @examples
#' \dontrun{
#' summariseBreakpointConsistency(computeBreakpointConsistency(facets_settings))
#' }
#'
#' @export
summariseBreakpointConsistency <- function(consistency_df,
                                           group_cols = character(),
                                           high_threshold = 0.3,
                                           prefix = "") {
  if (!is.data.frame(consistency_df)) {
    stop("consistency_df must be a data.frame")
  }
  if (!("consistency" %in% colnames(consistency_df))) {
    stop("consistency_df must include column 'consistency'")
  }
  missing_group <- setdiff(group_cols, colnames(consistency_df))
  if (length(missing_group) > 0) {
    stop(sprintf(
      "consistency_df is missing grouping columns: %s",
      paste(missing_group, collapse = ", ")
    ))
  }

  if (nrow(consistency_df) == 0) {
    out <- .emptyRow(c(
      "mean_consistency", "median_consistency", "frac_high",
      "frac_unanimous", "n_breakpoints"
    ))
    out$n_breakpoints <- 0L
    out$high_threshold <- high_threshold
    return(out)
  }

  grouped <- if (length(group_cols) > 0) {
    dplyr::group_by(consistency_df, dplyr::across(dplyr::all_of(group_cols)))
  } else {
    consistency_df
  }

  out <- dplyr::summarise(
    grouped,
    mean_consistency = mean(consistency, na.rm = TRUE),
    median_consistency = stats::median(consistency, na.rm = TRUE),
    frac_high = mean(consistency > high_threshold, na.rm = TRUE),
    frac_unanimous = mean(consistency >= 1, na.rm = TRUE),
    n_breakpoints = dplyr::n(),
    .groups = "drop"
  )
  out <- as.data.frame(out)
  out$high_threshold <- high_threshold

  if (nzchar(prefix)) {
    metric_cols <- setdiff(colnames(out), group_cols)
    colnames(out)[match(metric_cols, colnames(out))] <- paste0(prefix, metric_cols)
  }
  out
}


#' All genome-wide comparison metrics for one pair of segmentations
#'
#' @description
#' Runs the section 1 and section 2 metrics for a single query/subject pair and
#' flattens the results into one row.
#'
#' This is the function analysis code should normally call. Computing the
#' metrics one at a time is not wrong, but it repeats the bin assignment several
#' times and, more importantly, it leaves each caller to decide what to do when
#' a metric fails -- which is how a sweep ends up with a result table whose
#' columns depend on which settings happened to work.
#'
#' @section The column-set contract:
#' The returned column set does not depend on the input. Every requested metric
#' contributes its columns whether it succeeded, failed or was skipped; a
#' failure logs a warning and contributes `NA`. An empty segmentation, a
#' single-segment one and a real pair all return the same one-row shape, so
#' `dplyr::bind_rows()` over a whole sweep can never produce a ragged frame.
#' This is the property the downstream pipeline depends on, and it is the main
#' reason to prefer this function over assembling the calls by hand.
#'
#' @section What is not here:
#' The per-segment annotators of section 3 and the truth-free metrics of section
#' 4 are excluded. The first return a segmentation rather than a row, and the
#' second take a family of segmentations rather than a pair. Call them
#' separately.
#'
#' @section Bins:
#' A single `bins` argument serves both kinds of metric. When it is `NULL` the
#' state metrics -- which require a real grid -- are skipped with a warning, and
#' only the geometry metrics run on an auto-tiled grid. Pass the capture targets
#' to get everything.
#'
#' @param query A `cnvstrml` object or `GRanges` -- the segmentation scored.
#' @param subject A `cnvstrml` object or `GRanges` -- the reference.
#' @param bins A `GRanges`, a BED file path, or `NULL`. See the *Bins* section.
#' @param bin_width Numeric width of auto-tiled bins, used only where `bins` is
#'   `NULL`. Default `1e6`.
#' @param tolerance_bp Numeric. Breakpoint matching tolerance. Default `3e5`.
#' @param statecol Character. Name of the state column.
#' @param simplify,infer_as_neutral Logical. State canonicalisation; see
#'   [computeStateKappa()].
#' @param genome Character genome identifier for auto-tiling.
#' @param metrics Character vector selecting which metric groups to run. Any of
#'   `"kappa"`, `"recovery"`, `"boundary"`, `"ari"`, `"vi"`, `"fragmentation"`.
#'   Unselected groups still contribute their columns, as `NA`.
#'
#' @return A one-row `data.frame`. ARI and VI appear twice each, as
#'   `ari_segment` / `ari_state` and `vi_segment` / `vi_state`, never as a bare
#'   `ari` or `vi`: the two partitions answer different questions and a column
#'   name that does not say which is a trap for whoever reads the table next.
#'
#' @seealso [computeStateKappa()], [computeEventRecovery()],
#'   [computeBoundaryConcordance()], [computeSegmentationARI()],
#'   [computeSegmentationVI()], [computeFragmentation()]
#'
#' @examples
#' \dontrun{
#' rows <- lapply(names(facets_settings), function(m) {
#'   out <- segmentationMetrics(facets_settings[[m]], truth, bins = targets)
#'   cbind(method = m, out)
#' })
#' dplyr::bind_rows(rows) # guaranteed rectangular
#' }
#'
#' @export
segmentationMetrics <- function(query,
                                subject,
                                bins = NULL,
                                bin_width = 1e6,
                                tolerance_bp = 3e5,
                                statecol = "cnv.annotation",
                                simplify = FALSE,
                                infer_as_neutral = TRUE,
                                genome = NULL,
                                metrics = c(
                                  "kappa", "recovery", "boundary",
                                  "ari", "vi", "fragmentation"
                                )) {
  metrics <- match.arg(metrics, several.ok = TRUE)

  out_cols <- c(
    "kappa", "kappa_p_o", "kappa_p_e", "n_bins_scored", "n_bins_unscored",
    "recall_altered", "precision_altered", "sorensen",
    "n_truth_altered", "n_query_altered",
    "boundary_f1", "boundary_precision", "boundary_recall",
    "n_query_breakpoints", "n_subject_breakpoints",
    "ari_segment", "ri_segment", "ari_state",
    "vi_segment", "nvi_segment", "vi_state", "nvi_state",
    "split_ratio", "merge_ratio",
    "mean_query_per_subject", "mean_subject_per_query",
    "n_query_segments", "n_subject_segments"
  )
  out <- .emptyRow(out_cols)

  # A metric that blows up must not take the rest of the row with it: the
  # column-set contract is worth more than the individual value.
  attempt <- function(label, expr) {
    tryCatch(expr, error = function(e) {
      logger::log_warn("segmentationMetrics: {label} failed: {conditionMessage(e)}")
      NULL
    })
  }

  bins_available <- !is.null(bins)
  if (!bins_available && any(c("kappa", "recovery") %in% metrics)) {
    logger::log_warn(
      "segmentationMetrics: bins is NULL, skipping the state metrics that require a real grid"
    )
  }

  if ("kappa" %in% metrics && bins_available) {
    k <- attempt("kappa", computeStateKappa(
      query, subject, bins, statecol, simplify, infer_as_neutral
    ))
    if (!is.null(k)) {
      out$kappa <- k$kappa
      out$kappa_p_o <- k$p_o
      out$kappa_p_e <- k$p_e
      out$n_bins_scored <- k$n_bins_scored
      out$n_bins_unscored <- k$n_bins_unscored
    }
  }

  if ("recovery" %in% metrics && bins_available) {
    r <- attempt("event recovery", computeEventRecovery(
      query, subject, bins,
      statecol = statecol, simplify = simplify
    ))
    if (!is.null(r)) {
      out$recall_altered <- r$recall_altered
      out$precision_altered <- r$precision_altered
      out$sorensen <- r$sorensen
      out$n_truth_altered <- r$n_truth_altered
      out$n_query_altered <- r$n_query_altered
    }
  }

  if ("boundary" %in% metrics) {
    b <- attempt("boundary concordance", computeBoundaryConcordance(
      query, subject,
      tolerance_bp = tolerance_bp
    ))
    if (!is.null(b)) {
      out$boundary_f1 <- b$boundary_f1
      out$boundary_precision <- b$boundary_precision
      out$boundary_recall <- b$boundary_recall
      out$n_query_breakpoints <- b$n_query_breakpoints
      out$n_subject_breakpoints <- b$n_subject_breakpoints
    }
  }

  if ("ari" %in% metrics) {
    a_seg <- attempt("segment ARI", computeSegmentationARI(
      query, subject, bins, bin_width, "segment",
      statecol, simplify, infer_as_neutral, genome
    ))
    if (!is.null(a_seg)) {
      out$ari_segment <- a_seg$ari
      out$ri_segment <- a_seg$ri
    }
    a_state <- attempt("state ARI", computeSegmentationARI(
      query, subject, bins, bin_width, "state",
      statecol, simplify, infer_as_neutral, genome
    ))
    if (!is.null(a_state)) {
      out$ari_state <- a_state$ari
    }
  }

  if ("vi" %in% metrics) {
    v_seg <- attempt("segment VI", computeSegmentationVI(
      query, subject, bins, bin_width, "segment",
      statecol, simplify, infer_as_neutral, genome
    ))
    if (!is.null(v_seg)) {
      out$vi_segment <- v_seg$vi
      out$nvi_segment <- v_seg$nvi
    }
    v_state <- attempt("state VI", computeSegmentationVI(
      query, subject, bins, bin_width, "state",
      statecol, simplify, infer_as_neutral, genome
    ))
    if (!is.null(v_state)) {
      out$vi_state <- v_state$vi
      out$nvi_state <- v_state$nvi
    }
  }

  if ("fragmentation" %in% metrics) {
    f <- attempt("fragmentation", computeFragmentation(
      query, subject, bins, bin_width, genome
    ))
    if (!is.null(f)) {
      out$split_ratio <- f$split_ratio
      out$merge_ratio <- f$merge_ratio
      out$mean_query_per_subject <- f$mean_query_per_subject
      out$mean_subject_per_query <- f$mean_subject_per_query
      out$n_query_segments <- f$n_query_segments
      out$n_subject_segments <- f$n_subject_segments
    }
  }

  out
}
