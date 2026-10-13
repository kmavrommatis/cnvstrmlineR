# Metrics comparing individual segments between two segmentation solutions.



#' Annotate segments with per-segment state concordance against a reference
#'
#' @description
#' For each segment, the fraction of its bins whose copy-number state matches
#' the state a reference segmentation assigns to the same bins.
#'
#' This is the state-aware, size-neutral per-segment score. A short segment
#' whose bins are all correctly stated scores 1; a segment that overlaps its
#' counterpart perfectly but calls the wrong copy number scores near 0. It is
#' the complement of [annotateBinJaccard()], which rewards positional overlap
#' and is indifferent to state.
#'
#' The counters, not just the fraction, are attached, because they are what
#' aggregates correctly: pooling `state.n_matched` over `state.n_total` across
#' segments gives the exon-weighted accuracy of the whole segmentation, whereas
#' averaging the per-segment fractions weights a two-bin segment as heavily as a
#' two-thousand-bin one. Both are useful and they are not the same number; see
#' [summariseStateConcordance()], which reports each.
#'
#' @section Direction:
#' The function is directional and both directions are worth computing.
#' Annotating the prediction against the truth gives something that behaves like
#' precision -- of what this segment claims, how much is right. Annotating the
#' truth against the prediction gives recall -- of what was really there, how
#' much was found. The pipeline that produced this project's
#' `*_segments_jaccard_predicted.tsv` and `*_segments_jaccard_truth.tsv` files
#' runs both.
#'
#' @inheritSection computeStateKappa Choosing the bins
#'
#' @section NA handling:
#' A segment with no bin carrying a defined reference state gets
#' `state.concordance = NA` and zero counters -- there was nothing to score it
#' against. That is distinct from a concordance of 0, which means the segment
#' had bins and got all of them wrong, and the two must not be conflated when
#' aggregating: see [summariseStateConcordance()].
#'
#' @param segments A `cnvstrml` object or `GRanges` -- the segmentation to
#'   annotate.
#' @param reference A `cnvstrml` object or `GRanges` providing the states to
#'   compare against.
#' @param bins A `GRanges` of bin-level regions, or a BED file path. Required.
#' @param statecol Character. Name of the state column, in both arguments.
#' @param simplify,infer_as_neutral Logical. State canonicalisation; see
#'   [computeStateKappa()].
#' @param reference_state Optional character vector of length `length(bins)`
#'   giving the reference state already assigned to each bin, as returned by
#'   [assignBinState()]. Purely an optimisation: when the same reference is
#'   scored against many segmentations the projection can be computed once and
#'   reused, which is the difference between one bin assignment per sweep and
#'   one per setting.
#'
#' @return The segments as a `GRanges` with three added `mcols`:
#' \describe{
#'   \item{state.concordance}{Fraction of the segment's scorable bins whose
#'     state matches the reference, or `NA`.}
#'   \item{state.n_matched}{Matching bins.}
#'   \item{state.n_total}{Scorable bins assigned to the segment -- the
#'     denominator.}
#' }
#'
#' @seealso [assignBinState()], [annotateBinJaccard()],
#'   [summariseStateConcordance()], [computeEventRecovery()]
#'
#' @examples
#' \dontrun{
#' # Hoist the reference projection across a sweep.
#' truth_state <- assignBinState(targets, truth)
#' scored <- lapply(facets_settings, function(fo) {
#'   annotateBinStateConcordance(fo, truth, targets, reference_state = truth_state)
#' })
#' }
#'
#' @export
annotateBinStateConcordance <- function(segments,
                                        reference,
                                        bins,
                                        statecol = "cnv.annotation",
                                        simplify = FALSE,
                                        infer_as_neutral = TRUE,
                                        reference_state = NULL) {
  seg_gr <- .asSegments(segments, "segments", statecol)
  ref_gr <- .asSegments(reference, "reference", statecol)
  bins <- .resolveBins(bins, allow_auto = FALSE)

  if (!is.null(reference_state) && length(reference_state) != length(bins)) {
    stop("reference_state must have one entry per bin")
  }

  n_seg <- length(seg_gr)
  concordance <- rep(NA_real_, n_seg)
  n_matched <- rep(0L, n_seg)
  n_total <- rep(0L, n_seg)

  if (n_seg > 0 && length(bins) > 0) {
    bin_seg <- .bestBinSegment(bins, seg_gr)
    ref_state <- if (is.null(reference_state)) {
      .binStateVector(bins, ref_gr, statecol, simplify, infer_as_neutral)
    } else {
      .canonicalStates(reference_state, simplify = simplify,
                       infer_as_neutral = infer_as_neutral)
    }

    seg_state_all <- .canonicalStates(
      S4Vectors::mcols(seg_gr)[[statecol]],
      simplify = simplify, infer_as_neutral = infer_as_neutral
    )

    pred_state <- rep(NA_character_, length(bins))
    has_seg <- !is.na(bin_seg)
    pred_state[has_seg] <- seg_state_all[bin_seg[has_seg]]

    scorable <- has_seg & !is.na(ref_state)
    if (any(scorable)) {
      seg_of_bin <- bin_seg[scorable]
      is_match <- pred_state[scorable] == ref_state[scorable]

      total_by_seg <- tapply(is_match, seg_of_bin, length)
      match_by_seg <- tapply(is_match, seg_of_bin, sum)
      seg_ids <- as.integer(names(total_by_seg))

      n_total[seg_ids] <- as.integer(total_by_seg)
      n_matched[seg_ids] <- as.integer(match_by_seg)
      concordance[seg_ids] <- as.numeric(match_by_seg) / as.numeric(total_by_seg)
    }
  }

  S4Vectors::mcols(seg_gr)$state.concordance <- concordance
  S4Vectors::mcols(seg_gr)$state.n_matched <- n_matched
  S4Vectors::mcols(seg_gr)$state.n_total <- n_total
  seg_gr
}


#' Annotate segments with their best bin-set Jaccard against another segmentation
#'
#' @description
#' For every segment on both sides, the Jaccard overlap -- measured in bins, not
#' base pairs -- with its best-matching counterpart, together with the state of
#' that counterpart and the three bin counts the score is built from.
#'
#' Measuring in bins rather than base pairs is the point. On capture data the
#' base-pair extent of a segment is dominated by the unsequenced gaps it spans,
#' so a base-pair Jaccard behaves largely as a segment-length metric and rewards
#' a caller for placing boundaries in deserts where no evidence exists. Counting
#' bins measures the evidence that was actually available.
#'
#' @section Reading the bin counts:
#' `bin.jaccard.index` alone cannot distinguish the two ways of scoring low, and
#' they call for opposite remedies. The three counters resolve it:
#'
#' - `n_self` much greater than `n_match`: this segment swallowed several
#'   counterparts -- a merge, wanting more sensitivity.
#' - `n_self` much less than `n_match`: this segment is one shard of a
#'   counterpart -- a split, wanting less.
#' - both similar with a small `n_shared`: the segments are simply in different
#'   places.
#'
#' @section Differences from the project-local compareSegmentOverlapRegions:
#' The bin counts here come from single assignment -- each bin belongs to
#' exactly one segment, the one it overlaps most -- so `n_self + n_match -
#' n_shared` is a true set union. The project-local implementation counts with
#' `GenomicRanges::countOverlaps()`, which counts a bin straddling a segment
#' boundary in both segments and inflates the union, deflating the score. The
#' two agree closely in practice, since straddling bins are a small minority,
#' and this implementation is never the larger of the two. Do not expect exact
#' equality when diffing them.
#'
#' @inheritSection computeStateKappa Choosing the bins
#'
#' @section NA handling:
#' A segment with no overlapping counterpart gets `bin.jaccard.index = NA` and
#' `bin.jaccard.status = "NO_OVERLAP"`. For per-segment inspection that `NA` is
#' correct -- there is no score. For an average over segments it is a genuine
#' zero, since nothing was matched, and [summariseSegmentJaccard()] treats it as
#' one while also reporting how often it happened.
#'
#' @param query A `cnvstrml` object or `GRanges`.
#' @param subject A `cnvstrml` object or `GRanges`.
#' @param bins A `GRanges` of bin-level regions, or a BED file path. Required.
#' @param statecol Character. Name of the state column, used to report which
#'   state the winning counterpart carried.
#' @param query_bin_counts,subject_bin_counts Optional integer vectors of bins
#'   per segment, one entry per segment on that side. An optimisation for sweeps
#'   where one side is fixed; recomputed when absent.
#'
#' @return A named `list` with elements `query` and `subject`, each the input
#'   `GRanges` with five added `mcols`: `bin.jaccard.index`,
#'   `bin.jaccard.status`, `bin.jaccard.n_self`, `bin.jaccard.n_match` and
#'   `bin.jaccard.n_shared`.
#'
#' @seealso [annotateBinStateConcordance()] for the state-aware companion,
#'   [cnvstrmliner::compareSegmentOverlap()] for the base-pair version,
#'   [summariseSegmentJaccard()] for the rollup.
#'
#' @examples
#' \dontrun{
#' cmp <- annotateBinJaccard(facets_obj, truth, bins = targets)
#' # Segments the caller invented.
#' sum(is.na(S4Vectors::mcols(cmp$query)$bin.jaccard.index))
#' }
#'
#' @export
annotateBinJaccard <- function(query,
                               subject,
                               bins,
                               statecol = "cnv.annotation",
                               query_bin_counts = NULL,
                               subject_bin_counts = NULL) {
  query_gr <- .asSegments(query, "query", statecol)
  subject_gr <- .asSegments(subject, "subject", statecol)
  bins <- .resolveBins(bins, allow_auto = FALSE)

  if (!is.null(query_bin_counts) && length(query_bin_counts) != length(query_gr)) {
    stop("query_bin_counts must have one entry per query segment")
  }
  if (!is.null(subject_bin_counts) && length(subject_bin_counts) != length(subject_gr)) {
    stop("subject_bin_counts must have one entry per subject segment")
  }

  # Single assignment on each side, so the per-segment bin sets are disjoint.
  query_assign <- .bestBinSegment(bins, query_gr)
  subject_assign <- .bestBinSegment(bins, subject_gr)

  if (is.null(query_bin_counts)) {
    query_bin_counts <- tabulate(query_assign, nbins = length(query_gr))
  }
  if (is.null(subject_bin_counts)) {
    subject_bin_counts <- tabulate(subject_assign, nbins = length(subject_gr))
  }

  # Bins assigned on both sides give the pairwise intersections directly.
  paired <- !is.na(query_assign) & !is.na(subject_assign)
  shared <- if (any(paired)) {
    as.data.frame(
      table(
        query_idx = query_assign[paired],
        subject_idx = subject_assign[paired]
      ),
      responseName = "n_shared", stringsAsFactors = FALSE
    )
  } else {
    data.frame(query_idx = character(), subject_idx = character(),
               n_shared = numeric(), stringsAsFactors = FALSE)
  }
  shared <- shared[shared$n_shared > 0, , drop = FALSE]
  shared$query_idx <- as.integer(shared$query_idx)
  shared$subject_idx <- as.integer(shared$subject_idx)

  query_state <- as.character(S4Vectors::mcols(query_gr)[[statecol]])
  subject_state <- as.character(S4Vectors::mcols(subject_gr)[[statecol]])

  query_best <- .bestBinJaccard(
    shared$query_idx, shared$subject_idx, shared$n_shared,
    query_bin_counts, subject_bin_counts, subject_state, length(query_gr)
  )
  subject_best <- .bestBinJaccard(
    shared$subject_idx, shared$query_idx, shared$n_shared,
    subject_bin_counts, query_bin_counts, query_state, length(subject_gr)
  )

  attach_cols <- function(gr, best) {
    S4Vectors::mcols(gr)$bin.jaccard.index <- best$index
    S4Vectors::mcols(gr)$bin.jaccard.status <- best$status
    S4Vectors::mcols(gr)$bin.jaccard.n_self <- best$n_self
    S4Vectors::mcols(gr)$bin.jaccard.n_match <- best$n_match
    S4Vectors::mcols(gr)$bin.jaccard.n_shared <- best$n_shared
    gr
  }

  list(
    query = attach_cols(query_gr, query_best),
    subject = attach_cols(subject_gr, subject_best)
  )
}


#' Annotate segments with the recurrence of their boundaries
#'
#' @description
#' Attaches, to each segment, how consistently its two boundaries are recovered
#' across a family of segmentations -- typically the settings of a sensitivity
#' sweep, or a panel of callers.
#'
#' The reason to want this without a truth: a boundary that almost every setting
#' finds is a feature of the data, and one that a single setting finds is
#' usually a feature of that setting. A segmentation built mostly out of
#' idiosyncratic boundaries is over-fitted to its parameters, and its mean
#' consistency says so before any truth is consulted.
#'
#' @section Caveat:
#' Recurrence is agreement, not correctness. If most settings share a failure
#' mode -- and settings of the same algorithm frequently do -- a wrong boundary
#' can be perfectly consistent. Treat a high score as the absence of one
#' specific defect rather than as evidence of quality, and read it alongside a
#' truth-based score where one exists.
#'
#' @param segments A `cnvstrml` object or `GRanges` to annotate.
#' @param consistency A `data.frame` from [computeBreakpointConsistency()].
#' @param tolerance_bp Numeric. Half-width of the window within which a segment
#'   boundary is considered to be at a consistency cluster. Should match the
#'   value used to build `consistency`.
#' @param high_threshold Numeric. Consistency above which a boundary counts as
#'   well supported, for the `frac_high` summary. Default `0.3`.
#'
#' @return A named `list`:
#' \describe{
#'   \item{segments}{The input as `GRanges` with `breakpoint.consistency.start`
#'     and `breakpoint.consistency.end` added. `NA` where a boundary matches no
#'     cluster, which happens only for boundaries excluded when the clusters
#'     were built -- ordinarily the sequence ends.}
#'   \item{summary}{A one-row `data.frame` with `mean_consistency`, `frac_high`,
#'     `n_breakpoints` and `high_threshold`.}
#' }
#'
#' @seealso [computeBreakpointConsistency()], [computeBoundaryConcordance()]
#'
#' @examples
#' \dontrun{
#' consistency <- computeBreakpointConsistency(facets_settings)
#' ann <- annotateBreakpointConsistency(facets_settings[[1]], consistency)
#' ann$summary$mean_consistency
#' }
#'
#' @export
annotateBreakpointConsistency <- function(segments,
                                          consistency,
                                          tolerance_bp = 3e5,
                                          high_threshold = 0.3) {
  seg_gr <- .asSegments(segments, "segments")

  empty_summary <- data.frame(
    mean_consistency = NA_real_, frac_high = NA_real_,
    n_breakpoints = 0L, high_threshold = high_threshold,
    stringsAsFactors = FALSE
  )

  if (length(seg_gr) == 0) {
    return(list(segments = seg_gr, summary = empty_summary))
  }

  if (!is.data.frame(consistency)) {
    stop("consistency must be a data.frame from computeBreakpointConsistency()")
  }

  if (nrow(consistency) == 0) {
    S4Vectors::mcols(seg_gr)$breakpoint.consistency.start <- NA_real_
    S4Vectors::mcols(seg_gr)$breakpoint.consistency.end <- NA_real_
    return(list(segments = seg_gr, summary = empty_summary))
  }

  missing_cols <- setdiff(c("seqnames", "pos", "consistency"), colnames(consistency))
  if (length(missing_cols) > 0) {
    stop(sprintf(
      "consistency must include column%s %s",
      if (length(missing_cols) > 1) "s" else "",
      paste0("'", missing_cols, "'", collapse = ", ")
    ))
  }

  cluster_gr <- GenomicRanges::GRanges(
    seqnames = consistency$seqnames,
    ranges = IRanges::IRanges(
      start = pmax(consistency$pos - tolerance_bp, 1),
      end = consistency$pos + tolerance_bp
    )
  )

  score_at <- function(pos) {
    edge_gr <- GenomicRanges::GRanges(
      seqnames = GenomicRanges::seqnames(seg_gr),
      ranges = IRanges::IRanges(start = pos, width = 1)
    )
    out <- rep(NA_real_, length(seg_gr))
    hits <- GenomicRanges::findOverlaps(edge_gr, cluster_gr, ignore.strand = TRUE)
    if (length(hits) == 0) {
      return(out)
    }
    # A boundary can sit inside more than one cluster window; the best-supported
    # cluster is the one that describes it.
    qh <- S4Vectors::queryHits(hits)
    value <- consistency$consistency[S4Vectors::subjectHits(hits)]
    best <- tapply(value, qh, max)
    out[as.integer(names(best))] <- as.numeric(best)
    out
  }

  consistency_start <- score_at(GenomicRanges::start(seg_gr))
  consistency_end <- score_at(GenomicRanges::end(seg_gr))

  S4Vectors::mcols(seg_gr)$breakpoint.consistency.start <- consistency_start
  S4Vectors::mcols(seg_gr)$breakpoint.consistency.end <- consistency_end

  all_values <- c(consistency_start, consistency_end)
  all_values <- all_values[!is.na(all_values)]

  summary_df <- data.frame(
    mean_consistency = if (length(all_values) > 0) mean(all_values) else NA_real_,
    frac_high = if (length(all_values) > 0) mean(all_values > high_threshold) else NA_real_,
    n_breakpoints = length(all_values),
    high_threshold = high_threshold,
    stringsAsFactors = FALSE
  )

  list(segments = seg_gr, summary = summary_df)
}
#' Compare segment status between two CNV profiles
#'
#' @description
#' Compares overlapping CNV segments between a query and a subject \code{cnvstrml} object.
#' Overlaps are computed bidirectionally and intersected so only true shared regions are
#' evaluated. For each overlapping region, the function assesses agreement at three levels:
#' overall CNV annotation (e.g., AMP, NEUT, DEL), total copy number, and clonal fraction (CCF),
#' with CCF agreement controlled by the \code{wiggle} tolerance.
#'
#' @param query A \code{cnvstrml} object providing the CNV segmentation to compare.
#' @param subject A \code{cnvstrml} object to serve as the reference or baseline.
#' @param wiggle Numeric tolerance for CCF agreement. Two CCF values are considered the same
#'   if their absolute difference is less than \code{wiggle}. For example, with \code{wiggle = 0.1},
#'   CCF 0.30 and 0.39 are considered equal.
#'
#' @return A named list with five elements describing per-segment comparisons and summaries:
#' \describe{
#'   \item{table}{Data frame of overlapping regions with per-region comparison results.
#'     Columns include: \code{chrom}, \code{start}, \code{end}, \code{query} (annotation),
#'     \code{subject} (annotation), \code{query_CNV}, \code{subject_CNV}, \code{query_CCF},
#'     \code{subject_CCF}, and three comparison labels:
#'     \describe{
#'       \item{status}{Agreement on CNV annotation. \strong{TP} if annotations match;
#'         otherwise \strong{FP} if \code{query == "NEUT"} (query neutral, subject aberrant),
#'         or \strong{FN} when query is aberrant and subject differs.}
#'       \item{cnv}{Agreement on total copy number when \code{status == "TP"}:
#'         \strong{TP} if copy numbers match; \strong{FP} if \code{query_CNV < subject_CNV};
#'         \strong{FN} if \code{query_CNV > subject_CNV}. When \code{status != "TP"},
#'         returns \strong{SAME}/\strong{DIFF} for equality check without FP/FN assignment.}
#'       \item{ccf}{Agreement on clonal fraction using \code{wiggle} when \code{status == "TP"}:
#'         \strong{TP} if \code{|query_CCF - subject_CCF| < wiggle}; \strong{FP} if \code{query_CCF < subject_CCF};
#'         \strong{FN} if \code{query_CCF > subject_CCF}. When \code{status != "TP"},
#'         returns \strong{SAME}/\strong{DIFF} based on the tolerance.}
#'     }}
#'
#'   \item{summary}{Genome-level aggregation excluding sex chromosomes (chrX, chrY).
#'     Groups by \code{query_CNV}, \code{subject_CNV}, and \code{status}, reporting
#'     \code{segments} (count), \code{total_bp} (bases), and \code{genome_fraction}
#'     (fraction of total bases across all overlapping regions). Sorted by \code{total_bp}.}
#'
#'   \item{stats}{Genome-level totals per \code{status} (TP/FP/FN), including \code{segments},
#'     \code{total_bp}, and \code{genome_fraction}. Useful for overall agreement metrics.}
#'
#'   \item{summary_by_chr}{Chromosome-level aggregation (excluding chrX/chrY), grouped by
#'     \code{chrom}, \code{query_CNV}, \code{subject_CNV}, and \code{status}. Adds
#'     \code{chromosome_fraction} = \code{total_bp / chromosome_length} using package
#'     reference lengths. Sorted by chromosome and \code{total_bp}.}
#'
#'   \item{stats_by_chr}{Chromosome-level totals per \code{status}, reporting \code{segments},
#'     \code{total_bp}, and \code{chromosome_fraction}. Provides per-chromosome agreement.}
#' }
#'
#' @details
#' \strong{How to use the outputs in CNV analysis}
#' \itemize{
#'   \item \emph{Annotation agreement (\code{status})}: Quantifies whether broad CNV calls
#'     (AMP/DEL/NEUT) match, supporting validation of pipelines or cross-method comparisons.
#'   \item \emph{Copy-number agreement (\code{cnv})}: Tests granularity beyond annotation;
#'     distinguish cases where both methods call AMP but disagree on total copies.
#'   \item \emph{Clonality agreement (\code{ccf})}: Evaluates consistency of clonality estimates
#'     under a tolerance \code{wiggle}, helpful for assessing subclonal differences.
#'   \item \emph{Genome-wide metrics (\code{summary}, \code{stats})}: Summarize burden and agreement
#'     across the whole genome; \code{genome_fraction} highlights dominant agreement/disagreement.
#'   \item \emph{Per-chromosome metrics (\code{summary_by_chr}, \code{stats_by_chr})}: Reveal
#'     chromosome-specific concordance patterns (e.g., consistent disagreement on chr7 could indicate
#'     pipeline-specific segmentation or normalization issues).
#' }
#'
#' \strong{Typical applications}
#' \itemize{
#'   \item Cross-caller benchmarking (e.g., FACETS vs CNVkit) using TP/FP/FN at annotation, CNV, and CCF levels.
#'   \item QC for over-segmentation: high \code{segments} with low \code{total_bp} suggests many focal differences.
#'   \item Investigating aneuploidy: per-chromosome fractions highlight broad events with disagreement.
#'   \item Clonality validation: use \code{ccf} TP rate under \code{wiggle} to judge consistency of subclonal calls.
#' }
#'
#' @examples
#' \dontrun{
#' # Compare two CNV profiles with a 0.1 CCF tolerance
#' cmp <- compareSegmentStatus(query = cnv_q, subject = cnv_s, wiggle = 0.1)
#'
#' # Overall agreement by status (TP/FP/FN)
#' cmp$stats
#'
#' # Genome-level summary by copy number and status
#' head(cmp$summary)
#'
#' # Per-chromosome agreement fractions
#' head(cmp$stats_by_chr)
#'
#' # Inspect detailed per-region comparisons
#' head(cmp$table)
#' }
#'
#' @seealso
#' \code{\link{compareSegments}} for base-range overlaps,
#' \code{\link{compareSegmentStatus}} for status-level comparisons,
#' and \code{\link{cnvMetrics}} for sample-wise CNV comparison.
#'
#' @export

compareSegmentStatus <- function(query, subject, wiggle = 0.1) {
  query.gr <- query@segments
  subject.gr <- subject@segments

  queryOv <- GenomicRanges::pintersect(IRanges::findOverlapPairs(query.gr, subject.gr))
  subjectOv <- GenomicRanges::pintersect(IRanges::findOverlapPairs(subject.gr, query.gr))

  keep <- intersect(which(subjectOv$hit == TRUE), which(queryOv$hit == TRUE))

  queryOv <- queryOv[keep]

  subjectOv <- subjectOv[keep]
  # all the ranges
  all <- length(subjectOv)
  all_bp <- sum(BiocGenerics::width(subjectOv))

  # find the common ranges (same cnv.flag)
  tt <- cbind(
    chrom = GenomeInfoDb::seqnames(queryOv) %>% as.character(),
    start = BiocGenerics::start(queryOv),
    end = BiocGenerics::end(queryOv),
    query = queryOv$cnv.annotation,
    subject = subjectOv$cnv.annotation,
    query_CNV = queryOv$total.allele.copies,
    subject_CNV = subjectOv$total.allele.copies,
    query_CCF = queryOv$cnv.cf,
    subject_CCF = subjectOv$cnv.cf
  ) %>% as.data.frame()
  tt <- tt %>%
    dplyr::mutate(across(c(start, end, query_CNV, subject_CNV, query_CCF, subject_CCF), as.numeric)) %>%
    dplyr::mutate(status = ifelse(query == subject, "TP", ifelse(query == "NEUT", "FP", "FN"))) %>%
    dplyr::mutate(cnv = ifelse(status == "TP", ifelse(query_CNV == subject_CNV, "TP", ifelse(query_CNV < subject_CNV, "FP", "FN")), ifelse(query_CNV == subject_CNV, "SAME", "DIFF"))) %>%
    dplyr::mutate(ccf = ifelse(status == "TP", ifelse(abs(query_CCF - subject_CCF) < wiggle, "TP", ifelse(query_CCF < subject_CCF, "FP", "FN")), ifelse(abs(query_CCF - subject_CCF) < wiggle, "SAME", "DIFF")))


  tt.summary <- tt %>%
    dplyr::filter(!chrom %in% c("chrX", "chrY")) %>%
    dplyr::group_by(query_CNV, subject_CNV, status) %>%
    dplyr::summarise(segments = dplyr::n(), total_bp = sum(end - start + 1), .groups = "keep") %>%
    dplyr::ungroup() %>%
    dplyr::mutate(genome_fraction = total_bp / sum(total_bp)) %>%
    dplyr::arrange(desc(total_bp))

  tt.stats <- tt.summary %>%
    dplyr::group_by(status) %>%
    dplyr::summarise(segments = sum(segments), total_bp = sum(total_bp)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(genome_fraction = total_bp / sum(total_bp))


  arms_lengths <- cnvstrmliner::chrom_arms[[query@genome]] %>% as.data.frame()
  chr_lengths <- arms_lengths %>%
    dplyr::group_by(seqnames) %>%
    dplyr::summarise(width = sum(width), .groups = c("keep")) %>%
    dplyr::ungroup()

  tt.summary_chr <- tt %>%
    dplyr::filter(!chrom %in% c("chrX", "chrY")) %>%
    dplyr::group_by(chrom, query_CNV, subject_CNV, status) %>%
    dplyr::summarise(segments = dplyr::n(), total_bp = sum(end - start + 1), .groups = "keep") %>%
    dplyr::ungroup() %>%
    dplyr::left_join(chr_lengths, by = c("chrom" = "seqnames")) %>%
    dplyr::mutate(chromosome_fraction = total_bp / width) %>%
    dplyr::arrange(chrom, desc(total_bp))

  tt.stats_chr <- tt.summary_chr %>%
    dplyr::group_by(chrom, status, width) %>%
    dplyr::summarise(segments = sum(segments), total_bp = sum(total_bp), .groups = "keep") %>%
    dplyr::ungroup() %>%
    dplyr::mutate(chromosome_fraction = total_bp / width) %>%
    dplyr::arrange(chrom, desc(total_bp)) %>%
    dplyr::select(-width)

  return(list(
    table = tt,
    summary = tt.summary,
    stats = tt.stats,
    summary_by_chr = tt.summary_chr,
    stats_by_chr = tt.stats_chr
  ))
}


#' Compare segment overlaps using per-state Jaccard similarity
#'
#' @description
#' Compares CNV segments between a query and a subject \code{cnvstrml} object and
#' annotates each segment with the best overlap match based on Jaccard similarity.
#' For each segment, all overlaps against the other object are evaluated, grouped by
#' overlapping \code{cnv.annotation} state, and the maximum Jaccard score across states
#' is retained.
#'
#' The Jaccard similarity for two segments \eqn{A} and \eqn{B} is:
#' \deqn{J(A, B) = \frac{|A \cap B|}{|A \cup B|}}
#'
#' where \eqn{|A \cap B|} is the intersecting width in base pairs and
#' \eqn{|A \cup B|} is the union width in base pairs.
#'
#' @param query A \code{cnvstrml} object to compare.
#' @param subject A \code{cnvstrml} object to compare against.
#'
#' @return A named list with two \code{GRanges} objects:
#' \describe{
#'   \item{query}{The query segments with two added metadata columns:
#'   \code{jaccard.index} (numeric best Jaccard score) and
#'   \code{jaccard.status} (state with best score, or \code{"NO_OVERLAP"}).}
#'   \item{subject}{The subject segments with the same added metadata columns,
#'   computed against query segments.}
#' }
#'
#' @export
compareSegmentOverlap <- function(query, subject) {
  if (!inherits(query, "cnvstrml")) {
    stop("query must be a cnvstrml object")
  }
  if (!inherits(subject, "cnvstrml")) {
    stop("subject must be a cnvstrml object")
  }

  if (!("cnv.annotation" %in% colnames(S4Vectors::mcols(query@segments)))) {
    stop("query@segments must include column 'cnv.annotation'")
  }
  if (!("cnv.annotation" %in% colnames(S4Vectors::mcols(subject@segments)))) {
    stop("subject@segments must include column 'cnv.annotation'")
  }

  .best_overlap_by_state <- function(source_gr, target_gr) {
    overlaps <- GenomicRanges::findOverlaps(source_gr, target_gr, ignore.strand = TRUE)

    best_index <- rep(NA_real_, length(source_gr))
    best_status <- rep("NO_OVERLAP", length(source_gr))

    if (length(overlaps) == 0) {
      return(list(index = best_index, status = best_status))
    }

    qh <- S4Vectors::queryHits(overlaps)
    sh <- S4Vectors::subjectHits(overlaps)

    intersection_gr <- GenomicRanges::pintersect(source_gr[qh], target_gr[sh], ignore.strand = TRUE)
    intersection_bp <- as.numeric(BiocGenerics::width(intersection_gr))
    union_bp <- as.numeric(BiocGenerics::width(source_gr[qh])) +
      as.numeric(BiocGenerics::width(target_gr[sh])) -
      intersection_bp

    overlap_tbl <- data.frame(
      source_idx = qh,
      state = as.character(S4Vectors::mcols(target_gr)$cnv.annotation[sh]),
      jaccard = intersection_bp / union_bp,
      stringsAsFactors = FALSE
    )

    best_by_state <- stats::aggregate(
      jaccard ~ source_idx + state,
      data = overlap_tbl,
      FUN = max
    )

    best_by_state <- best_by_state[order(best_by_state$source_idx, -best_by_state$jaccard, best_by_state$state), ]
    best_per_segment <- best_by_state[!duplicated(best_by_state$source_idx), ]

    best_index[best_per_segment$source_idx] <- best_per_segment$jaccard
    best_status[best_per_segment$source_idx] <- best_per_segment$state

    list(index = best_index, status = best_status)
  }

  query_gr <- query@segments
  subject_gr <- subject@segments

  query_best <- .best_overlap_by_state(query_gr, subject_gr)
  subject_best <- .best_overlap_by_state(subject_gr, query_gr)

  S4Vectors::mcols(query_gr)$jaccard.index <- query_best$index
  S4Vectors::mcols(query_gr)$jaccard.status <- query_best$status

  S4Vectors::mcols(subject_gr)$jaccard.index <- subject_best$index
  S4Vectors::mcols(subject_gr)$jaccard.status <- subject_best$status

  return(list(
    query = query_gr,
    subject = subject_gr
  ))
}
