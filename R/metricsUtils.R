# Shared internal helpers used by the metrics files. Not part of the public API.

# =============================================================================
# SEGMENTATION COMPARISON METRICS
#
# Every metric in this file answers the same shape of question: given two CNV
# segmentations of the same genome, how much do they agree? Some of the pairs
# are prediction-versus-truth (the simulation cohort), some are
# caller-versus-caller (real samples, where no truth exists). The file is
# organised by what a metric is *about* rather than by which of those two uses
# it serves, because most of them serve both.
#
# The file is written to be dropped verbatim into the R/ directory of the
# cnvstrmliner package. Consequences of that, all deliberate:
#
#   * Every external call is namespace qualified. Nothing here depends on
#     scripts/R/functions.R or scripts/R/alternative_metrics.R, and sourcing
#     this file alone into a fresh session must work.
#   * Functions are lowerCamelCase; arguments and locals are snake_case. The
#     package's own two metric functions (cnvMetrics, cnvEntropy) are
#     Prefix_snake, which is a minority form among its ~70 exports. That
#     inconsistency is pre-existing and is not resolved here.
#   * New per-segment mcols follow the package's dotted convention
#     (cnv.annotation, jaccard.index): bin.jaccard.*, state.*,
#     breakpoint.consistency.*.
#   * Roxygen is markdown = TRUE. Markdown is used for prose; Rd macros only
#     where markdown has no equivalent (\deqn, \describe, \dontrun).
#   * No dependency outside the package's existing Imports. In particular
#     infotheo is not used: variation of information is three lines on a
#     contingency table, and infotheo discretises, which is wrong here.
#
# Deliberate divergences from the project-local implementations these functions
# were ported from are documented per function in an "@section Differences from
# ..." block. There are three, and they all matter:
#
#   1. Auto-tiled bin grids depend only on (genome, bin_width, seqlevels), not
#      on the inputs. The reference computeARI/computeVI tile from the min/max
#      coordinate of the two inputs, which puts every comparison on a different
#      grid and makes the numbers non-comparable across callers.
#   2. Bin-set Jaccard assigns each bin to exactly one segment, so the union is
#      a true set union. The reference counts with countOverlaps, which counts
#      a boundary-straddling bin twice.
#   3. The Rand index is the full pair-agreement rate. The reference reports
#      only the same-cluster-in-both term.
#
# The names of the reference implementations are avoided throughout, because
# the analysis notebooks source() those into the global environment and a
# clashing name would let source() order decide silently which maths runs.
# =============================================================================


# ===== 0. INTERNAL HELPERS ===================================================
#
# Two families. The first normalises inputs: segmentations, bin grids and CNV
# state vocabularies all arrive in more than one shape and every exported
# function needs them in exactly one. The second does the arithmetic that the
# information-theoretic metrics share -- kappa, ARI and VI are each a handful
# of lines once a contingency table between two per-bin label vectors exists,
# and building that table once is the main reason this file is shorter than the
# sum of the implementations it replaces.
# =============================================================================


#' Coerce a segmentation argument to GRanges
#'
#' @description
#' Internal. Accepts a `cnvstrml` object or a bare `GRanges` and returns the
#' segments as `GRanges`, raising a contract error otherwise. Optionally checks
#' that named `mcols` columns are present.
#'
#' @param x A `cnvstrml` object or a `GRanges`.
#' @param what Character. Name of the argument, used in error messages.
#' @param require_cols Character vector of `mcols` column names that must exist.
#'
#' @return A `GRanges` object.
#'
#' @keywords internal
.asSegments <- function(x, what = "x", require_cols = character()) {
  if (inherits(x, "cnvstrml")) {
    x <- x@segments
  }
  if (!inherits(x, "GRanges")) {
    stop(sprintf("%s must be a cnvstrml object or a GRanges object", what))
  }
  if (length(require_cols) > 0) {
    .requireMcols(x, require_cols, what)
  }
  x
}


#' Check that a GRanges carries the required metadata columns
#'
#' @description
#' Internal. Mirrors the error wording used by
#' [cnvstrmliner::compareSegmentOverlap()] so that a missing column reports the
#' same way wherever it is detected.
#'
#' @param gr A `GRanges` object.
#' @param cols Character vector of required `mcols` column names.
#' @param what Character. Name of the argument, used in error messages.
#'
#' @return `TRUE`, invisibly. Raises an error if a column is absent.
#'
#' @keywords internal
.requireMcols <- function(gr, cols, what = "x") {
  present <- colnames(S4Vectors::mcols(gr))
  missing_cols <- setdiff(cols, present)
  if (length(missing_cols) > 0) {
    stop(sprintf(
      "%s must include column%s %s",
      what,
      if (length(missing_cols) > 1) "s" else "",
      paste0("'", missing_cols, "'", collapse = ", ")
    ))
  }
  invisible(TRUE)
}


#' Canonicalise CNV state labels
#'
#' @description
#' Internal. Applies the two state collapses the package recognises, in the
#' order [cnvstrmliner::cnvMetrics()] applies them. Factored out so that no two
#' metrics in this file can drift apart on what counts as the same state.
#'
#' @param x Character vector (or anything coercible) of state labels.
#' @param simplify Logical. Collapse `AMP-LOH` to `AMP`, `NEUT-LOH` to `NEUT`,
#'   and both `HETDEL` and `HOMDEL` to `DEL`.
#' @param infer_as_neutral Logical. Map `INFER` to `NEUT`.
#'
#' @return A character vector the same length as `x`. `NA` is preserved.
#'
#' @keywords internal
.canonicalStates <- function(x, simplify = FALSE, infer_as_neutral = TRUE) {
  x <- as.character(x)
  if (infer_as_neutral) {
    x[!is.na(x) & x == "INFER"] <- "NEUT"
  }
  if (simplify) {
    x[!is.na(x) & x == "AMP-LOH"] <- "AMP"
    x[!is.na(x) & x == "NEUT-LOH"] <- "NEUT"
    x[!is.na(x) & x %in% c("HETDEL", "HOMDEL")] <- "DEL"
  }
  x
}


#' Divide, returning NA rather than NaN on an empty denominator
#'
#' @description
#' Internal. A zero denominator in these metrics means "the question does not
#' apply to this input" -- a caller with no altered segments has no precision --
#' and that is `NA`, not zero. Zero would assert that it got everything wrong.
#'
#' @param num,den Numeric vectors.
#' @param zero Value to return where `den` is zero or not finite.
#'
#' @return A numeric vector.
#'
#' @keywords internal
.safeDiv <- function(num, den, zero = NA_real_) {
  out <- rep(zero, length(num))
  ok <- is.finite(den) & den != 0
  out[ok] <- num[ok] / den[ok]
  out
}


#' One all-NA row with a fixed column set
#'
#' @description
#' Internal. Used by the rollups and by [segmentationMetrics()] so that a
#' degenerate input still returns the same columns as a successful one, and
#' `dplyr::bind_rows()` over a sweep never produces a ragged frame.
#'
#' @param cols Character vector of column names.
#'
#' @return A one-row `data.frame` of `NA_real_`.
#'
#' @keywords internal
.emptyRow <- function(cols) {
  out <- as.data.frame(
    stats::setNames(rep(list(NA_real_), length(cols)), cols),
    stringsAsFactors = FALSE
  )
  out
}


#' Resolve a bin grid argument
#'
#' @description
#' Internal. Turns the several accepted spellings of a bin grid into one
#' `GRanges`: an existing `GRanges` passes through, a character path is read as
#' BED, and `NULL` auto-tiles the genome when the calling metric permits it.
#'
#' The auto-tiled grid is built from the genome's sequence lengths, restricted
#' to the sequence levels actually present in `segments_list`. It therefore
#' depends only on `(genome, bin_width, seqlevels)` and **not** on the segment
#' coordinates, so two callers scored against the same reference land on the
#' same grid and their metrics are comparable. Deriving the grid from the inputs
#' instead -- as the project-local `computeARI()` and `computeVI()` do -- breaks
#' exactly the comparison a sensitivity sweep exists to make.
#'
#' @param bins A `GRanges`, a character path to a BED file, or `NULL`.
#' @param segments_list Optional list of `GRanges` used to restrict the
#'   auto-tiled grid to the sequence levels in play.
#' @param bin_width Numeric width of the auto-tiled bins.
#' @param genome Character genome identifier, e.g. `"hg38"`. Required when
#'   auto-tiling.
#' @param what Character. Name of the argument, used in error messages.
#' @param allow_auto Logical. Whether `NULL` is acceptable. Metrics in which a
#'   bin is a *unit of evidence* set this `FALSE`; metrics in which a bin is
#'   merely a discretisation of the coordinate axis set it `TRUE`.
#'
#' @return A `GRanges` object.
#'
#' @keywords internal
.resolveBins <- function(bins,
                         segments_list = NULL,
                         bin_width = 1e6,
                         genome = NULL,
                         what = "bins",
                         allow_auto = TRUE) {
  if (inherits(bins, "GRanges")) {
    return(bins)
  }

  if (is.character(bins)) {
    if (length(bins) != 1L) {
      stop(sprintf("%s must be a single file path when given as a character", what))
    }
    if (!file.exists(bins)) {
      stop(sprintf("%s file does not exist: %s", what, bins))
    }
    if (!requireNamespace("rtracklayer", quietly = TRUE)) {
      stop(sprintf(
        "reading %s from a BED file requires the rtracklayer package", what
      ))
    }
    return(rtracklayer::import.bed(bins))
  }

  if (!is.null(bins)) {
    stop(sprintf("%s must be a GRanges object, a BED file path, or NULL", what))
  }

  if (!allow_auto) {
    stop(sprintf(
      paste(
        "%s is required for this metric and has no default. A bin here is a",
        "unit of evidence, not a discretisation of the coordinate axis: on",
        "whole-exome data a uniform base-pair tiling would score agreement",
        "mostly over unsequenced territory. Pass the capture target regions."
      ),
      what
    ))
  }

  if (is.null(genome) || !nzchar(genome)) {
    stop(sprintf(
      "cannot auto-tile %s: no genome is available. Pass bins explicitly, or a genome identifier.",
      what
    ))
  }

  seq_lengths <- .genomeSeqLengths(genome)
  keep <- names(seq_lengths)
  if (!is.null(segments_list)) {
    present <- unique(unlist(lapply(
      segments_list,
      function(g) as.character(GenomeInfoDb::seqlevelsInUse(g))
    )))
    present <- present[!is.na(present)]
    if (length(present) > 0) {
      keep <- intersect(keep, present)
    }
  }
  if (length(keep) == 0) {
    stop(sprintf("cannot auto-tile %s: no sequence level is shared with the genome", what))
  }

  tiles <- GenomicRanges::tileGenome(
    seq_lengths[keep],
    tilewidth = bin_width,
    cut.last.tile.in.chrom = TRUE
  )
  logger::log_debug(
    "auto-tiled {what}: {length(tiles)} bins of {bin_width} bp over {length(keep)} sequences ({genome})"
  )
  tiles
}


#' Sequence lengths for a genome from the package's bundled data
#'
#' @description
#' Internal. Prefers `cnvstrmliner::seqinfo`, which is a `Seqinfo` per genome,
#' and falls back to `cnvstrmliner::chrom_size`, a two-column data.frame per
#' genome.
#'
#' @param genome Character genome identifier, e.g. `"hg38"`.
#'
#' @return A named integer vector of sequence lengths.
#'
#' @keywords internal
.genomeSeqLengths <- function(genome) {
  si <- try(cnvstrmliner::seqinfo, silent = TRUE)
  if (!inherits(si, "try-error") && is.list(si) && genome %in% names(si)) {
    lens <- GenomeInfoDb::seqlengths(si[[genome]])
    lens <- lens[!is.na(lens)]
    if (length(lens) > 0) {
      return(lens)
    }
  }
  cs <- try(cnvstrmliner::chrom_size, silent = TRUE)
  if (!inherits(cs, "try-error") && is.list(cs) && genome %in% names(cs)) {
    tab <- cs[[genome]]
    lens <- stats::setNames(as.numeric(tab[[2]]), as.character(tab[[1]]))
    lens <- lens[!is.na(lens)]
    if (length(lens) > 0) {
      return(lens)
    }
  }
  stop(sprintf("no sequence lengths are available for genome '%s'", genome))
}


#' Genome identifier of a segmentation argument
#'
#' @description
#' Internal. Reads the `genome` slot of a `cnvstrml`, or the `genome` field of a
#' `GRanges` `seqinfo`, returning `NULL` when neither is populated.
#'
#' @param x A `cnvstrml` object or a `GRanges`.
#'
#' @return A single character genome identifier, or `NULL`.
#'
#' @keywords internal
.segmentsGenome <- function(x) {
  if (inherits(x, "cnvstrml")) {
    g <- x@genome
    if (length(g) == 1 && !is.na(g) && nzchar(g)) {
      return(g)
    }
    x <- x@segments
  }
  if (inherits(x, "GRanges")) {
    g <- unique(GenomeInfoDb::genome(x))
    g <- g[!is.na(g) & nzchar(g)]
    if (length(g) == 1) {
      return(g)
    }
  }
  NULL
}


#' Assign each bin to the single segment it overlaps most
#'
#' @description
#' Internal. Returns, for every bin, the index of the segment with which it
#' shares the most base pairs, breaking ties by the first segment, and `NA` for
#' bins overlapping no segment.
#'
#' The single-assignment guarantee is what makes bin *sets* well defined: a bin
#' straddling a segment boundary belongs to one segment, so a union of two bin
#' sets is a true set union rather than a double count.
#'
#' @param bins A `GRanges` of bins.
#' @param seg_gr A `GRanges` of segments.
#'
#' @return An integer vector of length `length(bins)`.
#'
#' @keywords internal
.bestBinSegment <- function(bins, seg_gr) {
  idx <- rep(NA_integer_, length(bins))
  if (length(bins) == 0 || length(seg_gr) == 0) {
    return(idx)
  }
  ov <- GenomicRanges::findOverlaps(bins, seg_gr, ignore.strand = TRUE)
  if (length(ov) == 0) {
    return(idx)
  }
  qh <- S4Vectors::queryHits(ov)
  sh <- S4Vectors::subjectHits(ov)
  shared_width <- as.numeric(BiocGenerics::width(
    GenomicRanges::pintersect(bins[qh], seg_gr[sh], ignore.strand = TRUE)
  ))
  ord <- order(qh, -shared_width, sh)
  keep <- !duplicated(qh[ord])
  idx[qh[ord][keep]] <- sh[ord][keep]
  idx
}


#' Per-bin CNV state vector
#'
#' @description
#' Internal. [.bestBinSegment()] followed by a state lookup and
#' [.canonicalStates()].
#'
#' @param bins A `GRanges` of bins.
#' @param seg_gr A `GRanges` of segments.
#' @param statecol Character. Name of the state column.
#' @param simplify,infer_as_neutral Passed to [.canonicalStates()].
#'
#' @return A character vector of length `length(bins)`, `NA` where unassigned.
#'
#' @keywords internal
.binStateVector <- function(bins, seg_gr, statecol = "cnv.annotation",
                            simplify = FALSE, infer_as_neutral = TRUE) {
  out <- rep(NA_character_, length(bins))
  idx <- .bestBinSegment(bins, seg_gr)
  assigned <- !is.na(idx)
  if (any(assigned)) {
    out[assigned] <- as.character(S4Vectors::mcols(seg_gr)[[statecol]][idx[assigned]])
  }
  .canonicalStates(out, simplify = simplify, infer_as_neutral = infer_as_neutral)
}


#' Per-bin segment identity vector
#'
#' @description
#' Internal. The segment *index* each bin belongs to, which is a different
#' partition of the bins from the one [.binStateVector()] produces: identity
#' distinguishes two adjacent segments that happen to share a state, state does
#' not.
#'
#' @param bins A `GRanges` of bins.
#' @param seg_gr A `GRanges` of segments.
#'
#' @return An integer vector of length `length(bins)`, `NA` where unassigned.
#'
#' @keywords internal
.binSegmentVector <- function(bins, seg_gr) {
  .bestBinSegment(bins, seg_gr)
}


#' Contingency table between two aligned per-bin label vectors
#'
#' @description
#' Internal. The shared substrate of Cohen's kappa, the adjusted Rand index and
#' variation of information. Bins that are `NA` on either side are dropped and
#' counted, never imputed.
#'
#' @param a,b Label vectors of equal length.
#' @param w Optional numeric weights, one per bin. `NULL` weights every bin
#'   equally.
#' @param common_levels Logical. Force both margins onto the union of the two
#'   label sets, giving a square table. Kappa needs this because it compares
#'   like with like; ARI and VI do not, because they compare partitions whose
#'   labels are arbitrary.
#'
#' @return A named list with `tab` (the contingency matrix, weighted if `w` was
#'   supplied), `n` (total weight retained), `n_scored` (bins retained) and
#'   `n_dropped` (bins discarded for an `NA` on either side).
#'
#' @keywords internal
.contingency <- function(a, b, w = NULL, common_levels = FALSE) {
  stopifnot(length(a) == length(b))
  keep <- !is.na(a) & !is.na(b)
  n_dropped <- sum(!keep)

  if (!any(keep)) {
    return(list(
      tab = matrix(numeric(0), nrow = 0, ncol = 0),
      n = 0, n_scored = 0L, n_dropped = as.integer(n_dropped)
    ))
  }

  a_keep <- a[keep]
  b_keep <- b[keep]

  if (common_levels) {
    lev <- sort(unique(c(as.character(a_keep), as.character(b_keep))))
    a_keep <- factor(as.character(a_keep), levels = lev)
    b_keep <- factor(as.character(b_keep), levels = lev)
  } else {
    a_keep <- factor(a_keep)
    b_keep <- factor(b_keep)
  }

  if (is.null(w)) {
    tab <- table(a_keep, b_keep)
    tab <- matrix(as.numeric(tab), nrow = nrow(tab), ncol = ncol(tab),
                  dimnames = dimnames(tab))
  } else {
    w_keep <- as.numeric(w)[keep]
    tab <- tapply(w_keep, list(a_keep, b_keep), sum)
    tab[is.na(tab)] <- 0
    tab <- matrix(as.numeric(tab), nrow = nrow(tab), ncol = ncol(tab),
                  dimnames = list(levels(a_keep), levels(b_keep)))
  }

  list(
    tab = tab,
    n = sum(tab),
    n_scored = as.integer(sum(keep)),
    n_dropped = as.integer(n_dropped)
  )
}


#' Cohen's kappa from a square contingency table
#'
#' @description
#' Internal. Expects a table built with `common_levels = TRUE`, so that the
#' diagonal is the agreement.
#'
#' @param tab A square numeric contingency matrix.
#'
#' @return A named list with `p_o`, `p_e`, `kappa` and `marginals`.
#'
#' @keywords internal
.kappaFromTable <- function(tab) {
  n <- sum(tab)
  if (n == 0 || nrow(tab) == 0 || nrow(tab) != ncol(tab)) {
    return(list(p_o = NA_real_, p_e = NA_real_, kappa = NA_real_,
                marginals = data.frame(
                  state = character(), p_a = numeric(), p_b = numeric(),
                  stringsAsFactors = FALSE
                )))
  }

  p_o <- sum(diag(tab)) / n
  p_a <- rowSums(tab) / n
  p_b <- colSums(tab) / n
  p_e <- sum(p_a * p_b)

  # p_e == 1 means both sides are constant and equal: every bin the same single
  # state on both. Chance already explains all of the agreement, so there is no
  # skill to credit. Returning (p_o - p_e) / 0 here would give NaN on exactly
  # the all-neutral-versus-all-neutral case that scoreNullCaller() exists to
  # check, so the conventional kappa = 0 is used instead.
  kappa_value <- if (isTRUE(all.equal(p_e, 1))) {
    logger::log_warn("degenerate kappa: expected agreement is 1, returning kappa = 0")
    0
  } else {
    (p_o - p_e) / (1 - p_e)
  }

  list(
    p_o = p_o,
    p_e = p_e,
    kappa = kappa_value,
    marginals = data.frame(
      state = rownames(tab),
      p_a = as.numeric(p_a),
      p_b = as.numeric(p_b),
      stringsAsFactors = FALSE
    )
  )
}


#' Rand and adjusted Rand index from a contingency table
#'
#' @description
#' Internal. Pair-counting agreement between two partitions of the same bins.
#'
#' @param tab A numeric contingency matrix of bin counts.
#'
#' @return A named list with `ari` and `ri`.
#'
#' @keywords internal
.ariFromTable <- function(tab) {
  n <- sum(tab)
  if (n < 2) {
    return(list(ari = NA_real_, ri = NA_real_))
  }

  comb2 <- function(x) x * (x - 1) / 2

  sum_ij <- sum(comb2(as.numeric(tab)))
  sum_i <- sum(comb2(rowSums(tab)))
  sum_j <- sum(comb2(colSums(tab)))
  total_pairs <- comb2(n)

  # The full Rand index: pairs the two partitions agree about, whether they are
  # together in both or apart in both. The project-local computeARI() reports
  # only the together-in-both term, which understates agreement badly when most
  # pairs are apart -- as they always are with many segments.
  ri <- (total_pairs + 2 * sum_ij - sum_i - sum_j) / total_pairs

  expected_index <- sum_i * sum_j / total_pairs
  max_index <- (sum_i + sum_j) / 2
  ari <- if (isTRUE(all.equal(max_index, expected_index))) {
    0
  } else {
    (sum_ij - expected_index) / (max_index - expected_index)
  }

  list(ari = ari, ri = ri)
}


#' Entropies, mutual information and variation of information from a table
#'
#' @description
#' Internal. All quantities are in nats.
#'
#' @param tab A numeric contingency matrix of bin counts.
#'
#' @return A named list with `vi`, `nvi`, `h_a`, `h_b` and `mi`.
#'
#' @keywords internal
.viFromTable <- function(tab) {
  n <- sum(tab)
  if (n == 0) {
    return(list(vi = NA_real_, nvi = NA_real_, h_a = NA_real_,
                h_b = NA_real_, mi = NA_real_))
  }

  p_ij <- tab / n
  p_a <- rowSums(p_ij)
  p_b <- colSums(p_ij)

  entropy <- function(p) {
    p <- p[p > 0]
    if (length(p) == 0) return(0)
    -sum(p * log(p))
  }

  h_a <- entropy(p_a)
  h_b <- entropy(p_b)

  p_outer <- outer(p_a, p_b, "*")
  nz <- p_ij > 0 & p_outer > 0
  mi <- if (any(nz)) sum(p_ij[nz] * log(p_ij[nz] / p_outer[nz])) else 0

  vi <- h_a + h_b - 2 * mi
  # Floating point can put VI a hair below zero when the two partitions agree
  # exactly; it is a distance and must not be reported negative.
  vi <- max(vi, 0)

  # Raw VI is unbounded and grows with the number of bins, so it cannot be
  # compared across samples with different grids. log(n) is the maximum possible
  # value, giving a [0, 1] companion.
  nvi <- if (n > 1) vi / log(n) else NA_real_

  list(vi = vi, nvi = nvi, h_a = h_a, h_b = h_b, mi = mi)
}


#' De-duplicated segment boundaries as a data frame
#'
#' @description
#' Internal. Both ends of every segment, de-duplicated, optionally dropping the
#' first and last position on each sequence -- those are artefacts of where the
#' segmentation happens to begin and end rather than detected events.
#'
#' @param seg_gr A `GRanges` of segments.
#' @param exclude_chr_ends Logical.
#'
#' @return A `data.frame` with columns `seqnames` and `pos`.
#'
#' @keywords internal
.breakpointTable <- function(seg_gr, exclude_chr_ends = TRUE) {
  empty <- data.frame(seqnames = character(), pos = numeric(),
                      stringsAsFactors = FALSE)
  if (length(seg_gr) == 0) {
    return(empty)
  }

  chrs <- as.character(GenomicRanges::seqnames(seg_gr))
  starts <- GenomicRanges::start(seg_gr)
  ends <- GenomicRanges::end(seg_gr)

  bp <- data.frame(
    seqnames = c(chrs, chrs),
    pos = c(starts, ends),
    stringsAsFactors = FALSE
  )
  bp <- unique(bp)

  if (exclude_chr_ends && nrow(bp) > 0) {
    first_pos <- stats::aggregate(list(pos = starts), by = list(seqnames = chrs), FUN = min)
    last_pos <- stats::aggregate(list(pos = ends), by = list(seqnames = chrs), FUN = max)
    edge_keys <- c(
      paste(first_pos$seqnames, first_pos$pos, sep = ":"),
      paste(last_pos$seqnames, last_pos$pos, sep = ":")
    )
    bp <- bp[!paste(bp$seqnames, bp$pos, sep = ":") %in% edge_keys, , drop = FALSE]
  }

  bp <- bp[order(bp$seqnames, bp$pos), , drop = FALSE]
  rownames(bp) <- NULL
  bp
}


#' Tolerance matching between two breakpoint sets
#'
#' @description
#' Internal. Counts how many breakpoints on each side have a counterpart on the
#' other within `tolerance_bp`. Matching is done in both directions rather than
#' as a one-to-one assignment, so the two counts can differ: three query
#' breakpoints clustered inside one tolerance window all match, but they match
#' the same single truth breakpoint.
#'
#' @param query_bp,subject_bp `GRanges` of width-1 breakpoints.
#' @param tolerance_bp Numeric half-width of the matching window.
#'
#' @return A named list with `n_matched_query` and `n_matched_subject`.
#'
#' @keywords internal
.matchBreakpoints <- function(query_bp, subject_bp, tolerance_bp) {
  if (length(query_bp) == 0 || length(subject_bp) == 0) {
    return(list(n_matched_query = 0L, n_matched_subject = 0L))
  }

  common_chrs <- intersect(
    GenomeInfoDb::seqlevels(query_bp),
    GenomeInfoDb::seqlevels(subject_bp)
  )
  if (length(common_chrs) > 0) {
    GenomeInfoDb::seqlevels(query_bp, pruning.mode = "coarse") <- common_chrs
    GenomeInfoDb::seqlevels(subject_bp, pruning.mode = "coarse") <- common_chrs
  }
  if (length(query_bp) == 0 || length(subject_bp) == 0) {
    return(list(n_matched_query = 0L, n_matched_subject = 0L))
  }

  query_window <- suppressWarnings(query_bp + tolerance_bp)
  subject_window <- suppressWarnings(subject_bp + tolerance_bp)

  query_hits <- GenomicRanges::findOverlaps(query_window, subject_bp, ignore.strand = TRUE)
  subject_hits <- GenomicRanges::findOverlaps(subject_window, query_bp, ignore.strand = TRUE)

  list(
    n_matched_query = length(unique(S4Vectors::queryHits(query_hits))),
    n_matched_subject = length(unique(S4Vectors::queryHits(subject_hits)))
  )
}


#' Cluster sorted positions with a running-median tolerance
#'
#' @description
#' Internal. Walks a sorted vector of positions and opens a new cluster whenever
#' the next position is further than `tolerance_bp` from the running median of
#' the open one. Vectorised port of the loop inside the project-local
#' `computeBreakpointStability()`; the running median, rather than the last
#' position, is what stops a chain of closely spaced breakpoints from merging
#' into one arbitrarily wide cluster.
#'
#' @param pos Numeric vector of positions, assumed sorted ascending.
#' @param tolerance_bp Numeric.
#'
#' @return An integer vector of cluster identifiers, one per position.
#'
#' @keywords internal
.clusterPositions <- function(pos, tolerance_bp) {
  n <- length(pos)
  if (n == 0) {
    return(integer(0))
  }
  if (n == 1) {
    return(1L)
  }

  cluster_id <- integer(n)
  cluster_id[1] <- 1L
  current <- 1L
  members <- pos[1]

  for (i in 2:n) {
    if (pos[i] - stats::median(members) <= tolerance_bp) {
      members <- c(members, pos[i])
    } else {
      current <- current + 1L
      members <- pos[i]
    }
    cluster_id[i] <- current
  }

  cluster_id
}


#' Bins-by-methods logical matrix of altered calls
#'
#' @description
#' Internal. The substrate of the truth-free agreement metrics: one column per
#' segmentation, `TRUE` where that segmentation calls the bin altered. Bins the
#' segmentation does not cover are `FALSE` -- for an altered/not question, "no
#' call here" and "not altered here" are the same answer from the point of view
#' of a peer looking for support.
#'
#' @param cnv_list Named list of `cnvstrml` objects or `GRanges`.
#' @param bins A `GRanges` of bins.
#' @param statecol Character. Name of the state column.
#' @param neutral_states Character vector of states counted as not altered.
#'
#' @return A logical matrix with `length(bins)` rows and `length(cnv_list)`
#'   columns, column names taken from `cnv_list`.
#'
#' @keywords internal
.alteredBinMatrix <- function(cnv_list, bins, statecol = "cnv.annotation",
                              neutral_states = c("NEUT", "INFER")) {
  method_names <- names(cnv_list)
  if (is.null(method_names)) {
    method_names <- paste0("method_", seq_along(cnv_list))
  }

  cols <- lapply(seq_along(cnv_list), function(i) {
    seg_gr <- .asSegments(cnv_list[[i]], what = method_names[i], require_cols = statecol)
    state <- .binStateVector(bins, seg_gr, statecol = statecol,
                             simplify = FALSE, infer_as_neutral = FALSE)
    !is.na(state) & !(state %in% neutral_states)
  })

  out <- matrix(unlist(cols), nrow = length(bins), ncol = length(cols))
  colnames(out) <- method_names
  out
}


#' Jaccard between two logical vectors
#'
#' @description
#' Internal. Set overlap of the `TRUE` positions. Two empty sets are perfectly
#' concordant about the absence of events, so the empty/empty case is 1 rather
#' than `NA`.
#'
#' @param a,b Logical vectors of equal length.
#'
#' @return A single numeric in `[0, 1]`.
#'
#' @keywords internal
.jaccardLogical <- function(a, b) {
  n_shared <- sum(a & b)
  n_union <- sum(a | b)
  if (n_union == 0) {
    return(1)
  }
  n_shared / n_union
}


#' Bin partition contingency table shared by ARI, VI and fragmentation
#'
#' @description
#' Internal. Resolves the grid, labels every bin on both sides according to
#' `partition_by`, and cross-tabulates.
#'
#' @param query,subject Segmentation arguments.
#' @param bins,bin_width,genome Grid arguments; see [.resolveBins()].
#' @param partition_by `"segment"` or `"state"`.
#' @param statecol,simplify,infer_as_neutral State arguments, used only in
#'   `"state"` mode.
#'
#' @return The list returned by [.contingency()].
#'
#' @keywords internal
.partitionBins <- function(query, subject, bins, bin_width, partition_by,
                           statecol, simplify, infer_as_neutral, genome) {
  require_cols <- if (partition_by == "state") statecol else character()
  query_gr <- .asSegments(query, "query", require_cols)
  subject_gr <- .asSegments(subject, "subject", require_cols)

  if (is.null(genome)) {
    genome <- .segmentsGenome(query)
  }
  if (is.null(genome)) {
    genome <- .segmentsGenome(subject)
  }

  bins <- .resolveBins(
    bins,
    segments_list = list(query_gr, subject_gr),
    bin_width = bin_width, genome = genome, allow_auto = TRUE
  )

  if (length(bins) == 0 || length(query_gr) == 0 || length(subject_gr) == 0) {
    return(list(tab = matrix(numeric(0), nrow = 0, ncol = 0),
                n = 0, n_scored = 0L, n_dropped = length(bins)))
  }

  if (partition_by == "state") {
    query_label <- .binStateVector(bins, query_gr, statecol, simplify, infer_as_neutral)
    subject_label <- .binStateVector(bins, subject_gr, statecol, simplify, infer_as_neutral)
  } else {
    query_label <- .binSegmentVector(bins, query_gr)
    subject_label <- .binSegmentVector(bins, subject_gr)
  }

  .contingency(query_label, subject_label, common_levels = FALSE)
}


# ===== 3. PER-SEGMENT ANNOTATION =============================================
#
# Sections 1 and 2 collapse a whole segmentation to a handful of numbers. That
# is what a sweep needs to rank settings, and it is useless for the question
# that follows: given the setting chosen, which of its segments can be
# believed? These functions answer that by attaching the comparison to each
# segment, so a segmentation can be filtered rather than accepted or rejected
# whole.
#
# They all return their input as GRanges with new mcols, following
# cnvstrmliner::compareSegmentOverlap, and they all use single bin assignment,
# so a bin belongs to exactly one segment on each side.
# =============================================================================


#' Copy-number state of each bin under a segmentation
#'
#' @description
#' Labels every bin with the state of the segment it overlaps most. Bins
#' overlapping no segment are `NA`.
#'
#' This is the projection that turns a segmentation into a per-bin label vector,
#' and it is the building block of every bin-level metric in this file. It is
#' exported because it is also the hoisting hook: when one reference is scored
#' against many segmentations -- forty-two settings of a parameter sweep, say --
#' the reference projection is computed once here and passed to
#' [annotateBinStateConcordance()] through its `reference_state` argument.
#'
#' @param bins A `GRanges` of bin-level regions.
#' @param segments A `cnvstrml` object or a `GRanges` of segments.
#' @param statecol Character. Name of the `mcols` column holding the state.
#' @param simplify,infer_as_neutral Logical. State canonicalisation; see
#'   [computeStateKappa()].
#'
#' @return A character vector of length `length(bins)`.
#'
#' @seealso [annotateBinStateConcordance()], [computeStateKappa()]
#'
#' @examples
#' \dontrun{
#' truth_state <- assignBinState(targets, truth)
#' table(truth_state, useNA = "ifany")
#' }
#'
#' @export
assignBinState <- function(bins,
                           segments,
                           statecol = "cnv.annotation",
                           simplify = FALSE,
                           infer_as_neutral = TRUE) {
  seg_gr <- .asSegments(segments, "segments", statecol)
  if (!inherits(bins, "GRanges")) {
    stop("bins must be a GRanges object")
  }
  .binStateVector(bins, seg_gr, statecol, simplify, infer_as_neutral)
}


#' Best bin-set Jaccard partner for each segment on one side
#'
#' @description
#' Internal. Given the shared-bin counts of every overlapping segment pair,
#' picks the partner maximising the Jaccard for each source segment. Ties break
#' on the partner's state label and then its index, so the result does not
#' depend on input ordering.
#'
#' @param source_idx,target_idx Integer vectors of paired segment indices.
#' @param n_shared Numeric vector of shared bin counts, same length.
#' @param source_bins,target_bins Integer vectors of bins per segment on each
#'   side.
#' @param target_state Character vector of states, one per target segment.
#' @param n_source Integer. Number of segments on the source side.
#'
#' @return A named list of five vectors, each of length `n_source`.
#'
#' @keywords internal
.bestBinJaccard <- function(source_idx, target_idx, n_shared,
                            source_bins, target_bins, target_state, n_source) {
  out <- list(
    index = rep(NA_real_, n_source),
    status = rep("NO_OVERLAP", n_source),
    n_self = rep(NA_integer_, n_source),
    n_match = rep(NA_integer_, n_source),
    n_shared = rep(NA_integer_, n_source)
  )

  if (length(source_idx) == 0) {
    return(out)
  }

  n_self <- source_bins[source_idx]
  n_match <- target_bins[target_idx]
  n_union <- n_self + n_match - n_shared
  jaccard <- ifelse(n_union > 0, n_shared / n_union, 0)
  state <- target_state[target_idx]

  ord <- order(source_idx, -jaccard, state, target_idx)
  keep <- !duplicated(source_idx[ord])
  win <- ord[keep]

  out$index[source_idx[win]] <- jaccard[win]
  out$status[source_idx[win]] <- state[win]
  out$n_self[source_idx[win]] <- as.integer(n_self[win])
  out$n_match[source_idx[win]] <- as.integer(n_match[win])
  out$n_shared[source_idx[win]] <- as.integer(n_shared[win])
  out
}


#' Resolve the grid and build the altered-bin matrix for a caller family
#'
#' @description
#' Internal. Shared entry point of [computeCallerAgreement()] and
#' [computeConsensusAgreement()], so the two cannot disagree about what a grid
#' or an altered call is.
#'
#' @param cnv_list Named list of segmentations.
#' @param bins,bin_width,genome Grid arguments.
#' @param statecol,neutral_states State arguments.
#'
#' @return A named list with `altered` (the logical matrix) and `bins`.
#'
#' @keywords internal
.prepareAlteredMatrix <- function(cnv_list, bins, bin_width, statecol,
                                  neutral_states, genome) {
  if (!is.list(cnv_list)) {
    stop("cnv_list must be a list of cnvstrml or GRanges objects")
  }

  method_names <- names(cnv_list)
  if (is.null(method_names)) {
    method_names <- paste0("method_", seq_along(cnv_list))
    names(cnv_list) <- method_names
  }

  if (length(cnv_list) == 0) {
    return(list(
      altered = matrix(logical(0), nrow = 0, ncol = 0,
                       dimnames = list(NULL, character())),
      bins = GenomicRanges::GRanges()
    ))
  }

  seg_list <- lapply(seq_along(cnv_list), function(i) {
    .asSegments(cnv_list[[i]], method_names[i], statecol)
  })

  if (is.null(genome)) {
    for (x in cnv_list) {
      genome <- .segmentsGenome(x)
      if (!is.null(genome)) break
    }
  }

  bins <- .resolveBins(bins, segments_list = seg_list, bin_width = bin_width,
                       genome = genome, allow_auto = TRUE)

  altered <- .alteredBinMatrix(
    stats::setNames(seg_list, method_names), bins,
    statecol = statecol, neutral_states = neutral_states
  )

  list(altered = altered, bins = bins)
}
