#' Filter Low-Quality Segments
#'
#' @description
#' Removes segments with insufficient evidence (low marker count or small size)
#' to prevent them from inflating the complexity penalty in BIC calculations.
#' Low-quality segments increase k without adding real information, unfairly
#' penalizing otherwise good models.
#'
#' @param segments A \code{GRanges} object containing CNV segments with metadata.
#' @param min_markers Numeric. Minimum number of markers (cnv.markers) required.
#'   Default is 5. Use NULL to disable this filter.
#' @param min_size Numeric. Minimum segment size in base pairs. Default is 50000 (50kb).
#'   Use NULL to disable this filter.
#' @param marker_col Character. Name of the metadata column containing marker counts.
#'   Default is "cnv.markers".
#'
#' @return A filtered \code{GRanges} object with low-quality segments removed.
#'   Returns all segments if no quality columns are found.
#'
#' @details
#' This function addresses a critical issue in BIC-based model selection:
#' spurious segments inflate the complexity penalty (k) without contributing
#' real signal, making good models appear worse than they are.
#'
#' Default thresholds are conservative:
#' \itemize{
#'   \item min_markers = 5: Segments with < 5 supporting markers are likely noise
#'   \item min_size = 50kb: Very small segments are often artifacts or low-confidence calls
#' }
#'
#' The filtering is applied with AND logic: segments must satisfy BOTH criteria
#' (or whichever filters are enabled).
#'
#' For methods without marker information, the function gracefully handles:
#' \itemize{
#'   \item Missing marker column: skips marker filtering
#'   \item All NA values: skips marker filtering (assumes all segments pass)
#'   \item All 0 values: skips marker filtering (assumes all segments pass)
#' }
#' This ensures methods that don't report marker counts aren't unfairly penalized.
#'
#' @examples
#' \dontrun{
#' # Default filtering (5 markers, 50kb)
#' filtered_segs <- filterLowQualitySegments(facets.cnv@segments)
#'
#' # More stringent filtering
#' filtered_segs <- filterLowQualitySegments(
#'   segments = cfreec.cnv@segments,
#'   min_markers = 10,
#'   min_size = 100000
#' )
#'
#' # Only size filtering (for methods without marker counts)
#' filtered_segs <- filterLowQualitySegments(
#'   segments = truth.cnv@segments,
#'   min_markers = NULL,
#'   min_size = 50000
#' )
#' }
#'
#' @seealso
#' \code{\link{computeBIC}}
#'
#' @export
filterLowQualitySegments <- function(segments,
                                     min_markers = 5,
                                     min_size = 50000,
                                     marker_col = "cnv.markers") {
  if (length(segments) == 0) {
    logger::log_info("No segments to filter")
    return(segments)
  }

  original_n <- length(segments)
  keep_mask <- rep(TRUE, original_n)

  # Apply marker filter if requested and column exists
  if (!is.null(min_markers) && marker_col %in% names(S4Vectors::mcols(segments))) {
    marker_counts <- S4Vectors::mcols(segments)[[marker_col]]

    # Check if marker data is actually available (not all NA or all 0)
    all_na <- all(is.na(marker_counts))
    all_zero <- !all_na && all(marker_counts == 0, na.rm = TRUE)

    if (all_na || all_zero) {
      # Method doesn't provide marker counts - skip marker filter
      logger::log_info(sprintf(
        "Marker filter: all values are %s - assuming all segments pass marker filter",
        ifelse(all_na, "NA", "0")
      ))
    } else {
      # Apply marker filter
      marker_filter <- !is.na(marker_counts) & marker_counts >= min_markers
      keep_mask <- keep_mask & marker_filter

      n_removed_markers <- sum(!marker_filter)
      logger::log_info(sprintf(
        "Marker filter (>= %d): removed %d segments (%.1f%%)",
        min_markers,
        n_removed_markers,
        100 * n_removed_markers / original_n
      ))
    }
  } else if (!is.null(min_markers)) {
    logger::log_info(sprintf("Column '%s' not found - skipping marker filter", marker_col))
  }

  # Apply size filter if requested
  if (!is.null(min_size)) {
    size_filter <- GenomicRanges::width(segments) >= min_size
    keep_mask <- keep_mask & size_filter

    n_removed_size <- sum(!size_filter & rep(TRUE, original_n))
    logger::log_info(sprintf(
      "Size filter (>= %.0f kb): removed %d segments (%.1f%%)",
      min_size / 1000,
      n_removed_size,
      100 * n_removed_size / original_n
    ))
  }

  # Apply combined filter
  filtered_segments <- segments[keep_mask]
  n_removed_total <- original_n - length(filtered_segments)

  logger::log_info(sprintf(
    "Total: kept %d / %d segments (removed %d, %.1f%%)",
    length(filtered_segments),
    original_n,
    n_removed_total,
    100 * n_removed_total / original_n
  ))

  return(filtered_segments)
}
