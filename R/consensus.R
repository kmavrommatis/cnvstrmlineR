#' Find Consensus CNV Segments Across Multiple Samples
#'
#' @description
#' This function identifies consensus copy number variant (CNV) segments across multiple
#' samples by finding common genomic regions and aggregating copy number information.
#' The algorithm splits and aligns segments from different CNV calls to identify regions
#' that are consistently detected across samples, enabling robust multi-sample CNV analysis.
#'
#'
#' Each consensus segment contains the copy number values from all input samples,
#' allowing for statistical analysis and identification of recurrent alterations.
#'
#' @param cnvList A named list of \code{cnvstrml} objects representing CNV segments
#'   from different samples or callers. If named, these names will be used as column
#'   headers in the output; otherwise, sample names from the objects will be used.
#'   All objects must have the same genome build.
#' @param ignore.samples Logical. If \code{FALSE} (default), all \code{cnvstrml}
#'   objects in the list must have the same sample name, which is useful when
#'   comparing multiple CNV calling methods on the same sample. If \code{TRUE},
#'   different sample names are allowed, enabling multi-sample consensus analysis.
#' @param method Character. Method for aggregating consensus values. Options:
#'   \itemize{
#'     \item \code{"quantitative"} - For numeric columns (copy numbers, log ratios).
#'       Automatically selected when \code{qual.col} is numeric.
#'     \item \code{"qualitative"} - For categorical columns (e.g., CNV annotations
#'       like "AMP", "DEL", "NEUT"). Automatically selected when \code{qual.col}
#'       is "cnv.annotation".
#'   }
#'   The method is automatically inferred based on \code{qual.col} if not specified.
#' @param qual.col Character. The column name from the CNV segments to use for
#'   consensus calculation. Options include:
#'   \itemize{
#'     \item \code{"total.allele.copies"} - Total copy number (default)
#'     \item \code{"minor.allele.copies"} - Minor allele copy number
#'     \item \code{"cnv.lrr"} - Log R ratio values
#'     \item \code{"cnv.cf"} - Cellular fraction
#'     \item \code{"cnv.ccf"} - Cancer cell fraction
#'     \item \code{"cnv.markers"} - Number of informative markers
#'     \item \code{"cnv.totalmarkers"} - Total number of markers
#'     \item \code{"cnv.totalmarkers_density"} - Marker density (total)
#'     \item \code{"cnv.markers_density"} - Marker density (informative)
#'     \item \code{"cnv.annotation"} - Categorical CNV type (AMP/DEL/NEUT/LOH)
#'     \item \code{"cnv.cnaqc.pass"} - CNAqc quality control pass/fail status
#'   }
#'   If \code{NULL}, all available columns are processed.
#' @param fill.gaps Logical. Currently not implemented. Reserved for future
#'   functionality to fill gaps between segments with neutral regions.
#' @param segment.annotation Character. Currently not implemented. Reserved for
#'   annotation of filled gap segments.
#' @param fn Character. Function for computing consensus summary statistics.
#'   Options:
#'   \itemize{
#'     \item \code{"majority"} - Reports the most common value across samples
#'       (median for quantitative, mode for qualitative)
#'     \item \code{"stats"} - Reports mean, median, and standard deviation
#'       (quantitative method only)
#'   }
#'   Currently not fully implemented.
#'
#' @return A \code{cnvstrml} object containing consensus segments with the following
#'   characteristics:
#'   \describe{
#'     \item{\code{@segments}}{GRanges object where each segment represents a
#'       consensus region. Metadata columns include:}
#'     \item{Sample columns}{One column per input sample containing the value of
#'       \code{qual.col} for that sample. Values are \code{NA} where the sample
#'       has no segment covering that region.}
#'     \item{\code{stats.mean}}{Mean value across samples (quantitative method only)}
#'     \item{\code{stats.median}}{Median value across samples (quantitative method only)}
#'     \item{\code{stats.sd}}{Standard deviation across samples (quantitative method only)}
#'     \item{\code{stats.majority}}{Most common value (median for quantitative,
#'       mode for qualitative)}
#'     \item{\code{stats.majority.count}}{Number of samples supporting the majority value}
#'     \item{\code{stats.count}}{Number of samples with non-NA values for each segment}
#'     \item{\code{@purity}}{Median purity across all input samples}
#'     \item{\code{@ploidy}}{Median ploidy across all input samples}
#'     \item{\code{@method}}{Comma-separated list of CNV calling methods used}
#'     \item{\code{@sample}}{Sample name (first sample if ignore.samples=FALSE)}
#'     \item{\code{@sex}}{Sex (prioritizes "F" if any female samples present)}
#'     \item{\code{@quality}}{List containing consensus metadata including
#'       \code{qual.col} and \code{method} used}
#'   }
#'
#' @details
#' \strong{Algorithm Overview:}
#'
#' The function performs the following steps:
#' \enumerate{
#'   \item Validates that all input objects have compatible genomes and sample names
#'   \item Identifies all unique breakpoints across samples using \code{\link{findCommonGR}}
#'   \item Creates consensus segments at the intersection of all breakpoints
#'   \item Annotates each consensus segment with values from each input sample
#'   \item Computes summary statistics (majority, mean, median, etc.)
#' }
#'
#' The consensus segments represent the intersection of all input segments, subdivided
#' at breakpoint boundaries to create a unified genomic coordinate system:
#'
#' \preformatted{
#' sample 1 : < --- seg1 ---- >< -------- seg 2 ------- >
#' sample 2 : < -- seg1 -- >< --- seg 2 --->< - seg 3 - >
#' consensus: < ---------- ><-><------------><---------- >
#' }
#' \strong{Coverage Considerations:}
#'
#' \emph{Important:} If not all samples cover the complete genome, some genomic
#' regions may be excluded from the consensus. For example, if gene X is covered
#' by segments in n-1 samples but completely missing in one sample, gene X will
#' not appear in the consensus output.
#'
#' To ensure complete genome coverage, consider using \code{\link{fillGaps}} on
#' individual samples before calling \code{findConsensus}.
#'
#' \strong{Qualitative vs Quantitative Analysis:}
#'
#' \itemize{
#'   \item \strong{Quantitative:} Used for numeric copy number data. The "majority"
#'     is defined as the median value. Statistical summaries (mean, median, SD)
#'     are computed across samples.
#'   \item \strong{Qualitative:} Used for categorical annotations (AMP, DEL, NEUT, LOH).
#'     The "majority" is the most frequent annotation. Segments annotated as "INFER"
#'     (from \code{fillGaps}) are excluded from majority calculations.
#' }
#'
#' \strong{Purity and Ploidy:}
#'
#' The consensus object's purity and ploidy are set to the median values across
#' all input samples, providing a representative estimate for the cohort.
#'
#' @note
#' \itemize{
#'   \item All elements in \code{cnvList} must share the same genome build (e.g., hg38)
#'   \item When \code{ignore.samples = FALSE}, all objects must have identical sample names
#'   \item For multi-sample analysis, set \code{ignore.samples = TRUE}
#'   \item Use \code{\link{applyConsensus}} to map consensus segments back to individual samples
#'   \item Consider running \code{\link{fillGaps}} before consensus to ensure complete coverage
#' }
#'
#' @seealso
#' \code{\link{fillGaps}} for filling gaps in individual samples before consensus
#'
#' \code{\link{applyConsensus}} for mapping consensus segments back to samples
#'
#' \code{\link{findCommonGR}} for the underlying algorithm that identifies common regions
#'
#' \code{\link{annotateCommonGR}} for annotating common regions with sample-specific values
#'
#' @examples
#' \dontrun{
#' # Compare multiple CNV callers on the same sample
#' cnv_facets <- parseFACETS("sample1_facets.vcf.gz")
#' cnv_purecn <- parsePureCN("sample1_purecn/")
#' cnv_sequenza <- parseSCLUST("sample1_sequenza.seg")
#'
#' # Create consensus using total copy number
#' consensus <- findConsensus(
#'   cnvList = list(
#'     FACETS = cnv_facets,
#'     PureCN = cnv_purecn,
#'     Sequenza = cnv_sequenza
#'   ),
#'   ignore.samples = FALSE,
#'   qual.col = "total.allele.copies"
#' )
#'
#' # Check consensus statistics
#' consensus@segments$stats.majority # Most common copy number
#' consensus@segments$stats.mean # Mean copy number across callers
#'
#' # Multi-sample consensus analysis
#' cnv_list <- list(
#'   sample1 = parseFACETS("sample1.vcf.gz"),
#'   sample2 = parseFACETS("sample2.vcf.gz"),
#'   sample3 = parseFACETS("sample3.vcf.gz")
#' )
#'
#' # Find recurrent alterations across cohort
#' cohort_consensus <- findConsensus(
#'   cnvList = cnv_list,
#'   ignore.samples = TRUE,
#'   qual.col = "cnv.annotation",
#'   method = "qualitative"
#' )
#'
#' # Identify regions with high agreement
#' recurrent <- cohort_consensus@segments[
#'   cohort_consensus@segments$stats.majority.count >= 2
#' ]
#' }
#'
#' @export

findConsensus <- function(cnvList,
                          ignore.samples = FALSE,
                          method = c("quantitative", "qualitative"), # not used
                          qual.col = c(
                            "total.allele.copies",
                            "minor.allele.copies",
                            "cnv.lrr",
                            "cnv.cf",
                            "cnv.ccf",
                            "cnv.markers",
                            "cnv.totalmarkers",
                            "cnv.totalmarkers_density",
                            "cnv.markers_density",
                            "cnv.annotation",
                            "cnv.cnaqc.pass"
                          ),
                          fill.gaps = FALSE, # not used
                          segment.annotation = "INFER", # not used
                          fn = c("majority", "stats") # not used
) {
  method <- match.arg(method)
  if (!is.null(qual.col)) {
    qual.col <- match.arg(qual.col)
  }
  if (is.null(qual.col)) {
    method <- "unprocessed"
    # we will annotate all the columns
    qual.col <- c(
      "total.allele.copies",
      "minor.allele.copies",
      "cnv.lrr",
      "cnv.cf", "cnv.ccf",
      "cnv.markers",
      "cnv.totalmarkers",
      "cnv.totalmarkers_density",
      "cnv.markers_density",
      "cnv.annotation"
    )
    aa <- sapply(cnvList, function(X) {
      "cnv.cnaqc.pass" %in% colnames(S4Vectors::mcols(X@segments))
    })
    if (all(aa == TRUE)) {
      qual.col <- c(qual.col, "cnv.cnaqc.pass")
    }
  } else if (qual.col == "cnv.annotation") {
    method <- "qualitative"
  } else {
    method <- "quantitative"
  }
  logger::log_info("Finding consensus using the column {qual.col} with method {method}", namespace = "cnvstrmliner")
  samples <- sapply(cnvList, function(x) {
    x@sample
  }) %>% unique()
  genomes <- sapply(cnvList, function(x) {
    x@genome
  }) %>% unique()
  sex <- sapply(cnvList, function(x) {
    x@sex
  }) %>% unique()
  if ("F" %in% sex) {
    sex <- "F"
  } else {
    sex <- "M"
  }
  errormesg <- c()
  error <- FALSE
  # we expect to have all objects coming from the same sample.
  if (ignore.samples == FALSE & length(samples) > 1) {
    msg <- "The samples in this list have different names. All the elements of the list should have the same sample name."
    errormesg <- c(errormesg, msg)
    error <- TRUE
  }
  # all samples must be on the same genome
  if (length(genomes) > 1) {
    msg <- "The samples in this list have different genomes All the elements of the list should have the same genome"
    errormesg <- c(errormesg, msg)
    error <- TRUE
  }
  if (error == TRUE) {
    logger::log_error("Find consensus found problems and cannot proceed")
    for (e in errormesg) {
      logger::log_error("{e}")
    }
    stop(paste(errormesg, "\n"))
  }

  # find all the common regions
  segsList <- lapply(cnvList, function(X) {
    X@segments
  })
  consensus <- findCommonGR(segsList)

  # annotate the consensus with the column of interest
  consensus <- annotateCommonGR(segsList, consensus, qual.col)

  # if we annotate a single column we may need to add the proper names to the columns
  # add names to the columns
  # if we expect to have the same sample in all the segments
  # we name the columns by the names of the elements of the segmentsList
  if (method %in% c("quantitative", "qualitative")) {
    if (ignore.samples == FALSE & length(names(cnvList) == length(cnvList))) {
      colnames <- names(cnvList) %>% make.unique()
    }
    # otherwise we name the columns with the samples included in each of the segments object
    if (ignore.samples == TRUE) {
      colnames <- sapply(cnvList, function(x) {
        x@sample
      }) %>% make.unique()
    }
    colnames(S4Vectors::mcols(consensus)) <- colnames
  }
  # add the majority and majority.count columns along with the other statistics

  if (method == "quantitative") {
    count <- rowSums(!is.na(S4Vectors::mcols(consensus)))
    majority <- S4Vectors::mcols(consensus) %>%
      apply(., 1, function(x) {
        median(x, na.rm = TRUE)
      })
    majority.count <- S4Vectors::mcols(consensus) %>%
      apply(., 1, function(x) {
        length(which(x == median(x, na.rm = TRUE)))
      })
    mean <- S4Vectors::mcols(consensus) %>% apply(., 1, function(x) {
      mean(x, na.rm = TRUE)
    })
    median <- S4Vectors::mcols(consensus) %>% apply(., 1, function(x) {
      median(x, na.rm = TRUE)
    })
    sd <- S4Vectors::mcols(consensus) %>% apply(., 1, function(x) {
      sd(x, na.rm = TRUE)
    })
    S4Vectors::mcols(consensus)["stats.mean"] <- mean
    S4Vectors::mcols(consensus)["stats.median"] <- median
    S4Vectors::mcols(consensus)["stats.sd"] <- sd
    S4Vectors::mcols(consensus)["stats.majority"] <- majority
    S4Vectors::mcols(consensus)["stats.majority.count"] <- majority.count
    S4Vectors::mcols(consensus)["stats.count"] <- count
  }

  abundance <- function(x) {
    if (all(is.na(x))) {
      return(list("stats.majority.count" = 0, "stats.majority" = NA))
    }
    x <- x[!is.na(x)]
    x <- factor(x)
    ux <- unique(x)
    tab <- tabulate(match(x, ux))
    m <- max(tab)
    l <- ux[which.max(tab)]
    return(list("stats.majority.count" = m, "stats.majority" = as.character(l)))
  }

  # add the majority and count columns
  if (method == "qualitative") {
    count <- rowSums(!is.na(S4Vectors::mcols(consensus)))
    tl <- S4Vectors::mcols(consensus) %>%
      as.data.frame() %>%
      apply(., 1, function(x) {
        abundance(x)
      }) %>%
      bind_rows()
    # majority=S4Vectors::mcols(consensus)  %>%
    #   apply(. , 1, function(x){ x[max(tabulate(factor(x)))]}  ) %>%unlist()
    # majority.count=S4Vectors::mcols(consensus)  %>%
    #   apply(. , 1, function(x){tx=table(x);  aa=sort(tx,decreasing = TRUE); aa[1]})%>%unlist()
    # majority.count[which(majority<0)]=0
    # majority[which(majority<0)]=NA
    S4Vectors::mcols(consensus)["stats.majority"] <- tl$stats.majority
    S4Vectors::mcols(consensus)["stats.majority.count"] <- tl$stats.majority.count
    S4Vectors::mcols(consensus)["stats.count"] <- count
  }

  purities <- sapply(cnvList, function(x) {
    x@purity
  }) %>% median(., na.rm = TRUE)
  ploidies <- sapply(cnvList, function(x) {
    x@ploidy
  }) %>% median(., na.rm = TRUE)
  methods <- sapply(cnvList, function(x) {
    x@method
  }) %>%
    unique() %>%
    paste(., collapse = ",")


  obj <- cnvstrml(
    segments = consensus,
    purity = purities,
    ploidy = ploidies,
    method = methods,
    sample = samples,
    quality = list("consensus column" = qual.col, "consensus method" = method),
    sex = sex,
    segtype = paste("consensus", method),
    genome = genomes
  )
  obj
}


#' Fill Gaps Between CNV Segments with Neutral Copy Number Regions
#'
#' @description
#' This S3 generic function identifies gaps between called CNV segments and fills
#' them with neutral copy number segments to create complete chromosome coverage.
#' This is essential for ensuring that the entire genome is represented in the
#' segmentation, which is particularly important for consensus analysis and
#' visualization.
#'
#' The function automatically handles sex chromosomes, telomeric and centromeric
#' regions, and adjusts marker counts based on observed marker density in the
#' existing segments.
#'
#' @param x An object containing CNV segments. Can be either a \code{cnvstrml}
#'   object or a \code{GRanges} object.
#' @param ... Additional arguments passed to specific methods. See
#'   \code{\link{fillGaps.cnvstrml}} or \code{\link{fillGaps.GRanges}} for
#'   method-specific parameters.
#'
#' @return An updated object of the same class as the input, with gaps filled
#'   by neutral segments. Filled segments are flagged with \code{cnv.fill = TRUE}.
#'
#' @details
#' This is a generic function that dispatches to specific methods based on the
#' class of the input object. The function is designed to prepare CNV data for
#' downstream analysis by ensuring complete genomic coverage.
#'
#' @seealso
#' \code{\link{fillGaps.cnvstrml}} for the cnvstrml method
#'
#' \code{\link{fillGaps.GRanges}} for the GRanges method
#'
#' \code{\link{findConsensus}} for consensus analysis that benefits from filled gaps
#'
#' \code{\link{applyConsensus}} typically used before filling gaps
#'
#' @export
fillGaps <- function(x, ...) {
  UseMethod("fillGaps")
}


#' Fill Gaps Between CNV Segments (cnvstrml method)
#'
#' @description
#' This method identifies gaps between called CNV segments in a \code{cnvstrml}
#' object and fills them with neutral copy number segments to ensure complete
#' chromosome coverage. This is particularly useful for consensus analysis and
#' ensuring that all genomic regions are represented in downstream analyses.
#'
#' The function automatically removes telomeric and centromeric regions, adjusts
#' marker counts based on observed density, and handles sex chromosomes appropriately.
#' All parameters are automatically extracted from the \code{cnvstrml} object unless
#' explicitly provided.
#'
#' @param x A \code{cnvstrml} object containing CNV segments to be gap-filled.
#' @param ploidy Integer. The expected total copy number for neutral regions.
#'   If \code{NULL} (default), uses the effective ploidy from the \code{cnvstrml}
#'   object (\code{x@effective_ploidy}). For diploid genomes, this is typically 2.
#' @param major.allele.copies Integer. The expected number of copies of the major
#'   allele in neutral regions. If \code{NULL} (default), automatically set to
#'   \code{ceiling(ploidy/2)} to create a balanced heterozygous state. For diploid
#'   neutral regions, this would be 1 (creating 1:1 major:minor).
#' @param gapsize Integer. Minimum gap size in base pairs to fill. Only gaps
#'   larger than this threshold will be filled with neutral segments. Default is
#'   30,000,000 (30 Mb) to avoid creating segments for small inter-segment spaces.
#' @param segment.annotation Character. The annotation label to assign to newly
#'   created segments in the \code{cnv.annotation} column. Options:
#'   \itemize{
#'     \item \code{"NEUT"} - Marks segments as neutral copy number (default)
#'     \item \code{"INFER"} - Marks segments as inferred/imputed
#'   }
#' @param ... Additional arguments passed through for S3 generic consistency.
#'
#' @return A \code{cnvstrml} object with filled gaps. Filled segments have:
#'   \code{cnv.fill = TRUE}, appropriate copy numbers, and estimated marker counts.
#'
#' @details
#' The function performs gap filling, excludes telomeric/centromeric regions,
#' handles sex chromosomes appropriately, and estimates marker counts based on
#' average density. It is recommended to run this \emph{after} \code{\link{applyConsensus}}.
#'
#' @seealso
#' \code{\link{fillGaps.GRanges}} for the GRanges method
#'
#' \code{\link{findConsensus}} for consensus analysis
#'
#' \code{\link{applyConsensus}} typically used before filling gaps
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' cnv <- parseFACETS("sample.vcf.gz")
#' cnv_filled <- fillGaps(cnv)
#'
#' # With custom parameters
#' cnv_filled <- fillGaps(cnv,
#'   gapsize = 10000000,
#'   segment.annotation = "INFER"
#' )
#' }
#'
#' @export
fillGaps.cnvstrml <- function(x,
                              ploidy = NULL,
                              major.allele.copies = NULL,
                              gapsize = 30000000,
                              segment.annotation = c("NEUT", "INFER"),
                              ...) {
  cnvObj <- x
  if (!is(cnvObj, "cnvstrml")) {
    stop("The input to the function fillGaps should be a cnvstrml object")
  }
  logger::log_warn("Filling in gaps. This will modify the segments by introducing segments in the gaps larger than {gapsize} annotated as {segment.annotation}.", namespace = "cnvstrmliner")
  segment.annotation <- match.arg(segment.annotation)
  segs <- cnvObj@segments
  # find the ploidy to apply
  if (is.null(ploidy)) {
    ploidy <- cnvObj@effective_ploidy
  }

  # do we need to update the markers based on the density of markers?
  new_segments <- fillGaps(
    segs,
    ploidy = ploidy,
    major.allele.copies = major.allele.copies,
    gapsize = gapsize,
    genome = cnvObj@genome,
    sex = cnvObj@sex,
    purity = cnvObj@purity,
    segment.annotation = segment.annotation
  )
  GenomeInfoDb::seqlevels(segs) <- GenomeInfoDb::seqlevels(new_segments)
  GenomeInfoDb::seqinfo(new_segments) <- GenomeInfoDb::seqinfo(segs)
  cnvObj@segments <- new_segments
  return(cnvObj)
}


#' Fill Gaps Between CNV Segments (GRanges method)
#'
#' @description
#' This method identifies gaps between segments in a \code{GRanges} object and
#' fills them with neutral copy number segments. This is the low-level implementation
#' used by \code{\link{fillGaps.cnvstrml}}, providing direct control over all
#' parameters including genome build, sex, and purity.
#'
#' @param x A \code{GRanges} object containing CNV segments with appropriate
#'   metadata columns. If marker density columns exist, they will be used to
#'   estimate marker counts for filled segments.
#' @param ploidy Integer. Expected total copy number for neutral regions. Must be
#'   specified explicitly. Typical values: 2 (diploid), 3 (triploid).
#' @param major.allele.copies Integer. Expected major allele copies in neutral
#'   regions. If \code{NULL} (default), set to \code{ceiling(ploidy/2)}.
#' @param gapsize Integer. Minimum gap size in base pairs to fill. Default is
#'   30,000,000 (30 Mb).
#' @param genome Character. Genome build identifier (e.g., \code{"hg38"},
#'   \code{"hg19"}). Used to load telomere, centromere, and short arm coordinates.
#'   Default is \code{"hg38"}.
#' @param sex Character. Sample sex for sex chromosome handling. \code{"F"} (default)
#'   or \code{"M"}.
#' @param purity Numeric. Sample purity (tumor fraction) between 0 and 1. Default
#'   is 1.
#' @param segment.annotation Character. Annotation label for new segments:
#'   \code{"NEUT"} (default) or \code{"INFER"}.
#' @param ... Additional arguments passed through for S3 generic consistency.
#'
#' @return A \code{GRanges} object with filled gaps, sorted by genomic coordinates.
#'   Filled segments have \code{cnv.fill = TRUE} and estimated marker counts.
#'
#' @details
#' This low-level method requires explicit specification of all parameters.
#' It excludes telomeric/centromeric regions and handles sex chromosomes
#' appropriately. Marker counts are estimated from average density if density
#' columns are present.
#'
#' @seealso
#' \code{\link{fillGaps.cnvstrml}} for the higher-level cnvstrml method
#'
#' \code{\link{fillGaps}} for the generic function
#'
#' @examples
#' \dontrun{
#' library(GenomicRanges)
#' segs <- GRanges(
#'   seqnames = c("chr1", "chr1"),
#'   ranges = IRanges(
#'     start = c(1000000, 100000000),
#'     end = c(50000000, 150000000)
#'   )
#' )
#' mcols(segs)$total.allele.copies <- c(2, 3)
#' mcols(segs)$major.allele.copies <- c(1, 2)
#' mcols(segs)$minor.allele.copies <- c(1, 1)
#'
#' filled <- fillGaps(segs,
#'   ploidy = 2,
#'   genome = "hg38",
#'   sex = "F"
#' )
#' }
#'
#' @export
fillGaps.GRanges <- function(x,
                             ploidy = NULL,
                             major.allele.copies = NULL,
                             gapsize = 30000000,
                             genome = "hg38",
                             sex = "F",
                             purity = 1,
                             segment.annotation = c("NEUT", "INFER"),
                             ...) {
  segs <- x
  segment.annotation <- match.arg(segment.annotation, several.ok = FALSE)
  if ("cnv.markers_density" %in% colnames(S4Vectors::mcols(segs)) &
    "cnv.totalmarkers_density" %in% colnames(S4Vectors::mcols(segs))) {
    density <- TRUE
    average_totaldensity <- mean(S4Vectors::mcols(segs)[, "cnv.totalmarkers_density"], na.rm = TRUE)
    average_density <- mean(S4Vectors::mcols(segs)[, "cnv.markers_density"], na.rm = TRUE)
  } else {
    density <- FALSE
    average_density <- 0
    average_totaldensity <- 0
  }
  logger::log_debug("fillgaps will update the markers and marker density for the segments wwith average density :{average_density}/{average_total_density} for markers/total markers.", namespace = "cnvstrmliner")

  logger::log_debug("fillgaps will add the new segments with ploidy: {ploidy}", namespace = "cnvstrmliner")
  if (is.null(major.allele.copies)) {
    major.allele.copies <- ceiling(ploidy / 2)
  }
  logger::log_debug("fillgaps will add the new segments with major allele copy number: {major.allele.copies}", namespace = "cnvstrmliner")


  # find the gaps we need to fill in
  # this covers the full genome outside the segments
  # which include centromeres and telomeres
  gaps <- GenomicRanges::gaps(segs)
  gaps <- gaps[BiocGenerics::strand(gaps) == "*"]
  gaps <- BiocGenerics::subset(gaps, width > gapsize)
  if (length(gaps) == 0) {
    return(segs)
  }

  # add the mcols we need in the gaps
  S4Vectors::mcols(gaps) <- gaps %>%
    data.frame() %>%
    dplyr::mutate(
      total.allele.copies = dplyr::case_when(
        (sex == "M" & seqnames %in% c("Y", "chrY")) ~ 1,
        (sex == "M" & seqnames %in% c("X", "chrX")) ~ ploidy - 1,
        .default = ploidy
      ),
      major.allele.copies = major.allele.copies,
      minor.allele.copies = total.allele.copies - major.allele.copies,
      cnv.totalmarkers = ifelse(density == TRUE, round((end - start + 1) * average_totaldensity), 0),
      cnv.markers = ifelse(density == TRUE, round((end - start + 1) * average_density), 0),
      cnv.cf = purity,
      cnv.ccf = 1,
      cnv.annotation = "INFER",
      cnv.loh = FALSE,
      cnv.lrr = 0,
      cnv.fill = TRUE
    ) %>%
    dplyr::select(-c(seqnames, start, end, width, strand))

  if (density == TRUE) {
    S4Vectors::mcols(gaps) <- S4Vectors::mcols(gaps) %>%
      data.frame() %>%
      dplyr::mutate(
        cnv.markers_density = average_density,
        cnv.totalmarkers_density = average_totaldensity
      )
  }
  S4Vectors::mcols(segs)["cnv.fill"] <- FALSE

  # let's create now a combined GRanges with the segments and the gaps
  # this contains the known segments and the inferred segments (including centromeres and telomeres)
  new_segments <- c(segs, gaps) %>%
    GenomeInfoDb::sortSeqlevels() %>%
    sort()


  # remove known gaps
  # first find the known gaps (centromere, telomere, short arms)
  # and then reverse the gaps to find the all_seq = all the regions that are covered by sequencing
  logger::log_info("Filling in {length(gaps)} gaps larger than {gapsize} bp for {genome}", namespace = "cnvstrmliner")
  gaps.centromere <- chrom_gaps[[genome]][["centromere"]]
  gaps.telomere <- chrom_gaps[[genome]][["telomere"]]
  gaps.short_arms <- chrom_gaps[[genome]][["short_arms"]]
  logger::log_info("Found {length(gaps.centromere)} centromeric regions, {length(gaps.telomere)} telomeric regions and {length(gaps.short_arms)} short arm regions to exclude from the filled segments.", namespace = "cnvstrmliner")
  # we reverse the gaps to find the all_seq = all the regions that are covered by sequencing
  all_gaps <- c(
    GenomicRanges::trim(gaps.centromere),
    GenomicRanges::trim(gaps.telomere),
    GenomicRanges::trim(gaps.short_arms)
  ) %>%
    GenomeInfoDb::sortSeqlevels() %>%
    sort()
  all_seq <- GenomicRanges::gaps(all_gaps)
  all_seq <- all_seq[BiocGenerics::strand(all_seq) == "*"]


  # now we can overlap the new_segments with the all_seq
  # this will give us all the segments & previously identified gaps
  # in sequencable regions of the genome, and exclude the known gaps (centromeres, telomeres, short arms)
  ov <- GenomicRanges::findOverlaps(query = all_seq, subject = new_segments)
  # find the unique regions in teh all_seqs (the regions that are covered by sequencing except the known gaps)
  new_segments <- GenomicRanges::pintersect(new_segments[S4Vectors::subjectHits(ov)], all_seq[S4Vectors::queryHits(ov)])
  rm(list = c("gaps.centromere", "gaps.telomere", "gaps.short_arms"))
  rm_idx <- which(BiocGenerics::width(new_segments) < gapsize & S4Vectors::mcols(new_segments)[, "cnv.annotation"] == "INFER")

  if (length(rm_idx) > 0) {
    new_segments <- new_segments[-rm_idx]
  }
  new_segments$cnv.annotation[which(new_segments$cnv.annotation == "INFER")] <- segment.annotation
  return(new_segments)
}


#' Map Consensus Segments Back to Individual Sample
#'
#' @description
#' This function maps consensus segments (generated by \code{\link{findConsensus}})
#' back to an individual sample's CNV segments. It splits the sample's segments
#' according to consensus breakpoints and annotates which regions are part of the
#' consensus versus sample-specific (unique) alterations.
#'
#' The function creates a refined segmentation where each segment is classified as
#' either shared with the consensus (annotated with the consensus majority value) or
#' unique to the sample. This enables identification of recurrent alterations across
#' a cohort and distinguishes them from sample-specific events.
#'
#' @param cnvObj A \code{cnvstrml} object containing CNV segments from an individual
#'   sample or CNV caller. This object will be updated with consensus annotations.
#'   All original segment metadata is preserved.
#' @param consensus A \code{cnvstrml} object containing consensus segments, typically
#'   generated by \code{\link{findConsensus}}. The consensus object must contain
#'   a \code{stats.majority} column (for quantitative data) or \code{majority} column
#'   (for qualitative data) in its segment metadata.
#'
#' @return A \code{cnvstrml} object with updated segments. The returned object contains:
#'   \describe{
#'     \item{Refined segmentation}{Segments are split at consensus breakpoints,
#'       creating alignment between sample and consensus coordinate systems}
#'     \item{\code{consensus} metadata column}{Added to \code{@segments}:
#'       \itemize{
#'         \item For segments overlapping consensus: contains the consensus
#'           \code{stats.majority} or \code{majority} value
#'         \item For sample-specific segments: \code{NA}
#'       }}
#'     \item{Preserved metadata}{All original segment annotations remain intact}
#'     \item{Sorted coordinates}{Segments sorted by genomic position}
#'   }
#'
#' @details
#' \strong{Algorithm:}
#'
#' \enumerate{
#'   \item Identifies overlaps between sample and consensus segments
#'   \item Splits sample segments at consensus breakpoints using \code{pintersect}
#'   \item Annotates consensus regions with majority values
#'   \item Identifies sample-specific regions using \code{setdiff}
#'   \item Combines and sorts all segments
#' }
#'
#' \strong{Segment Types:}
#'
#' \itemize{
#'   \item \strong{Consensus segments}: Overlap with consensus, have \code{consensus}
#'     column populated (represent recurrent/validated alterations)
#'   \item \strong{Unique segments}: Not in consensus, \code{consensus = NA}
#'     (represent sample-specific events)
#' }
#'
#' \strong{Use Cases:}
#'
#' \itemize{
#'   \item \strong{Method comparison}: Identify concordant vs caller-specific calls
#'   \item \strong{Cohort analysis}: Distinguish recurrent from private alterations
#'   \item \strong{Quality filtering}: Prioritize consensus-supported regions
#' }
#'
#' \strong{Workflow:}
#'
#' Recommended order:
#' \preformatted{
#' consensus <- findConsensus(cnv_list)
#' cnv <- applyConsensus(cnv, consensus)  # Before fillGaps
#' cnv <- fillGaps(cnv)
#' }
#'
#' Run \emph{before} \code{\link{fillGaps}} to avoid splitting filled neutral segments.
#'
#' @note
#' \itemize{
#'   \item Requires \code{stats.majority} or \code{majority} column in consensus
#'   \item Preserves all original metadata during splitting
#'   \item Modifies segment boundaries but not copy number values
#'   \item Sample-specific regions have \code{consensus = NA}
#' }
#'
#' @seealso
#' \code{\link{findConsensus}} for generating consensus from multiple samples
#'
#' \code{\link{fillGaps}} typically used after applying consensus
#'
#' \code{\link{findCommonGR}} for identifying common regions
#'
#' @examples
#' \dontrun{
#' # Method comparison
#' facets <- parseFACETS("sample.vcf.gz")
#' purecn <- parsePureCN("sample_purecn/")
#'
#' consensus <- findConsensus(
#'   list(FACETS = facets, PureCN = purecn),
#'   ignore.samples = FALSE
#' )
#'
#' facets <- applyConsensus(facets, consensus)
#' purecn <- applyConsensus(purecn, consensus)
#'
#' # Check consensus coverage
#' table(is.na(facets@segments$consensus))
#'
#' # Cohort analysis
#' cnv_list <- list(
#'   s1 = parseFACETS("s1.vcf.gz"),
#'   s2 = parseFACETS("s2.vcf.gz")
#' )
#' cohort_consensus <- findConsensus(cnv_list, ignore.samples = TRUE)
#' cnv_list <- lapply(cnv_list, applyConsensus, consensus = cohort_consensus)
#'
#' # Find sample-specific alterations
#' unique_segs <- cnv_list$s1@segments[is.na(cnv_list$s1@segments$consensus)]
#' }
#'
#' @export
applyConsensus <- function(cnvObj, consensus) {
  # map them back to each element of the segmentList to get the
  S4Vectors::mcols(cnvObj@segments)["consensus"] <- NA

  # find the segments that are shared between the CNV segments and the consensus
  ov <- GenomicRanges::findOverlaps(query = cnvObj@segments, subject = consensus@segments)
  common_segments <- GenomicRanges::pintersect(cnvObj@segments[S4Vectors::queryHits(ov)], consensus@segments[S4Vectors::subjectHits(ov)])

  # add the annotation of the consensus segments
  ov3 <- GenomicRanges::findOverlaps(query = common_segments, subject = consensus@segments)
  S4Vectors::mcols(common_segments)[S4Vectors::queryHits(ov3), "consensus"] <- S4Vectors::mcols(consensus@segments)[S4Vectors::subjectHits(ov3), "majority"]

  # find the unique segments, i.e. the regions that are not covered by the consensus
  unique_segments <- GenomicRanges::setdiff(cnvObj@segments, consensus@segments)
  ov2 <- GenomicRanges::findOverlaps(query = unique_segments, subject = cnvObj@segments)
  unique_segments <- GenomicRanges::pintersect(
    cnvObj@segments[S4Vectors::subjectHits(ov2)],
    unique_segments[S4Vectors::queryHits(ov2)]
  )


  S4Vectors::mcols(unique_segments) <- mcols(cnvObj@segments)[S4Vectors::subjectHits(ov2), ]
  segs <- c(common_segments, unique_segments) %>%
    GenomeInfoDb::sortSeqlevels() %>%
    sort()
  GenomeInfoDb::seqinfo(segs) <- GenomeInfoDb::seqinfo(cnvObj@segments)
  cnvObj@segments <- segs

  return(cnvObj)
}


#' get the consensus sequence from a list of GRanges objects that represent the CNV segments
#'
#' The consensus is a new Granges that contains only the shared segments
#' The metadata of the consensus have the cnv.annotations for each of the elements of the segments
#'
#' @param segmentsList A list of GRanges objects.
#' @param method Character value of either 'qualitative' or 'quantitative'.
#' @param qual.col Character metadata column used for quantitative consensus.
#' @param fn Character function name to apply in the quantitative consensus case.
#' 'qualitative' it returns an object with qualitative information about the segments (e.g NEUT, AMP etc)
#' If it is set to 'quantitative' it retuns an object with the actual number of copy numbers.
#' @returns a Granges with all the segments that are shared across all samples.


findConsensusGR <- function(segmentsList, method = c("qualitative", "quantitative"), qual.col = "total.allele.copies", fn = c("majority", "stats")) {
  # find the regions that are shared across all predictions
  method <- match.arg(method)
  fn <- match.arg(fn)
  # get the shared ranges
  consensus <- segmentsList %>% Reduce(function(x, y) {
    ov <- GenomicRanges::findOverlaps(query = x, subject = y)
    gr <- GenomicRanges::pintersect(x[S4Vectors::queryHits(ov)], y[S4Vectors::subjectHits(ov)])
    keep <- which(BiocGenerics::width(gr) > 1)
    gr[keep]
  }, .)
  S4Vectors::mcols(consensus) <- NULL
  GenomeInfoDb::seqinfo(consensus) <- GenomeInfoDb::seqinfo(segmentsList[[1]])
  # add qualitative information to the consensus
  if (method == "qualitative") {
    metadata <- segmentsList %>%
      lapply(function(X) {
        d <- rep(NA, length(consensus))
        ov <- GenomicRanges::findOverlaps(query = consensus, subject = X)
        d[S4Vectors::queryHits(ov)] <- X[S4Vectors::subjectHits(ov)]$cnv.annotation
        d
      }) %>%
      dplyr::bind_cols()

    majority <- metadata %>% apply(., 1, function(x) {
      tx <- table(x)
      tx["INFER"] <- -1
      aa <- sort(tx, decreasing = TRUE)
      names(aa)[1]
    })
    majority.count <- metadata %>% apply(., 1, function(x) {
      tx <- table(x)
      tx["INFER"] <- -1
      aa <- sort(tx, decreasing = TRUE)
      aa[1]
    })
  }

  # add quantitative information to the consensus
  if (method == "quantitative") {
    metadata <- segmentsList %>%
      lapply(function(X) {
        d <- rep(NA, length(consensus))
        ov <- GenomicRanges::findOverlaps(query = consensus, subject = X)
        d[S4Vectors::queryHits(ov)] <- S4Vectors::mcols(X)[S4Vectors::subjectHits(ov), qual.col]
        infers <- S4Vectors::mcols(X)[S4Vectors::subjectHits(ov), "cnv.annotation"]
        bad.infers <- which(infers == "INFER")
        d[bad.infers] <- NA
        d
      }) %>%
      dplyr::bind_cols()
    if (fn == "majority") {
      majority <- metadata %>% apply(., 1, function(x) {
        median(x, na.rm = TRUE)
      })
      majority.count <- metadata %>% apply(., 1, function(x) {
        length(which(x == median(x, na.rm = TRUE)))
      })
    }
    if (fn == "stats") {
      mean <- metadata %>% apply(., 1, function(x) {
        mean(x, na.rm = TRUE)
      })
      median <- metadata %>% apply(., 1, function(x) {
        median(x, na.rm = TRUE)
      })
      sd <- metadata %>% apply(., 1, function(x) {
        sd(x, na.rm = TRUE)
      })
    }
  }

  if (!is.null(names(segmentsList))) {
    colnames(metadata) <- names(segmentsList)
  }
  S4Vectors::mcols(consensus) <- metadata
  if (fn == "majority") {
    S4Vectors::mcols(consensus)$majority <- majority
    S4Vectors::mcols(consensus)$majority.count <- majority.count
  }
  if (fn == "stats") {
    S4Vectors::mcols(consensus)$mean <- mean
    S4Vectors::mcols(consensus)$median <- median
    S4Vectors::mcols(consensus)$sd <- sd
  }

  consensus
}


#' get the common segments from a list of GRanges objects that represent the CNV segments
#'
#' The common is a new Granges that contains all the shared segments.
#' When the segments are processed the function gaps is used to fill regions of the genome
#' that are not included in individual samples.
#' The metadata of the consensus are set to NULL
#'
#' @param segsList A list of GRanges objects representing segments.
#' @returns a Granges with all the segments that are shared across all samples.
#'
findCommonGR <- function(segsList) {
  segmentsList.gaps <- lapply(segsList, function(X) {
    g <- GenomicRanges::gaps(X)
    d <- c(X, g) %>%
      GenomeInfoDb::sortSeqlevels() %>%
      BiocGenerics::sort()
    d[BiocGenerics::strand(d) == "*"]
  })
  # now we find the common GRanges
  common <- segmentsList.gaps %>% Reduce(function(x, y) {
    ov <- GenomicRanges::findOverlaps(query = x, subject = y)
    gr <- GenomicRanges::pintersect(x[S4Vectors::queryHits(ov)], y[S4Vectors::subjectHits(ov)])
    keep <- which(BiocGenerics::width(gr) > 1)
    gr[keep]
  }, .)
  S4Vectors::mcols(common) <- NULL
  names(common) <- NULL
  GenomeInfoDb::seqinfo(common) <- GenomeInfoDb::seqinfo(segmentsList.gaps[[1]])
  if (length(common) == 0) {
    warning("There are no common segments for these samples. This can cause major problems in some downstream analyses.")
  }
  common
}


#' annotate the common regions with columns from the constituent segments
#'
#' @param segmentsList A list of GRanges objects, the same used to generate the common regions.
#' @param common A GRanges object that includes the segments to annotate.
#' @param columns Character vector of columns that exist in segmentsList and will be added in the common object.
#' @returns an updated common object with the values of the column in the individual samples.
annotateCommonGR <- function(segmentsList, common, columns) {
  counter <- 1
  metadata <- lapply(1:length(segmentsList), function(n) {
    X <- segmentsList[[n]]
    name <- names(segmentsList)[n]
    label <- paste0(name, "||", columns)
    if (is.null(name)) {
      name <- paste0("sample_", counter)
    }
    d <- matrix(nrow = length(common), ncol = length(columns)) %>% as.data.frame() # rows = the length of common. #cols is the columns we get from each segList
    ov <- GenomicRanges::findOverlaps(query = common, subject = X)
    d[S4Vectors::queryHits(ov), ] <- S4Vectors::mcols(X)[S4Vectors::subjectHits(ov), columns] %>%
      as.data.frame()
    d <- magrittr::set_colnames(d, label)
    counter <<- counter + 1
    d
  }) %>% dplyr::bind_cols()
  # if(! is.null(names(segmentsList))){
  #   colnames( metadata)=names(segmentsList)
  # }
  S4Vectors::mcols(common) <- metadata
  common
}


#' Remove regions that are known gaps e.g. telomeres and centromeres
#'
#' @param cnvObj (cnvstremlinr) an object with segments
#' @returns an updated object without the telo/centromeric regions

removeKnownGaps <- function(cnvObj) {
  if (!is(cnvObj, "cnvstrml")) {
    stop("The input to the function removeKnownGaps should be a cnvstrml object")
  }
  logger::log_warn("Removing regions mapping on Centromeres and Telomeres. This will change the segments of the object.", namespace = "cnvstrmliner")
  segs <- cnvObj@segments
  # do we need to update the markers based on the density of markers?

  # load known gaps of the sequence assembly
  # these are the objects
  # gaps.centromere, gaps.short_arms, gaps.telomeres
  # load(file=file.path(
  #   system.file('extdata', package = 'cnvstrmliner'),
  #   'seqinfo',
  #   paste0( cnvObj@genome, ".gaps.rds"))
  # )

  # make continuous gaps for the centromere
  gaps.telomere <- chrom_gaps[[cnvObj@genome]][["telomere"]]
  gaps.centromere <- chrom_gaps[[cnvObj@genome]][["centromere"]]
  gaps.centromere <- GenomicRanges::reduce(gaps.centromere, min.gapwidth = 1e6)
  gaps.short_arms <- chrom_gaps[[cnvObj@genome]][["short_arms"]]

  # remove known gaps
  # we reverse the gaps to find the all_seq = all the regions that are covered by sequencing
  all_gaps <- c(
    GenomicRanges::trim(gaps.centromere),
    GenomicRanges::trim(gaps.telomere),
    GenomicRanges::trim(gaps.short_arms)
  ) %>%
    GenomeInfoDb::sortSeqlevels() %>%
    sort()
  all_seq <- GenomicRanges::gaps(all_gaps)
  all_seq <- all_seq[BiocGenerics::strand(all_seq) == "*"]


  ov <- GenomicRanges::findOverlaps(query = all_seq, subject = new_segments)
  new_segments <- GenomicRanges::pintersect(new_segments[S4Vectors::subjectHits(ov)], all_seq[S4Vectors::queryHits(ov)])
  rm(list = c("gaps.centromere", "gaps.telomere", "gaps.short_arms"))

  cnvObj@segments <- new_segments
  GenomeInfoDb::seqinfo(cnvObj@segments) <- GenomeInfoDb::seqinfo(cnvObj@segments)
  return(cnvObj)
}
