# Metrics computed on a single cnvstrml object (or each object of a list) independently.



#' Calculate the entropy of copy number aberration patterns across samples
#'
#' @description
#' This function computes the conditional entropy for each sample's copy number profile
#' as a quality metric that quantifies the complexity and fragmentation of the CNV landscape.
#' Entropy, borrowed from information theory, measures the uncertainty or randomness in
#' the distribution of copy number states across genomic segments.
#'
#' **What is Entropy?**
#'
#' Entropy is a fundamental concept from information theory that quantifies the amount of
#' uncertainty or information content in a probability distribution. In the context of
#' copy number analysis, it measures how "spread out" or unpredictable the copy number
#' states are across the genome. Mathematically, the conditional entropy H(Y|X) is computed as:
#'
#' \deqn{H(Y|X) = -\sum_{x,y} p(x,y) \log p(y|x)}
#'
#' where p(x,y) is the joint probability distribution of the discretized copy number values,
#' and p(y|x) is the conditional probability.
#'
#' **How is it computed?**
#'
#' For each sample in the consensus object, the computation follows these steps:
#' \enumerate{
#'   \item Copy number values across consensus segments are extracted
#'   \item Values are discretized into bins using equal-frequency discretization
#'   \item The conditional entropy is calculated using the discretized distribution
#'   \item The result quantifies the information content of the copy number pattern
#' }
#'
#' **Interpreting Entropy Values in CNV Context:**
#'
#' \itemize{
#'   \item **Low entropy (close to 0):** Indicates a simple, uniform CNV profile with few
#'         distinct copy number states. This typically represents:
#'         \itemize{
#'           \item High-quality samples with clean signals
#'           \item Genomes with minimal copy number aberrations (e.g., diploid baseline)
#'           \item Consistent, large-scale CNV events (e.g., whole-chromosome gains/losses)
#'         }
#'   \item **Moderate entropy:** Suggests a moderate level of genomic complexity with
#'         several distinct copy number states distributed across the genome. This is
#'         typical for:
#'         \itemize{
#'           \item Cancer samples with focal amplifications and deletions
#'           \item Good quality samples with genuine biological complexity
#'         }
#'   \item **High entropy:** Indicates a highly fragmented or noisy CNV profile with many
#'         different copy number states. This can signal:
#'         \itemize{
#'           \item Over-segmentation due to technical noise
#'           \item Poor quality samples with unreliable signals
#'           \item Highly unstable genomes (e.g., chromothripsis)
#'           \item Potential artifacts from sequencing or computational methods
#'         }
#' }
#'
#' **Using Entropy as a Quality Metric:**
#'
#' Entropy serves as an important QC metric in CNV analysis:
#' \itemize{
#'   \item **Quality filtering:** Samples with abnormally high entropy may indicate technical
#'         problems (noise, over-segmentation) and should be flagged for review
#'   \item **Method comparison:** When comparing different CNV calling methods, entropy helps
#'         identify which methods produce more stable, biologically plausible segmentations
#'   \item **Cohort analysis:** Outlier samples with extreme entropy values relative to the
#'         cohort distribution may require special attention
#'   \item **Complementary to other QC metrics:** Should be interpreted alongside metrics like
#'         segment count, segment size distribution, and purity/ploidy estimates
#'   \item **Batch effect detection:** Systematic differences in entropy across batches may
#'         indicate technical variability rather than biological differences
#' }
#'
#' Note: The input must be a consensus object (created by \code{\link{findConsensus}}) to
#' ensure all samples are evaluated on a common set of genomic segments, making entropy
#' The cnvObject must have \code{segtype}
#'   starting with "consensus".
#' If a single sample \code{cnvstrml} object is provided, it is advised to first intersect
#' the segments of the object with a reference set of segments to ensure meaningful
#' entropy calculation and enable comparisons across samples, this can be achieved by
#' providing the targets argument with the reference set of segments (e.g. exome coordinates).
#'
#' @param cnvObj Either a \code{cnvstrml} object of consensus type, created by applying
#'   \code{\link{findConsensus}} to a list of CNV objects,
#'   or a \code{cnvstrml} of a single sample. In this case the entropy will be calculated
#'   only for that sample.
#' @param targets (optional) Either \code{GRanges} object specifying genomic regions to focus
#'   the entropy calculation on or \code{character} pointing to a bed file with the genomic
#'   regions - in this case the bed file will be parsed with the rtracklayer::import.bed.
#'   If provided, only segments overlapping these regions will
#'   be considered. This is useful for restricting the analysis to exonic regions or other
#'   areas of interest.
#' @return A named numeric vector containing the conditional entropy value for each sample
#'   in the consensus object. Higher values indicate more complex/fragmented CNV patterns.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Create consensus from a list of CNV objects
#' consensus_obj <- findConsensus(cnvList, method = "quantitative", ignore.samples = TRUE)
#'
#' # Calculate entropy for all samples
#' entropy_values <- cnvEntropy(consensus_obj)
#'
#' # Identify potential quality issues
#' cutoff <- median(entropy_values) + 2 * sd(entropy_values)
#' high_entropy_samples <- names(entropy_values[entropy_values > cutoff])
#'
#'
#' # # Calculate entropy for a single sample CNV object
#' bed.file <- system.file("extdata/example_CNV/BRCA_XX_01/S07604514_Regions_v6_hg38_Targets.bed.gz", package = "cnvstrmliner")
#'
#' cnvEntropy(cnv, bed.file)
#' }
#'
#' @seealso
#' \code{\link{CNV_entropy_sample}} for the underlying calculation on individual samples
#' \code{\link{findConsensus}} for creating consensus objects
#' \code{\link{cnvMetrics}} for other CNV comparison metrics
#'
cnvEntropy <- function(cnvObj, targets = NULL) {
  # Validate cnvObj
  if (is.list(cnvObj)) {
    # Check if all elements in the list are cnvstrml objects
    if (!all(sapply(cnvObj, function(x) inherits(x, "cnvstrml")))) {
      stop("All elements in cnvObj list must be cnvstrml objects")
    }
  } else {
    # Check if single object is cnvstrml
    if (!inherits(cnvObj, "cnvstrml")) {
      stop("cnvObj must be a cnvstrml object or a list of cnvstrml objects")
    }
  }

  # Validate targets if provided
  if (!is.null(targets)) {
    if (!is.character(targets) && !inherits(targets, "GRanges")) {
      stop("targets must be a GRanges object or a character string (file path)")
    }
  }

  # if(!is.list( cnvList )){
  #   cnvList=list( cnvList )
  # }
  # # create a dataframe with shared segments across all samples of the cohort
  # consensus.quant=suppressMessages( findConsensus( cnvList, method="quantitative", ignore.samples = TRUE) )

  if (!startsWith(cnvObj@segtype, prefix = "consensus")) {
    if (is.null(targets)) {
      logger::log_warn("cnvEntropy is designed to work on consensus CNV objects or
      a single object with known targets. ")
      return(NA)
    }
    # load the bed file first if targets is a character
    if (is.character(targets)) {
      if (!file.exists(targets)) {
        logger::log_error("The provided targets file does not exist: {targets}")
        return(NA)
      }
      targets <- rtracklayer::import.bed(targets)
    }
    ov <- GenomicRanges::findOverlaps(query = cnvObj@segments, subject = targets)
    S4Vectors::mcols(targets) <- S4Vectors::DataFrame(matrix(NA,
      nrow = length(targets),
      ncol = ncol(S4Vectors::mcols(cnvObj@segments))
    ))
    colnames(S4Vectors::mcols(targets)) <- colnames(S4Vectors::mcols(cnvObj@segments))
    S4Vectors::mcols(targets)[S4Vectors::subjectHits(ov), ] <- S4Vectors::mcols(cnvObj@segments)[S4Vectors::queryHits(ov), ]

    df <- S4Vectors::mcols(targets)["total.allele.copies"] %>%
      as.data.frame() %>%
      tidyr::drop_na() %>%
      magrittr::set_colnames(cnvObj@sample)
    return(CNV_entropy_sample(df[[1]]))
  }

  dat <- cnvObj@segments %>%
    as.data.frame(row.names = NULL) %>%
    magrittr::set_colnames(c("seqnames", "start", "end", "width", "strand", colnames(S4Vectors::mcols(cnvObj@segments))))


  entropy <- dat %>%
    dplyr::select(-c("seqnames", "start", "end", "width", "strand", tidyselect::starts_with("stats."))) %>%
    apply(., 2, CNV_entropy_sample)

  return(entropy)
}

#' Calculate the conditional entropy for a list of segments
#'
#' @param X (vector) vector of int  . They correspond to the copy number of a sample
#' @returns the conditional entropy of a sample
#' @description
#'
#' Returns the conditional entropy of a sample.
#' The input is expected to be a vector of numbers that indicate the segment copy number.
#' In order to make sense it is better if all the samples are compared on the same set
#' of segments, thus the function cnvEntropy first finds the consensus of a set of
#' samples.


CNV_entropy_sample <- function(X) {
  tryCatch(
    {
      X %>%
        infotheo::discretize() %>%
        infotheo::condentropy()
    },
    error = function(x) {
      logger::log_error("Error calculating entropy for this sample ", e)
      return(NA)
    }
  )
}


#' Calculate the Gini Coefficient for Coverage Distribution
#'
#' @description
#' Computes the Gini coefficient as a measure of coverage uniformity across
#' genomic targets or bins. The Gini coefficient quantifies inequality in the
#' distribution of sequencing depth, where values closer to 0 indicate uniform
#' coverage and values closer to 1 indicate highly unequal coverage distribution.
#'
#' **Note:** This function is exported but requires external data in a specific format.
#' It does **not** accept \code{cnvstrmliner} objects. Instead, it expects a data frame
#' with pre-processed coverage data (e.g., from CNVkit .cnr files or similar formats).
#'
#' @param cnr_data A data frame containing coverage information with the following
#'   required column:
#'   \describe{
#'     \item{depth}{Numeric vector of read depth or coverage values for each
#'       genomic bin or target region}
#'   }
#'   Rows with missing (NA) or zero depth values are automatically excluded from
#'   the calculation.
#' @param use_log Logical. If \code{TRUE}, the function will compute the Gini coefficient using the log2-transformed coverage values.
#' This is not recommended for general use, but may be useful in specific contexts where the log-transformed distribution is more informative.
#' If \code{TRUE}, the function will filter out bins with log2 values <= -10 to exclude unmappable or poorly covered regions.
#' @return A numeric value representing the Gini coefficient (0 to 1), where:
#'   \itemize{
#'     \item Values near 0 indicate uniform coverage distribution (good quality)
#'     \item Values near 1 indicate highly unequal coverage (potential quality issues)
#'   }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Example with CNVkit .cnr data
#' cnr_data <- read.table("sample.cnr", header = TRUE, sep = "\t")
#' gini_score <- computeGini(cnr_data)
#' print(paste("Gini coefficient:", gini_score))
#'
#' # Interpret results
#' if (gini_score < 0.2) {
#'   message("Good coverage uniformity")
#' } else if (gini_score > 0.5) {
#'   message("Poor coverage uniformity - potential QC issue")
#' }
#' }
#'
#' @seealso
#' \code{\link{computeMAPD}} for calculating the MAPD metric
#'
#' \code{\link[ineq]{Gini}} for the underlying Gini coefficient calculation
#'
computeGini <- function(cnr_data, use_log = FALSE) {
  if (use_log) {
    valid_coverage <- cnr_data %>% dplyr::filter(log2 > -10)
    gini_score <- ineq::Gini(valid_coverage$log2)
  } else {
    valid_coverage <- cnr_data %>% dplyr::filter(!is.na(depth) & depth > 0)
    gini_score <- ineq::Gini(valid_coverage$depth)
  }
  return(gini_score)
}
#' Calculate MAPD (Median Absolute Pairwise Difference) for Log-Ratio Data
#'
#' @description
#' Computes the MAPD (Median Absolute Pairwise Difference) metric as a measure of
#' noise and variability in copy number log-ratio data. MAPD quantifies the median
#' absolute difference between adjacent genomic bins within each chromosome, providing
#' an assessment of signal smoothness and data quality.
#'
#' Lower MAPD values indicate smoother, higher quality copy number profiles, while
#' higher values suggest noisy data with excessive variability between consecutive bins.
#'
#' **Note:** This function is exported but requires external data in a specific format.
#' It does **not** accept \code{cnvstrmliner} objects. Instead, it expects a data frame
#' with pre-processed log-ratio data (e.g., from CNVkit .cnr files or similar formats).
#' The function uses \code{data.table} syntax for efficient computation by chromosome.
#'
#' @param cnr_data A data frame (or data.table) containing log-ratio information with
#'   the following required columns:
#'   \describe{
#'     \item{chromosome}{Character or factor indicating the chromosome for each bin}
#'     \item{log2}{Numeric vector of log2 copy number ratios for each genomic bin.
#'       Bins with log2 <= -10 are automatically filtered out as they typically
#'       represent unmappable or poorly covered regions}
#'   }
#' @param use_log Logical. If \code{TRUE} (default), the function will compute MAPD using the log2-transformed values.
#' If \code{FALSE}, it will compute MAPD using the raw depth.
#' This option is provided for flexibility, but using log2 values is generally recommended for assessing copy number data quality.
#' If \code{use_log = TRUE}, the function will filter out bins with log2 values <= -10 to exclude unmappable or poorly covered regions.
#' @return A numeric value representing the MAPD score. Typical interpretation:
#'   \itemize{
#'     \item MAPD < 0.2: Excellent quality, smooth signal
#'     \item MAPD 0.2-0.3: Good quality
#'     \item MAPD 0.3-0.5: Moderate quality, some noise
#'     \item MAPD > 0.5: Poor quality, high noise levels
#'   }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Example with CNVkit .cnr data
#' library(data.table)
#' cnr_data <- fread("sample.cnr")
#' mapd_score <- computeMAPD(cnr_data)
#' print(paste("MAPD score:", round(mapd_score, 4)))
#'
#' # Quality assessment
#' if (mapd_score < 0.3) {
#'   message("High quality copy number data")
#' } else {
#'   message("Noisy data - consider additional filtering or QC")
#' }
#' }
#'
#' @seealso
#' \code{\link{computeGini}} for calculating the Gini coefficient
#'
computeMAPD <- function(cnr_data, use_log = TRUE) {
  # MAPD (Median Absolute Pairwise Difference)
  # filter out bins that are not sequenced
  if (use_log) {
    cnr_data <- cnr_data %>% dplyr::filter(log2 > -10) # or apply a more lenient filter
    cnr_data <- cnr_data %>%
      dplyr::group_by(chromosome) %>%
      dplyr::mutate(logR_diff = abs(log2 - dplyr::lag(log2))) %>%
      dplyr::ungroup()
    # Calculate the MAPD (Median of Absolute Pairwise Differences)
    mapd_score <- median(cnr_data$logR_diff, na.rm = TRUE)
  } else {
    cnr_data <- cnr_data %>% dplyr::filter(!is.na(depth) & depth > 0)
    cnr_data <- cnr_data %>%
      dplyr::group_by(chromosome) %>%
      dplyr::mutate(depth_diff = abs(depth - dplyr::lag(depth))) %>%
      dplyr::ungroup()
    # Calculate the MAPD (Median of Absolute Pairwise Differences)
    mapd_score <- median(cnr_data$depth_diff, na.rm = TRUE)
  }
  # Calculate the difference between adjacent LogR values within each chromosome


  return(mapd_score)
}


#' Compute Coverage per Target from BAM and BED
#'
#' @description
#' Computes per-target coverage statistics from an aligned sequencing file (BAM)
#' over a set of target regions (BED/GRanges). For each target, it derives the
#' median depth and a genome-normalized log2 ratio, returning a GRanges object
#' with these metrics stored in the metadata columns.
#'
#' **Note:** This function is exported but requires external data inputs. It does
#' not accept or operate on \code{cnvstrmliner} objects. Provide a BAM file and
#' target regions (e.g., CNVkit or panel BED) prepared outside of cnvstrmliner.
#'
#' @param bamfile Character path to an indexed BAM file containing aligned reads.
#' @param bedfile Either a \code{GRanges} of target regions or a character path
#'   to a BED file describing targets. In typical workflows, the BED should be
#'   parsed into a \code{GRanges} prior to calling this function.
#'
#' @return A \code{GRanges} object for the provided target regions with the
#'   following metadata columns added:
#'   \describe{
#'     \item{depth}{Median read depth per target}
#'     \item{log2}{Genome-normalized log2 ratio: \eqn{\log_2(depth / genome\_coverage)}}
#'   }
#'
#' @section Implementation details:
#' Autosomal targets are used to estimate the genome-wide coverage baseline.
#' Reads flagged as secondary, QC failed, duplicates, or supplementary are
#' excluded when computing coverage.
#'
#' @examples
#' \dontrun{
#' # Load target regions from a BED file
#' bed <- rtracklayer::import.bed("targets.bed")
#'
#' # Compute coverage metrics from BAM over targets
#' cov_gr <- computeCoverage(bamfile = "tumor.bam", bedfile = bed)
#'
#' # Inspect added metrics
#' head(S4Vectors::mcols(cov_gr))
#' }
#'
#' @seealso
#' \code{\link{computeMAPD}} for assessing log-ratio noise
#' \code{\link{computeGini}} for coverage uniformity
#'
#' @export
computeCoverage <- function(bamfile, bedfile) {
  coverage_matrix <- bamsignals::bamCoverage(bamfile, bed, filteredFlag = 3840) # exclude secondary, QC fail, duplicates, supplementary
  autosomic <- which(!seqnames(bed) %in% c("chrX", "chrY", "chrMT", "chrM", "X", "Y", "MT", "M"))
  median_coverage <- lapply(coverage_matrix, function(x) {
    list(depth = median(x), sum = sum(x), n = length(x))
  }) %>% bind_rows()
  genome_coverage <- sum(median_coverage[autosomic, "sum"]) / sum(median_coverage[autosomic, "n"])
  median_coverage["log2"] <- log2(median_coverage$depth / genome_coverage)
  df <- S4Vectors::DataFrame(median_coverage[c("depth", "log2")])
  S4Vectors::mcols(bed) <- df
  return(bed)
}

#' Compute genomic sizes for each CNV state
#'
#' @description
#' Helper function to compute the total genomic size (in base pairs) for each
#' CNV state (NEUTRAL, AMP, HETDEL, HOMDEL) in a cnvstrml object.
#'
#' @param cnv_obj A cnvstrml object with \code{@segments} (GRanges) containing
#'   \code{total.allele.copies} column.
#'
#' @return A named numeric vector with genomic sizes (in bp) for each state:
#'   \code{NEUTRAL}, \code{AMP}, \code{HETDEL}, \code{HOMDEL}.
#'
#' @keywords internal
#' @export
#'
computeStateSizes <- function(cnv_obj) {
  if (!inherits(cnv_obj, "cnvstrml")) {
    stop("cnv_obj must be a cnvstrml object")
  }

  segments <- cnv_obj@segments

  # Define states based on copy number
  neutral_size <- sum(width(segments[segments$total.allele.copies == 2]))
  amp_size <- sum(width(segments[segments$total.allele.copies > 2]))
  hetdel_size <- sum(width(segments[segments$total.allele.copies == 1]))
  homdel_size <- sum(width(segments[segments$total.allele.copies == 0]))

  return(c(
    NEUTRAL = neutral_size,
    AMP = amp_size,
    HETDEL = hetdel_size,
    HOMDEL = homdel_size
  ))
}


#' Compute size-based weights for CNV state combination
#'
#' @description
#' Helper function to compute weights for combining state-specific metrics.
#' Weight for each state is the average of genomic sizes from both samples.
#'
#' @param truth_sizes Named numeric vector of genomic sizes for truth sample
#'   (output of \code{computeStateSizes}).
#' @param method_sizes Named numeric vector of genomic sizes for method sample
#'   (output of \code{computeStateSizes}).
#'
#' @return A named numeric vector with weights for each state.
#'
#' @keywords internal
#' @export
computeStateWeights <- function(truth_sizes, method_sizes) {
  states <- union(names(truth_sizes), names(method_sizes))

  weights <- sapply(states, function(state) {
    truth_val <- ifelse(state %in% names(truth_sizes), truth_sizes[state], 0)
    method_val <- ifelse(state %in% names(method_sizes), method_sizes[state], 0)
    return((truth_val + method_val) / 2)
  })

  return(weights)
}


# ===== 2. BOUNDARY AND PARTITION =============================================
#
# Where the state metrics of section 1 ask "did it get the copy number right",
# these ask "did it cut the genome in the right places". The distinction is
# worth keeping sharp: a caller can place every boundary correctly and assign
# every state wrongly, or the reverse, and a metric that conflates the two
# cannot say which happened. In a sensitivity sweep that is the whole
# diagnostic question -- over-segmentation and state error respond to opposite
# adjustments.
#
# All five functions here are state-blind, with the single exception of the
# partition_by = "state" mode of ARI and VI, which exists so the same
# pair-counting arithmetic can be pointed at either question.
# =============================================================================


#' Segment boundaries as a GRanges of positions
#'
#' @description
#' Every start and end coordinate in a segmentation, de-duplicated, as width-1
#' ranges. Adjacent segments share a boundary and contribute it once.
#'
#' @param segments A `cnvstrml` object or a `GRanges` of segments.
#' @param exclude_chr_ends Logical. Drop the first and last position on each
#'   sequence. Default `TRUE`: those positions record where the segmentation
#'   happens to begin and end -- usually the first and last covered target --
#'   rather than a detected transition, and counting them inflates agreement
#'   between any two segmentations of the same assay.
#'
#' @return A `GRanges` of width-1 ranges, sorted by sequence and position. An
#'   empty `GRanges` if there are no boundaries left after exclusion.
#'
#' @seealso [computeBoundaryConcordance()], [computeBreakpointConsistency()]
#'
#' @examples
#' \dontrun{
#' bp <- getSegmentBreakpoints(truth)
#' length(bp)
#' }
#'
#' @export
getSegmentBreakpoints <- function(segments, exclude_chr_ends = TRUE) {
  seg_gr <- .asSegments(segments, "segments")
  bp <- .breakpointTable(seg_gr, exclude_chr_ends = exclude_chr_ends)

  if (nrow(bp) == 0) {
    return(GenomicRanges::GRanges())
  }

  GenomicRanges::GRanges(
    seqnames = bp$seqnames,
    ranges = IRanges::IRanges(start = bp$pos, width = 1),
    seqinfo = GenomeInfoDb::seqinfo(seg_gr)
  )
}
#' Compute an integrated weighted BIC score for CNV segmentations
#'
#' @description
#' Computes a Bayesian Information Criterion (BIC)-style score that combines
#' evidence from segment-level fit of log-ratio values and allele-count
#' consistency. This metric is designed to compare alternative copy-number
#' segmentations on the same sample, where lower values indicate a preferred
#' balance between goodness-of-fit and model complexity.
#'
#' @param bins A \code{GRanges} object containing bin-level observations
#'   (e.g., CNVkit \code{.cnr} bins). Must include a numeric column specified by
#'   \code{logR_col}. If \code{weight_col} is present, it is used as per-bin
#'   reliability weight.
#' @param segments A \code{GRanges} object containing segmentation intervals
#'   (e.g., CNVkit \code{.cns} segments). For the B-allele-frequency term, this
#'   object is expected to contain metadata columns
#'   \code{major.allele.counts} and \code{minor.allele.counts}.
#' @param logR_col Character scalar. Name of the numeric metadata column in
#'   \code{bins} containing log-ratio values. Default is \code{"logR"}.
#' @param weight_col Character scalar. Name of the numeric metadata column in
#'   \code{bins} containing bin weights. If missing, all bins are assigned weight 1.
#'   Default is \code{"weight"}.
#' @param noisy Logical. If \code{TRUE}, applies a stronger complexity penalty,
#'   useful for noisy profiles that can otherwise over-segment. If \code{FALSE}
#'   (default), uses the standard piece-wise model penalty.
#' @param strict_bic Logical. If \code{TRUE}, computes a stricter BIC using the
#'   full Gaussian log-likelihood for the log-ratio term and a canonical parameter
#'   count in the penalty. If \code{FALSE} (default), keeps the current
#'   cnvstrmliner-style BIC approximation and optional \code{noisy} penalty.
#' @param winsorize Logical. If \code{TRUE} (default), bin-level log-ratios are
#'   capped at the 5th and 95th percentiles before the fit is computed.
#' @param mergeSegments Numeric or \code{NULL}. If numeric, adjacent segments closer
#'   than this distance (in base pairs) are merged with \code{\link{mergeSegments}}
#'   before scoring, to avoid over-counting near-duplicate breakpoints. Default is
#'   \code{NULL} (no merging).
#' @param baf Logical. If \code{TRUE} (default), includes the binomial
#'   B-allele-frequency log-likelihood term based on the allele-copy columns of
#'   \code{segments}. If \code{FALSE}, that term is set to zero.
#' @param size_weighted_penalty Logical. If \code{TRUE}, each segment contributes to
#'   the complexity penalty in proportion to the number of bins it contains (capped
#'   at 1), so short segments are penalised fractionally. Default is \code{FALSE}.
#' @param size_weight_L0_bins Numeric. Number of bins at which a segment reaches full
#'   weight when \code{size_weighted_penalty = TRUE}. Default is 15.
#' @param chromosome Character or \code{NA}. If supplied, restricts bins and segments
#'   to this chromosome before scoring. Default is \code{NA} (whole genome).
#'
#' @return A named list with the BIC-like score and its components: \code{BIC},
#'   \code{AIC}, \code{AICc}, the penalty terms \code{k} and \code{k_effective},
#'   the counts \code{n} (bins) and \code{m} (segments), \code{RSS}, the
#'   log-likelihood contributions \code{logLik_LR} and \code{logLik_BAF}, and the
#'   settings used (\code{noisy}, \code{strict_bic}, \code{winsorize},
#'   \code{size_weighted_penalty}, \code{size_weight_L0_bins}). Lower \code{BIC}
#'   values indicate better support for the segmentation under the model assumptions.
#'
#' @details
#' The score is built from two likelihood components:
#' \enumerate{
#'   \item \strong{Weighted logR fit:} Bin-level log-ratios are winsorized at the
#'         5th and 95th percentiles to reduce outlier influence, then compared to
#'         segment-level weighted means to compute a weighted residual sum of squares.
#'   \item \strong{Allele-count consistency:} A binomial log-likelihood term is
#'         computed from \code{major.allele.counts} and \code{minor.allele.counts}
#'         in \code{segments}.
#' }
#'
#' The final criterion follows the usual BIC form:
#' \deqn{\mathrm{BIC} = -2 \times \log L + k \times \log(n)}
#' where \eqn{n} is the number of bins and \eqn{k} is a penalty term derived from
#' the number of segments.
#'
#' @examples
#' \dontrun{
#' # bins_gr: GRanges with columns "logR" and optionally "weight"
#' # seg_gr:  GRanges with segment calls and allele count columns
#' bic <- computeBIC(
#'   bins = bins_gr,
#'   segments = seg_gr,
#'   logR_col = "logR",
#'   weight_col = "weight",
#'   noisy = FALSE
#' )
#'
#' # Compare with stronger penalty in noisy samples
#' bic_noisy <- computeBIC(bins_gr, seg_gr, noisy = TRUE)
#' }
#'
#' @seealso
#' \code{\link{computeMAPD}} for profile noise assessment
#' \code{\link{computeGini}} for coverage uniformity
#'
#' @export
computeBIC <- function(bins, segments, 
                      logR_col = "log2", 
                      weight_col = "weight", 
                      noisy = FALSE,
                      strict_bic = FALSE,
                      winsorize = TRUE,
                      mergeSegments=NULL,
                      baf=TRUE,
                      size_weighted_penalty = FALSE,
                      size_weight_L0_bins = 15,
                      chromosome=NA) {
  
  # 1. Setup Weights (Requirement: set to 1 if column is missing)
  if (!is.null(weight_col) && weight_col %in% colnames(S4Vectors::mcols(bins))) {
    bin_weights <- S4Vectors::mcols(bins)[[weight_col]]
  } else {
    bin_weights <- rep(1, length(bins))
  }
  
  # Run the BIC only on the selected chromosome 
  if(!is.na(chromosome)){
    bins <- bins[GenomicRanges::seqnames(bins) == chromosome]
    logger::log_info("Selecting {length(bins)} bins overlapping with chromosome {chromosome} for BIC calculation")
    segments <- segments[GenomicRanges::seqnames(segments) == chromosome]
    logger::log_info("Selecting {length(segments)} segments overlapping with chromosome {chromosome} for BIC calculation")
  }else{
    logger::log_info("Using all {length(bins)} bins and {length(segments)} segments for BIC calculation")
  }
  # remove the centromeres and telomeres if present
  segments <- segments[! S4Vectors::mcols(segments)$arm %in% c("Centromere","Telomere")]
  logger::log_info("After removing centromeres and telomeres, {length(segments)} segments remain for BIC calculation")
  # Merge the segments to avoid overcounting segments that are very close to each other 
  if(!is.null(mergeSegments)){
    segments <-cnvstrmliner::mergeSegments(segments, minDistance=mergeSegments)
  }
  logger::log_info("After merging close segments, {length(segments)} segments remain for BIC calculation")
  logR_vec <- S4Vectors::mcols(bins)[[logR_col]]
  # 2. Winsorize logR Data (Requirement a)
  # Mitigates outlier impact by capping at 5th/95th percentiles (S_R12, S_R30)
  
  if (winsorize) {
    lims <- quantile(logR_vec, probs = c(0.05, 0.95), na.rm = TRUE)
    logR_vec[logR_vec < lims[1]] <- lims[1]
    logR_vec[logR_vec > lims[2]] <- lims[2]
  }
  
  # 3. Map Bins to Segments
  hits <- GenomicRanges::findOverlaps(bins, segments)
  n_obs <- length(unique(S4Vectors::queryHits(hits)))
  m_segs <- length(segments)

  if (n_obs == 0) {
    stop("No overlapping bins between 'bins' and 'segments'.")
  }
  
  # 4. Calculate Weighted Log-Likelihood for logR
  # WRSS is the sum of weight * squared residuals (S_S20, S_R36)
  bin_data <- data.frame(
    val = logR_vec[S4Vectors::queryHits(hits)],
    wt = bin_weights[S4Vectors::queryHits(hits)],
    seg = S4Vectors::subjectHits(hits)
  )
  
  # Calculate weighted means per segment
  seg_stats <- aggregate(cbind(val * wt, wt) ~ seg, data = bin_data, FUN = sum)
  seg_stats$weighted_mean <- seg_stats$V1 / seg_stats$wt
  
  # Map segment weighted means back to each bin to calculate WRSS
  predicted_logR <- seg_stats$weighted_mean[match(bin_data$seg, seg_stats$seg)]
  wrss <- sum(bin_data$wt * (bin_data$val - predicted_logR)^2, na.rm = TRUE)

  # Compute effective segment complexity using segment size in bins.
  # Short segments contribute fractionally to the penalty when enabled.
  if (size_weighted_penalty) {
    seg_bin_counts <- as.numeric(table(factor(bin_data$seg, levels = seq_len(m_segs))))
    seg_weights <- pmin(1, seg_bin_counts / size_weight_L0_bins)
    k_complexity <- sum(seg_weights)
  } else {
    k_complexity <- m_segs
  }
  k_complexity <- max(1, k_complexity)
  
  # Gaussian log-likelihood approximation
  wrss <- max(wrss, .Machine$double.eps)
  if (strict_bic) {
    ll_logR <- -0.5 * n_obs * (log(2 * pi) + 1 + log(wrss / n_obs))
  } else {
    ll_logR <- -0.5 * n_obs * log(wrss / n_obs)
  }
  
  # 5. Calculate Log-Likelihood for BAF (Requirement c)
  # Uses 'major.allele.counts' (A) and 'minor.allele.counts' (B) (S_R20, S_R111)
  if( baf){
    major <- S4Vectors::mcols(segments)$major.allele.copies
    minor <- S4Vectors::mcols(segments)$minor.allele.copies
    total_counts <- S4Vectors::mcols(segments)$total.allele.copies
    
    # Binomial LL: B ~ Binom(Total, prob = minor/total)
    # This favors segmentations with internally consistent allele frequencies.
    valid_baf <- total_counts > 0
    obs_p <- minor[valid_baf] / total_counts[valid_baf]
    
    ll_baf <- sum(dbinom(minor[valid_baf], 
                        size = total_counts[valid_baf], 
                        prob = obs_p, 
                        log = TRUE))
  } else {
    ll_baf <- 0
  }
  
  # 6. Define Penalty Factor (Requirement b)
  
  if (strict_bic) {
    # k_penalty <- 2 * m_segs
    k_penalty <- k_complexity + 1
  } else if (noisy) {
    # Requested formula: k= m+3(m-1)
    k_penalty <- k_complexity + 3 * (k_complexity - 1) # Stronger penalty for noisy profiles (S_R124, S_R156)
  } else {
    # Standard BIC penalty for piece-wise constant models (S_R124, S_R156)
    k_penalty <- 2 * k_complexity - 1
  }
  
  # 7. Final BIC Calculation
  # BIC = -2 * (LL_total) + k * ln(n)
  
  total_ll <- ll_logR/n_obs + ll_baf/m_segs
  
  bic_value <- -2 * total_ll + k_penalty * log(n_obs)
  aic_value <- -2 * total_ll + 2 * k_penalty
  aicc_value <- aic_value + (2 * k_penalty * (k_penalty + 1)) / (n_obs - k_penalty - 1)
  
  return(list("BIC"=bic_value,
              "AIC"=aic_value,
              "AICc"=aicc_value,
              "k"=k_penalty,
              "k_effective"=k_complexity,
              "n"=n_obs,
              "m"=m_segs,
              "RSS"=wrss,
              "logLik_LR"=ll_logR/n_obs,
              "logLik_BAF"=ll_baf/m_segs,
              "noisy"=noisy,
              "strict_bic"=strict_bic,
              "winsorize"=winsorize,
              "size_weighted_penalty"=size_weighted_penalty,
              "size_weight_L0_bins"=size_weight_L0_bins))
}




#' Compute Modified BIC with Segment-Aware Penalty (mBIC_seg)
#'
#' @description
#' A variant of BIC that adjusts the complexity penalty based on segment properties.
#' Reduces penalty for segments with strong signals and long lengths, increases
#' penalty for segments in noisy regions.
#'
#' @param bins A GRanges object containing bin-level observations
#' @param segments A GRanges object containing segmentation intervals
#' @param logR_col Character. Name of logR column in bins
#' @param weight_col Character. Name of weight column in bins
#' @param signal_threshold Numeric. Threshold for "strong signal" (default 0.2)
#' @param min_segment_width Numeric. Minimum segment width in bp for penalty reduction
#' @param chromosome Character. Chromosome to analyze (NA = all)
#' @param winsorize Logical. Whether to winsorize logR values
#' @param baf Logical. Whether to include BAF component
#'
#' @return List with mBIC_seg value and diagnostic information
#'
#' @export
computeMBICseg <- function(bins, segments,
                            logR_col = "log2",
                            weight_col = "weight",
                            signal_threshold = 0.2,
                            min_segment_width = 1e6,
                            chromosome = NA,
                            winsorize = TRUE,
                            baf = TRUE) {
  # First compute standard BIC components
  bic_result <- computeBIC(bins, segments, logR_col, weight_col,
    noisy = FALSE, strict_bic = FALSE,
    winsorize = winsorize, baf = baf,
    chromosome = chromosome
  )

  # Calculate segment-specific penalties
  if (!is.na(chromosome)) {
    segments <- segments[GenomicRanges::seqnames(segments) == chromosome]
  }

  n_obs <- bic_result$n
  m_segs <- length(segments)

  if (m_segs == 0) {
    return(list(mBIC_seg = Inf, BIC = Inf, penalty_adjustment = 0))
  }

  # Map bins to segments to calculate segment properties
  hits <- GenomicRanges::findOverlaps(bins, segments)
  logR_vec <- S4Vectors::mcols(bins)[[logR_col]]

  # Work directly with vectors (option #3 - avoid intermediate data.frame)
  query_idx <- S4Vectors::queryHits(hits)
  subject_idx <- S4Vectors::subjectHits(hits)
  bin_values <- logR_vec[query_idx]

  # Calculate per-segment statistics using data.table (memory-efficient)
  # Column names must match downstream code: val.mean, val.sd, val.n
  bin_dt <- data.table::data.table(
    val = bin_values,
    seg = subject_idx
  )

  seg_stats <- bin_dt[, .(
    val.mean = mean(val, na.rm = TRUE),
    val.sd = sd(val, na.rm = TRUE),
    val.n = .N
  ), by = seg]

  # Convert to data.frame for compatibility with downstream code
  seg_stats <- as.data.frame(seg_stats)

  # Clean up intermediate objects
  rm(bin_dt, query_idx, subject_idx, bin_values)

  # Calculate global noise level
  global_noise <- sd(logR_vec, na.rm = TRUE)

  # Calculate segment widths
  seg_widths <- GenomicRanges::width(segments)
  median_width <- median(seg_widths)

  # Compute penalty adjustments for each segment
  penalty_adjustments <- sapply(seq_len(m_segs), function(i) {
    if (i %in% seg_stats$seg) {
      idx <- which(seg_stats$seg == i)
      mean_logR <- abs(seg_stats$val.mean[idx])
      local_noise <- seg_stats$val.sd[idx]
      seg_width <- seg_widths[i]

      # Base penalty
      base_penalty <- log(n_obs)

      # Signal factor: reduce penalty for strong signals
      signal_factor <- ifelse(mean_logR > signal_threshold,
        exp(-mean_logR / signal_threshold),
        1.0
      )

      # Length factor: reduce penalty for long segments
      length_factor <- ifelse(seg_width > min_segment_width,
        exp(-seg_width / median_width),
        1.0
      )

      # Noise factor: increase penalty in noisy regions
      noise_factor <- ifelse(!is.na(local_noise) && local_noise > 0,
        local_noise / global_noise,
        1.0
      )

      return(base_penalty * signal_factor * length_factor * noise_factor)
    } else {
      return(log(n_obs))
    }
  })

  # Calculate modified BIC
  total_penalty <- sum(penalty_adjustments)
  total_ll <- bic_result$logLik_LR + bic_result$logLik_BAF

  mBIC_seg_value <- -2 * total_ll + total_penalty

  return(list(
    mBIC_seg = mBIC_seg_value,
    BIC = bic_result$BIC,
    penalty_adjustment = total_penalty - bic_result$k * log(n_obs),
    segment_penalties = penalty_adjustments,
    n = n_obs,
    m = m_segs
  ))
}


#' Compute Cross-Validation Information Criterion (CVIC)
#'
#' @description
#' Uses k-fold cross-validation to assess segmentation quality. Trains on k-1 folds
#' and predicts on held-out fold, measuring prediction error.
#'
#' @param bins A GRanges object containing bin-level observations
#' @param segments A GRanges object containing segmentation intervals
#' @param logR_col Character. Name of logR column in bins
#' @param weight_col Character. Name of weight column in bins
#' @param n_folds Integer. Number of cross-validation folds (default 5)
#' @param chromosome Character. Chromosome to analyze (NA = all)
#' @param seed Integer. Random seed for reproducibility
#'
#' @return List with CVIC value and cross-validation diagnostics
#'
#' @export
computeCVIC <- function(bins, segments,
                        logR_col = "log2",
                        weight_col = "weight",
                        n_folds = 5,
                        chromosome = NA,
                        seed = 42) {
  set.seed(seed)

  # Filter to chromosome if specified
  if (!is.na(chromosome)) {
    bins <- bins[GenomicRanges::seqnames(bins) == chromosome]
    segments <- segments[GenomicRanges::seqnames(segments) == chromosome]
  }

  # Map bins to segments
  hits <- GenomicRanges::findOverlaps(bins, segments)
  if (length(hits) == 0) {
    return(list(CVIC = Inf, cv_error = Inf, complexity_penalty = Inf))
  }

  # Work with vectors directly (memory-efficient)
  query_idx <- S4Vectors::queryHits(hits)
  subject_idx <- S4Vectors::subjectHits(hits)

  logR_vec <- S4Vectors::mcols(bins)[[logR_col]]
  logR_values <- logR_vec[query_idx]

  # Setup weights
  if (!is.null(weight_col) && weight_col %in% colnames(S4Vectors::mcols(bins))) {
    bin_weights <- S4Vectors::mcols(bins)[[weight_col]]
    weight_values <- bin_weights[query_idx]
  } else {
    weight_values <- rep(1, length(query_idx))
  }

  # Create bin-segment mapping using data.table (memory-efficient)
  bin_dt <- data.table::data.table(
    bin_idx = query_idx,
    seg_idx = subject_idx,
    logR = logR_values,
    weight = weight_values
  )

  # Clean up vectors
  rm(query_idx, subject_idx, logR_values, weight_values)

  # Assign bins to folds randomly
  unique_bins <- unique(bin_dt$bin_idx)
  fold_assignment <- sample(rep(1:n_folds, length.out = length(unique_bins)))
  bin_dt[, fold := fold_assignment[match(bin_idx, unique_bins)]]

  # Perform k-fold cross-validation using data.table
  cv_errors <- numeric(n_folds)

  for (fold in 1:n_folds) {
    # Use data.table subsetting (efficient, no copies)
    train_idx <- bin_dt$fold != fold
    test_idx <- bin_dt$fold == fold

    if (sum(test_idx) == 0) next

    # Calculate segment means from training data using data.table
    seg_means <- bin_dt[train_idx, .(mean_logR = mean(logR, na.rm = TRUE)), by = seg_idx]

    # Predict on test data (efficient join)
    test_dt <- bin_dt[test_idx]
    test_dt <- seg_means[test_dt, on = "seg_idx"]

    # Calculate weighted prediction error
    cv_errors[fold] <- sum(test_dt$weight * (test_dt$logR - test_dt$mean_logR)^2,
      na.rm = TRUE
    ) / sum(test_dt$weight, na.rm = TRUE)

    # Clean up
    rm(test_dt, seg_means)
  }

  mean_cv_error <- mean(cv_errors[cv_errors > 0], na.rm = TRUE)

  # Complexity penalty based on number of segments
  n_obs <- nrow(bin_dt)
  m_segs <- length(segments)
  complexity_penalty <- m_segs * log(n_obs)

  # CVIC = CV error + complexity penalty
  cvic_value <- mean_cv_error + complexity_penalty

  # Clean up
  rm(bin_dt)

  return(list(
    CVIC = cvic_value,
    cv_error = mean_cv_error,
    complexity_penalty = complexity_penalty,
    cv_errors_by_fold = cv_errors,
    n = n_obs,
    m = m_segs
  ))
}


#' Compute Minimum Segment Quality BIC (MSQBIC)
#'
#' @description
#' Adds quality-based penalties to standard BIC for segments that don't meet
#' quality thresholds (too small, weak signal, high internal variance).
#'
#' @param bins A GRanges object containing bin-level observations
#' @param segments A GRanges object containing segmentation intervals
#' @param logR_col Character. Name of logR column in bins
#' @param weight_col Character. Name of weight column in bins
#' @param min_width_threshold Numeric. Minimum acceptable segment width (bp)
#' @param signal_threshold Numeric. Minimum acceptable signal strength
#' @param small_segment_penalty Numeric. Penalty for segments below width threshold
#' @param weak_signal_penalty Numeric. Penalty for weak signal segments
#' @param homogeneity_bonus Numeric. Bonus for homogeneous segments
#' @param chromosome Character. Chromosome to analyze (NA = all)
#' @param winsorize Logical. Whether to winsorize logR values
#' @param baf Logical. Whether to include BAF component
#'
#' @return List with MSQBIC value and quality diagnostics
#'
#' @export
computeMSQBIC <- function(bins, segments,
                          logR_col = "log2",
                          weight_col = "weight",
                          min_width_threshold = 1e6,
                          signal_threshold = 0.15,
                          small_segment_penalty = 2,
                          weak_signal_penalty = 1,
                          homogeneity_bonus = 0.5,
                          chromosome = NA,
                          winsorize = TRUE,
                          baf = TRUE) {
  # Compute standard BIC
  bic_result <- computeBIC(bins, segments, logR_col, weight_col,
    noisy = FALSE, strict_bic = FALSE,
    winsorize = winsorize, baf = baf,
    chromosome = chromosome
  )

  # Filter to chromosome if specified
  if (!is.na(chromosome)) {
    segments <- segments[GenomicRanges::seqnames(segments) == chromosome]
  }

  if (length(segments) == 0) {
    return(list(MSQBIC = Inf, BIC = Inf, quality_penalties = 0))
  }

  # Map bins to segments
  hits <- GenomicRanges::findOverlaps(bins, segments)
  logR_vec <- S4Vectors::mcols(bins)[[logR_col]]

  # Work directly with vectors (memory-efficient)
  query_idx <- S4Vectors::queryHits(hits)
  subject_idx <- S4Vectors::subjectHits(hits)
  bin_values <- logR_vec[query_idx]

  # Calculate per-segment statistics using data.table
  bin_dt <- data.table::data.table(
    val = bin_values,
    seg = subject_idx
  )

  seg_stats <- bin_dt[, .(
    val.mean = mean(abs(val), na.rm = TRUE),
    val.var = var(val, na.rm = TRUE),
    val.n = .N
  ), by = seg]

  seg_stats <- as.data.frame(seg_stats)
  rm(bin_dt, query_idx, subject_idx, bin_values)

  # Calculate segment widths
  seg_widths <- GenomicRanges::width(segments)

  # Calculate quality penalties using sapply (memory-efficient)
  segment_quality_penalties <- sapply(seq_along(segments), function(i) {
    penalty <- 0

    # Penalty for small segments
    if (seg_widths[i] < min_width_threshold) {
      penalty <- penalty + small_segment_penalty
    }

    # Penalty for weak signal segments
    if (i %in% seg_stats$seg) {
      idx <- which(seg_stats$seg == i)
      mean_signal <- seg_stats$val.mean[idx]

      if (mean_signal < signal_threshold) {
        penalty <- penalty + weak_signal_penalty
      }

      # Bonus for homogeneous segments (low internal variance)
      seg_var <- seg_stats$val.var[idx]
      if (!is.na(seg_var) && seg_var < 0.01) { # Very low variance
        penalty <- penalty - homogeneity_bonus
      }
    }

    return(penalty)
  })

  total_penalties <- sum(segment_quality_penalties)

  # MSQBIC = BIC + quality penalties
  msqbic_value <- bic_result$BIC + total_penalties * log(bic_result$n)

  return(list(
    MSQBIC = msqbic_value,
    BIC = bic_result$BIC,
    total_quality_penalty = total_penalties * log(bic_result$n),
    segment_quality_penalties = segment_quality_penalties,
    n = bic_result$n,
    m = bic_result$m
  ))
}


#' Compute Adaptive BIC with Data-Driven Complexity Penalty
#'
#' @description
#' Adjusts the BIC complexity penalty based on local data characteristics.
#' Uses lower penalty in noisy regions (allowing more segments to explain variation)
#' and standard penalty in clean regions.
#'
#' @param bins A GRanges object containing bin-level observations
#' @param segments A GRanges object containing segmentation intervals
#' @param logR_col Character. Name of logR column in bins
#' @param weight_col Character. Name of weight column in bins
#' @param window_size Integer. Number of segments to use for local noise estimation
#' @param penalty_low Numeric. Penalty multiplier for noisy regions (default 0.5)
#' @param penalty_standard Numeric. Penalty multiplier for clean regions (default 1.0)
#' @param chromosome Character. Chromosome to analyze (NA = all)
#' @param winsorize Logical. Whether to winsorize logR values
#' @param baf Logical. Whether to include BAF component
#'
#' @return List with adaptive_BIC value and diagnostics
#'
#' @export
computeAdaptiveBIC <- function(bins, segments,
                               logR_col = "log2",
                               weight_col = "weight",
                               window_size = 10,
                               penalty_low = 0.5,
                               penalty_standard = 1.0,
                               chromosome = NA,
                               winsorize = TRUE,
                               baf = TRUE) {
  # Compute standard BIC components
  bic_result <- computeBIC(bins, segments, logR_col, weight_col,
    noisy = FALSE, strict_bic = FALSE,
    winsorize = winsorize, baf = baf,
    chromosome = chromosome
  )

  # Filter to chromosome if specified
  if (!is.na(chromosome)) {
    bins <- bins[GenomicRanges::seqnames(bins) == chromosome]
    segments <- segments[GenomicRanges::seqnames(segments) == chromosome]
  }

  if (length(segments) == 0) {
    return(list(adaptive_BIC = Inf, BIC = Inf, adaptive_penalty = 0))
  }

  # Map bins to segments
  hits <- GenomicRanges::findOverlaps(bins, segments)
  logR_vec <- S4Vectors::mcols(bins)[[logR_col]]

  # Calculate global noise level
  global_noise <- sd(logR_vec, na.rm = TRUE)

  # Work directly with vectors (memory-efficient)
  query_idx <- S4Vectors::queryHits(hits)
  subject_idx <- S4Vectors::subjectHits(hits)
  bin_values <- logR_vec[query_idx]

  # Calculate per-segment statistics using data.table
  bin_dt <- data.table::data.table(
    val = bin_values,
    seg = subject_idx
  )

  seg_stats <- bin_dt[, .(
    val.sd = sd(val, na.rm = TRUE),
    val.n = .N
  ), by = seg]

  seg_stats <- as.data.frame(seg_stats)
  rm(bin_dt, query_idx, subject_idx, bin_values)

  # Calculate adaptive penalties for each segment using sapply (memory-efficient)
  n_obs <- bic_result$n

  adaptive_penalties <- sapply(seq_along(segments), function(i) {
    # Estimate local noise from nearby segments
    start_idx <- max(1, i - window_size %/% 2)
    end_idx <- min(length(segments), i + window_size %/% 2)
    local_segs <- start_idx:end_idx

    local_noise <- mean(seg_stats$val.sd[seg_stats$seg %in% local_segs], na.rm = TRUE)

    # Determine penalty multiplier based on local vs global noise
    if (!is.na(local_noise) && local_noise > global_noise) {
      # Noisy region: use lower penalty
      penalty_mult <- penalty_low
    } else {
      # Clean region: use standard penalty
      penalty_mult <- penalty_standard
    }

    return(penalty_mult * log(n_obs))
  })

  # Calculate adaptive BIC
  total_ll <- bic_result$logLik_LR + bic_result$logLik_BAF
  total_adaptive_penalty <- sum(adaptive_penalties)

  adaptive_bic_value <- -2 * total_ll + total_adaptive_penalty

  return(list(
    adaptive_BIC = adaptive_bic_value,
    BIC = bic_result$BIC,
    adaptive_penalty = total_adaptive_penalty,
    segment_penalties = adaptive_penalties,
    penalty_reduction = bic_result$k * log(n_obs) - total_adaptive_penalty,
    n = n_obs,
    m = length(segments)
  ))
}
