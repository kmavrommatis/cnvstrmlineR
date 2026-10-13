#' Create consensus CNV calls via ensemble voting across multiple callers
#'
#' @description
#' Generates a consensus \code{cnvstrml} object by combining CNV calls from multiple input
#' samples (e.g., from different callers or different versions of the same caller).
#' Uses majority voting at the annotation level to retain only segments with consistent
#' CNV status (AMP, DEL, NEUT, etc.) across at least half + 1 samples.
#' Adjusts per-sample copy numbers to a consensus ploidy, computes median purity/ploidy,
#' and aggregates marker density and logR across samples before merging adjacent segments.
#'
#' @param cnvList List of \code{cnvstrml} objects to combine. All must have the same genome build.
#' @param sampleName Character sample name for the output. If \code{NULL} and all inputs share
#'   a name, that name is used with a warning. If samples differ, the first sample name is used.
#' @param methodName Character method description for the output. If \code{NULL}, concatenates
#'   input method names with " ensemble" suffix.
#' @param info Logical. Currently unused; reserved for detailed logging.
#' @param fill.gaps Logical. Deprecated and ignored; retained for backwards compatibility.
#' @param min.size Integer minimum segment width in base pairs to retain after consensus
#'   (default: 200). Filters short, potentially spurious segments.
#' @param min.majority Integer minimum number of supporting samples for a consensus segment
#'   (default: 1, i.e., all segments kept). Increase to require stricter consensus.
#' @param .mergeSegments Logical. If TRUE, adjacent compatible segments are merged after consensus is built.
#'
#' @return A \code{cnvstrml} object with consensus segments and metadata:
#' \describe{
#'   \item{@segments}{A \code{GRanges} with qualitative consensus (annotation) and quantitative
#'     fields aggregated via majority vote or median:
#'     \describe{
#'       \item{total.allele.copies}{Median of copy numbers across supporting samples, adjusted to consensus ploidy.}
#'       \item{major.allele.copies}{Calculated from total and minor; derived via median.}
#'       \item{minor.allele.copies}{Median minor allele copies across supporting samples.}
#'       \item{cnv.annotation}{CNV status (AMP, DEL, NEUT) assigned via \code{annotateCNV()}.}
#'       \item{cnv.loh}{Logical; \code{TRUE} if minor allele copies = 0.}
#'       \item{cnv.lrr}{Median log-ratio.}
#'       \item{cnv.cf}{Median purity (cancer cell fraction proxy).}
#'       \item{cnv.ccf}{Clonal cancer cell fraction; mean across supporting samples.}
#'       \item{cnv.totalmarkers}{Total markers, derived from marker density and segment width.}
#'       \item{cnv.markers}{Used markers, derived from marker density.}
#'       \item{support}{Count of samples supporting each consensus segment.}
#'     }}
#'   \item{@sample}{Output sample name.}
#'   \item{@method}{Method description (concatenated input methods + " ensemble").}
#'   \item{@purity}{Median purity of input samples.}
#'   \item{@ploidy}{Median ploidy of input samples.}
#'   \item{@sex}{Sex of samples; all must agree.}
#'   \item{@genome}{Genome build; all must agree.}
#' }
#'
#' @details
#' \strong{Consensus and voting strategy}
#' \itemize{
#'   \item First pass (qualitative): Finds segments where CNV annotation (NEUT/AMP/DEL/etc.)
#'     matches in > 50% of samples; discards segments marked as INFER from gap-filling.
#'   \item Second pass (quantitative): For retained segments, computes median copy numbers,
#'     logR, and marker densities. Clonal fraction (CCF) is averaged; other numeric fields use median.
#'   \item Ploidy adjustment: All input samples are adjusted to the consensus (rounded median) ploidy
#'     before aggregating copy-number values, ensuring consistent interpretation.
#' }
#'
#' \strong{Filtering and refinement}
#' \itemize{
#'   \item \code{min.size}: Removes segments narrower than this threshold (default 200 bp).
#'   \item \code{min.majority}: Requires at least this many samples to support a segment (default 1).
#'   \item Segment merging: Adjacent segments with the same annotation are merged unless \code{.mergeSegments = FALSE}.
#' }
#'
#' \strong{Input requirements}
#' \itemize{
#'   \item All inputs must be \code{cnvstrml} objects.
#'   \item All must have the same \code{genome} build and \code{sex}.
#'   \item Each sample must have a unique name in the list.
#' }
#'
#' \strong{Use cases in multi-caller benchmarking}
#' \itemize{
#'   \item Identify robust CNV calls supported by multiple methods.
#'   \item Estimate consensus copy numbers and clonality.
#'   \item Reduce caller-specific false positives via voting.
#'   \item Create a "gold-standard" reference for method comparison.
#' }
#'
#' @examples
#' \dontrun{
#' # Ensemble three callers
#' ens <- ensemblCall(
#'   cnvList = list(facets_call, cnvkit_call, cfreec_call),
#'   sampleName = "S1_consensus",
#'   min.size = 200,
#'   min.majority = 2 # Require at least 2 callers
#' )
#'
#' # Inspect purity and ploidy
#' cat("Median purity:", ens@purity, "\nMedian ploidy:", ens@ploidy, "\n")
#'
#' # Check segment support
#' table(S4Vectors::mcols(ens@segments)$support)
#' }
#'
#' @seealso \code{\link{findConsensus}}, \code{\link{compareSegments}}, \code{\link{addCCF}}
#'
#' @export
#'

ensemblCall <- function(cnvList,
                        sampleName = NULL,
                        methodName = NULL,
                        info = FALSE,
                        fill.gaps = TRUE, # ignored
                        min.size = 200,
                        min.majority = 1,
                        .mergeSegments = TRUE) {
  # checks to confirm the proper inputs
  if (!is(cnvList, "list") | sapply(cnvList, function(X) {
    class(X)
  }) %>% unique() != "cnvstrml") {
    stop("ensemblCall can be used only with lists of cnvstrml objects", )
  }
  # check that all the samples have the same genome as baseline.
  genomes <- sapply(cnvList, function(x) {
    x@genome
  }) %>% unique()
  if (length(genomes) != 1) {
    stop("The genome version used in this list ", paste(genomes, collapse = ","), " is not consistent")
  }

  # create the method name for the ensembl sample by concatenating all the methods.
  method <- sapply(cnvList, function(x) {
    x@method
  }) %>% unique()
  if (length(method) > 1) {
    method <- paste(method, sep = ",")
  }

  # get the sex of the ensembl sample
  # all constituent samples should have the same sex
  sex <- sapply(cnvList, function(x) {
    x@sex
  }) %>% unique()
  if (length(sex) != 1) {
    stop("The samples of this list do not have consistent sex", paste(sex, collapse = ","))
  }

  # generate the sample name
  samples <- sapply(cnvList, function(x) {
    x@sample
  }) %>% unique()
  if (is.null(sampleName)) {
    sampleName <- samples[1]
    if (length(samples) != 1) {
      warning("The samples of this list are not consistent ", paste(samples, collapse = ","), ". We will use ", sampleName)
    }
  }
  if (length(samples) == 1) {
    stop("This function requires that each sample has a unique sampleName")
  }
  names(cnvList) <- samples

  # calculate the values for the purity and ploidy for the new object
  purities <- sapply(cnvList, function(x) {
    x@purity
  })
  newPurity <- purities %>% median(., na.rm = TRUE)
  logger::log_info("Calculated purity as {newPurity} from the original ones. Mean:{round(mean(purities,na.rm=TRUE), digits=3)}, Median:{round(median(purities,na.rm=TRUE), digits=3)}, Stand Dev:{round(sd(purities,na.rm=TRUE), digits=3)  }", namespace = "cnvstrmliner")

  ploidies <- sapply(cnvList, function(x) {
    x@ploidy
  })
  newPloidy <- ploidies %>% median(., na.rm = TRUE)
  rounded_ploidy <- round(newPloidy)
  logger::log_info("Calculated ploidy as { newPloidy}, from the original ones. Mean: {round(mean(ploidies,na.rm=TRUE), digits=3)} Median:{round(median(ploidies,na.rm=TRUE) , digits=3)} Stand Dev: {round(sd(ploidies,na.rm=TRUE),digits=3)}", namespace = "cnvstrmliner")

  # find the consensus among the comprising samples (in qualitative mode)
  # i.e. mark NEUT/AMP etc
  # we add gaps and mark them as INFER
  consensus <- suppressMessages({
    cnvstrmliner::findConsensus(cnvList,
      ignore.samples = TRUE,
      qual.col = "cnv.annotation"
    )
  })
  # keep only the segments where the majority.count is greater than half of the samples
  # but we will ignore the ones that are not containing any information

  S4Vectors::mcols(consensus@segments)["cutoff"] <- S4Vectors::mcols(consensus@segments)[, "stats.count"] / 2
  segs <- consensus@segments[consensus@segments$stats.majority.count > consensus@segments$cutoff] %>%
    magrittr::set_names(NULL)
  logger::log_info("{round(length(segs)/length(consensus@segments) *100, digits=2)}% of segments retained after removing segments without majority", namespace = "cnvstrmliner")
  segs <- segs[which(BiocGenerics::width(segs) >= min.size & segs$stats.majority.count >= min.majority)]
  logger::log_info("{round(length(segs)/length(consensus@segments) *100, digits=2)} % of segments retained after removing short segments", namespace = "cnvstrmliner")
  # create a mask that keeps only the samples that were consistent with the majority call for NEUTR/AMP/DEL etc
  mask <- createMask(segs, samples)

  # we need to do some updating to the values based on the final ploidy we want to use

  # iterate through the samples and check how we need to modify their counts
  # for each object we multiple the ratio of total.allele copies with the new ploidy
  # add a 0.01 at the end to be able to round it to the higher value if it is exactly .5
  # eg 4.5 will be rounded to 5
  cnvListAdjusted <- lapply(cnvList, function(C) {
    adjustPloidy(C, rounded_ploidy)
  })


  # consensus.segs is producing a new object that has all the common regions between
  # the samples, along with the values of the columns of interest
  consensus.segs <- suppressMessages({
    cnvstrmliner::findConsensus(cnvListAdjusted,
      ignore.samples = TRUE,
      qual.col = NULL
    )
  })

  ov <- GenomicRanges::findOverlaps(query = segs, subject = consensus.segs@segments)
  consensus.segs@segments <- consensus.segs@segments[S4Vectors::subjectHits(ov)]

  ### -- produce an object with the values of the column we want
  # column_of_interest='total.allele.copies'
  getValue <- function(column_of_interest) {
    consensus.qtotal <- consensus.segs
    cols_to_keep <- grep(column_of_interest, colnames(S4Vectors::mcols(consensus.segs@segments)))
    if (length(cols_to_keep) == 0) {
      majority <- rep(NA, nrow(S4Vectors::mcols(consensus.segs@segments)))
      return(as.data.frame(majority) %>% magrittr::set_colnames(column_of_interest))
    }
    S4Vectors::mcols(consensus.qtotal@segments) <- S4Vectors::mcols(consensus.qtotal@segments)[, cols_to_keep]
    new_colnames <- colnames(S4Vectors::mcols(consensus.qtotal@segments)) %>%
      gsub(pattern = paste0("\\|\\|", column_of_interest), replacement = "", .)
    colnames(S4Vectors::mcols(consensus.qtotal@segments)) <- new_colnames
    # find the major alleles from the consensus object
    consensus.qtotal <- suppressMessages(applyMask(consensus.qtotal, mask))

    if (startsWith(column_of_interest, "cnv")) {
      majority <- S4Vectors::mcols(consensus.qtotal@segments)[, samples] %>%
        as.matrix() %>%
        rowMeans(., na.rm = TRUE)
    } else {
      majority <- S4Vectors::mcols(consensus.qtotal@segments)[, samples] %>%
        as.matrix() %>%
        matrixStats::rowMedians(., na.rm = TRUE)
    }
    return(as.data.frame(majority) %>% magrittr::set_colnames(column_of_interest))
  }

  new_metadata <- future.apply::future_lapply(
    c(
      "total.allele.copies",
      "minor.allele.copies",
      "cnv.lrr",
      "cnv.cf",
      "cnv.ccf",
      "cnv.totalmarkers_density",
      "cnv.markers_density",
      "cnv.cnaqc.pass"
    ),
    function(X) {
      getValue(X)
    }
  ) %>% dplyr::bind_cols()
  S4Vectors::mcols(consensus.segs@segments) <- new_metadata

  consensus.segs@segments$cnv.totalmarkers <- round(consensus.segs@segments$cnv.totalmarkers_density * BiocGenerics::width(consensus.segs@segments))
  consensus.segs@segments$cnv.markers <- round(consensus.segs@segments$cnv.markers_density * BiocGenerics::width(consensus.segs@segments))

  # # add the number of samples that provide support for the inference of the ensembl value

  ov <- GenomicRanges::findOverlaps(query = segs, subject = consensus.segs@segments)
  S4Vectors::mcols(consensus.segs@segments)[S4Vectors::subjectHits(ov), "support"] <- segs$stats.majority.count[S4Vectors::queryHits(ov)]


  mcols_new <- S4Vectors::mcols(consensus.segs@segments) %>% as.data.frame()
  flags <- apply(mcols_new, 1, function(K) {
    cnvstrmliner::annotateCNV(
      minor.cn = K["minor.allele.copies"], total.cn = K["total.allele.copies"],
      ploidy = round(newPloidy)
    )
  })
  loh <- rep(FALSE, nrow(mcols_new))
  loh[which(mcols_new$minor.allele.copies == 0)] <- TRUE


  S4Vectors::mcols(consensus.segs@segments) <- mcols_new %>%
    dplyr::mutate(
      major.allele.copies = ifelse(is.na(minor.allele.copies),
        total.allele.copies,
        total.allele.copies - minor.allele.copies
      ),
      cnv.loh = loh,
      cnv.annotation = flags
    ) %>%
    dplyr::select(
      total.allele.copies, major.allele.copies, minor.allele.copies, cnv.totalmarkers,
      cnv.markers, cnv.cf, cnv.ccf, cnv.loh, cnv.annotation, cnv.lrr, cnv.cnaqc.pass, support
    )
  # merge segments that can be merged
  if (.mergeSegments == TRUE) {
    if (!"arm" %in% colnames(S4Vectors::mcols(consensus.segs@segments))) {
      consensus.segs <- cnvstrmliner::addChromosomeAnnotation(consensus.segs)
    }
    consensus.segs@segments <- tryCatch(
      {
        cnvstrmliner::mergeSegments(consensus.segs@segments, minDistance = min.size)
      },
      error = function(e) {
        stop("Unable to merge segments [", e, "]")
      },
      warning = function(w) {
        warning("While merging segments [", w, "]")
      }
    )
  }

  consensus.segs@sex <- sex
  consensus.segs@sample <- sampleName
  consensus.segs@method <- ifelse(is.null(methodName), paste(method, "ensemble"), methodName)
  consensus.segs@purity <- newPurity
  consensus.segs@ploidy <- newPloidy
  consensus.segs@effective_ploidy <- rounded_ploidy
  consensus.segs <- cnvstrmliner::addCCF(consensus.segs)
  return(consensus.segs)
}


#' apply a mask of TRUE/FALSE on the mcols of the consensus object
#'
#' @param cnvObj (cnvstrmlinr object) a consensus object
#' @param mask (GRanges) an object with TRUE/FALSE in the mcols
#'
#' @description
#' The mask should contain the  columns of the gr object we want to filter
#'
applyMask <- function(cnvObj, mask) {
  if (!is(cnvObj, "cnvstrml")) {
    stop("This function works only with cnvstrml objects. You provided a ", class(cnvObj))
  }
  if (grepl("consensus", cnvObj@segtype) == FALSE) {
    message("This function can work only with consensus containing objects")
    return(NULL)
  }
  mcols_mask <- S4Vectors::mcols(mask)
  mcols_gr <- S4Vectors::mcols(cnvObj@segments)[colnames(mcols_mask)]
  masked_values <- lapply(colnames(mcols_mask), function(S) {
    selection <- which(mcols_mask[, S] == FALSE)
    if (length(selection) == 0) {
      return(mcols_gr[, S])
    }
    mcols_gr[selection, S] <- NA
    mcols_gr[, S] %>% as.vector()
  }) %>%
    dplyr::bind_cols() %>%
    magrittr::set_colnames(colnames(mcols_mask))
  S4Vectors::mcols(cnvObj@segments)[colnames(mcols_mask)] <- masked_values[colnames(mcols_mask)]
  cnvObj
}


#' create a mask
#'
#' @param segs A GRanges object with multiple columns. It is mandatory to have the column 'majority'.
#' @param cols Character vector of column names used to create the mask. This should not include the 'majority' column.
#' @description
#'
#' create a GRanges object where the mcols have values TRUE/FALSE
#' where samples that their individual call matches the 'majority' call
#'
createMask <- function(segs, cols) {
  aa <- segs %>% as.data.frame(row.names = NULL, )
  aa_mcols <- lapply(make.names(cols), function(S) {
    id <- which(aa[S] == aa["stats.majority"])
    aa[S] <- FALSE
    aa[id, S] <- TRUE
    aa[S]
  }) %>% dplyr::bind_cols()
  mask <- segs
  S4Vectors::mcols(mask) <- aa_mcols
  colnames(S4Vectors::mcols(mask)) <- cols
  mask
}


#' adjust copy counts based on the ploidy
#' @param cnvObj A cnvstrml object whose copy numbers will be adjusted.
#' @param newploidy Integer ploidy value to apply to the object.
#' @description
#' Each segment of the object is first converted to a ratio with respect to the ploidy
#' and then multiplied with the new ploidy.
#' For segments where the new (adjusted) ploidy is a decimal we
#' add 0.01 if it is above the ploidy and subtract 0.01 if below the ploidy
#' this way if a segment gets a value e.g. 2.5 and the new ploidy is 2, the segment will be rounded to 3
#' conversely if the new ploidy is 1.5 it will be rounded to 1.
#' The new ploidy is stored in the effectivePloidy slot of the derived object
#'
adjustPloidy <- function(cnvObj, newploidy) {
  newcopies <- S4Vectors::mcols(cnvObj@segments)[, "total.allele.copies"] / cnvObj@effective_ploidy * newploidy
  newcopies[which(newcopies > newploidy)] <- newcopies[which(newcopies > newploidy)] + 0.01
  newcopies[which(newcopies < newploidy)] <- newcopies[which(newcopies < newploidy)] - 0.01
  S4Vectors::mcols(cnvObj@segments)[, "total.allele.copies"] <- round(newcopies)
  newcopies <- S4Vectors::mcols(cnvObj@segments)[, "minor.allele.copies"] / cnvObj@effective_ploidy * newploidy
  newcopies[which(newcopies < 1 & newcopies > 0)] <- 1
  S4Vectors::mcols(cnvObj@segments)[, "minor.allele.copies"] <- round(newcopies)

  S4Vectors::mcols(cnvObj@segments)[, "major.allele.copies"] <- S4Vectors::mcols(cnvObj@segments)[, "total.allele.copies"] -
    S4Vectors::mcols(cnvObj@segments)[, "minor.allele.copies"]
  cnvObj@effective_ploidy <- newploidy
  cnvObj
}


#' Merge tandem segments
#'
#'
#' Merge adjacent segments with consistent copy number and annotation
#'
#' @description
#' Merges tandem CNV segments when they are separated by a small distance and share
#' identical copy number states, annotation, and LOH status across a chromosomal arm.
#' Merged segments are represented by a single interval with aggregated metadata
#' (median copy numbers, weighted-average purity/clonality, summed marker counts).
#'
#' @param segments A \code{GRanges} object representing CNV segments. Must contain
#'   chromosome arm annotation (added via \code{addChromosomeAnnotation()}).
#'   Required metadata columns: \code{total.allele.copies}, \code{major.allele.copies},
#'   \code{minor.allele.copies}, \code{cnv.annotation}, \code{cnv.loh}, \code{cnv.lrr},
#'   \code{cnv.cf}, \code{cnv.ccf}, \code{cnv.totalmarkers}, \code{cnv.markers}.
#'
#' @param minDistance Integer; maximum gap (in base pairs) between two segments
#'   to be considered adjacent and mergeable (default: 1).
#' @param minlrrDiffPct Numeric. Currently unused; reserved for future implementation
#'   of log-ratio similarity thresholds.
#' @param minccfDiffPct Numeric. Currently unused; reserved for future implementation
#'   of clonal fraction similarity thresholds.
#'
#' @return A \code{GRanges} object with merged segments. Metadata columns are aggregated as follows:
#' \describe{
#'   \item{Chromosomal position}{Merged into single interval spanning merged segments.}
#'   \item{arm}{Preserved from constituent segments.}
#'   \item{total.allele.copies, major.allele.copies, minor.allele.copies}{Median across merged segments.}
#'   \item{cnv.annotation, cnv.loh}{Preserved (must be identical for merging to occur).}
#'   \item{cnv.lrr}{Weighted average by segment width.}
#'   \item{cnv.cf, cnv.ccf}{Weighted average by segment width.}
#'   \item{cnv.totalmarkers, cnv.markers}{Sum across merged segments.}
#'   \item{cnv.cnaqc.pass}{Weighted average (if present).}
#'   \item{support}{Minimum count across merged segments (if present).}
#' }
#'
#' @details
#' \strong{Merging criteria}
#' Two adjacent segments merge if:
#' \itemize{
#'   \item Distance between end of first and start of second < \code{minDistance} bp.
#'   \item Same \code{cnv.annotation}.
#'   \item Same \code{total.allele.copies}, \code{major.allele.copies}, \code{minor.allele.copies}.
#'   \item Same \code{cnv.loh} status.
#'   \item Belong to the same chromosome and arm.
#' }
#'
#' \strong{Aggregation strategy}
#' \itemize{
#'   \item \emph{Copy numbers}: Median ensures robust estimates when multiple segments collapse.
#'   \item \emph{Purity/clonality}: Weighted by segment width to respect larger regions.
#'   \item \emph{Markers}: Summed to reflect total evidence across the merged region.
#'   \item \emph{Support}: Minimum preserved, indicating the weakest supporting sample(s).
#' }
#'
#' \strong{Workflow context}
#' \code{mergeSegments()} is typically called after consensus building (e.g., in \code{ensemblCall})
#' to simplify final segment boundaries and reduce over-fragmentation. Arm-level annotation
#' is required; call \code{addChromosomeAnnotation()} first if not present.
#'
#' @examples
#' \dontrun{
#' # Merge segments within 200 bp, using median ploidy of 2
#' merged <- mergeSegments(cnv_object@segments, minDistance = 200)
#'
#' # Compare before and after
#' cat("Original segments:", length(cnv_object@segments), "\n")
#' cat("Merged segments:", length(merged), "\n")
#' }
#'
#' @seealso \code{\link{addChromosomeAnnotation}}, \code{\link{ensemblCall}}, \code{\link{compareSegments}}
#'
#' @export
mergeSegments <- function(segments, minDistance = 1, minlrrDiffPct = 0.05, minccfDiffPct = 0.05) {
  # split each segments by chromosome (to avoid merging segments from different chromosomes)
  # for any segment with TRUE value, we merge it with teh following segment

  starting_segments <- length(segments)
  if (!"arm" %in% colnames(S4Vectors::mcols(segments))) {
    stop("The segments should contain the annotation with arms. Please run the 'addChromosomeAnnotation' function. ")
  }

  segments <- segments %>%
    split(., f = paste(GenomeInfoDb::seqnames(.), S4Vectors::mcols(.)[, "arm"])) %>%
    lapply(., function(x) {
      flags <- abs(GenomicRanges::end(x) - dplyr::lead(GenomicRanges::start(x))) < minDistance &
        x$total.allele.copies == dplyr::lead(x$total.allele.copies) &
        x$minor.allele.copies == dplyr::lead(x$minor.allele.copies) &
        x$cnv.loh == dplyr::lead(x$cnv.loh) &
        x$cnv.annotation == dplyr::lead(x$cnv.annotation)
      flags[is.na(flags)] <- FALSE

      S4Vectors::mcols(x)["merge"] <- flags

      # merge_index <- which(flags == TRUE)
      # GenomicRanges::end(x)[merge_index] <- GenomicRanges::start(x)[merge_index + 1]
      x
    }) %>%
    as(., "GRangesList") %>%
    unlist()

  # Create merge groups after unlisting to ensure unique group IDs across all chromosomes/arms
  # A new group starts when the previous segment has merge=FALSE
  merge_flags <- S4Vectors::mcols(segments)[, "merge"]
  merge_group <- c(1, cumsum(!merge_flags[-length(merge_flags)]) + 1)
  S4Vectors::mcols(segments)["merge_group"] <- merge_group

  # Ensure segments in the same merge_group touch each other
  # If two consecutive segments have the same merge_group, extend the first to touch the second
  for (i in seq_len(length(segments) - 1)) {
    if (S4Vectors::mcols(segments)[i, "merge_group"] == S4Vectors::mcols(segments)[i + 1, "merge_group"]) {
      if (GenomicRanges::end(segments)[i] != GenomicRanges::start(segments)[i + 1]) {
        GenomicRanges::end(segments)[i] <- GenomicRanges::start(segments)[i + 1]
      }
    }
    if (S4Vectors::mcols(segments)[i, "merge_group"] != S4Vectors::mcols(segments)[i + 1, "merge_group"]) {
      if (abs(GenomicRanges::end(segments)[i] - GenomicRanges::start(segments)[i + 1]) < 2) {
        GenomicRanges::end(segments)[i] <- GenomicRanges::start(segments)[i + 1] - 2
      }
    }
  }


  # find which segments we need to merge
  # segmentsToMerge=sapply(2:length(segments),function(idx){
  #   cnvstrmliner:::compareSegments( segA=segments[idx-1, ],
  #                    segB=segments[idx ,],
  #                    minDistance=minDistance,
  #                    minlrrDiffPct=minlrrDiffPct,
  #                    minccfDiffPct=minccfDiffPct)
  # })
  # extend the mergable segments by few bp so we can merge them using reduce and a distance of 0
  # merge_index=which(segmentsToMerge==TRUE)
  # GenomicRanges::end(segments)[merge_index] = GenomicRanges::start(segments)[merge_index+1]
  segmentsReduced <- GenomicRanges::reduce(segments, min.gapwidth = 0, with.revmap = TRUE)
  # build the mcols for the reduced object
  revmap <- segmentsReduced$revmap
  mcols_dt <- as.data.frame(S4Vectors::mcols(segments))
  widths <- BiocGenerics::width(segments) # Precompute widths
  idx <- 0
  mcols_new <- lapply(revmap, function(X) {
    # Subset once
    subset_dt <- mcols_dt[X, , drop = FALSE]
    widths_sub <- widths[X]
    idx <<- idx + 1
    list(
      revmap = idx,
      arm = tryCatch(unique(subset_dt[["arm"]]), error = function(e) NA),
      total.allele.copies = median(subset_dt[["total.allele.copies"]], na.rm = TRUE),
      major.allele.copies = median(subset_dt[["major.allele.copies"]], na.rm = TRUE),
      minor.allele.copies = median(subset_dt[["minor.allele.copies"]], na.rm = TRUE),
      cnv.totalmarkers = sum(subset_dt[["cnv.totalmarkers"]], na.rm = TRUE),
      cnv.markers = sum(subset_dt[["cnv.markers"]], na.rm = TRUE),
      cnv.cf = weighted.mean(subset_dt[["cnv.cf"]], widths_sub, na.rm = TRUE),
      cnv.ccf = weighted.mean(subset_dt[["cnv.ccf"]], widths_sub, na.rm = TRUE),
      cnv.annotation = unique(subset_dt[["cnv.annotation"]]),
      cnv.loh = unique(subset_dt[["cnv.loh"]]),
      cnv.lrr = weighted.mean(subset_dt[["cnv.lrr"]], widths_sub, na.rm = TRUE),
      cnv.cnaqc.pass = tryCatch(weighted.mean(subset_dt[["cnv.cnaqc.pass"]], widths_sub, na.rm = TRUE), error = function(e) NA),
      support = tryCatch(min(subset_dt[["support"]], na.rm = TRUE), error = function(e) {
        NA
      }, warning = function(w) {})
    )
  }) %>% dplyr::bind_rows()

  S4Vectors::mcols(segmentsReduced) <- mcols_new
  final_segments <- length(segmentsReduced)
  logger::log_warn("Merging segments with distance < {minDistance}bp. This will modify the segments of the sample from {starting_segments} to {final_segments}.", namespace = "cnvstrmliner")
  segmentsReduced
}

#' compare two segments and decide if they can be merged
#' @param segA A single GRanges segment.
#' @param segB A single GRanges segment.
#' @param minDistance Integer minimum distance between two segments to consider them for merging.
#' @param minlrrDiffPct Numeric percentage difference between the two lrr values to accept for merging.
#' @param minccfDiffPct Numeric percentage difference between the two ccf values to accept for merging.
#' @return TRUE if the segments can be merged.
compareSegments <- function(segA, segB, minDistance = 1, minlrrDiffPct = 0.05, minccfDiffPct = 0.05) {
  # Precompute values to avoid redundant calls
  seqA <- as.character(GenomeInfoDb::seqnames(segA))
  seqB <- as.character(GenomeInfoDb::seqnames(segB))

  # Check chromosome match early
  if (seqA != seqB) {
    return(FALSE)
  }

  # Extract metadata only once
  mcolsA <- S4Vectors::mcols(segA)
  mcolsB <- S4Vectors::mcols(segB)

  # Check if 'arm' exists and compare
  if (!is.null(mcolsA$arm) && !is.null(mcolsB$arm) && mcolsA$arm != mcolsB$arm) {
    return(FALSE)
  }

  # Distance filter
  if (BiocGenerics::start(segB) - BiocGenerics::end(segA) > minDistance) {
    return(FALSE)
  }

  # Avoid division by zero and ensure meaningful percentage differences
  safeDiv <- function(x, y) ifelse(y == 0, Inf, x / y)

  lrrDiffPct <- safeDiv(abs(mcolsA$cnv.lrr - mcolsB$cnv.lrr), max(abs(mcolsA$cnv.lrr), abs(mcolsB$cnv.lrr), 1e-6))
  if (!is.nan(lrrDiffPct) && lrrDiffPct > minlrrDiffPct) {
    return(FALSE)
  }

  ccfDiffPct <- safeDiv(abs(mcolsA$cnv.cf - mcolsB$cnv.cf), max(abs(mcolsA$cnv.cf), abs(mcolsB$cnv.cf), 1e-6))
  if (!is.nan(ccfDiffPct) && ccfDiffPct > minccfDiffPct) {
    return(FALSE)
  }

  # Categorical column checks
  if (mcolsA$cnv.loh != mcolsB$cnv.loh) {
    return(FALSE)
  }
  if (mcolsA$cnv.annotation != mcolsB$cnv.annotation) {
    return(FALSE)
  }

  # Check total allele copy and major allele copy equality
  return(mcolsA$total.allele.copies == mcolsB$total.allele.copies &&
    mcolsA$major.allele.copies == mcolsB$major.allele.copies)
}
