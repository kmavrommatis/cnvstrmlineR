#' Parse Control-FREEC CNV and allelic imbalance output
#'
#' @description
#' High-level parser for Control-FREEC output files. Control-FREEC is a tool for detecting somatic
#' copy-number changes and allelic imbalances (including LOH) from sequencing data (WGS, WES, targeted panels).
#' It automatically computes, normalizes, and segments copy-number and BAF profiles, then calls CNV alterations and LOH.
#' The matched normal is optional for WGS but mandatory for WES/targeted data.
#' This function aggregates multiple Control-FREEC output files (_info, _CNVs/_CNVs.pvalue, _ratio, _BAF)
#' into a single \code{cnvstrml} object with standardized fields, marker counts, and logR values.
#'
#' @param cfreec.fn Character path to the main Control-FREEC output file (typically \code{_CNVs} or \code{_CNVs.p.value.txt}).
#'   The function infers the file stem and searches for associated files.
#' @param filterLowUncertainty Logical; if \code{TRUE}, segments with uncertainty = -1 (low confidence)
#'   are excluded (default: \code{FALSE}, keep all).
#' @param pvalThreshold Numeric p-value cutoff (0–1) for filtering segments if p-values are present
#'   in the CNVs file (default: 1, no filtering).
#' @param genome Character genome build (e.g., \code{"hg38"}) used to set \code{Seqinfo} and UCSC seqlevels.
#' @param sample Character sample name to store in the object. If \code{NULL}, the file stem is used.
#'
#' @return A \code{cnvstrml} object with segments and metadata:
#' \describe{
#'   \item{@segments}{A \code{GRanges} with metadata:
#'     \describe{
#'       \item{seqnames, start, end}{Segment coordinates from _CNVs file.}
#'       \item{copy.number}{Predicted copy number (integer or allele-specific if genotype available).}
#'       \item{status}{Type of alteration (gain, loss, or neutral) from _CNVs.}
#'       \item{genotype}{Inferred or predicted genotype (if BAF data available); \code{NA} if not determined.}
#'       \item{uncertainty}{Uncertainty metric from Control-FREEC (confidence measure).}
#'       \item{cnv.lrr}{Median log-ratio from _ratio.txt (or \code{NA} if file absent).}
#'       \item{cnv.markers}{Count of heterozygous SNPs used in segment from _BAF.txt (or density-derived).}
#'       \item{cnv.totalmarkers}{Total SNPs in segment from _BAF.txt (or density-derived).}
#'       \item{cnv.markers_density, cnv.totalmarkers_density}{Marker density per base pair.}
#'       \item{cnv.cf}{Cancer cell fraction (purity) from _info.txt.}
#'       \item{cnv.ccf}{Clonal cancer cell fraction set to 1 (subclones not parsed).}
#'       \item{cnv.annotation}{CNV status (AMP, DEL, NEUT) assigned via \code{annotateCNV()}.}
#'     }}
#'   \item{@purity}{Purity estimate from _info.txt.}
#'   \item{@ploidy}{Ploidy estimate from _info.txt.}
#'   \item{@effective_ploidy}{Rounded ploidy.}
#'   \item{@sample}{Sample name.}
#'   \item{@method}{\code{"FREEC"}.}
#'   \item{@quality}{List containing:
#'     \describe{
#'       \item{file}{Named character vector of parsed files (information, segments, ratios, BAF).}
#'       \item{quality}{Quality notes from _info.txt.}
#'       \item{best.ploidy}{Optimal ploidy from Control-FREEC log file.}
#'       \item{alt.ploidy}{Alternative ploidy estimate.}
#'     }}
#' }
#'
#' @details
#' \strong{Input files and discovery}
#' Control-FREEC produces multiple output files. This parser searches for and uses:
#' \itemize{
#'   \item \code{_info.txt}: Purity, ploidy, run quality information.
#'   \item \code{_CNVs} or \code{_CNVs.p.value.txt}: Segment coordinates, copy number, status, genotype, uncertainty, p-values.
#'   \item \code{_ratio.txt}: Log-ratio and copy-number windows; parsed for median logR per segment.
#'   \item \code{_BAF.txt}: B-allele frequencies; parsed for marker counts per segment.
#'   \item \code{.log}: Control-FREEC log file; parsed for ploidy estimates.
#' }
#' Missing files are tolerated with warnings; if _info.txt is absent, defaults are applied (ploidy=2, purity=1).
#'
#' \strong{Gap filling}
#' Control-FREEC returns only annotated segments (loss/gain/neutral-LOH). This parser calls \code{fillGaps()}
#' to add regions between segments shorter than 1 Mb, annotating gaps as NEUT and computing marker counts
#' from _ratio.txt and _BAF.txt.
#'
#' \strong{Genotype and LOH}
#' When BAF data is available, Control-FREEC infers genotype at each segment; \code{NA} genotype indicates
#' insufficient SNP coverage. Allele-specific copy numbers are derived from genotype and total copy number.
#'
#' \strong{Subclone and clonality handling}
#' Control-FREEC outputs are currently ambiguous regarding subclone coding (as of v11.5).
#' This parser uses only the main-clone information; \code{cnv.ccf} is set to 1 (clonal).
#' See \url{https://github.com/BoevaLab/FREEC/issues/84}.
#'
#' @examples
#' \dontrun{
#' # Parse Control-FREEC output (auto-discovers associated files)
#' cnv <- parseCfreec(
#'   cfreec.fn = "sample_CNVs.p.value.txt",
#'   filterLowUncertainty = FALSE,
#'   pvalThreshold = 1,
#'   genome = "hg38",
#'   sample = "S1"
#' )
#'
#' # Inspect segments
#' head(cnv@segments)
#' cat("Purity:", cnv@purity, "Ploidy:", cnv@ploidy, "\n")
#'
#' # Check quality metadata
#' cnv@quality$best.ploidy
#' cnv@quality$file
#' }
#'
#' @seealso \code{\link{parseCfreec_file}}, \code{\link{parseCfreec_info}}, \code{\link{parseCfreec_ratio}},
#'   \code{\link{parseCfreec_BAF}}, \code{\link{fillGaps}}, \code{\link{parseCNV}}
#'
#' @export
#' @description
#' Control-FREEC is a tool for detection of copy-number changes and allelic imbalances (including LOH)
#' using deep-sequencing data.
#'
#'
#'
parseCfreec <- function(cfreec.fn,
                        filterLowUncertainty = FALSE,
                        pvalThreshold = 1,
                        genome = "hg38",
                        sample = NULL) {
  # get the proper stem for these files
  freec.name <- predictName(cfreec.fn)
  # load the information file
  info_file <- paste0(freec.name, "_info.txt")
  if (file.exists(info_file)) {
    logger::log_info("Loading info file { basename( info_file) }")
    cf.info <- parseCfreec_info(filename = info_file)
  } else {
    cf.info <- list(
      sampleName = basename(freec.name),
      ploidy = 2,
      purity = 1,
      quality = "Unable to find info file"
    )
    info_file <- NA
  }
  if (!is.null(sample)) {
    cf.info["sampleName"] <- sample
  }
  # load the file with the segments.
  segfile <- paste0(freec.name, "_CNVs")
  if (file.exists(paste0(freec.name, "_CNVs.p.value.txt"))) {
    segfile <- paste0(freec.name, "_CNVs.p.value.txt")
    logger::log_info("Loading file with p-values { basename(segfile) }")
  } else {
    logger::log_info("Loading segments without p-values { basename(segfile)}")
  }
  cf.segments <- parseCfreec_file(segfile, filterLowUncertainty, pvalThreshold, ploidy = round(cf.info$ploidy))

  # load the ratio file to add the median ratio
  ratio_file <- paste0(freec.name, "_ratio.txt")
  if (file.exists(ratio_file)) {
    logger::log_info("Loading ratios file { basename(ratio_file) }")
    cf.segments <- parseCfreec_ratio(filename = paste0(freec.name, "_ratio.txt"), cf.segments = cf.segments, genome = genome)
  } else {
    ratio_file <- NA
  }

  # load the Log File
  logfile <- "ControlFreec.log"
  if (file.exists(paste0(freec.name, ".log"))) {
    logfile <- paste0(freec.name, ".log")
  }
  logger::log_info("Loading log file { basename(logfile) }")
  cf.log <- tryCatch(
    {
      parseCfreec_LOG(logfile)
    },
    error = function(e) {
      logger::log_error("Cannot find the log file from ControlFREEC {logfile}")
      return(NULL)
    }
  )
  # load the BAF file to add the number of markers
  baf_file <- paste0(freec.name, "_BAF.txt")
  if (file.exists(baf_file)) {
    logger::log_info("Loading BAF file { basename(baf_file) }")
    cf.segments <- parseCfreec_BAF(filename = paste0(freec.name, "_BAF.txt"), cf.segments = cf.segments, genome = genome)
  } else {
    baf_file <- NA
  }
  S4Vectors::mcols(cf.segments)["cnv.cf"] <- cf.info$purity
  S4Vectors::mcols(cf.segments)["cnv.ccf"] <- 1 # we don't consider subclones
  cf.segments <- cnvstrmliner::fillGaps(
    cf.segments,
    ploidy = cf.info$ploidy,
    gapsize = 1e6,
    genome = genome,
    segment.annotation = "NEUT"
  )
  obj <- cnvstrml(
    segments = cf.segments,
    method = "FREEC",
    purity = cf.info$purity,
    ploidy = cf.info$ploidy,
    effective_ploidy = round(cf.info$ploidy),
    sample = cf.info$sampleName,
    quality = list(
      "file" = c("segments" = segfile, "info" = info_file, "ratios" = ratio_file, "BAF" = baf_file),
      "quality" = cf.info$quality,
      "best.ploidy" = ifelse(is.null(cf.log$best.ploidy), cf.info$ploidy, cf.log$best.ploidy),
      "alt.ploidy" = ifelse(is.null(cf.log$alt.ploidy), cf.info$ploidy, cf.log$alt.ploidy)
    ),
    genome = genome
  )

  # the object so far does not include markers for the neutral regions, although the information is available in the ratio file.
  S4Vectors::mcols(obj@segments)[, "cnv.markers"] <- round(BiocGenerics::width(obj@segments) * S4Vectors::mcols(obj@segments)[, "cnv.markers_density"])
  S4Vectors::mcols(obj@segments)[, "cnv.totalmarkers"] <- round(BiocGenerics::width(obj@segments) * S4Vectors::mcols(obj@segments)[, "cnv.totalmarkers_density"])


  obj@quality$file <- c("information" = info_file, "segments" = segfile, "ratios" = ratio_file, "BAF" = baf_file)
  obj
}

#' Generate a cnv Granges from a dataframe created by ControlFreec (file ending in .CNVs or contains pvalues.)
#'
#' When the genotype is - it means that the progrm did not find any SNPs to decide the correct genotype, and
#' the estimated copy number depends on the logRR
#' The output includes
#' * chromosome
#' * start position
#' * end position
#' * predicted copy number
#' * type of alterations
#' * genotype (if BAF option was used)
#' * uncertainty of the assignment (not the uncertainty of CNA)
#' * WilcoxonRankSumTestPvalue
#' * KolmogorovSmirnovPvalue
#' * status (somatic or germline)
#' * percentage of germline CNA/LOH
#' Note: at least for v 11.5 it is not possible to understand how subclones are coded in the output files.
#' @param cfreec.fn Character input filename with CNV segments.
#' @param filterLowUncertainty Logical. If TRUE, segments with uncertainty equal to -1 will be discarded.
#' @param pvalThreshold Numeric threshold of p-values to filter when p-values are included in the file.
#' @param genome Character reference genome version to use.
#' @param ploidy Numeric ploidy used when interpreting copy-number states.
#' @return A GRanges object with the segments predicted by ControlFreec

parseCfreec_file <- function(cfreec.fn,
                             filterLowUncertainty = TRUE,
                             pvalThreshold = 0.01,
                             genome = "hg38",
                             ploidy = 2) {
  controlfreec.df <- data.table::fread(cfreec.fn, header = FALSE, sep = "\t")
  if (controlfreec.df[1, 1] == "chr") {
    colnames(controlfreec.df) <- as.character(controlfreec.df[1, ])
    controlfreec.df <- controlfreec.df[-1, ]
  } else {
    colnames(controlfreec.df) <- c("chr", "start", "end", "copy number", "status", "genotype", "uncertainty")
  }
  controlfreec.df <- controlfreec.df %>%
    dplyr::mutate(
      `copy number` = as.numeric(as.character(`copy number`)),
      uncertainty = as.numeric(as.character(uncertainty))
    )
  if ("WilcoxonRankSumTestPvalue" %in% colnames(controlfreec.df)) {
    controlfreec.df <- controlfreec.df %>%
      dplyr::mutate(
        WilcoxonRankSumTestPvalue = as.numeric(as.character(WilcoxonRankSumTestPvalue)),
        KolmogorovSmirnovPvalue = as.numeric(as.character(KolmogorovSmirnovPvalue))
      )
  }
  seqinfo <- getSeqInfo(genome)
  controlfreec.gr <- GenomicRanges::makeGRangesFromDataFrame(controlfreec.df,
    keep.extra.columns = TRUE,
    ignore.strand = TRUE,
  ) %>%
    GenomeInfoDb::sortSeqlevels() %>%
    BiocGenerics::sort()

  GenomeInfoDb::seqlevelsStyle(controlfreec.gr) <- "UCSC"
  GenomeInfoDb::seqlevels(controlfreec.gr) <- GenomeInfoDb::seqlevels(seqinfo)
  GenomeInfoDb::seqinfo(controlfreec.gr) <- seqinfo

  if (any(is.na(S4Vectors::mcols(controlfreec.gr)["genotype"]))) {
    logger::log_warn("The genotype column contains NA values. This is not expected, I will use an estimate")
  }
  mcols_new <- S4Vectors::mcols(controlfreec.gr) %>%
    tibble::as_tibble() %>%
    dplyr::bind_rows(dplyr::tibble(WilcoxonRankSumTestPvalue = numeric(), KolmogorovSmirnovPvalue = numeric())) %>%
    dplyr::mutate(
      total.allele.copies = as.numeric(copy.number),
      major.allele.copies = ifelse((genotype == "-" | is.na(genotype)), ceiling(total.allele.copies / 2), stringr::str_count(genotype, "A")),
      minor.allele.copies = ifelse((genotype == "-" | is.na(genotype)), total.allele.copies - major.allele.copies, stringr::str_count(genotype, "B")),
      uncertainty = as.numeric(uncertainty),
      WilcoxonRankSumTestPvalue = ifelse(is.na(WilcoxonRankSumTestPvalue), 0, as.numeric(WilcoxonRankSumTestPvalue)),
      KolmogorovSmirnovPvalue = ifelse(is.na(KolmogorovSmirnovPvalue), 0, as.numeric(KolmogorovSmirnovPvalue))
    )
  # flags <- apply(mcols_new, 1, function(K) {
  #   cnvstrmliner::annotateCNV(major.cn = K["major.allele.copies"], minor.cn = K["minor.allele.copies"], ploidy = ploidy)
  # })
  # # loh=cnvstrmliner::annotateLOH( flags)
  # mcols_new <- mcols_new %>%
  #   dplyr::mutate(
  #     cnv.annotation = flags,
  #     # cnv.loh=loh,
  #     cnv.lrr = log2(total.allele.copies) - 1
  #   )

  S4Vectors::mcols(controlfreec.gr) <- mcols_new

  if (filterLowUncertainty == TRUE) {
    keep <- which(S4Vectors::mcols(controlfreec.gr)[, "uncertainty"] > -1)
    controlfreec.gr <- controlfreec.gr[keep]
  }
  # if we have the pvalues in the data frame we filter based on those
  if ("WilcoxonRankSumTestPvalue" %in% colnames(S4Vectors::mcols(controlfreec.gr))) {
    keep <- which(S4Vectors::mcols(controlfreec.gr)[, "WilcoxonRankSumTestPvalue"] < pvalThreshold &
      S4Vectors::mcols(controlfreec.gr)[, "KolmogorovSmirnovPvalue"] < pvalThreshold)
    controlfreec.gr <- controlfreec.gr[keep]
  }

  return(controlfreec.gr)
}

#' Parse the ControlFreec info file
#'
#' @param filename (character) the file to parse (extension _info.txt)
#' @description
#' File with information about the FREEC run
#' Program_Version
#' * Sample_Name
#' * Control_Used - True if Control-FREEC used a matching control file for data analysis
#' * CGcontent_Used - True if Control-FREEC used GC-content to normalize read counts
#' * Mappability_Used - True if Control-FREEC used read mappability data to exclude low mappability regions
#' * Looking_For_Subclones - True if Control-FREEC was looking for subclonal CNAs
#' * Breakpoint_Threshold
#' * Window - 0 if Control-FREEC did not use the window approach but counted reads per exon/amplicon
#' * Number_Of_Reads|Pairs_In_Sample
#' * Number_Of_Reads|Pairs_In_Control
#' * Output_Ploidy
#' * Sample_Purity - 1 if Control-FREEC did not find contamination or was nor looking for it
#' * Good_Polynomial_Fit - True if all EM converged
#'

parseCfreec_info <- function(filename) {
  df <- read.table(filename,
    header = FALSE,
    row.names = 1,
    sep = "\t"
  )

  return(
    list(
      sampleName = df["Sample_Name", ],
      ploidy = as.numeric(df["Output_Ploidy", ]),
      purity = as.numeric(df["Sample_Purity", ]),
      quality = ifelse(df["Good_Polynomial_Fit", ] == "True", "Good fit", "Not good fit")
    )
  )
}

#' Parse the ControlFreec _ratio.txt file
#' @param filename (character) the file to parse (extension _ratio.txt)
#' @param cf.segments GRanges object with the parsed Control-FREEC segments.
#'   Used to map window-level ratio medians back onto segment-level calls.
#' @param genome Character genome build used to load sequence metadata.
#'   Default is \code{"hg38"}.
#' @description
#'
#' file with ratios and predicted copy number alterations for each window. Information for each window:
#' * chromosome
#' * start position
#' * ratio
#' * median ratio for the whole fragment resulted from segmentation
#' * predicted copy number
#' * median(abs(BAF-0.5)) per window if `BAF` option is used
#' if the `BAF` option is used
#' * estimated BAF
#' * genotype
#' * percentage of uncertainty of the predicted genotype, max=100
#' if minimalSubclonePresence is used
#' * copy number in the subclone
#' * subclone prevalence
#'
#' The median ratio for a segment is addded to the segmetns under the logrr column
#'
parseCfreec_ratio <- function(filename, cf.segments, genome = "hg38") {
  seqinfo <- getSeqInfo(genome)

  # Read in chunks to avoid loading the full file into memory.
  # Deduplicate on-the-fly: collapse contiguous windows with the same
  # (Chromosome, MedianRatio) into interval ranges before building GRanges.
  chunk_size <- 10000L
  all_names <- names(data.table::fread(filename, nrows = 0L, header = TRUE))
  col_select <- c("Chromosome", "Start", "MedianRatio")
  col_idx <- match(col_select, all_names)

  result_chunks <- list()
  prev_tail <- NULL
  rows_read <- 0L

  repeat {
    chunk <- tryCatch(
      data.table::fread(
        filename,
        select = col_idx,
        header = FALSE,
        sep = "\t",
        skip = 1L + rows_read,
        nrows = chunk_size
      ),
      error = function(e) data.table::data.table()
    )
    if (nrow(chunk) == 0L) break
    data.table::setnames(chunk, col_select)
    rows_read <- rows_read + nrow(chunk)

    # Prepend the trailing run from the previous chunk so it can be
    # merged if it continues into this chunk.
    if (!is.null(prev_tail)) {
      chunk <- data.table::rbindlist(list(prev_tail, chunk))
    }

    chunk[, run_id := data.table::rleid(Chromosome, MedianRatio)]

    # Hold back the last run — it may continue into the next chunk.
    last_id <- chunk[.N, run_id]
    prev_tail <- chunk[run_id == last_id, ..col_select]
    chunk <- chunk[run_id != last_id]

    if (nrow(chunk) > 0L) {
      deduped <- chunk[, .(
        Chromosome = Chromosome[1L],
        Start = min(Start),
        End = max(Start),
        MedianRatio = MedianRatio[1L]
      ), by = run_id]
      deduped[, run_id := NULL]
      result_chunks <- c(result_chunks, list(deduped))
    }
  }

  # Flush the last held-back run.
  if (!is.null(prev_tail) && nrow(prev_tail) > 0L) {
    deduped <- prev_tail[, .(
      Chromosome = Chromosome[1L],
      Start = min(Start),
      End = max(Start),
      MedianRatio = MedianRatio[1L]
    )]
    result_chunks <- c(result_chunks, list(deduped))
  }

  dt <- data.table::rbindlist(result_chunks)
  rm(result_chunks, prev_tail)

  df <- dt %>%
    GenomicRanges::makeGRangesFromDataFrame(
      keep.extra.columns = TRUE,
      start.field = "Start",
      end.field = "End",
      ignore.strand = TRUE
    ) %>%
    GenomeInfoDb::sortSeqlevels() %>%
    BiocGenerics::sort()
  rm(dt)

  GenomeInfoDb::seqlevelsStyle(df) <- "UCSC"
  GenomeInfoDb::seqlevels(df) <- GenomeInfoDb::seqlevels(seqinfo)
  GenomeInfoDb::seqinfo(df) <- seqinfo

  ov <- GenomicRanges::findOverlaps(query = cf.segments, subject = df)

  lrr <-
    S4Vectors::mcols(df)[S4Vectors::subjectHits(ov), , drop = FALSE] %>%
    as.data.frame() %>%
    dplyr::mutate(MedianRatio = as.numeric(MedianRatio)) %>%
    dplyr::mutate(LogR = suppressWarnings(dplyr::case_when(
      is.na(MedianRatio) ~ NA,
      MedianRatio == 0 ~ NA,
      MedianRatio == -1 ~ -100,
      MedianRatio > 0 ~ log2(MedianRatio),
      .default = -13
    ))) %>%
    dplyr::pull(LogR)
  S4Vectors::mcols(cf.segments)[S4Vectors::queryHits(ov), "cnv.lrr"] <- lrr
  rm(df)
  gc()
  return(cf.segments)
}


#' parse the ControlFreec BAF file
#' @param filename Character filename to parse.
#' @param cf.segments A GRanges object with the segments from the CNV file.
#' @param genome Character genome version to use.
#' @return Updated cf.segments with the number of markers used.
#' @description
#' File B-allele frequencies for each possibly heterozygous SNP position. Information for each window:
#' * chromosome
#' * start position
#' * BAF
#' * fitted frequency of A-allele
#' * fitted frequency of B-allele
#' * inferred frequency of A-allele
#' * inferred frequency of B-allele
#' * percentage of uncertainty of the predicted genotype, max=100.
parseCfreec_BAF <- function(filename, cf.segments, genome = "hg38") {
  seqinfo <- getSeqInfo(genome)
  df <- data.table::fread(filename,
    header = TRUE,
    select = c("Chromosome", "Position", "uncertainty"),
    sep = "\t"
  ) %>%
    GenomicRanges::makeGRangesFromDataFrame(
      keep.extra.columns = TRUE,
      start.field = "Position",
      end.field = "Position",
      ignore.strand = TRUE
    ) %>%
    GenomeInfoDb::sortSeqlevels() %>%
    BiocGenerics::sort()

  GenomeInfoDb::seqlevelsStyle(df) <- "UCSC"
  GenomeInfoDb::seqlevels(df) <- GenomeInfoDb::seqlevels(seqinfo)
  GenomeInfoDb::seqinfo(df) <- seqinfo

  all.markers <- GenomicRanges::countOverlaps(query = cf.segments, subject = df)
  het.markers <- GenomicRanges::countOverlaps(query = cf.segments, subject = df[which(df$uncertainty > -1)])
  S4Vectors::mcols(cf.segments)["cnv.markers"] <- het.markers
  S4Vectors::mcols(cf.segments)["cnv.totalmarkers"] <- all.markers

  return(cf.segments)
}


#' Parse the ControlFREEC log file
#' @param filename Character filename of the log file produced by Control-FREEC.
#' @description
#' The log file that controlFREEC produces contains information about the quality of calls for each
#' ploidy. In particular it contains the RSS Score, Percentage of Genome explained and Contamination
#' Based on these values the algorith picks the best solution, however the second best may also be acceptable
#' e.g.
#' Ploidy  RSS score   Percentage of Genome Explained  Contamination
#' 2   0.112551    0.909453    0.280795
#' 3   0.082889    0.854484    0.229635
#' 4   150.726 0.90639 0
#' ..Best ploidy set to 3 according to the RSS score
#' ..Best ploidy could have been set to 2 according to the percentage of the copy number changes explained by a model with a given ploidy
#'
#' Based on these values we can determine if we should accept the solution proposed by ContolFREEC or consider another one.
#' For example if the best solution is the one with ploidy 3, but the ploidy 2 is also acceptable one may opt for the latter as it is more
#' parsimonious
#' @return list with results containing the details of the solution, best and alternative ploidy
parseCfreec_LOG <- function(filename) {
  if (!file.exists(filename)) {
    return(list(
      scores = NULL,
      best.ploidy = NULL,
      best.ploidy.reason = "Unparsable",
      alt.ploidy = NULL,
      alt.ploidy.reason = "Unparsable"
    ))
  }
  # find where to start parsing from:
  s <- grep("^Ploidy\tRSS score", readLines(filename))
  e <- grep("Best ploidy set to", readLines(filename))
  model.table <- read.table(filename, header = TRUE, skip = s - 1, nrows = e - s - 1, sep = "\t")

  best.ploidy.parse <- readLines(filename)[e] %>% stringr::str_match("Best ploidy set to (\\d+) according to the (.*)")
  best.ploidy <- best.ploidy.parse[2]
  best.ploidy.reason <- best.ploidy.parse[3]
  alternative.ploidy.parse <- readLines(filename)[e + 1] %>% stringr::str_match("Best ploidy could have been set to (\\d+) according to (.*) by")
  alternative.ploidy <- alternative.ploidy.parse[2]
  alternative.ploidy.reason <- alternative.ploidy.parse[3]

  return(list(
    scores = model.table,
    best.ploidy = as.numeric(best.ploidy),
    best.ploidy.reason = best.ploidy.reason,
    alt.ploidy = as.numeric(alternative.ploidy),
    alt.ploidy.reason = alternative.ploidy.reason
  ))
}
