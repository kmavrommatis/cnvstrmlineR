#' parse the output of PureCN
#'
#' this function is a high level parser that can handle
#' either data frames from PureCN or vcf files
#'
#' @param purecn.filename Character input file name root.
#'        PureCN produces a few files with this name root and different extensions.
#' @param sample Character sample name. If the input is a VCF file it will use the sample name found in the file.
#' @param genome Character genome version assumed for this dataset.
#' @param methodName Character method label stored in the resulting object.
#' @return a cnvstrml object
#' @description PureCN produces a number of files.
#' This function is composing the contents of the cnvstrml object by integrating the information found in these files.
#' It tries to find as many as possible
#' The files looking for are
#' <purecn.filename>.csv # contains general information about the purity, ploidy, contamination etc.
#' <purecn.filename>_loh.csv # contains the segments along with the copy number information.
#' additional files that are available but not parsed
#' <purecn.filename>_genes.csv
#' <purecn.filename>_variants.csv
#' If the estimated allele copy numbers are not integers a warning is issued and the values are rounded to the closest integer.
#' If a copy number change is present in only a fraction of the tumor cells (sub-clonality),
#' the average copy number for that segment will be a non-integer.
#' For example, if half the tumor cells have 2 copies and half have 3, the resulting value will be 2.5.
#' In other cases in samples with low purity or high noise (e.g., poor quality FFPE data),
#' the algorithm may struggle to find a perfect integer fit.
#' The effective ploidy is calculated from the segments as the most represented allele-specific copy number state across the genome.
#' It is used as the baseline for the copy number based CCF of subclonal segments (see \code{\link{assignSubclonalFraction}})
#' and is passed to the cnvstrml object, which updates the cnv.annotation field based on it.
#' The cnv.ccf field is estimated for segments with a subclonal copy number (CN.SUBCLONAL in the variants file)
#' and is 1 for all other segments.
#' @export
parsePureCN <- function(
  purecn.filename,
  sample = NA_character_,
  genome = "hg38",
  methodName = "PureCN"
) {
  logger::log_debug("PureCN root is {purecn.filename}", namespace = "parser")
  # get the proper stem for these files
  purecn.filename <- predictName(purecn.filename)
  logger::log_debug("after predictName PureCN root is {purecn.filename}", namespace = "parser")

  info_file <- paste0(purecn.filename, ".csv")
  logger::log_debug("Loading general information from {info_file}", namespace = "parser")
  info <- parsePureCN_info(info_file)

  segments_file <- paste0(purecn.filename, "_loh.csv")
  logger::log_debug("Loading segments from {segments_file}", namespace = "parser")
  segments <- parsePureCN_loh(segments_file)
  logger::log_info("Removing segments with NA copy number values", namespace = "parser")
  segments <- segments[!is.na(segments$M)]

  fractions_file <- paste0(purecn.filename, "_variants.csv")
  logger::log_info("Loading segment fractions for non clonal events", namespace = "parser")
  fractions <- parsePureCN_fractions(fractions_file, genome)
  genes_file <- paste0(purecn.filename, "_genes.csv")

  logger::log_debug("Processing segments and adding CNV annotations", namespace = "parser")

  if (any(S4Vectors::mcols(segments)$C != round(S4Vectors::mcols(segments)$C))) {
    logger::log_warn("Some CNV annotations could not be assigned precisely due to non-integer copy number values. Rounding to closest integer. Compare 'C' and 'M' values to 'total.allele.copies' and 'major.allele.copies' ", namespace = "parser")
  }

  mcols_new <- S4Vectors::mcols(segments) %>% as.data.frame()
  # we give the ploidy called by the method not effective ploidy

  S4Vectors::mcols(segments) <- mcols_new %>%
    dplyr::mutate(
      total.allele.copies = round(C),
      minor.allele.copies = ifelse(is.na(M), round(C), round(M)),
      major.allele.copies = ifelse(is.na(minor.allele.copies), total.allele.copies, total.allele.copies - minor.allele.copies),
      cnv.totalmarkers = num.mark,
      cnv.markers = num.snps,
      cnv.cf = as.numeric(info[1, "Purity"]),
      cnv.ccf = 1,
      # cnv.annotation = "NA",
      cnv.lrr = seg.mean
    )
  # calculate effective ploidy, the baseline for the copy number based CCF.
  # the cnvstrml constructor annotates the segments with this value
  effective_ploidy <- calculateEffectivePloidy(segments)
  logger::log_debug("Effective ploidy is {effective_ploidy}", namespace = "parser")

  logger::log_info("Estimating the cellular fraction of subclonal segments", namespace = "parser")
  segments <- assignSubclonalFraction(
    segments, fractions,
    purity = as.numeric(info[1, "Purity"]),
    ploidy = effective_ploidy
  )

  if (is.na(sample)) {
    sample <- as.character(info["Sampleid"])
  }
  logger::log_debug("Sample identified as {sample}", namespace = "parser")
  sex <- info[1, "Sex"]
  if (is.logical(sex)) {
    sex <- NA_character_
  } # purecn was not able to detect the sex
  logger::log_debug("Creating cnvstrml object with purity {info[1,'Purity']} and ploidy {info[1,'Ploidy']}", namespace = "parser")
  obj <- cnvstrml(
    segments = segments,
    purity = as.numeric(info[1, "Purity"]),
    ploidy = as.numeric(info[1, "Ploidy"]),
    effective_ploidy = effective_ploidy,
    sex = sex,
    method = methodName,
    sample = sample,
    quality = list(
      "comment" = info[1, "Comment"],
      "contamination" = info[1, "Contamination"],
      "file" = c("segments" = purecn.filename, "info" = info_file, "fractions" = fractions_file, "genes" = genes_file)
    ),
    genome = genome
  )
  obj
}



#' parse information
#'
#' Parse the file that has the information about purity/ploidy etc
#' @param filename (character) the filename with curated information to parse (extension .csv)
#' @return a dataframe with the parsed information.
#'
parsePureCN_info <- function(filename) {
  filename <- normalizePath(filename)
  if (!file.exists(filename)) {
    logger::log_error("Filename {filename} does not exist. Aborting !!!", namespace = "parser")
    stop("Filename ", filename, " does not exist. Aborting !!!")
  }

  logger::log_debug("Reading PureCN info file {filename}", namespace = "parser")
  df <- read.csv(filename,
    header = TRUE
  )

  quality <- df[1, "Comment"]

  quality <- quality %>%
    stringr::str_split(";", simplify = TRUE) %>%
    sapply(function(x) {
      if (grepl("POOR GOF", x)) {
        return("POOR GOF")
      } else {
        return(x)
      }
    }) %>%
    paste(., collapse = ";")
  df[1, "Comment"] <- quality
  logger::log_debug("Processed PureCN quality information", namespace = "parser")
  df
}


#' parse segments
#'
#'
#' Parse the _variants.csv file with the copy number
#' @param filename (character) the filename withs somatic status and cellular fraction to parse (extension _variants.csv)
#' @return a  dataframe with the parsed information
parsePureCN_fractions <- function(filename, genome) {
  filename <- normalizePath(filename)
  if (!file.exists(filename)) {
    logger::log_error("Filename {filename} does not exist. Aborting !!!", namespace = "parser")
    stop("Filename ", filename, " does not exist. Aborting !!!")
  }

  logger::log_debug("Reading PureCN fractions file {filename}", namespace = "parser")
  seqinfo <- getSeqInfo(genome)
  df <- read.table(filename,
    header = TRUE,
    sep = ","
  ) %>%
    GenomicRanges::makeGRangesFromDataFrame(.,
      keep.extra.columns = TRUE,
      ignore.strand = TRUE,
      seqinfo = seqinfo,
    )

  df
}


#' parse amplifications/deletions
#'
#'
#' Parse the _genes.csv with gene level information
#' @param filename (character) the filename withs somatic status and cellular fraction to parse (extension _variants.csv)
#' @return a  dataframe with the parsed information
parsePureCN_ampdel <- function(filename) {
  filename <- normalizePath(filename)
  if (!file.exists(filename)) {
    logger::log_error("Filename {filename} does not exist. Aborting !!!", namespace = "parser")
    stop("Filename ", filename, " does not exist. Aborting !!!")
  }
  logger::log_debug("Reading PureCN amplifications/deletions file {filename}", namespace = "parser")
  df <- read.table(filename,
    header = TRUE,
    sep = ","
  )

  df
}


#' parse LOH (segments)
#'
#'
#' Parse the _loh.csv with segment level information.
#' This is the main function for parsing the segment coordinates and data.
#' @param filename Character filename with segment-level LOH information to parse.
#' @param genome Character genome build assumed for the LOH file.
#' @return a  dataframe with the parsed information
parsePureCN_loh <- function(filename, genome = "hg38") {
  filename <- normalizePath(filename)
  if (!file.exists(filename)) {
    logger::log_error("Filename {filename} does not exist. Aborting !!!", namespace = "parser")
    stop("Filename ", filename, " does not exist. Aborting !!!")
  }

  logger::log_debug("Reading PureCN LOH segments file {filename} for genome {genome}", namespace = "parser")
  seqinfo <- getSeqInfo(genome)

  df <- data.table::fread(
    filename,
    header = TRUE,
    sep = ","
  ) %>%
    GenomicRanges::makeGRangesFromDataFrame(
      keep.extra.columns = TRUE,
      ignore.strand = TRUE,
      seqinfo = seqinfo
    )

  logger::log_debug("Created GRanges object with {length(df)} segments", namespace = "parser")
  df
}


#' estimate the clonal fraction of subclonal PureCN segments
#'
#' Intersects the segments with the variants returned by \code{parsePureCN_fractions}.
#' Segments that overlap at least one variant with \code{CN.SUBCLONAL} set to TRUE
#' are treated as subclonal and get a cellular fraction in \code{cnv.ccf}.
#' Two methods are used, in order of preference:
#'
#' \strong{1. Somatic variant based.}
#' For subclonal segments with at least \code{min.somatic} somatic variants (\code{ML.SOMATIC})
#' \deqn{CCF = AR * (p * C + (1 - p) * 2) / (p * M)}
#' where AR is the median allelic ratio of the somatic variants in the segment,
#' p is the purity, C is the fitted (non rounded) total copy number of the segment
#' and M is the multiplicity (default 1).
#'
#' \strong{2. Copy-number based.}
#' All other subclonal segments. The fitted copy number \code{C} is the cell-weighted average
#' of the baseline state \code{b} (the rounded \code{ploidy} argument) and an altered integer state \code{n}:
#' \deqn{C = f * n + (1 - f) * b \Rightarrow f = (C - b) / (n - b)}
#' with \code{n = floor(C)} for losses (\code{C < b}) and \code{ceiling(C)} for gains.
#' If \code{n == b} the fraction is undefined and the previous \code{cnv.ccf} (1) is kept.
#'
#' The estimates are bounded to the range 0 to 1. Segments that are not subclonal keep their
#' previous \code{cnv.ccf} (1 if absent).
#'
#' @param segments GRanges of segments with the \code{C} metadata column.
#' @param fractions GRanges of variants as returned by \code{parsePureCN_fractions}.
#'        Needs the columns \code{CN.SUBCLONAL}, \code{ML.SOMATIC} and the column named in \code{ar.column}.
#' @param purity Numeric purity of the sample (p).
#' @param ploidy Numeric baseline copy number, normally the effective ploidy of the sample
#'        (see \code{calculateEffectivePloidy}). It is rounded to an integer.
#' @param min.somatic Integer minimum number of somatic variants required for the somatic based estimate.
#' @param ar.column Character name of the column in \code{fractions} holding the allelic ratio.
#'        Defaults to \code{"AR"}. In PureCN output \code{ALLELIC.IMBALANCE} is not a ratio
#'        (it holds negative log-scale values).
#' @param M Numeric multiplicity of the mutation. Default 1.
#' @return the segments with the \code{cnv.ccf} column updated for the subclonal segments and the
#'   additional column \code{cnv.ccf.method} (\code{"somatic"}, \code{"copynumber"}
#'   or NA when no estimate was possible).
#' @export
assignSubclonalFraction <- function(segments, fractions, purity, ploidy, min.somatic = 5, ar.column = "AR", M = 1) {
  stopifnot(
    methods::is(segments, "GRanges"), methods::is(fractions, "GRanges"),
    is.numeric(purity), length(purity) == 1, !is.na(purity), purity > 0,
    is.numeric(ploidy), length(ploidy) == 1, !is.na(ploidy)
  )
  fr <- S4Vectors::mcols(fractions)
  missing.cols <- setdiff(c("CN.SUBCLONAL", "ML.SOMATIC", ar.column), colnames(fr))
  if (length(missing.cols) > 0) {
    stop("fractions is missing the column(s): ", paste(missing.cols, collapse = ", "))
  }
  seg.cols <- colnames(S4Vectors::mcols(segments))
  if (!"C" %in% seg.cols) {
    stop("segments needs the C column")
  }
  n <- length(segments)
  C <- S4Vectors::mcols(segments)$C

  # subclonal segments
  sub.hits <- GenomicRanges::findOverlaps(segments, fractions[as.logical(fr$CN.SUBCLONAL) %in% TRUE], ignore.strand = TRUE)
  is.sub <- seq_len(n) %in% S4Vectors::queryHits(sub.hits)

  # copy number based estimate, using the original (non rounded) C
  baseline <- round(ploidy)
  alt <- ifelse(C < baseline, floor(C), ceiling(C))
  ccf.cn <- ifelse(alt == baseline, NA_real_, (C - baseline) / (alt - baseline))

  # somatic variant based estimate, using the original (non rounded) C
  somatic <- as.logical(fr$ML.SOMATIC) %in% TRUE
  som.hits <- GenomicRanges::findOverlaps(segments, fractions[somatic], ignore.strand = TRUE)
  som.ar <- as.numeric(fr[[ar.column]])[somatic][S4Vectors::subjectHits(som.hits)]
  n.somatic <- tabulate(S4Vectors::queryHits(som.hits), nbins = n)
  med.ar <- rep(NA_real_, n)
  if (length(som.ar) > 0) {
    m <- tapply(som.ar, S4Vectors::queryHits(som.hits), stats::median, na.rm = TRUE)
    med.ar[as.integer(names(m))] <- m
  }
  use.somatic <- is.sub & n.somatic >= min.somatic & !is.na(med.ar)
  med.ar[!use.somatic] <- NA_real_
  ccf.som <- med.ar * (purity * C + (1 - purity) * 2) / (purity * M)

  # somatic where possible, copy number for the rest of the subclonal segments
  use.cn <- is.sub & !use.somatic & !is.na(ccf.cn)
  ccf <- rep(NA_real_, n)
  ccf[use.somatic] <- ccf.som[use.somatic]
  ccf[use.cn] <- ccf.cn[use.cn]
  ccf <- pmin(pmax(ccf, 0), 1)
  method <- rep(NA_character_, n)
  method[use.somatic] <- "somatic"
  method[use.cn] <- "copynumber"

  if (!"cnv.ccf" %in% seg.cols) {
    S4Vectors::mcols(segments)$cnv.ccf <- 1
  }
  update <- use.somatic | use.cn
  S4Vectors::mcols(segments)$cnv.ccf[update] <- ccf[update]
  S4Vectors::mcols(segments)$cnv.ccf.method <- method
  logger::log_debug(
    "Estimated CCF for {sum(is.sub)} subclonal segments: {sum(use.somatic)} somatic based, {sum(use.cn)} copy-number based",
    namespace = "parser"
  )
  segments
}
