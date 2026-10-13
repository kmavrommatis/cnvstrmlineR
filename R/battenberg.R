#' parse the output of Battenberg
#' @description
#' this function is a high level parser that can handle
#' the data produced by Battenberg by combining the information found in various files.
#'
#' @param battenberg.name Character input file name root.
#'        Battenberg produces a few files with this name root and different extensions.
#' @param sample Character sample name. If the input is a VCF file it will use the sample name found in the file.
#' @param genome Character genome version assumed for this dataset.
#' @param methodName Character method label stored in the resulting object.
#' @param pvalThreshold Numeric threshold of the p-value used to decide if a segment is subclonal.
#' @return a cnvstrml object
#' @details Battenberg produces a number of files.
#' This function is composing the contents of the cnvstrml object by integrating the information found in these files.\cr
#' The files looking for are \cr
#' \itemize{
#' \item<Battenberg.name>_subclones.txt # contains the segments with their CNV annotation.
#' \item<Battenberg.name>_cellularity_ploidy.txt# contains the information about purity and ploidy or alternative <Battenberg.name>._rho_and_psi.txt
#' }
#'
#' The Battenberg output contains multiple solutions with potentially multiple clones, but typically 1 or 2 per solution.
#' This parser keeps the first solution and assigns the segment with the highest fraction (i.e. most abundant) to the Clone1,
#' and if there is a second clone in the same solution for the same segment, that is assigned to Clone2.
#' As a result the 'Clone1' and 'Clone2' are arbitrarily populated with segments. \cr
#' The default object from this parser is populated by Clone1.
#'
#' @export
parseBattenberg <- function(
  battenberg.name,
  sample = NA_character_,
  genome = "hg38",
  methodName = "Battenberg",
  pvalThreshold = 0.05
) {
  # get the proper stem for these files
  battenberg.name <- predictName(battenberg.name)

  # dispatch to the right function

  info_file <- paste0(battenberg.name, "_cellularity_ploidy.txt")
  if (!file.exists(info_file)) {
    info_file <- findFile(paste0(battenberg.name, "_rho_and_psi*.txt"))
  }
  logger::log_debug("Parsing file ", info_file, " for getting ploidy, purity etc")
  info.orig <- parseBattenberg_info(info_file)
  if (nrow(info.orig) > 1) {
    info <- info.orig[which(info.orig$is.best == TRUE), ]
    info["cellularity"] <- info["rho"]
  } else {
    info <- info.orig
  }

  segments_file <- findFile(paste0(battenberg.name, "_subclones*.txt"))
  segments <- parseBattenberg_segments(segments_file,
    pval_threshold = pvalThreshold,
    genome = genome
  )


  mcols_new <- S4Vectors::mcols(segments) %>%
    as.data.frame() %>%
    rowwise() %>%
    dplyr::mutate(
      Clone1_A = ifelse(pval > pvalThreshold, nMaj1_A, ifelse(frac1_A > frac2_A, nMaj1_A, nMaj2_A)),
      Clone1_B = ifelse(pval > pvalThreshold, nMin1_A, ifelse(frac1_A > frac2_A, nMin1_A, nMin2_A)),
      Clone1_ccf = max(frac1_A, frac2_A, na.rm = TRUE),
      Clone2_A = ifelse(pval > pvalThreshold, NA, ifelse(frac1_A < frac2_A, nMaj1_A, nMaj2_A)),
      Clone2_B = ifelse(pval > pvalThreshold, NA, ifelse(frac1_A < frac2_A, nMin1_A, nMin2_A)),
      Clone2_ccf = ifelse(pval > pvalThreshold, NA, min(frac2_A, frac2_A, na.rm = TRUE))
    ) %>%
    dplyr::ungroup()

  suppressWarnings({
    S4Vectors::mcols(segments) <- mcols_new %>%
      dplyr::mutate(
        total.allele.copies = Clone1_A + Clone1_B,
        major.allele.copies = Clone1_A,
        minor.allele.copies = Clone1_B,
        cnv.markers = NA,
        cnv.ccf = Clone1_ccf,
        cnv.cf = Clone1_ccf * as.numeric(info["cellularity"]),
        # cnv.annotation = NA,
        cnv.lrr = LogR
      )
  })
  # effective_ploidy =cnvstrmliner:::calculateEffectivePloidy( segments )
  # flags=apply(mcols_new, 1, function(K){
  #   cnvstrmliner::annotateCNV(major.cn = K['Clone1_A'],
  #                             minor.cn = K['Clone1_B'] ,
  #                             ploidy=effective_ploidy)})
  # S4Vectors::mcols( segments )['cnv.annotation']=flags
  # loh=cnvstrmliner::annotateLOH( flags )

  if (is.na(sample)) {
    sample <- "not provided"
  }
  # sex=toupper(info[1,'sex_estimated'])
  sex <- NA_character_
  obj <- cnvstrml(
    segments = segments,
    purity = as.numeric(info["cellularity"]),
    ploidy = as.numeric(info["ploidy"]),
    # effective_ploidy = effective_ploidy,
    sex = sex,
    method = methodName,
    sample = sample,
    quality = list(bb_info = info.orig, info_file = basename(info_file)),
    genome = genome
  )
  obj
}


#' parse segments
#'
#' @description
#' Parse the _variants.csv file with the copy number
#' @param filename Character filename with Battenberg segments to parse.
#' @param genome Character genome build assumed for the segment file.
#' @param pval_threshold Numeric threshold of the p-value used to determine if a subclone is present.
#' @return a  dataframe with the parsed information
#' @details
#' This function parses the segments from the file generated by Battenberg \cr
#' In general Battenberg adds the following columns (https://github.com/Wedge-lab/battenberg) \cr
#' \itemize{
#' \item chr	The chromosome of the segment
#' \item startpos	Start position on the chromosome
#' \item endpos	End position on the chromosome
#' \item BAF	The B-allele frequency of the segment
#' \item pval	P-value that is obtained when testing whether this segment should be represented by one or two states. A low p-value will result in the fitting of a second copy number state
#' \item LogR	The log ratio of normalised tumour coverage versus its matched normal sequencing sample
#' \item ntot	An internal total copy number value used to determine the priority of solutions. NOTE: This is not the total copy number of this segment!
#' \item nMaj1_A	The major allele copy number of state 1 from solution A
#' \item nMin1_A	The minor allele copy number of state 1 from solution A
#' \item frac1_A	Fraction of tumour cells carrying state 1 in solution A
#' \item nMaj2_A	The major allele copy number of state 2 from solution A. This value can be NA
#' \item nMin2_A	The minor allele copy number of state 2 from solution A. This value can be NA
#' \item frac2_A	Fraction of tumour cells carrying state 2 in solution A. This value can be NA
#' }
#'
#' Based on comments found in the Battenberg parser from phyloWGS
#' The p-values correspond "to whether a
#  segment should be clonal or subclonal copynumber. We first fit a
#  clonal copynumber profile for the whole sample and then perform a
#  simple two-sided t-test twhere the null hypothesis is: A particular
#  segment is clonal. And the alternative: It is subclonal."
#
#  Thus: if t-test falls below significance threshold, we push cnv1 to
#  clonal frequency.

parseBattenberg_segments <- function(filename,
                                     genome = "hg38",
                                     pval_threshold = 0.05) {
  filename <- normalizePath(filename)
  if (!file.exists(filename)) {
    stop("Filename ", filename, " does not exist. Aborting !!!")
  }
  # seqlevels(Battenberg_subclones.gr)=GenomeInfoDb::seqlevels( seqinfo )

  seqinfo <- getSeqInfo(genome)
  df <- read.table(filename,
    header = TRUE,
    sep = "\t"
  ) %>%
    dplyr::mutate(
      START_T = ifelse(startpos < endpos, startpos, endpos),
      END_T = ifelse(startpos < endpos, endpos, startpos),
      startpos = START_T,
      endpos = END_T,
      chr = ifelse(chr %in% c(seq(1, 22), "X", "Y", "MT"), paste0("chr", chr), chr),
      chr = ifelse(chr == "chrMT", "chrM", chr)
    ) %>%
    tidyr::drop_na(startpos, endpos, frac1_A, nMaj1_A, nMin1_A) %>%
    dplyr::select(-c(START_T, END_T)) %>%
    GenomicRanges::makeGRangesFromDataFrame(.,
      keep.extra.columns = TRUE,
      start.field = "startpos",
      end.field = "endpos",
      ignore.strand = TRUE,
      seqinfo = seqinfo
    ) %>%
    GenomicRanges::trim(use.names = TRUE)
  df
}


#' parse information
#'
#' Parse the file that has the information about purity/ploidy etc
#' @param filename (character) the filename with curated information to parse (extension .csv)
#' @return a dataframe with the parsed information.
#'
#' The text file has the format
#' cellularity	ploidy	psi
#' 0.97	2.00908702395302	2.00352916192448
parseBattenberg_info <- function(filename) {
  filename <- normalizePath(filename)
  if (!file.exists(filename)) {
    stop("Filename ", filename, " does not exist. Aborting !!!")
  }

  df <- read.table(filename,
    header = TRUE
  )

  df
}


#' update subclones
#'
#' Battenberg contains information about subclones
#' This function creates a new object with information coming only from the
#' requested clone
#' @param obj (cnvstrml) the cnvstrml object
#' @param clone (character) the name of the clone to upate with (e.g. Clone1). should match the name of the clone in the segments slot.
#' @return a new cnvstrml object with the information of the clone.
#' @export
#'
updateBattenbergClone <- function(obj, clone) {
  segments <- obj@segments

  mcols_new <- S4Vectors::mcols(segments) %>% as.data.frame()

  major.cn.column <- paste0(clone, "_A")
  minor.cn.column <- paste0(clone, "_B")
  ccf.column <- paste0(clone, "_ccf")

  # check if the clone exists
  if (any(grepl(major.cn.column, colnames(mcols_new))) &&
    any(grepl(minor.cn.column, colnames(mcols_new)))) {
    # we found the clone
  } else {
    # we don't have the clone so we need to show an error
    stop("Clone ", clone, " does not exist in this object")
  }

  mcols_new <- mcols_new %>%
    dplyr::mutate(
      A = ifelse(is.na(.[[major.cn.column]]), Clone1_A, .[[major.cn.column]]),
      B = ifelse(is.na(.[[minor.cn.column]]), Clone1_B, .[[minor.cn.column]]),
      CopyNr = A + B
    )

  S4Vectors::mcols(segments) <- mcols_new %>%
    dplyr::mutate(
      total.allele.copies = CopyNr,
      minor.allele.copies = B,
      major.allele.copies = A,
      cnv.cf = get(ccf.column) * as.numeric(obj@purity),
      # cnv.annotation = NA,
      cnv.loh = NA,
      cnv.lrr = cnv.lrr
    )
  # effective_ploidy <- cnvstrmliner:::calculateEffectivePloidy(segments)
  # flags <- apply(mcols_new, 1, function(K) {
  #   cnvstrmliner::annotateCNV(total.cn = K["CopyNr"], minor.cn = K["B"], ploidy = obj@effective_ploidy)
  # })
  # loh=cnvstrmliner::annotateLOH( flags )
  # S4Vectors::mcols(segments)["cnv.annotation"] <- flags
  # S4Vectors::mcols(segments)['cnv.loh']=loh

  obj <- cnvstrml(
    segments = segments,
    purity = as.numeric(obj@purity),
    ploidy = as.numeric(obj@ploidy),
    # effective_ploidy = obj@effective_ploidy,
    sex = obj@sex,
    method = obj@method,
    sample = paste(obj@sample, clone),
    quality = obj@quality,
    genome = obj@genome
  )
  obj
}
