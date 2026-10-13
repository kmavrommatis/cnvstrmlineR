#' parse the output of SCLUST
#' @details
#' This function is composing the contents of the cnvstrml object by integrating the information found in these files.\cr
#' It tries to find as many as possible in the same directory.\cr
#' The files looking for are: \cr
#' \itemize{
#' \item <sclust.name>_cn_summary.txt # contains general information about the purity, ploidy, contamination etc.\cr
#' \item <sclust.name>_subclonal_cn.txt# contains the segments along with the copy number information\cr
#' \item <sclust.name>_allelic_states.txt # contains the information of segments with allele copies\cr
#' }
#' The following files are not parsed..\cr
#' <sclust.name>_uncorr_cn.seg\cr
#' <sclust.name>_iCN.seg\cr
#' <sclust.name>_muts_expAF.txt\cr
#'
#' @param sclust.name Character input file name root.\cr
#'        SCLUST produces a few files with this name root and different extensions.
#' @param sample Character sample name. If the input is a VCF file it will use the sample name found in the file.
#' @param genome Character genome version assumed for this dataset.
#' @param methodName Character method label stored in the resulting object.
#' @return a cnvstrml object
#' @description SCLUST produces a number of files.
#' This function is a high level parser that can handle
#' the data produced by SCLUST by combining the information found in the various files it produces
#' to a single cnvstrmliner object
#'
#' @export
parseSCLUST <- function(
  sclust.name,
  sample = NA_character_,
  genome = "hg38",
  methodName = "SCLUST"
) {
  logger::log_debug("SCLUST root is ", sclust.name)
  # get the proper stem for these files
  sclust.name <- predictName(sclust.name)
  logger::log_debug("after predictName SCLUST root is ", sclust.name)
  # dispatch to the right function

  info_file <- findFile(paste0(sclust.name, "_cn_summary.txt"))

  info <- parseSCLUST_info(info_file)

  segments_file <- findFile(paste0(sclust.name, "_allelic_states.txt"))
  segments <- parseSCLUST_segments(segments_file)

  subclones_file <- findFile(paste0(sclust.name, "_subclonal_cn.txt"))
  subclones <- parseSCLUST_subclones(subclones_file)
  # merge the information from the two objects
  mo <- GenomicRanges::findOverlaps(subject = segments, query = subclones)
  for (i in setdiff(colnames(S4Vectors::mcols(subclones)), colnames(S4Vectors::mcols(segments)))) {
    S4Vectors::mcols(segments)[i] <- NA
    S4Vectors::mcols(segments)[S4Vectors::subjectHits(mo), i] <- S4Vectors::mcols(subclones)[S4Vectors::queryHits(mo), i]
  }
  mcols_new <- S4Vectors::mcols(segments) %>% as.data.frame()

  suppressWarnings({
    S4Vectors::mcols(segments) <- mcols_new %>%
      dplyr::mutate(
        total.allele.copies = CopyNr,
        minor.allele.copies = B,
        major.allele.copies = A,
        cnv.markers = n_SNPs,
        cnv.totalmarkers = n_SNPs,
        cnv.cf = as.numeric(info[1, "purity"]),
        cnv.ccf = 1, # since this is the totality of the clones
        # cnv.annotation = NA,
        # cnv.loh=loh,
        cnv.lrr = ifelse(Copy_Nr_Raw <= 0, -12, log2(Copy_Nr_Raw / as.numeric(info[1, "ploidy"])))
      )
  })
  # effective_ploidy = cnvstrmliner:::calculateEffectivePloidy( segments )
  # flags=apply(mcols_new, 1, function(K){
  #   cnvstrmliner::annotateCNV(total.cn = K['CopyNr'],
  #                             minor.cn = K['B'] ,
  #                             ploidy=effective_ploidy  )})
  # S4Vectors::mcols(segments)['cnv.annotation']=flags
  # loh=cnvstrmliner::annotateLOH( flags )

  if (is.na(sample)) {
    sample <- as.character(info["sample_name"])
  }
  sex <- toupper(info[1, "sex_estimated"])
  if (is.logical(sex)) {
    sex <- NA_character_
  }

  obj <- cnvstrml(
    segments = segments,
    purity = as.numeric(info[1, "purity"]),
    ploidy = as.numeric(info[1, "ploidy"]),
    sex = sex,
    method = methodName,
    # effective_ploidy = effective_ploidy,
    sample = sample,
    quality = list(
      "status" = info[1, "status"],
      "file" = c("segments" = segments_file, "subclones" = subclones_file, "info" = info_file)
    ),
    genome = genome
  )
  obj
}


#' parse segments
#'
#'
#' Parse the _variants.csv file with the copy number
#' @param filename Character filename of the SCLUST segments file.
#' @param genome Character genome build assumed for the segment file.
#' @return a  dataframe with the parsed information
parseSCLUST_segments <- function(filename, genome = "hg38") {
  filename <- normalizePath(filename)
  if (!file.exists(filename)) {
    stop("Filename ", filename, " does not exist. Aborting !!!")
  }
  # seqlevels(sclust_subclones.gr)=GenomeInfoDb::seqlevels( seqinfo )
  seqinfo <- getSeqInfo(genome)
  df <- read.table(filename,
    header = TRUE,
    sep = "\t"
  ) %>%
    dplyr::mutate(
      START_T = ifelse(Start < End, Start, End),
      END_T = ifelse(Start < End, End, Start),
      Start = START_T,
      End = END_T
    ) %>%
    dplyr::select(-c(START_T, END_T)) %>%
    GenomicRanges::makeGRangesFromDataFrame(.,
      keep.extra.columns = TRUE,
      ignore.strand = TRUE,
      seqinfo = seqinfo
    ) %>%
    GenomicRanges::trim(use.names = TRUE)
  df
}


#' parse subclones
#'
#' Parse the _subclonal_cn.txt file with the copy number
#' @param filename Character filename with subclones produced by SCLUST.
#' @param genome Character genome build assumed for the subclone file.
#' @return a  dataframe with the parsed information

parseSCLUST_subclones <- function(filename, genome = "hg38") {
  filename <- normalizePath(filename)
  if (!file.exists(filename)) {
    stop("Filename ", filename, " does not exist. Aborting !!!")
  }
  sclust_subclones <- read.table(filename, header = TRUE) %>%
    dplyr::mutate(
      START_T = ifelse(Start < End, Start, End),
      END_T = ifelse(Start < End, End, Start),
      Start = START_T,
      End = END_T
    ) %>%
    dplyr::select(-c(START_T, END_T))
  if (nrow(sclust_subclones) > 0) {
    seqinfo <- getSeqInfo(genome)
    sclust_subclones.gr <- GenomicRanges::makeGRangesFromDataFrame(sclust_subclones,
      keep.extra.columns = TRUE,
      ignore.strand = TRUE,
      seqinfo = seqinfo
    ) %>%
      GenomicRanges::trim(use.names = TRUE)
    # we now have two objects, one with the segments and one with subclonal information
  } else {
    return(GenomicRanges::GRanges())
  }
  return(sclust_subclones.gr)
}


#' parse information
#'
#' Parse the file that has the information about purity/ploidy etc
#' @param filename (character) the filename with curated information to parse (extension .csv)
#' @return a dataframe with the parsed information.
#'
parseSCLUST_info <- function(filename) {
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
#' Sclust contains information about subclones
#' This function creates a new object with information coming only from the
#' requested clone
#' @param obj (cnvstrml) the cnvstrml object
#' @param clone (character) the name of the clone to upate with (e.g. Clone1). should match the name of the clone in the segments slot.
#' @return a new cnvstrml object with the information of the clone.
#' @export
#'
updateSCLUSTClone <- function(obj, clone) {
  segments <- obj@segments

  mcols_new <- S4Vectors::mcols(segments) %>% as.data.frame()

  major.cn.column <- paste0(clone, "_A")
  minor.cn.column <- paste0(clone, "_B")
  cf.column <- paste0(clone, "_Fraction")

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
      A = ifelse(is.na(.[[major.cn.column]]), A, .[[major.cn.column]]),
      B = ifelse(is.na(.[[minor.cn.column]]), B, .[[minor.cn.column]]),
      CopyNr = A + B
    )


  S4Vectors::mcols(segments) <- mcols_new %>%
    dplyr::mutate(
      total.allele.copies = CopyNr,
      minor.allele.copies = B,
      major.allele.copies = A,
      cnv.ccf = get(cf.column),
      cnv.cf = cnv.ccf * as.numeric(obj@purity),
      # cnv.annotation = NA,
      cnv.loh = NA,
      cnv.lrr = Copy_Nr_Raw
    )
  S4Vectors::mcols(segments)[which(is.na(S4Vectors::mcols(segments)["cnv.ccf"])), "cnv.ccf"] <- 1
  S4Vectors::mcols(segments)[which(is.na(S4Vectors::mcols(segments)["cnv.cf"])), "cnv.cf"] <- obj@purity


  # effective_ploidy = cnvstrmliner:::calculateEffectivePloidy(segments)
  # flags=apply(mcols_new, 1, function(K){
  #   cnvstrmliner::annotateCNV(total.cn = K['CopyNr'], minor.cn = K['B'],ploidy = effective_ploidy )})
  # S4Vectors::mcols(segments)['cnv.annotation']=flags
  # loh=cnvstrmliner::annotateLOH( flags )
  # S4Vectors::mcols(segments)['cnv.loh']=loh
  obj <- cnvstrml(
    segments = segments,
    purity = as.numeric(obj@purity),
    ploidy = as.numeric(obj@ploidy),
    # effective_ploidy = effective_ploidy,
    sex = obj@sex,
    method = obj@method,
    sample = paste(obj@sample, clone),
    quality = obj@quality,
    genome = obj@genome
  )

  obj
}
