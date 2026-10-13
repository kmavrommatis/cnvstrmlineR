#' parse the output of Sequenza
#' @details
#' This function is composing the contents of the cnvstrml object by integrating the information found in these files.\cr
#' It tries to find as many as possible in the same directory.\cr
#' The files looking for are: \cr
#' \itemize{
#' \item <sequenzat.name>_cn_summary.txt # contains general information about the purity, ploidy, contamination etc.\cr
#' \item <sequenza.name>_subclonal_cn.txt# contains the segments along with the copy number information\cr
#' \item <sequenza.name>_allelic_states.txt # contains the information of segments with allele copies\cr
#' }
#' The following files are not parsed..\cr
#' <sequenza.name>_uncorr_cn.seg\cr
#' <sequenza.name>_iCN.seg\cr
#' <sequenza.name>_muts_expAF.txt\cr
#'
#' @param sequenza.name Character input file name root.\cr
#'        Sequenza produces a few files with this name root and different extensions.
#' @param sample Character sample name. If the input is a VCF file it will use the sample name found in the file.
#' @param genome Character genome version assumed for this dataset.
#' @param methodName Character method label stored in the resulting object.
#' @return a cnvstrml object
#' @description Sequenza produces a number of files.
#' This function is a high level parser that can handle
#' the data produced by Sequenza by combining the information found in the various files it produces
#' to a single cnvstrmliner object
#'
#' @export
parseSequenza <- function(
  sequenza.name,
  sample = NA_character_,
  genome = "hg38",
  methodName = "Sequenza"
) {
  logger::log_debug("Sequenza root is {sequenza.name}", namespace = "parser")
  # get the proper stem for these files
  sequenza.name <- predictName(sequenza.name)
  logger::log_debug("after predictName Sequenza root is {sequenza.name}", namespace = "parser")
  # dispatch to the right function
  general_file <- findFile(paste0(sequenza.name, "_seqzinfo.txt"))
  logger::log_debug("Loading general information from a *seqzinfo.txt", namespace = "parser")
  general <- tryCatch(
    {
      parseSequenza_info(general_file)
    },
    error = function(e) {
      logger::log_warn("There is not general log info file {general_file}, continuing without it.", namespace = "parser")
      return(NA)
    }
  )


  info_file <- findFile(paste0(sequenza.name, "_confints_CP.txt"))
  if (length(info_file) > 0) {
    logger::log_debug("Loading information from {info_file}", namespace = "parser")

    info <- parseSequenza_info(info_file)
    purity <- as.numeric(info[2, "cellularity"])
    purity <- min(1, purity)
    ploidy <- as.numeric(info[2, "ploidy.estimate"])
  } else {
    info <- NULL
    purity <- NA
    ploidy <- NA
  }
  alt_file <- findFile(paste0(sequenza.name, "_alternative_solutions.txt"))
  if (length(alt_file) > 0) {
    alt <- parseSequenza_info(alt_file)
  } else {
    alt <- NULL
  }
  segments_file <- findFile(paste0(sequenza.name, "_segments.txt"))
  logger::log_debug("Loading segments from {segments_file}", namespace = "parser")
  segments <- parseSequenza_segments(segments_file, genome = genome)


  mcols_new <- S4Vectors::mcols(segments) %>% as.data.frame()
  # chromosome	start.pos	end.pos	Bf	N.BAF	sd.BAF	depth.ratio	N.ratio	sd.ratio	CNt	A	B	LPP
  suppressWarnings({
    S4Vectors::mcols(segments) <- mcols_new %>%
      dplyr::mutate(
        total.allele.copies = CNt,
        minor.allele.copies = B,
        major.allele.copies = A,
        cnv.markers = N.BAF,
        cnv.totalmarkers = N.BAF,
        cnv.cf = as.numeric(purity),
        cnv.ccf = 1, # since this is the totality of the clones
        # cnv.annotation = NA,
        # cnv.loh=loh,
        cnv.lrr = ifelse(CNt <= 0, -12, log2(depth.ratio))
      )
  })
  # effective_ploidy = cnvstrmliner:::calculateEffectivePloidy( segments )
  # flags=apply(mcols_new, 1, function(K){
  #   cnvstrmliner::annotateCNV(total.cn = K['CNt'],
  #                             minor.cn = K['B'] ,
  #                             ploidy=effective_ploidy  )})
  # S4Vectors::mcols(segments)['cnv.annotation']=flags
  # loh=cnvstrmliner::annotateLOH( flags )

  if (is.na(sample)) {
    sample <- as.character(general["sample"])
  }

  obj <- cnvstrml(
    segments = segments,
    purity = as.numeric(purity),
    ploidy = as.numeric(ploidy),
    method = methodName,
    # effective_ploidy = effective_ploidy,
    sample = sample,
    quality = list(info = info, alternative.solutions = alt, runinfo = general),
    genome = genome
  )
  obj@quality$file <- c("runinfo" = general_file, "information" = info_file, "alternative solutions" = alt_file, "segments" = segments_file)
  obj
}


#' parse segments
#'
#'
#' Parse the _variants.csv file with the copy number
#' @param filename Character filename with segment-level Sequenza output to parse.
#' @param genome Character genome build assumed for the segment file.
#' @return a  dataframe with the parsed information
parseSequenza_segments <- function(filename, genome = "hg38") {
  filename <- normalizePath(filename)
  if (!file.exists(filename)) {
    stop("Filename ", filename, " does not exist. Aborting !!!")
  }
  # seqlevels(sequenza_subclones.gr)=GenomeInfoDb::seqlevels( seqinfo )

  # chromosome	start.pos	end.pos	Bf	N.BAF	sd.BAF	depth.ratio	N.ratio	sd.ratio	CNt	A	B	LPP
  seqinfo <- getSeqInfo(genome)
  df <- read.table(filename,
    header = TRUE,
    sep = "\t"
  ) %>%
    dplyr::mutate(
      START_T = ifelse(start.pos < end.pos, start.pos, end.pos),
      END_T = ifelse(start.pos < end.pos, end.pos, start.pos),
      start.pos = START_T,
      end.pos = END_T
    ) %>%
    tidyr::drop_na(c("CNt", "A", "B")) %>%
    dplyr::select(-c(START_T, END_T)) %>%
    GenomicRanges::makeGRangesFromDataFrame(.,
      seqnames.field = "chromosome",
      start.field = "start.pos",
      end.field = "end.pos",
      keep.extra.columns = TRUE,
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
parseSequenza_info <- function(filename) {
  filename <- normalizePath(filename)
  if (!file.exists(filename)) {
    warning("Filename ", filename, " with info information does not exist. Aborting !!!")
    return(NULL)
  }

  df <- read.table(filename,
    header = TRUE
  )
  df
}
