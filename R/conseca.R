#' create a cnvstrml object from a dataframe with Conseca
#'
#' @param df Data frame produced by Conseca.
#' @param sample Character sample name.
#' @param constituent Character method to extract prediction from (gatk or freec).
#' @param genome Character genome version assumed for this dataset.
#' @param methodName Character method label stored in the resulting object.
#' @return a cnvstrml object
#'
parseConseca_df <- function(
  df,
  sample = NA_character_,
  constituent = "freec",
  genome = "hg38",
  methodName = "Conseca"
) {
  # convert the input to a dataframe
  if (!is.data.frame(df)) {
    df <- as.data.frame(df)
  }
  # filter by the program used to call the segments
  df <-
    df %>%
    dplyr::filter(grepl(constituent, support)) %>%
    dplyr::select("chromosome", "start", "end", "status", "confidence", starts_with(constituent)) %>%
    dplyr::filter(status != "neutral")


  # convert to GRanges
  # and normalize the chromosome name to UCSC
  conseca.gr <- GenomicRanges::makeGRangesFromDataFrame(df,
    keep.extra.columns = TRUE,
    ignore.strand = TRUE
  ) %>%
    GenomeInfoDb::sortSeqlevels() %>%
    BiocGenerics::sort()
  GenomeInfoDb::seqlevelsStyle(conseca.gr) <- "UCSC"


  # Add the flags
  if (constituent == "gatk") {
    mcols_new <- suppressWarnings(S4Vectors::mcols(conseca.gr) %>%
      tibble::as_tibble() %>%
      dplyr::mutate(
        gatk_cn = as.numeric(gatk_cn),
        total.allele.copies = ifelse(gatk_call == "0", 2, gatk_cn),
        major.allele.copies = total.allele.copies - 1,
        minor.allele.copies = 1,
        lr = gatk_mean_log2_copy_ratio
      ))
  }
  if (constituent == "freec") {
    mcols_new <- suppressWarnings(S4Vectors::mcols(conseca.gr) %>%
      tibble::as_tibble() %>%
      dplyr::mutate(
        freec_cn = as.numeric(freec_cn),
        freec_cn1 = as.numeric(freec_cn1),
        freec_cn2 = as.numeric(freec_cn2),
        freec_cn1 = ifelse(is.na(freec_cn1) & is.na(freec_cn2), round(freec_cn / 2), freec_cn1),
        freec_cn2 = ifelse(is.na(freec_cn2), freec_cn - freec_cn1, freec_cn2),
        total.allele.copies = freec_cn,
        major.allele.copies = freec_cn1,
        minor.allele.copies = freec_cn2,
        lr = 0
      ))
  }


  S4Vectors::mcols(conseca.gr) <- mcols_new

  purity <- 1
  ploidy <- 2
  obj <- cnvstrml(
    segments = conseca.gr,
    purity = as.numeric(purity),
    ploidy = as.numeric(ploidy),
    method = paste(constituent, methodName),
    sample = sample,
    genome = genome
  )
  obj
}


#' parse the output of Conseca
#'
#' this function is a high level parser that can handle
#' data from Conseca.
#' Conseca is an aggregate tool that contains the information fromo ControlFreec and GATK
#'
#'
#' @param conseca.data Character path or data frame containing Conseca data.
#'        if the type of the argument is character it is assumed to be a file.
#'        if the type of the argument is data.frame it is parsed as such.
#' @param sample Character sample name. If the input is a VCF file it will use the sample name found in the file.
#' @param genome Character genome version assumed for this dataset.
#' @return a cnvstrml object
#' @export
parseConseca <- function(
  conseca.data,
  sample = NA_character_,
  genome = "hg38"
) {
  # dispatch to the right function
  if (is.character(conseca.data) && file.exists(conseca.data)) {
    # message("Preparing the segments from the vcf file ", basename( conseca.data),". Purity and Ploidy are expected to be found in the header of the file.")
    conseca.data <- read.table(conseca.data, header = TRUE, sep = "\t")
  }
  if (is.data.frame(conseca.data)) {
    # message("Preparing the segments from a data frame. Purity and Ploidy have to be set independently.")
    obj <- parseConseca_df(conseca.data, sample = sample, genome = genome)
  }


  mcols_new <- S4Vectors::mcols(obj@segments) %>% as.data.frame()
  flags <- apply(mcols_new, 1, function(K) {
    cnvstrmliner::annotateCNV(major.cn = K["major.allele.copies"], minor.cn = K["minor.allele.copies"])
  })

  S4Vectors::mcols(obj@segments) <- mcols_new %>%
    dplyr::mutate(
      cnv.annotation = flags,
      cnv.loh = ifelse(grepl("LOH", cnv.annotation), TRUE, FALSE),
      cnv.lrr = lr
    )

  obj
}
