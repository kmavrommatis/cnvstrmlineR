#' parse the output of SEG file
#'
#' this function is a high level parser that can handle
#' the data produced various methods in the IGV Seg format
#' https://igv.org/doc/desktop/#FileFormats/DataTracks/#seg
#' @param sample.file Character input file name.
#' @param sample Character sample name. If the input is a VCF file it will use the sample name found in the file.
#' @param genome Character genome version assumed for this dataset.
#' @param methodName Character method label stored in the resulting object.
#' @param header Logical indicating whether the file contains a header.
#' @param last.column Character value of either 'CopyNr' or 'lrr' describing the meaning of the last column.
#' @param ploidy Integer expected ploidy. Default is 2.
#' @return a cnvstrml object
#' @description a SEG file (typically .seg or .cbs)
#' is a tab-delimited text file that lists loci and associated numeric values.
#' The first row contains column headings, and each subsequent row contains
#' a locus and an associated numeric value. IGV ignores the column headings.
#' It reads the first four columns as track name, chromosome, start location, and end location.
#' It reads the last column as the numeric value for that locus
#' (if the value is non-numeric, it throws an error).
#'
#' If the last column is set to be a logR value, it will arbitrarily set teh copy numbers to a value based on the assumption that
#'
#' lrr ~ 0  equals ploidy (defaults to 2)
#'
#' lrr < -0.41 equals loss of one copy
#'
#' lrr < -10 equals total loss of all copies
#'
#' lrr > 0.32  equals gain of one copy
#'
#' lrr > 1 equals two or more gains
#'
#'
#' @export
parseSEG <- function(
  sample.file,
  sample = NA_character_,
  genome = "hg38",
  methodName = "SEG",
  header = FALSE,
  last.column = c("CopyNr", "lrr"),
  ploidy = 2
) {
  segments_file <- sample.file

  last.column <- match.arg(last.column)
  segments <- parseSEG_segments(segments_file, header = header)

  if (last.column == "CopyNr") {
    S4Vectors::mcols(segments)["CopyNr"] <- S4Vectors::mcols(segments)[, "Values"]
    S4Vectors::mcols(segments)["cnv.lrr"] <- ifelse(S4Vectors::mcols(segments)[, "CopyNr"] == 0, -12, log2(S4Vectors::mcols(segments)[, "CopyNr"] / ploidy))
  }
  if (last.column == "lrr") {
    S4Vectors::mcols(segments)["CopyNr"] <- cut(S4Vectors::mcols(segments)[, "Values"],
      breaks = c(-10e3, -10, -0.41, 0.32, 1, 10e3),
      labels = c("0", as.character(ploidy - 1), as.character(ploidy), as.character(ploidy + 1), as.character(ploidy + 2))
    ) %>%
      as.character() %>%
      as.numeric()
    S4Vectors::mcols(segments)["cnv.lrr"] <- S4Vectors::mcols(segments)[, "Values"]
  }

  mcols_new <- S4Vectors::mcols(segments) %>% as.data.frame()

  s <- unique(mcols_new$sample)
  if (length(s) > 1) {
    stop("The file ", sample.file, " contains multiple samples.", paste(s), " This is not valid for this parser")
  }
  # merge the information from the two objects
  # flags <- apply(mcols_new, 1, function(K) {
  #   cnvstrmliner::annotateCNV(total.cn = K["CopyNr"], minor.cn = K["minor"])
  # })
  # loh=cnvstrmliner::annotateLOH( flags )

  S4Vectors::mcols(segments) <- mcols_new %>%
    dplyr::select(-Values) %>%
    dplyr::mutate(
      total.allele.copies = CopyNr,
      minor.allele.copies = NA,
      major.allele.copies = NA,
      cnv.markers = n_SNPs,
      cnv.totalmarkers = n_SNPs,
      cnv.cf = 1,
      cnv.loh = NA
    )

  if (is.na(sample)) {
    sample <- as.character(s)
  }
  sex <- NA_character_
  obj <- cnvstrml(
    segments = segments,
    sex = sex,
    method = methodName,
    sample = sample,
    genome = genome
  )
  obj
}


#' parse segments
#'
#'
#' Parse the _variants.csv file with the copy number
#' @param filename Character filename of the SEG file to parse.
#' @param genome Character genome build assumed for this SEG file.
#' @param header Logical indicating whether the file contains a header row.
#' @return a  dataframe with the parsed information
parseSEG_segments <- function(filename, genome = "hg38", header = FALSE) {
  filename <- normalizePath(filename)
  if (!file.exists(filename)) {
    stop("Filename ", filename, " does not exist. Aborting !!!")
  }
  # seqlevels(sclust_subclones.gr)=GenomeInfoDb::seqlevels( seqinfo )
  seqinfo <- getSeqInfo(genome)
  df <- data.table::fread(filename,
    header = header,
    sep = "\t"
  )

  colnames(df)[1:4] <- c("sample", "chromosome", "start", "end")
  if (ncol(df) > 4) {
    colnames(df)[ncol(df)] <- "Values"
  }
  # get the number of markers used for this segment. This is TRUE for sclust, but I don't know if this is still true for other methods
  if (ncol(df) == 6) {
    colnames(df)[5] <- "n_SNPs"
  }
  df <- df %>%
    dplyr::select(., sample, chromosome, start, end, n_SNPs, Values) %>%
    GenomicRanges::makeGRangesFromDataFrame(.,
      keep.extra.columns = TRUE,
      ignore.strand = TRUE,
      seqinfo = seqinfo
    ) %>%
    GenomicRanges::trim(use.names = TRUE)
  df
}
