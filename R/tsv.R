#' parse a simple TSV file
#'
#' this function is a high level parser that can handle
#' the results in a tab delimited file
#' with two columns called cn1 and cn2
#'
#' @param tsv.filename (character) the input file name
#' @param genome (character)  the genome version assumed for this dataset
#' @param sample (sample)  the name of the sample.
#' @param gapsize (numeric)  the minimum gap size to fill between segments.
#' @param ploidy (numeric)  the ploidy of the sample, default is 2.
#' @return a cnvstrml object
#' @description This function parses a tab delimited file
#' It expects 5 or 6 columns:chromosome, start, end, cn1, cn2, snps
#' where
#' * cn1 and cn2 are the copy numbers of the different alleles
#' * snps is the number of heterozygous SNPs in that segment, this
#' value is not necessary.
#'
#' Assuming a normal ploidy of two it calculates the logR
#' If snps is not provided a value for markers that is derived from the
#' length of the segment.
#' A header with column names 'chromosome','start','end','cn1','cn2','snps'
#' is expected
#' @export
parseTSV <- function(
  tsv.filename,
  sample = NA_character_,
  genome = "hg38",
  gapsize = 3e+07,
  ploidy = 2
) {
  purity <- 1


  segments <- suppressWarnings(data.table::fread(
    file = tsv.filename,
    sep = "\t",
    header = TRUE
  ))


  expected_cols <- c("chromosome", "start", "end", "cn1", "cn2")
  if (!setequal(intersect(expected_cols, colnames(segments)), expected_cols)) {
    message <- paste("This input file ", tsv.filename, " does not include allele copies.")
    logger::log_error(message)
    stop(message)
  }
  if (!"snps" %in% colnames(segments)) {
    segments <- segments %>% dplyr::mutate(snps = round((end - start) / 10000))
  }


  segments <- segments %>%
    dplyr::mutate(
      total.allele.copies = cn1 + cn2,
      major.allele.copies = pmax(cn1, cn2),
      minor.allele.copies = pmin(cn1, cn2),
      cnv.markers = snps,
      cnv.totalmarkers = snps,
      cnv.cf = as.numeric(purity),
      cnv.ccf = 1L, # since this is the totality of the clones
      # cnv.annotation = NA,
      cnv.lrr = ifelse((cn1 + cn2) > 0, log2((cn1 + cn2) / 2), -12)
    )

  seqinfo <- getSeqInfo(genome)
  segments <- segments %>%
    GenomicRanges::makeGRangesFromDataFrame(.,
      keep.extra.columns = TRUE,
      ignore.strand = TRUE
    )
  GenomeInfoDb::seqlevelsStyle(segments) <- "UCSC"
  GenomeInfoDb::seqlevels(segments, pruning.mode = "coarse") <- GenomeInfoDb::seqlevels(seqinfo)
  GenomeInfoDb::seqinfo(segments, pruning.mode = "coarse") <- seqinfo
  segments <- cnvstrmliner::fillGaps(
    segments,
    ploidy = ploidy,
    major.allele.copies = ploidy - 1,
    segment.annotation = "NEUT",
    purity = purity,
    genome = genome,
    gapsize = gapsize
  )

  mcols_new <- S4Vectors::mcols(segments) %>% as.data.frame()

  # effective_ploidy <- cnvstrmliner:::calculateEffectivePloidy(segments)
  # #
  # flags <- apply(mcols_new, 1, function(K) {
  #   cnvstrmliner::annotateCNV(
  #     major.cn = K["major.allele.copies"],
  #     minor.cn = K["minor.allele.copies"],
  #     ploidy = effective_ploidy
  #   )
  # })
  # S4Vectors::mcols(segments)["cnv.annotation"] <- flags


  obj <- cnvstrml(
    segments = segments,
    quality = list(
      parser =
        paste("The major and minor allele copies are inferred from total copy number since they were not provided in the input file.",
          "CNVkit does not provide purity information. Arbitrarily set to 1.",
          sep = ";"
        )
    ),
    # effective_ploidy = effective_ploidy,
    method = "CNVkit call",
    sample = sample,
    genome = genome
  )
  obj@quality$file <- c("segments" = tsv.filename)
  obj
}
