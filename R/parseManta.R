#' Parse the manta ouput (vcf v4.1)
#'
#' @param sv.vcf (character) a file name (in vcf format)
#' @param genome (character) the genome version to use (hg38). This is passed to VariantAnnotation::readVcf
#' @param sample (character) the name of the sample
#' @return a cnvstrml object, Purity and Ploidy are not set
#' @export


parseManta <- function(sv.vcf, sample = NULL, genome = "hg38") {
  # load the file
  vcf <- VariantAnnotation::readVcf(sv.vcf)
  gr <- StructuralVariantAnnotation::breakpointRanges(vcf) %>%
    # StructuralVariantAnnotation::breakendRanges(vcf)
    GenomeInfoDb::sortSeqlevels() %>%
    S4Vectors::sort()
  if (is.null(sample) || length(sample) == 0) {
    sample <- VariantAnnotation::samples(VariantAnnotation::header(vcf))
    # message("Extracting sample name from the header of the VCF file :" , sample)
  }
  if (is.null(sample) || length(sample) == 0) {
    sample <- sv.vcf %>%
      basename() %>%
      stringi::stri_replace_all_fixed(., c(".vcf.gz", ".vcf.bgz", ".vcf"), replacement = "", vectorize_all = FALSE)
    # message("Extracting sample name from the filename :" , sample)
  }
  # keep only the svs that pass the filters
  gr <- gr[gr$FILTER == "PASS"]
  # obj=cnvstrml(
  #   segments=gr,
  #   purity = 1,
  #   ploidy = 1,
  #   method='manta',
  #   sample=sample,
  #   genome=genome
  # )
  # obj
  gr
}


#' Match CNV segments with SV breakpoints
#'
#' In order to match the SV breakpoints to the CNV segments
#' we expect to find a breakpoint (BND) in the vicinity of the
#' starting and ending coordinates of the CNV segment
#' Furthermore the two breakpoints need to be partners.
#'
#'
#' @param cnvsegments (cnvstrml) the segments predicted to be copy altered
#' @param svsegments (GRanges) the SV ranges from breakpointRanges
#' @param wiggle (numeric) distance between the CNV boundaries and sv breakpoints
#' @description This function is using the SV breakpoints that were called by manta
#' to confirm the segments that were called by the CNV calling algorithm
#' One would expect that the boundaries of a CNV segment will match the breakpoints
#' of SV events. Ideally the same SV event should have breakpoints near the boundaries
#' of the segment, however given how complex some SV events may be this is not always the case.
#' This function overlaps the CNV segments with the SV breakpoints +/- a wiggle distance to account for
#' the differences in the prediction accuracy of each approach.
#' This function adds the 'support' column which takes a numeric value
#' @export
#' @return an updates GRanges with segments and support from SV
#' The resulting segments have the additional columns
#' sv.support: TRUE if a BND is near the boundary of the CNV segment (+/- wiggle)
#' sv.support.chr: The BND of the SV events are on the same chromosome (i.e. not a intra chromosomal change)
#' sv.support.partner: TRUE if The parnter BND for the SV matches the BND on the other end of the CNV segment (highest congruence between the methods)
#' sv.support.class: TRUE if The type of the CNV (e.g. deletion/amplification) and the type of the SV (e.g. deletion) match
matchCNVSV <- function(cnvsegments, svsegments, wiggle = 1000) {
  cnv <- cnvsegments@segments
  gr <- svsegments
  GenomicRanges::end(cnv) <- GenomicRanges::start(cnv) + 1
  cnv <- GenomicRanges::flank(cnv, wiggle)
  ov.start <- GenomicRanges::findOverlaps(query = cnv, subject = gr)
  if (is.null(names(cnv))) {
    names(cnv) <- rownames(as.data.frame(cnv))
  }
  svs.start <- data.frame(
    "segment" = names(cnv)[S4Vectors::queryHits(ov.start)],
    "segment.chr" = seqnames(cnv)[S4Vectors::queryHits(ov.start)],
    "segment.start" = start(cnv)[S4Vectors::queryHits(ov.start)],
    "segment.end" = end(cnv)[S4Vectors::queryHits(ov.start)],
    "segment.width" = width(cnv)[S4Vectors::queryHits(ov.start)],
    "segment.type" = S4Vectors::mcols(cnv)[S4Vectors::queryHits(ov.start), "cnv.annotation"],
    "SV5" = names(gr)[S4Vectors::subjectHits(ov.start)],
    "SV5.chr" = seqnames(gr[S4Vectors::subjectHits(ov.start)]),
    "SV5.start" = start(gr[S4Vectors::subjectHits(ov.start)]),
    "SV5.end" = end(gr[S4Vectors::subjectHits(ov.start)]),
    "SV5.len" = S4Vectors::mcols(gr)[S4Vectors::subjectHits(ov.start), "svLen"],
    "SV5.type" = S4Vectors::mcols(gr)[S4Vectors::subjectHits(ov.start), "svtype"],
    "SV5.partner" = gr[S4Vectors::subjectHits(ov.start)]$partner
  ) %>%
    dplyr::mutate(
      SV5.partner.chr = seqnames(gr[SV5.partner]) %>% as.vector(),
      SV5.partner.start = start(gr[SV5.partner]),
      SV5.partner.end = end(gr[SV5.partner])
    )

  # svs.start=mergeByOverlaps( query=cnv, subject=gr)
  # svs.start=svs.start[,c('cnv','cnv.annotation', 'cnv.cf','gr', 'sourceId','partner','svtype' )]
  cnv <- cnvsegments@segments
  GenomicRanges::start(cnv) <- GenomicRanges::end(cnv) - 1
  cnv <- GenomicRanges::flank(cnv, wiggle)
  if (is.null(names(cnv))) {
    names(cnv) <- rownames(as.data.frame(cnv))
  }
  ov.end <- GenomicRanges::findOverlaps(query = cnv, subject = gr)
  svs.end <- data.frame(
    "segment" = names(cnv)[S4Vectors::queryHits(ov.end)],
    "segment.chr" = seqnames(cnv)[S4Vectors::queryHits(ov.end)],
    "segment.start" = start(cnv)[S4Vectors::queryHits(ov.end)],
    "segment.end" = end(cnv)[S4Vectors::queryHits(ov.end)],
    "segment.width" = width(cnv)[S4Vectors::queryHits(ov.end)],
    "segment.type" = S4Vectors::mcols(cnv)[S4Vectors::queryHits(ov.end), "cnv.annotation"],
    "SV3" = names(gr)[S4Vectors::subjectHits(ov.end)],
    "SV3.chr" = seqnames(gr[S4Vectors::subjectHits(ov.end)]),
    "SV3.start" = start(gr[S4Vectors::subjectHits(ov.end)]),
    "SV3.end" = end(gr[S4Vectors::subjectHits(ov.end)]),
    "SV3.len" = S4Vectors::mcols(gr)[S4Vectors::subjectHits(ov.end), "svLen"],
    "SV3.type" = S4Vectors::mcols(gr)[S4Vectors::subjectHits(ov.end), "svtype"],
    "SV3.partner" = gr[S4Vectors::subjectHits(ov.end)]$partner
  ) %>%
    dplyr::mutate(
      SV3.partner.chr = seqnames(gr[SV3.partner]) %>% as.vector(),
      SV3.partner.start = start(gr[SV3.partner]),
      SV3.partner.end = end(gr[SV3.partner])
    )
  #
  # svs.end=mergeByOverlaps( query=cnv, subject=gr)
  # svs.end=svs.end[,c('cnv','cnv.annotation', 'cnv.cf','gr', 'sourceId','partner','svtype' )]
  # segments where the start and end have a hit to SV and the other end of the segment to the parnter


  rr <- dplyr::inner_join(svs.start, svs.end, by = c("segment" = "segment"), multiple = "all", suffix = c(".start", ".end"), relationship = "many-to-many") %>%
    dplyr::select(segment,
      segment.chr = segment.chr.start,
      segment.start = segment.start.start,
      segment.end = segment.end.start,
      segment.width = segment.width.start,
      segment.type = segment.type.start,
      SV5,
      SV5.chr,
      SV5.start,
      SV5.end,
      SV5.type,
      SV5.partner,
      SV5.partner.chr,
      SV5.partner.start,
      SV5.partner.end,
      SV5.len,
      SV3,
      SV3.chr,
      SV3.start,
      SV3.end,
      SV3.type,
      SV3.partner,
      SV3.partner.chr,
      SV3.partner.start,
      SV3.partner.end,
      SV3.len
    ) %>%
    dplyr::mutate(
      sv.support = TRUE,
      sv.support.chr = ifelse(SV3.chr == SV5.chr, TRUE, FALSE),
      sv.support.partner = ifelse(SV3.partner == SV5, TRUE, FALSE),
      sv.support.class = TRUE,
      sv.support.class = ifelse(segment.type == "AMP" & (SV5.type == "DEL" | SV3.type == "DEL"), FALSE, TRUE),
      sv.support.class = ifelse(segment.type == "HOMDEL" & (SV5.type == "DUP" | SV3.type == "DUP"), FALSE, TRUE),
      sv.support.class = ifelse(segment.type == "HETDEL" & (SV5.type == "DUP" | SV3.type == "DUP"), FALSE, TRUE)
    ) %>%
    dplyr::filter(sv.support == TRUE, sv.support.class == TRUE)
  # dplyr::inner_join( svs.start, svs.end, by=dplyr::join_by( segment.start==segment.end , SV.start==SV.end.partner, SV.start.partner==SV.end)  , multiple='all')

  S4Vectors::mcols(cnvsegments@segments)[, c("sv.support", "sv.support.chr", "sv.support.partner", "sv.support.class")] <- FALSE

  S4Vectors::mcols(cnvsegments@segments)[as.numeric(rr$segment), c("sv.support", "sv.support.chr", "sv.support.partner", "sv.support.class")] <- rr[, c("sv.support", "sv.support.chr", "sv.support.partner", "sv.support.class")]
  return(cnvsegments)
}
