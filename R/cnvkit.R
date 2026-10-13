#' Parse CNVkit call output (.cns) into a cnvstrml object
#'
#' @description
#' High-level parser for CNVkit "call" output files (typically with extension \code{.cns}).
#' It reads segments, harmonizes column names, drops heavy gene annotations, and constructs a
#' \code{cnvstrml} object with standardized fields. If allele-specific columns are missing,
#' the function infers them from total copy number and records this decision in the object quality.
#'
#' @param cnvkit.filename Character path to a CNVkit \code{.cns} file produced by \code{cnvkit call}.
#' @param sample Character sample identifier to store in the resulting object.
#' @param genome Character genome build (e.g., \code{"hg38"}) used to set \code{Seqinfo} and UCSC \code{seqlevels}.
#'
#' @return A \code{cnvstrml} object whose \code{@segments} is a \code{GRanges} with standardized metadata:
#' \describe{
#'   \item{total.allele.copies}{From CNVkit \code{cn}.}
#'   \item{major.allele.copies}{\code{max(cn1, cn2)} or inferred.}
#'   \item{minor.allele.copies}{\code{min(cn1, cn2)} or inferred.}
#'   \item{cnv.markers}{From \code{probes}.}
#'   \item{cnv.totalmarkers}{Same as \code{probes}.}
#'   \item{cnv.cf}{Purity proxy set to 1 (CNVkit does not provide purity).}
#'   \item{cnv.ccf}{Clonal cancer cell fraction set to 1 (all clones).}
#'   \item{cnv.lrr}{Log-ratio from \code{log2}.}
#'   \item{cnv.annotation}{Status assigned via \code{annotateCNV()} using effective ploidy.}
#' }
#'
#' @details
#' \strong{Input columns and flexible mappings}
#' \itemize{
#'   \item Chromosome: \code{chrom} or \code{chromosome}.
#'   \item Start: \code{svstart} or \code{start} (1-based).
#'   \item End: \code{svend} or \code{end} (inclusive).
#'   \item Required CNVkit fields: \code{cn} (total copies); optionally \code{cn1}, \code{cn2}.
#'   \item Additional fields: \code{probes} (markers), \code{log2} (segment log-ratio).
#' }
#'
#' \strong{Allele inference when \code{cn1}/\code{cn2} are missing}
#' If \code{cn1} and \code{cn2} are absent, alleles are inferred as
#' \code{cn1 = ifelse(cn > 1, cn - 1, cn)} and \code{cn2 = ifelse(cn > 1, 1, 0)}.
#' The object's \code{quality$parser} records that allele copies were inferred and purity set to 1.
#'
#' \strong{Genome harmonization and annotation}
#' The parser converts \code{seqlevelsStyle} to UCSC, prunes to the requested \code{genome}'s
#' \code{Seqinfo}, estimates effective ploidy via \code{calculateEffectivePloidy()}, and assigns
#' \code{cnv.annotation} with \code{annotateCNV(total.cn, minor.cn, ploidy)}.
#'
#' \strong{Performance consideration}
#' The CNVkit \code{gene} column is not read to reduce memory footprint; only required fields are loaded.
#'
#' @examples
#' \dontrun{
#' cnv <- parseCNVkit(
#'   cnvkit.filename = "sample.cns",
#'   sample = "S1",
#'   genome = "hg38"
#' )
#'
#' # Inspect segments and annotations
#' cnv@segments
#' table(S4Vectors::mcols(cnv@segments)$cnv.annotation)
#' }
#'
#' @seealso \code{\link{parseCNV}}, \code{\link{annotateCNV}}, \code{\link{calculateEffectivePloidy}}
#'
#' @export
parseCNVkit <- function(
  cnvkit.filename,
  sample = NA_character_,
  genome = "hg38"
) {
  purity <- 1

  segments <- suppressWarnings(data.table::fread(
    file = cnvkit.filename,
    sep = "\t",
    header = TRUE,
    drop = "gene"
  ))


  if ("chrom" %in% colnames(segments)) {
    segments <- segments %>%
      dplyr::rename(chromosome = chrom)
  }
  if ("svstart" %in% colnames(segments)) {
    segments <- segments %>%
      dplyr::rename(start = svstart)
  }
  if ("svend" %in% colnames(segments)) {
    segments <- segments %>%
      dplyr::rename(end = svend)
  }

  full_format <- TRUE
  expected_cols <- c("chromosome", "start", "end", "cn", "cn1", "cn2")
  if (!setequal(intersect(expected_cols, colnames(segments)), expected_cols)) {
    message <- paste("This input file ", cnvkit.filename, " does not include allele copies. We expected columns cn, cn1, cn2. We will arbitrarily set them to 1/1")
    full_format <- FALSE
    expected_cols <- c("chromosome", "start", "end", "cn")
    if (!setequal(intersect(expected_cols, colnames(segments)), expected_cols)) {
      logger::log_error("This input file ", cnvkit.filename, " is not the proper format. We expect at least the column cn.")
      stop()
    } else {
      logger::log_warn(message)
    }
  }

  seqinfo <- getSeqInfo(genome)
  segments <- segments %>%
    GenomicRanges::makeGRangesFromDataFrame(.,
      keep.extra.columns = TRUE,
      ignore.strand = TRUE
    )
  GenomeInfoDb::seqlevelsStyle(segments) <- "UCSC"
  GenomeInfoDb::seqlevels(segments, pruning.mode = "coarse") <- GenomeInfoDb::seqlevels(seqinfo)
  GenomeInfoDb::seqinfo(segments, pruning.mode = "coarse") <- seqinfo

  mcols_new <- S4Vectors::mcols(segments) %>% as.data.frame()

  if (!setequal(intersect(c("cn1", "cn2"), colnames(mcols_new)), c("cn1", "cn2"))) {
    mcols_new <- mcols_new %>%
      dplyr::mutate(
        cn1 = ifelse(cn > 1, cn - 1, cn),
        cn2 = ifelse(cn > 1, 1, 0)
      )
  }


  suppressWarnings({
    S4Vectors::mcols(segments) <- mcols_new %>%
      dplyr::mutate(
        total.allele.copies = cn,
        major.allele.copies = pmax(cn1, cn2),
        minor.allele.copies = pmin(cn1, cn2),
        cnv.markers = probes,
        cnv.totalmarkers = probes,
        cnv.cf = as.numeric(purity),
        cnv.ccf = 1, # since this is the totality of the clones
        # cnv.annotation = NA,
        # cnv.loh=loh,
        cnv.lrr = log2
      )
  })
  # effective_ploidy <- cnvstrmliner:::calculateEffectivePloidy(segments)
  # flags <- apply(mcols_new, 1, function(K) {
  #   cnvstrmliner::annotateCNV(
  #     total.cn = K["cn"],
  #     minor.cn = K["cn2"],
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
    method = "CNVkit call",
    sample = sample,
    genome = genome
  )
  obj@quality$file <- c("segments" = cnvkit.filename)
  obj
}
