#' Parse GATK CopyNumber segment output
#'
#' @description
#' Parser for GATK CopyNumber (\code{gatk-tools ModelSegments}) output files, typically with
#' extension \code{.seg}. GATK CopyNumber performs segmentation of copy-ratio and allele-fraction
#' profiles to identify regions of somatic copy-number alteration. This parser reads tab-delimited
#' GATK segments, estimates copy numbers from log-ratios and CALL values using pragmatic cutoffs,
#' and constructs a \code{cnvstrml} object. Note: GATK does not provide allele-specific copy numbers
#' or LOH calls; minor alleles are estimated heuristically.
#'
#' @param sample.file Character path to a GATK CopyNumber output file (typically \code{.seg}).
#'   Must be a tab-delimited text file with standard GATK column names.
#' @param sample Character sample identifier. If \code{NA} (default), the filename (without extension)
#'   is used as the sample name.
#' @param genome Character genome build (e.g., \code{"hg38"}) used to set \code{Seqinfo} and UCSC seqlevels.
#' @param ploidy Numeric baseline ploidy for the sample (default: 2). Used as the baseline for
#'   inferring total copy numbers from log-ratio cutoffs.
#'
#' @return A \code{cnvstrml} object with:
#' \describe{
#'   \item{@segments}{A \code{GRanges} with metadata:
#'     \describe{
#'       \item{CONTIG}{Chromosome/contig name.}
#'       \item{START, END}{Segment start and end coordinates.}
#'       \item{NUM_POINTS_COPY_RATIO}{Number of log-ratio data points in segment.}
#'       \item{MEAN_LOG2_COPY_RATIO}{Mean log2 copy-ratio.}
#'       \item{CALL}{GATK call (0, +, -; 0=neutral, +=gain, -=loss).}
#'       \item{total.allele.copies}{Inferred total copy number (see Details).}
#'       \item{minor.allele.copies}{\code{NA}; not provided by GATK. Set to \code{floor(total/2)}.}
#'       \item{major.allele.copies}{\code{NA}; not provided by GATK.}
#'       \item{cnv.markers}{Equal to NUM_POINTS_COPY_RATIO.}
#'       \item{cnv.lrr}{Log2 copy-ratio.}
#'       \item{cnv.cf}{Cancer cell fraction set to 1 (purity not provided).}
#'       \item{cnv.annotation}{CNV status (AMP, DEL, NEUT) assigned via \code{annotateCNV()}.}
#'       \item{cnv.loh}{\code{NA}; not provided by GATK.}
#'     }}
#'   \item{@sample}{Sample name (from file or user-supplied).}
#'   \item{@method}{\code{"GATK CopyNumber"}.}
#'   \item{@genome}{Genome build.}
#'   \item{@sex}{\code{NA}; not inferred by this parser.}
#' }
#'
#' @details
#' \strong{Input file format}
#' GATK CopyNumber output is tab-delimited with header row. Key columns:
#' \itemize{
#'   \item \code{CONTIG}: Chromosome/contig name (e.g., chr1, 1, chrX).
#'   \item \code{START}, \code{END}: Segment boundaries (1-based).
#'   \item \code{NUM_POINTS_COPY_RATIO}: Count of copy-ratio probes/bins in segment.
#'   \item \code{MEAN_LOG2_COPY_RATIO}: Mean log2(sample/diploid-reference) per segment.
#'   \item \code{CALL}: Categorical call from GATK (0=neutral, +=gain, -=loss).
#' }
#' Comments (lines starting with \code{@}) are skipped.
#'
#' \strong{Copy number inference}
#' GATK provides log-ratio and categorical call, but not allele-specific copy numbers.
#' This parser estimates total copy number using log-ratio cutoffs and the supplied \code{ploidy}:
#' \itemize{
#'   \item \code{CALL == "0"}: Total copies = \code{ploidy} (neutral).
#'   \item \code{MEAN_LOG2_COPY_RATIO < -2} and \code{CALL == "-"}: Total copies = 0 (homozygous deletion).
#'   \item \code{-2 <= MEAN_LOG2_COPY_RATIO < 0} and \code{CALL == "-"}: Total copies = \code{ploidy - 1}.
#'   \item \code{0 < MEAN_LOG2_COPY_RATIO <= 1} and \code{CALL == "+"}: Total copies = \code{ploidy + 1}.
#'   \item \code{MEAN_LOG2_COPY_RATIO > 1} and \code{CALL == "+"}: Total copies = \code{ploidy + 2}.
#' }
#' These are \emph{indicative} estimates; actual values depend on sample-specific callability and normalization.
#'
#' \strong{Allele-specific copies and LOH}
#' \itemize{
#'   \item GATK does not provide minor/major allele copy numbers.
#'   \item Minor alleles are set to \code{floor(total / 2)} as a placeholder.
#'   \item LOH information is not provided; \code{cnv.loh} set to \code{NA}.
#' }
#'
#' \strong{Limitations and considerations}
#' \itemize{
#'   \item No purity information; \code{cnv.cf} set to 1 (assumed clonal, pure sample).
#'   \item No subclonal or clonal fraction data.
#'   \item Copy-number cutoffs are heuristic; calibration may differ across cohorts and assays.
#'   \item Sex determination not implemented; use \code{getSex()} post-parsing if needed.
#' }
#'
#' @examples
#' \dontrun{
#' # Parse GATK CopyNumber output
#' cnv <- parseGATK(
#'   sample.file = "sample.seg",
#'   sample = "S1",
#'   genome = "hg38",
#'   ploidy = 2
#' )
#'
#' # Inspect segments and calls
#' head(cnv@segments)
#' table(S4Vectors::mcols(cnv@segments)$CALL)
#' table(S4Vectors::mcols(cnv@segments)$cnv.annotation)
#' }
#'
#' @seealso \code{\link{parseCNV}}, \code{\link{annotateCNV}}, \code{\link{getSex}}
#'
#' @export
parseGATK <- function(
  sample.file,
  sample = NA_character_,
  genome = "hg38",
  ploidy = 2
) {
  segments_file <- sample.file
  segments_file <- normalizePath(segments_file)
  if (!file.exists(segments_file)) {
    stop("Filename ", segments_file, " does not exist. Aborting !!!")
  }
  # seqlevels(sclust_subclones.gr)=GenomeInfoDb::seqlevels( seqinfo )
  seqinfo <- getSeqInfo(genome)
  df <- read.table(segments_file,
    comment.char = "@",
    header = TRUE,
    sep = "\t"
  )

  df <- df %>%
    GenomicRanges::makeGRangesFromDataFrame(.,
      seqnames.field = "CONTIG",
      start.field = "START",
      end.field = "END",
      keep.extra.columns = TRUE,
      ignore.strand = TRUE,
      seqinfo = seqinfo
    ) %>%
    GenomicRanges::trim(use.names = TRUE)

  mcols_new <- S4Vectors::mcols(df) %>%
    as.data.frame() %>%
    dplyr::mutate(
      CopyNr = ifelse(CALL == "0", ploidy,
        ifelse(MEAN_LOG2_COPY_RATIO < -2, 0,
          ifelse(MEAN_LOG2_COPY_RATIO < 0, ploidy - 1,
            ifelse(MEAN_LOG2_COPY_RATIO > 1, ploidy + 2,
              ifelse(MEAN_LOG2_COPY_RATIO > 0, ploidy + 1, ploidy)
            )
          )
        )
      ),
      minor = floor(CopyNr / 2)
    )


  # merge the information from the two objects
  # flags <- apply(mcols_new, 1, function(K) {
  # cnvstrmliner::annotateCNV(total.cn = K["CopyNr"], minor.cn = K["minor"])
  # })
  # loh=cnvstrmliner::annotateLOH( flags )

  S4Vectors::mcols(df) <- mcols_new %>%
    dplyr::mutate(
      total.allele.copies = CopyNr,
      minor.allele.copies = NA,
      major.allele.copies = NA,
      cnv.markers = NUM_POINTS_COPY_RATIO,
      cnv.cf = 1,
      cnv.annotation = flags,
      cnv.loh = NA,
      cnv.lrr = MEAN_LOG2_COPY_RATIO
    )

  if (is.na(sample)) {
    s <- basename(sample.file) %>% stringr::str_split(., pattern = "\\.", simplify = TRUE)
    sample <- as.character(s[1])
  }
  sex <- NA_character_
  obj <- cnvstrml(
    segments = df,
    sex = sex,
    method = "GATK CopyNumber",
    sample = sample,
    genome = genome
  )
  obj
}
