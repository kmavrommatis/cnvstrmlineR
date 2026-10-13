#' Parse a VCF file in BND format with CNV data
#'
#' This function reads and parses a VCF file containing copy number variation (CNV) data
#' in BND (breakend) format. Such files are typically created using the \code{\link{saveVCF}}
#' function from the cnvstrmliner package. The function extracts segment information,
#' copy number states, purity, ploidy, and quality metrics from the VCF file and returns

#' a cnvstrml object.
#'
#' @param vcf A character string specifying the path to the input VCF file.
#'   The file can be compressed (`.vcf.gz` or `.vcf.bgz`) or uncompressed (`.vcf`).
#' @param sample A character string specifying the sample name. If \code{NA} (default),
#'   the sample name is first extracted from the VCF header. If not available in the header,
#'   the sample name is derived from the input filename by removing the `.vcf`, `.vcf.gz`,
#'   or `.vcf.bgz` extension.
#' @param genome A character string specifying the genome assembly version.
#'   Default is \code{"hg38"}. Common values include \code{"hg19"} and \code{"hg38"}.
#'
#' @return A \code{\link{cnvstrml}} object containing:
#' \describe{
#'   \item{segments}{A \code{GRanges} object with CNV segments and associated metadata
#'     including total, major, and minor allele copy numbers, log-ratio values,
#'     cellular fraction, cancer cell fraction, CNV annotations, and LOH status.}
#'   \item{purity}{Numeric value representing tumor purity extracted from VCF header.}
#'   \item{ploidy}{Numeric value representing tumor ploidy extracted from VCF header.}
#'   \item{effective_ploidy}{Calculated effective ploidy based on segment data.}
#'   \item{method}{The CNV calling method used, extracted from VCF header.}
#'   \item{sample}{The sample identifier.}
#'   \item{quality}{A list of quality metrics extracted from VCF header fields
#'     prefixed with "quality".}
#'   \item{genome}{The genome assembly version.}
#' }
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Reads the VCF file using \code{VariantAnnotation::readVcf}.
#'   \item Extracts metadata from the VCF header including purity, ploidy, sex,
#'     CNV method, and quality metrics.
#'   \item Filters segments to retain only those with FILTER status "neutral" or "PASS".
#'   \item Extracts copy number information from the INFO fields including:
#'     \itemize{
#'       \item TOTAL_CN: Total copy number
#'       \item MAJOR_CN: Major allele copy number
#'       \item MINOR_CN: Minor allele copy number
#'       \item NUM_MARK: Total number of markers
#'       \item NUM_HET_MARK: Number of heterozygous markers
#'       \item CN_LR_MEDIAN: Median log-ratio
#'       \item CF: Cellular fraction
#'       \item CCF: Cancer cell fraction
#'       \item CNV_ANN: CNV annotation
#'       \item LOH: Loss of heterozygosity status
#'     }
#'   \item Adjusts segment coordinates based on SVLEN and END fields.
#'   \item Creates and returns a cnvstrml object.
#' }
#'
#' @seealso
#' \code{\link{saveVCF}} for creating VCF files from cnvstrml objects.
#' \code{\link{cnvstrml}} for the cnvstrml class definition.
#' \code{\link{parseFACETS}}, \code{\link{parsePureCN}}, \code{\link{parseSequenza}}
#' for parsing output from specific CNV callers.
#'
#' @examples
#' \dontrun{
#' # Parse a VCF file with default parameters
#' cnv_obj <- parseVCF("sample_cnv.vcf.gz")
#'
#' # Parse with explicit sample name and genome version
#' cnv_obj <- parseVCF(
#'   vcf = "sample_cnv.vcf.gz",
#'   sample = "SAMPLE_001",
#'   genome = "hg38"
#' )
#'
#' # Access the parsed data
#' segments <- cnv_obj@segments
#' purity <- cnv_obj@purity
#' ploidy <- cnv_obj@ploidy
#' }
#'
#' @importFrom VariantAnnotation readVcf meta header samples info
#' @importFrom SummarizedExperiment rowRanges
#' @importFrom S4Vectors mcols
#' @importFrom GenomicRanges end start width
#' @importFrom tibble as_tibble
#' @importFrom tidyr replace_na
#' @importFrom dplyr mutate select
#' @importFrom purrr map_chr map_lgl
#' @importFrom stringi stri_replace_all_fixed
#' @importFrom logger log_debug
#'
#' @export
parseVCF <- function(
  vcf,
  sample = NA_character_,
  genome = "hg38"
) {
  logger::log_debug("Loading VCF file {vcf}", namespace = "parser")
  ff.vcf <- VariantAnnotation::readVcf(vcf, genome = genome)
  logger::log_debug("Parsing metadata from VCF header", namespace = "parser")
  meta <- VariantAnnotation::meta(VariantAnnotation::header(ff.vcf))
  if ("PURITY" %in% names(meta)) {
    purity <- meta$PURITY[1, 1]
  } else if ("purity" %in% names(meta)) {
    purity <- meta$purity[1, 1]
  } else {
    purity <- NA
  }
  purity <- as.numeric(purity)

  if ("PLOIDY" %in% names(meta)) {
    ploidy <- meta$PLOIDY[1, 1]
  } else if ("ploidy" %in% names(meta)) {
    ploidy <- meta$ploidy[1, 1]
  } else {
    ploidy <- NA
  }
  ploidy <- as.numeric(ploidy)

  sex <- meta$sex[1, 1] %>% as.character()
  cnvmethod <- meta$cnvmethod[, 1]

  quality_fields <- names(meta)[grepl("^quality", names(meta), ignore.case = FALSE)]


  # Convert ranges to a dataframe
  logger::log_debug("Processing VCF ranges and filtering segments", namespace = "parser")
  ff.gr <- SummarizedExperiment::rowRanges(ff.vcf)
  keep <- which(ff.gr$FILTER %in% c("neutral", "PASS"))
  ff.vcf <- ff.vcf[keep]
  ff.gr <- SummarizedExperiment::rowRanges(ff.vcf)

  logger::log_debug("Adding copy number and allelic information to segments", namespace = "parser")
  S4Vectors::mcols(ff.gr) <- cbind(S4Vectors::mcols(ff.gr), VariantAnnotation::info(ff.vcf))
  mcols_new <- S4Vectors::mcols(ff.gr) %>%
    tibble::as_tibble() %>%
    tidyr::replace_na(list(MINOR_CN = 0, CF = 0, CCF = 1)) %>%
    dplyr::mutate(
      total.allele.copies = TOTAL_CN,
      major.allele.copies = MAJOR_CN,
      minor.allele.copies = MINOR_CN,
      cnv.totalmarkers = NUM_MARK,
      cnv.markers = NUM_HET_MARK,
      cnv.lrr = CN_LR_MEDIAN,
      cnv.cf = CF,
      cnv.ccf = CCF,
      cnv.annotation = purrr::map_chr(CNV_ANN, 1),
      cnv.loh = purrr::map_lgl(LOH, ~ .x[1] == 1)
    ) %>%
    dplyr::select(
      total.allele.copies,
      major.allele.copies,
      minor.allele.copies,
      cnv.totalmarkers,
      cnv.markers,
      cnv.lrr,
      cnv.cf,
      cnv.ccf,
      cnv.annotation,
      cnv.loh
    )

  GenomicRanges::end(ff.gr) <- ff.gr$END

  idx <- which((ff.gr$SVLEN == GenomicRanges::width(ff.gr)) == FALSE)
  if (length(idx) > 0) {
    GenomicRanges::start(ff.gr)[idx] <- ff.gr$END[idx] - ff.gr$SVLEN[idx] + 1
  }
  S4Vectors::mcols(ff.gr) <- mcols_new
  if (is.na(sample) || length(sample) == 0) {
    sample <- VariantAnnotation::samples(VariantAnnotation::header(ff.vcf))
    logger::log_debug("Extracting sample name from VCF header: {sample}", namespace = "parser")
    # message("Extracting sample name from the header of the VCF file :" , sample)
  }
  if (is.na(sample) || length(sample) == 0) {
    sample <- vcf %>%
      basename() %>%
      stringi::stri_replace_all_fixed(., c(".vcf.gz", ".vcf.bgz", ".vcf"), replacement = "", vectorize_all = FALSE)
    logger::log_debug("Extracting sample name from filename: {sample}", namespace = "parser")
    # message("Extracting sample name from the filename :" , sample)
  }

  qual.df <- list()
  for (field in quality_fields) {
    # Remove 'quality' prefix from the field name
    clean_name <- sub("^quality", "", field)
    qual.df[[clean_name]] <- meta[[field]][1, 1]
  }
  logger::log_debug("Creating cnvstrml object for sample {sample} with purity {purity} and ploidy {ploidy}", namespace = "parser")
  obj <- cnvstrml(
    segments = ff.gr,
    purity = as.numeric(purity),
    ploidy = as.numeric(ploidy),
    effective_ploidy = calculateEffectivePloidy(ff.gr),
    method = cnvmethod,
    sample = sample,
    quality = qual.df,
    genome = genome
  )
  obj
}
