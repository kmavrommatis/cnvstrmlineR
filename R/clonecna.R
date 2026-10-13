#' Parse CloneCNA copy number output
#'
#' @description
#' Parser for CloneCNA (version 2.0) output files. CloneCNA is a tool for detecting somatic copy
#' number alterations (CNAs) and inferring clonal structure from sequencing data. It estimates
#' tumor purity, ploidy, and identifies subclonal populations with their cellularity.
#' This function extracts metadata from the file header and segment-level copy number information,
#' constructing a \code{cnvstrml} object with standardized fields and CNV annotations.
#'
#' @param clonecna.file Character path to a CloneCNA output file. The file contains a header section
#'   (typically the first ~10 lines) with summary statistics, followed by a tab-delimited table
#'   of copy number segments with columns: \code{Chr}, \code{StartPos}, \code{EndPos}, \code{CN},
#'   \code{mCN}, \code{Aberration}, \code{Cellularity}, \code{Score}.
#' @param genome Character genome build (e.g., \code{"hg38"}, \code{"hg19"}) used to set
#'   \code{Seqinfo} and UCSC seqlevels style. Default: \code{"hg38"}.
#' @param methodName Character method identifier to store in the \code{cnvstrml} object.
#'   Default: \code{"cloneCNA"}.
#' @param sample Character sample identifier. If \code{NA} (default), the sample name is
#'   extracted from the filename by removing file extensions and CloneCNA-specific suffixes.
#'
#' @return A \code{cnvstrml} object with:
#' \describe{
#'   \item{@segments}{A \code{GRanges} with metadata columns:
#'     \describe{
#'       \item{Chr, StartPos, EndPos}{Original chromosome coordinates (preserved).}
#'       \item{CN}{Total copy number from CloneCNA.}
#'       \item{mCN}{Minor allele copy number from CloneCNA.}
#'       \item{Aberration}{Aberration type from CloneCNA (e.g., AHET, NHET, NLOH, ALOH).}
#'       \item{Cellularity}{Clone-specific cellularity for this segment.}
#'       \item{Score}{CloneCNA quality score for the segment.}
#'       \item{total.allele.copies}{Standardized total CN (equal to \code{CN}).}
#'       \item{major.allele.copies}{Calculated as \code{CN - mCN}.}
#'       \item{minor.allele.copies}{Standardized minor CN (equal to \code{mCN}).}
#'       \item{cnv.cf}{Segment-specific cellularity (equal to \code{Cellularity}).}
#'       \item{cnv.ccf}{Cancer cell fraction, calculated as \code{Cellularity / purity}.}
#'       \item{cnv.markers}{Set to \code{NA_integer_} (not provided by CloneCNA).}
#'       \item{cnv.totalmarkers}{Set to \code{NA_integer_} (not provided by CloneCNA).}
#'       \item{cnv.lrr}{Set to \code{NA_real_} (not provided by CloneCNA output).}
#'     }}
#'   \item{@purity}{Tumor purity extracted from "Proportion of abnormal cells in the sample"
#'     in the file header.}
#'   \item{@ploidy}{Average copy number (ploidy) extracted from "Average copy number" in the header.}
#'   \item{@effective_ploidy}{Calculated ploidy via \code{calculateEffectivePloidy()}.}
#'   \item{@sample}{Sample identifier (from filename or user-supplied).}
#'   \item{@method}{Method identifier (default: \code{"cloneCNA"}).}
#'   \item{@quality}{List containing CloneCNA-specific metadata:
#'     \describe{
#'       \item{lcr_baseline_shift}{LCR baseline shift parameter.}
#'       \item{sigma_lcr}{Sigma of t-distribution for LCR data.}
#'       \item{sigma_baf}{Sigma of t-distribution for BAF data.}
#'       \item{df_lcr}{Degrees of freedom for LCR t-distribution.}
#'       \item{df_baf}{Degrees of freedom for BAF t-distribution.}
#'       \item{file}{Full path to the CloneCNA input file.}
#'     }}
#' }
#'
#' @details
#' \strong{File format}
#'
#' CloneCNA output files have a two-part structure:
#' \enumerate{
#'   \item \strong{Header section}: Contains summary statistics and model parameters, including:
#'     \itemize{
#'       \item Proportion of abnormal cells (purity)
#'       \item Average copy number (ploidy)
#'       \item LCR baseline shift
#'       \item Sigma and degrees of freedom for t-distributions
#'       \item Cellularity of subclonal populations
#'     }
#'   \item \strong{Segment table}: Tab-delimited data starting after a line beginning with
#'     \code{"Chr  StartPos  EndPos"}. Each row describes one genomic segment with total and
#'     minor allele copy numbers, aberration type, clone-specific cellularity, and quality score.
#' }
#'
#' \strong{Chromosome naming}
#'
#' CloneCNA typically outputs chromosomes without the "chr" prefix (e.g., "1", "2", "X").
#' This parser automatically prepends "chr" to ensure UCSC-style naming (e.g., "chr1", "chr2").
#'
#' \strong{Subclonal CNAs}
#'
#' Unlike some callers that report a single clonal cellularity, CloneCNA assigns segment-specific
#' cellularity values, allowing inference of subclonal CNAs. The parser stores these in both
#' \code{cnv.cf} (raw cellularity) and \code{cnv.ccf} (cancer cell fraction adjusted by purity).
#'
#' \strong{Missing data}
#'
#' CloneCNA output does not include per-segment marker counts or log-ratio values. These fields
#' (\code{cnv.markers}, \code{cnv.totalmarkers}, \code{cnv.lrr}) are set to \code{NA}.
#'
#' @examples
#' \dontrun{
#' # Parse a CloneCNA output file
#' cnv_obj <- parseCloneCNA(
#'   clonecna.file = "sample_clonecna_output.txt",
#'   genome = "hg38",
#'   sample = "SAMPLE_001"
#' )
#'
#' # View purity and ploidy
#' cnv_obj@purity
#' cnv_obj@ploidy
#'
#' # Access segments
#' segments <- cnv_obj@segments
#'
#' # View quality metrics
#' cnv_obj@quality
#' }
#'
#' @seealso
#' \code{\link{parseFACETS}}, \code{\link{parsePureCN}}, \code{\link{parseSCLUST}}
#' for other CNV callers; \code{\link{cnvstrml-class}} for object structure;
#' \code{\link{annotateCNV}} for CNV annotation; \code{\link{calculateEffectivePloidy}}
#' for ploidy calculation.
#'
#' @export
parseCloneCNA <- function(clonecna.file, genome = "hg38", methodName = "cloneCNA", sample = NA) {
  logger::log_debug("Parsing CloneCNA file {clonecna.file}", namespace = "parser")
  # get the proper stem for these files
  clonecna.name <- predictName(basename(clonecna.file))
  # Read the file
  clonecna.file <- normalizePath(clonecna.file)
  if (!file.exists(clonecna.file)) {
    logger::log_error("Filename {clonecna.file} does not exist. Aborting !!!", namespace = "parser")
    stop("Filename ", clonecna.file, " does not exist. Aborting !!!")
  }

  # Read all lines
  lines <- readLines(clonecna.file)
  logger::log_debug("Read {length(lines)} lines from file", namespace = "parser")

  # Parse header information
  logger::log_debug("Parsing header metadata", namespace = "parser")
  purity <- NA_real_
  avg_cn <- NA_real_
  lcr_shift <- NA_real_
  sigma_lcr <- NA_real_
  sigma_baf <- NA_real_
  df_lcr <- NA_real_
  df_baf <- NA_real_
  cellularity <- NA_character_

  for (line in lines) {
    if (grepl("Proportion of abnormal cells in the sample:", line)) {
      purity <- as.numeric(sub(".*:\\s*([0-9.]+).*", "\\1", line))
      logger::log_debug("Extracted purity: {purity}", namespace = "parser")
    } else if (grepl("Average copy number:", line)) {
      avg_cn <- as.numeric(sub(".*:\\s*([0-9.]+).*", "\\1", line))
      logger::log_debug("Extracted average copy number (ploidy): {avg_cn}", namespace = "parser")
    } else if (grepl("LCR baseline shift:", line)) {
      lcr_shift <- as.numeric(sub(".*:\\s*([-0-9.]+).*", "\\1", line))
    } else if (grepl("Sigma of t-distributions for LCR and BAF data:", line)) {
      sigmas <- stringr::str_extract_all(line, "[0-9.]+")[[1]]
      if (length(sigmas) >= 2) {
        sigma_lcr <- as.numeric(sigmas[1])
        sigma_baf <- as.numeric(sigmas[2])
      }
    } else if (grepl("Degree of freedom of t-distributions for LCR and BAF data:", line)) {
      dfs <- stringr::str_extract_all(line, "[0-9]+")[[1]]
      if (length(dfs) >= 2) {
        df_lcr <- as.numeric(dfs[1])
        df_baf <- as.numeric(dfs[2])
      }
    } else if (grepl("Cellularity of subclonal populations in the sample:", line)) {
      cellularity <- sub(".*:\\s*(.*)$", "\\1", line)
    }
  }

  # Find where the segment data starts (after the header line with Chr, StartPos, etc)
  header_line <- which(grepl("^Chr\\s+StartPos\\s+EndPos", lines))
  if (length(header_line) == 0) {
    logger::log_error("Could not find segment header line in file", namespace = "parser")
    stop("Could not find segment header line starting with 'Chr StartPos EndPos' in file")
  }

  logger::log_debug("Found segment data starting at line {header_line + 1}", namespace = "parser")

  # Read the segment data as a table
  segments_df <- data.table::fread(
    file = clonecna.file,
    skip = header_line - 1,
    header = TRUE
  )

  logger::log_debug("Parsing {nrow(segments_df)} segments", namespace = "parser")

  # Get seqinfo for the genome
  seqinfo <- getSeqInfo(genome)

  # Prepare dataframe for GRanges conversion
  segments_df <- segments_df %>%
    dplyr::mutate(
      seqnames = paste0("chr", Chr),
      start = StartPos,
      end = EndPos
    )

  # Convert to GRanges
  segments.gr <- GenomicRanges::makeGRangesFromDataFrame(
    segments_df,
    seqnames.field = "seqnames",
    start.field = "start",
    end.field = "end",
    keep.extra.columns = TRUE,
    ignore.strand = TRUE,
    seqinfo = seqinfo
  ) %>%
    GenomeInfoDb::sortSeqlevels() %>%
    BiocGenerics::sort()

  # GenomeInfoDb::seqlevelsStyle(segments.gr) <- "UCSC"

  logger::log_debug("Adding copy number annotations", namespace = "parser")

  # Add required columns for cnvstrml object
  mcols_new <- S4Vectors::mcols(segments.gr) %>%
    tibble::as_tibble() %>%
    dplyr::mutate(
      total.allele.copies = CN,
      major.allele.copies = mCN,
      minor.allele.copies = CN - mCN,
      cnv.cf = ifelse(Cellularity > 0, Cellularity, 1),
      cnv.ccf = cnv.cf / purity,
      cnv.markers = NA_integer_, # not provided by CloneCNA
      cnv.totalmarkers = NA_integer_, # not provided by CloneCNA
      cnv.lrr = NA_real_ # not provided by CloneCNA
    )

  S4Vectors::mcols(segments.gr) <- mcols_new

  # Handle sample name
  if (is.na(sample)) {
    sample <- clonecna.name
    logger::log_debug("Extracted sample name from filename: {sample}", namespace = "parser")
  }

  # Create quality list with metadata
  quality_list <- list(
    "lcr_baseline_shift" = lcr_shift,
    "sigma_lcr" = sigma_lcr,
    "sigma_baf" = sigma_baf,
    "df_lcr" = df_lcr,
    "df_baf" = df_baf,
    "file" = clonecna.file
  )

  logger::log_debug("Creating cnvstrml object for sample {sample} with purity {purity} and ploidy {avg_cn}", namespace = "parser")

  # Create cnvstrml object
  obj <- cnvstrml(
    segments = segments.gr,
    purity = as.numeric(purity),
    ploidy = as.numeric(avg_cn),
    effective_ploidy = calculateEffectivePloidy(segments.gr),
    method = methodName,
    sample = sample,
    quality = quality_list,
    genome = genome
  )

  logger::log_info("Successfully parsed CloneCNA file with {length(segments.gr)} segments", namespace = "parser")
  obj
}
