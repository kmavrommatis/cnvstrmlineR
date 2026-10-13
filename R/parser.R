#' Parse CNV data from various CNV calling methods
#'
#' @description
#' A unified dispatcher function that parses copy number variation (CNV) output files
#' from multiple CNV calling tools. This function acts as a central entry point that
#' delegates to the appropriate method-specific parser based on the `cnv.method` argument.
#'
#' @details
#' This function provides a convenient wrapper around all supported CNV parsing functions
#' in the cnvstrmliner package. Instead of calling method-specific parsers directly
#' (e.g., `parseFACETS()`, `parseSCLUST()`), users can use this unified interface and
#' specify the method name.
#'
#' \subsection{Supported Methods}{
#' The following CNV calling methods are currently supported (case-insensitive):
#' \describe{
#'   \item{\strong{FACETS}}{Fraction and Allele-Specific Copy Number Estimates from Tumor Sequencing.
#'     Accepts data frames or VCF files. Delegates to \code{\link{parseFACETS}}.}
#'   \item{\strong{PureCN}}{Copy number calling and SNV classification using targeted short
#'     read sequencing. Delegates to \code{\link{parsePureCN}}.}
#'   \item{\strong{SCLUST}}{Subclonal copy number and mutation analysis. Delegates to
#'     \code{\link{parseSCLUST}}.}
#'   \item{\strong{Battenberg}}{Subclonal copy number analysis from whole genome sequencing.
#'     Delegates to \code{\link{parseBattenberg}}.}
#'   \item{\strong{Sequenza}}{Copy number estimation from tumor-normal sequencing.
#'     Delegates to \code{\link{parseSequenza}}.}
#'   \item{\strong{CNVkit}}{Copy number from targeted DNA sequencing. Delegates to
#'     \code{\link{parseCNVkit}}.}
#'   \item{\strong{GATK}}{GATK CNV calling pipeline output. Delegates to \code{\link{parseGATK}}.}
#'   \item{\strong{FREEC / Control-FREEC}}{Control-FREEC copy number and allelic imbalance
#'     detection. Delegates to \code{\link{parseCfreec}}.}
#'   \item{\strong{Conseca}}{Consensus copy number analysis. Delegates to \code{\link{parseConseca}}.}
#'   \item{\strong{VCF}}{Generic VCF file with CNV information in the header and INFO fields.
#'     Expects purity, ploidy, and copy number annotations. Delegates to \code{\link{parseVCF}}.}
#'   \item{\strong{TSV}}{Tab-separated values file with columns: chromosome, start, end, cn1, cn2.
#'     A simple format for allele-specific copy number data. Delegates to \code{\link{parseTSV}}.}
#'   \item{\strong{SEG}}{Standard SEG format file (6 columns: sample, chrom, start, end, markers, value).
#'     Compatible with IGV and GISTIC input. Delegates to \code{\link{parseSEG}}.}
#' }
#' }
#'
#' \subsection{Method Name Matching}{
#' Method names are matched case-insensitively and support common variations:
#' \itemize{
#'   \item "FACETS", "facets", "Facets" → \code{parseFACETS()}
#'   \item "PureCN", "purecn", "PURECN" → \code{parsePureCN()}
#'   \item "SCLUST", "sclust", "Sclust" → \code{parseSCLUST()}
#'   \item "Battenberg", "battenberg", "BATTENBERG" → \code{parseBattenberg()}
#'   \item "Sequenza", "sequenza", "SEQUENZA" → \code{parseSequenza()}
#'   \item "CNVkit", "cnvkit", "CNVKIT" → \code{parseCNVkit()}
#'   \item "GATK", "gatk", "Gatk" → \code{parseGATK()}
#'   \item "FREEC", "freec", "controlfreec", "Control-FREEC", "cfreec" → \code{parseCfreec()}
#'   \item "Conseca", "conseca", "CONSECA" → \code{parseConseca()}
#'   \item "VCF", "vcf" → \code{parseVCF()}
#'   \item "TSV", "tsv" → \code{parseTSV()}
#'   \item "SEG", "seg" → \code{parseSEG()}
#' }
#' }
#'
#' @param cnv.data The input CNV data. This can be:
#'   \itemize{
#'     \item A file path (character string) pointing to the CNV output file(s)
#'     \item A data frame containing CNV segments (for methods that support it, e.g., FACETS)
#'   }
#'   The exact format depends on the CNV method being used. See the documentation
#'   for individual parser functions for details.
#' @param cnv.method Character string specifying the CNV calling method. Must be one of
#'   the supported methods listed in Details. Matching is case-insensitive.
#' @param sample Character string specifying the sample name. If \code{NA} (default),
#'   the parser will attempt to extract the sample name from the input file.
#' @param genome Character string specifying the reference genome build. Default is
#'   \code{"hg38"}. Common values include \code{"hg19"}, \code{"hg38"}, \code{"mm9"},
#'   \code{"mm10"}.
#' @param ... Additional arguments passed to the method-specific parser function.
#'   See individual parser documentation for available options:
#'   \itemize{
#'     \item \code{parseFACETS}: \code{methodName}
#'     \item \code{parsePureCN}: \code{methodName}
#'     \item \code{parseSCLUST}: \code{methodName}
#'     \item \code{parseBattenberg}: \code{methodName}, \code{pvalThreshold}
#'     \item \code{parseSequenza}: \code{methodName}
#'     \item \code{parseCNVkit}: (no additional arguments)
#'     \item \code{parseGATK}: \code{ploidy}
#'     \item \code{parseCfreec}: \code{filterLowUncertainty}, \code{pvalThreshold}
#'     \item \code{parseConseca}: (no additional arguments)
#'     \item \code{parseVCF}: (no additional arguments)
#'     \item \code{parseTSV}: (no additional arguments)
#'     \item \code{parseSEG}: \code{methodName}, \code{header}, \code{last.column}, \code{ploidy}
#'   }
#'
#' @return A cnvstrmliner object containing the parsed CNV segments with standardized
#'   metadata fields including:
#'   \itemize{
#'     \item Sample information (name, sex)
#'     \item Purity and ploidy estimates
#'     \item CNV segments as a GRanges object with allele-specific copy numbers
#'     \item Method-specific quality metrics
#'   }
#'
#' @note
#' \itemize{
#'   \item This function is the recommended entry point for parsing CNV data as it
#'     provides a consistent interface across all supported methods.
#'   \item If an unsupported method is specified, the function will stop with an
#'     informative error message listing all supported methods.
#'   \item Some parsers may require additional files to be present in the same
#'     directory as the input file (e.g., info files, subclone files).
#' }
#'
#' @seealso
#' Individual parser functions:
#' \code{\link{parseFACETS}}, \code{\link{parsePureCN}}, \code{\link{parseSCLUST}},
#' \code{\link{parseBattenberg}}, \code{\link{parseSequenza}}, \code{\link{parseCNVkit}},
#' \code{\link{parseGATK}}, \code{\link{parseCfreec}}, \code{\link{parseConseca}},
#' \code{\link{parseVCF}}, \code{\link{parseTSV}}, \code{\link{parseSEG}}
#'
#' @examples
#' \dontrun{
#' # Parse FACETS output
#' cnv_obj <- parseCNV("path/to/facets_output.vcf.gz", cnv.method = "FACETS")
#'
#' # Parse PureCN output with custom sample name
#' cnv_obj <- parseCNV("path/to/purecn_output",
#'   cnv.method = "PureCN",
#'   sample = "Sample001"
#' )
#'
#' # Parse Battenberg with custom p-value threshold
#' cnv_obj <- parseCNV("path/to/battenberg_output",
#'   cnv.method = "Battenberg",
#'   pvalThreshold = 0.01
#' )
#'
#' # Parse CNVkit output for hg19 genome
#' cnv_obj <- parseCNV("path/to/cnvkit.cns",
#'   cnv.method = "CNVkit",
#'   genome = "hg19"
#' )
#'
#' # Parse Control-FREEC output
#' cnv_obj <- parseCNV("path/to/freec_output",
#'   cnv.method = "FREEC",
#'   filterLowUncertainty = TRUE
#' )
#'
#' # Parse from a data frame (FACETS)
#' facets_df <- read.table("facets_segments.tsv", header = TRUE)
#' cnv_obj <- parseCNV(facets_df, cnv.method = "FACETS", sample = "Sample001")
#'
#' # Parse a generic VCF file with CNV data
#' cnv_obj <- parseCNV("path/to/cnv_output.vcf.gz", cnv.method = "VCF")
#'
#' # Parse a TSV file with allele-specific copy numbers
#' cnv_obj <- parseCNV("path/to/segments.tsv",
#'   cnv.method = "TSV",
#'   sample = "Sample001"
#' )
#'
#' # Parse a SEG file (e.g., from IGV or for GISTIC input)
#' cnv_obj <- parseCNV("path/to/sample.seg",
#'   cnv.method = "SEG",
#'   last.column = "lrr", header = TRUE
#' )
#' }
#'
#' @export

parseCNV <- function(cnv.data,
                     cnv.method,
                     sample = NA_character_,
                     genome = "hg38",
                     ...) {
  # Validate that method is provided
  supported_methods <- cnvstrmliner::getSupportedMethods()
  if (missing(cnv.method) || is.null(cnv.method) || is.na(cnv.method)) {
    stop("cnv.method must be specified. Supported methods: ", paste(supported_methods, collapse = ", "))
  }

  # Normalize method name to lowercase for matching
  method_lower <- tolower(cnv.method)

  # Dispatch to the appropriate parser based on method
  cnv_obj <- switch(method_lower,
    # FACETS
    "facets" = cnvstrmliner::parseFACETS(facets.data = cnv.data, sample = sample, genome = genome, ...),

    # PureCN
    "purecn" = cnvstrmliner::parsePureCN(purecn.filename = cnv.data, sample = sample, genome = genome, ...),

    # SCLUST
    "sclust" = cnvstrmliner::parseSCLUST(sclust.name = cnv.data, sample = sample, genome = genome, ...),

    # Battenberg
    "battenberg" = cnvstrmliner::parseBattenberg(battenberg.name = cnv.data, sample = sample, genome = genome, ...),
    "bb" = cnvstrmliner::parseBattenberg(battenberg.name = cnv.data, sample = sample, genome = genome, ...),

    # Sequenza
    "sequenza" = cnvstrmliner::parseSequenza(sequenza.name = cnv.data, sample = sample, genome = genome, ...),

    # CNVkit
    "cnvkit" = cnvstrmliner::parseCNVkit(cnvkit.filename = cnv.data, sample = sample, genome = genome, ...),
    "cnv-kit" = cnvstrmliner::parseCNVkit(cnvkit.filename = cnv.data, sample = sample, genome = genome, ...),

    # GATK
    "gatk" = cnvstrmliner::parseGATK(sample.file = cnv.data, sample = sample, genome = genome, ...),

    # Control-FREEC (multiple aliases)
    "freec" = cnvstrmliner::parseCfreec(cfreec.fn = cnv.data, sample = sample, genome = genome, ...),
    "cfreec" = cnvstrmliner::parseCfreec(cfreec.fn = cnv.data, sample = sample, genome = genome, ...),
    "controlfreec" = cnvstrmliner::parseCfreec(cfreec.fn = cnv.data, sample = sample, genome = genome, ...),
    "control-freec" = cnvstrmliner::parseCfreec(cfreec.fn = cnv.data, sample = sample, genome = genome, ...),

    # CloneCNA
    "clonecna" = cnvstrmliner::parseCloneCNA(clonecna.file = cnv.data, sample = sample, genome = genome, ...),


    # Conseca
    "conseca" = cnvstrmliner::parseConseca(conseca.data = cnv.data, sample = sample, genome = genome, ...),

    # VCF
    "vcf" = cnvstrmliner::parseVCF(vcf = cnv.data, sample = sample, genome = genome, ...),

    # SEG
    "seg" = cnvstrmliner::parseSEG(sample.file = cnv.data, sample = sample, genome = genome, ...),

    # TSV
    "tsv" = cnvstrmliner::parseTSV(tsv.filename = cnv.data, sample = sample, genome = genome, ...),


    # Unknown method - provide helpful error
    {
      stop(sprintf(
        "Unknown CNV method '%s'. Supported methods are: %s",
        cnv.method, paste(supported_methods, collapse = ", ")
      ))
    }
  )

  return(cnv_obj)
}


#' @rdname parseCNV
#' @export
loadCNV <- parseCNV


#' @rdname parseCNV
#' @export
readCNV <- parseCNV
