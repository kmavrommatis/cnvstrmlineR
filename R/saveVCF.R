#' Create a VCF header to store the output from the tools
#'
#' @param cnvobj (cnvstrmliner object) the object containing the data to store
#' @returns a string with the header of the vcf file.

make_header <- function(cnvobj) {
  # extra: Named vector of additional information. E.g. c(purity=0.5, ploidy= 2.1)
  genome <- cnvobj@genome
  header <- c(
    "##fileformat=VCFv4.2",
    sprintf("##reference=%s", genome)
  )
  # Contigs
  # --------------------------

  seqinfo <- getSeqInfo(cnvobj@genome)
  # chrom_size<- NULL
  # if(gbuild == 'hg19') {chrom_size<- genomes$HG19}
  # if(gbuild == 'hg38') {chrom_size<- genomes$HG38}
  # if(gbuild == 'mm9') {chrom_size<- genomes$MM9}
  # if(gbuild == 'mm10') {chrom_size<- genomes$MM10}

  stopifnot(!is.null(seqinfo))
  size <- GenomeInfoDb::seqlengths(seqinfo)

  name <- GenomeInfoDb::seqlevels(seqinfo)
  header <- c(header, sprintf("##contig=<ID=%s,length=%s>", name, size))

  # --------------------------

  header <- c(header, '##FILTER=<ID=PASS,Description="All filters passed">')
  header <- c(header, '##FILTER=<ID=neutral,Description="Copy number neutral">')
  header <- c(header, '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">')
  header <- c(header, '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">')
  header <- c(header, '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">')
  header <- c(header, '##INFO=<ID=NUM_MARK,Number=1,Type=Integer,Description="Number of SNPs in the segment">')
  header <- c(header, '##INFO=<ID=NUM_HET_MARK,Number=1,Type=Integer,Description="Number of SNPs that are deemed heterozygous">')
  header <- c(header, '##INFO=<ID=CN_LR_MEDIAN,Number=1,Type=Float,Description="Median log-ratio (logR) of the segment. logR is defined by the log-ratio of total read depth in the tumor versus that in the normal">')
  header <- c(header, '##INFO=<ID=CF,Number=1,Type=Float,Description="Cellular fraction, fraction of DNA associated with the aberrant genotype. Set to 1 for normal diploid">')
  header <- c(header, '##INFO=<ID=CCF,Number=1,Type=Float,Description="Cancer cellular fraction, fraction of DNA associated with the aberrant genotype within the cancer cells. Set to 1 for clonal aberrations">')
  header <- c(header, '##INFO=<ID=TOTAL_CN,Number=1,Type=Integer,Description="Total copy number. 2 for normal diploid">')
  header <- c(header, '##INFO=<ID=MAJOR_CN,Number=1,Type=Integer,Description="Major copy number. 1 for normal diploid">')
  header <- c(header, '##INFO=<ID=MINOR_CN,Number=1,Type=Integer,Description="Lesser (minor) copy number. 1 for normal diploid">')
  header <- c(header, '##INFO=<ID=CNV_ANN,Number=.,Type=String,Description="Annotation features assigned to this CNV">')
  header <- c(header, '##INFO=<ID=LOH,Number=.,Type=Integer,Description="Is this segment predicted to be LOH (1) or not (0)">')
  header <- c(header, '##INFO=<ID=CNV_QC_PASS,Number=.,Type=Integer,Description="Is this segment passing QC (1) or not (0), -1 for NA">')
  header <- c(header, paste0("##cnvstrmliner=", packageVersion("cnvstrmliner")))


  # fields from the cnvstrmliner object
  # that contain standard information expected always in the
  # output VCF file
  extra <- c(
    "purity" = cnvobj@purity,
    "ploidy" = cnvobj@ploidy,
    "effectiveploidy" = cnvobj@effective_ploidy,
    "sex" = cnvobj@sex,
    "cnvmethod" = cnvobj@method
  )

  # additional fields that we capture in the quality slot
  # which may or may not be present
  extra2 <- tryCatch(
    {
      X <- lapply(cnvobj@quality, function(x) {
        if (is.numeric(x) | is.character(x) | is.logical(x)) {
          return(x)
        } else {
          return(NA)
        }
      })
      X <- X[!is.na(X)]

      names(X) <- paste("quality_", names(X), sep = "")
      names(X) <- gsub(pattern = " ", replacement = "", names(X)) %>%
        gsub(pattern = "_", replacement = "")
      X
    },
    warning = function(w) {
      message(w)
      NA
    },
    error = function(e) {
      message(e)
      NA
    }
  )
  if (!all(is.na(extra2))) {
    extra <- c(extra, extra2)
  }


  for (i in 1:length(extra)) {
    key <- names(extra)[i]
    value <- extra[i]
    stopifnot(key != "")
    header <- c(header, sprintf("##%s=%s", key, value))
  }
  header <- c(header, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
  return(paste(header, collapse = "\n"))
}

#' Convert a record (segment) to a VCF record
#'
#' @param segments (GRanges) the segments of this object


recordToVcf <- function(segments) {
  stopifnot(class(segments) == "GRanges")
  stopifnot(length(segments) > 0)
  # Convert the annotated facets record to a VCF record.
  info <- data.frame(
    "SVTYPE" = paste("SVTYPE",
      dplyr::case_when(
        S4Vectors::mcols(segments)[, "cnv.annotation"] %>% grepl("DEL", .) ~ "DEL",
        S4Vectors::mcols(segments)[, "cnv.annotation"] %>% grepl("AMP", .) ~ "DUP",
        TRUE ~ "CNV"
      ),
      sep = "="
    ),
    "SVLEN" = paste("SVLEN", BiocGenerics::end(segments) - BiocGenerics::start(segments), sep = "="),
    "END" = paste("END", BiocGenerics::end(segments), sep = "="),
    "NUM_MARK" = paste("NUM_MARK", S4Vectors::mcols(segments)[, "cnv.totalmarkers"], sep = "="),
    "NUM_HET_MARK" = paste("NUM_HET_MARK", S4Vectors::mcols(segments)[, "cnv.markers"], sep = "="),
    "CN_LR_MEDIAN" = paste("CN_LR_MEDIAN", round(S4Vectors::mcols(segments)[, "cnv.lrr"], 3), sep = "="),
    "CF" = paste("CF", S4Vectors::mcols(segments)[, "cnv.cf"], sep = "="),
    "CCF" = paste("CCF", S4Vectors::mcols(segments)[, "cnv.ccf"], sep = "="),
    "TOTAL_CN" = paste("TOTAL_CN", S4Vectors::mcols(segments)[, "total.allele.copies"], sep = "="),
    "MAJOR_CN" = paste("MAJOR_CN", S4Vectors::mcols(segments)[, "major.allele.copies"], sep = "="),
    "MINOR_CN" = paste("MINOR_CN", S4Vectors::mcols(segments)[, "minor.allele.copies"], sep = "="),
    "CNV_ANN" = paste("CNV_ANN", S4Vectors::mcols(segments)[, "cnv.annotation"], sep = "="),
    "LOH" = paste("LOH", as.integer(S4Vectors::mcols(segments)[, "cnv.loh"]), sep = "=")
  )

  if ("cnv.cnaqc.pass" %in% colnames(S4Vectors::mcols(segments))) {
    cbind(info, "CNA_QC_PASS" = ifelse(is.na(S4Vectors::mcols(segments)[, "cnv.cnaqc.pass"]),
      "CNV_QC_PASS=-1",
      paste("CNV_QC_PASS", as.integer(S4Vectors::mcols(segments)[, "cnv.cnaqc.pass"]), sep = "=")
    ))
  }
  # CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO'

  vcf <- data.frame(
    CHROM = GenomeInfoDb::seqnames(segments) %>% as.character(),
    POS = BiocGenerics::start(segments),
    REF = "N",
    ALT = "<CNV>",
    QUAL = ".",
    FILTER = ifelse(S4Vectors::mcols(segments)[, "cnv.annotation"] == "NEUT", "neutral", "PASS")
  ) %>%
    dplyr::mutate(ID = dplyr::row_number(), .after = "POS") %>%
    dplyr::mutate(INFO = info %>% apply(., 1, paste, collapse = ";"))


  return(vcf)
}

#' Write a VCF file with the contents of the CNV object
#'
#' @description
#' Exports copy number variation (CNV) segments from a cnvstrmliner object to a
#' VCF (Variant Call Format) file. The output file is automatically compressed
#' using bgzip and indexed with tabix for efficient random access and compatibility
#' with downstream bioinformatics tools.
#'
#' @details
#' This function performs the following operations:
#' \enumerate{
#'   \item Sanitizes the filename by removing any existing compression or VCF extensions
#'   \item Creates a temporary VCF file with a properly formatted header containing:
#'     \itemize{
#'       \item VCF version (4.2)
#'       \item Reference genome information
#'       \item Contig definitions with chromosome lengths
#'       \item INFO field definitions for CNV-specific attributes
#'       \item Sample metadata (purity, ploidy, sex, CNV calling method)
#'       \item Quality metrics from the cnvstrmliner object
#'     }
#'   \item Converts each CNV segment to a VCF record with the following INFO fields:
#'     \itemize{
#'       \item \strong{SVTYPE}: Type of structural variant (DEL, DUP, or CNV)
#'       \item \strong{SVLEN}: Length of the variant
#'       \item \strong{END}: End position of the segment
#'       \item \strong{NUM_MARK}: Total number of SNP markers in the segment
#'       \item \strong{NUM_HET_MARK}: Number of heterozygous markers
#'       \item \strong{CN_LR_MEDIAN}: Median log-ratio of the segment
#'       \item \strong{CF}: Cellular fraction
#'       \item \strong{CCF}: Cancer cell fraction
#'       \item \strong{TOTAL_CN}: Total copy number
#'       \item \strong{MAJOR_CN}: Major allele copy number
#'       \item \strong{MINOR_CN}: Minor allele copy number
#'       \item \strong{CNV_ANN}: CNV annotation (e.g., AMP, DEL, NEUT)
#'       \item \strong{LOH}: Loss of heterozygosity indicator (1 or 0)
#'     }
#'   \item Compresses the VCF file using bgzip (produces .vcf.gz)
#'   \item Creates a tabix index file (.vcf.gz.tbi) for random access
#' }
#'
#' The FILTER column is set to "PASS" for all CNV events except neutral segments,
#' which are marked with the "neutral" filter. Neutral segments are typically
#' ones that have total allele copies equal to ploidy.
#'
#' @param filename Character string specifying the output file path. The function
#'   automatically handles file extensions - any existing `.vcf`, `.vcf.gz`, `.bgz`,
#'   or `.bcf` extensions are stripped and the final output will have `.vcf.gz`
#'   extension. The tabix index will be created with `.vcf.gz.tbi` extension.
#' @param cnvObj A cnvstrmliner object containing the CNV segments to export.
#'   The object must have valid segments in the `@segments` slot (as a GRanges object)
#'   and populated metadata fields including purity, ploidy, sex, method, and genome.
#'
#' @return Character string containing the path to the output compressed VCF file
#'   (with `.vcf.gz` extension), or NULL if an error occurred during file writing
#'   or compression.
#'
#' @note
#' \itemize{
#'   \item The function requires the \pkg{Rsamtools} package for bgzip compression
#'     and tabix indexing.
#'   \item The output VCF follows the VCF 4.2 specification with structural variant
#'     conventions for representing CNV events.
#'   \item Temporary files are created in the same directory as the output file
#'     and are cleaned up after successful compression.
#'   \item If the destination file already exists, it will be overwritten.
#' }
#'
#' @seealso
#' \code{\link{make_header}} for VCF header generation,
#' \code{\link[Rsamtools]{bgzip}} for file compression,
#' \code{\link[Rsamtools]{indexTabix}} for tabix indexing
#'
#' @examples
#' \dontrun{
#' # Load or create a cnvstrmliner object
#' cnv_obj <- parseFACETS("path/to/facets_output.tsv")
#'
#' # Export to VCF (will create output.vcf.gz and output.vcf.gz.tbi)
#' output_path <- writeVCF("path/to/output", cnv_obj)
#'
#' # The function handles various input filename formats
#' writeVCF("results/sample1.vcf.gz", cnv_obj) # -> results/sample1.vcf.gz
#' writeVCF("results/sample1.vcf", cnv_obj) # -> results/sample1.vcf.gz
#' writeVCF("results/sample1", cnv_obj) # -> results/sample1.vcf.gz
#' }
#'
#' @export
#' @importFrom Rsamtools bgzip indexTabix

writeVCF <- function(filename, cnvObj) {
  filename <- filename %>%
    gsub(".gz$", "", .) %>%
    gsub(".bgz$", "", .) %>%
    gsub(".vcf$", "", .) %>%
    gsub(".bcf$", "", .)

  vcf_tmp <- tempfile(pattern = paste0(basename(filename), "."), tmpdir = dirname(filename), fileext = ".vcf")

  header <- make_header(cnvObj)
  body <- recordToVcf(cnvObj@segments)
  tryCatch(
    {
      write(header, vcf_tmp)
      write.table(body, vcf_tmp, append = TRUE, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
    },
    warning = function(w) {
      message("Warning: ", w)
    },
    error = function(e) {
      message("Error:", e)
    }
  )
  dest <- NULL
  dest <- tryCatch(
    {
      dest <- paste0(filename, ".vcf.gz")
      x_ <- Rsamtools::bgzip(vcf_tmp, dest = dest, overwrite = TRUE)
      x_ <- Rsamtools::indexTabix(paste0(filename, ".vcf.gz"), format = "vcf4")
      return(dest)
    },
    warning = function(w) {
      message("Warning: ", w)
    },
    error = function(e) {
      message("Error:", e)
    },
    finally = {
      unlink(vcf_tmp)
    }
  )
  return(dest)
}


#' @rdname writeVCF
#' @export
saveVCF <- writeVCF
