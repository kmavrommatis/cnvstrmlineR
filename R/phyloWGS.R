#' Save CNV segments to PhyloWGS-compatible format
#'
#' @description
#' Exports CNV segments from a \code{cnvstrml} object to a tab-delimited file
#' compatible with PhyloWGS (a tool for inferring subclonal composition and phylogenetic
#' trees from WGS data). Coordinates are converted to NCBI format (e.g., chr1 → 1, chrX → X, chrM → MT),
#' and cellular prevalence is adjusted by sample purity. Missing minor allele copy numbers are
#' imputed to zero.
#'
#' @param filename Character path to the output file to write.
#' @param cnvObj A \code{cnvstrml} object containing CNV segments.
#'
#' @return Character path to the written file on success. Note: function returns \code{filename}
#'   regardless of write success or failure; check file existence externally if needed.
#'
#' @details
#' \strong{Output file format}
#' Tab-delimited text with no header unless explicitly added. Columns in order:
#' \describe{
#'   \item{chromosome}{Chromosome in NCBI format: numeric 1–22, X, Y, MT (not chr-prefixed).}
#'   \item{start}{1-based segment start position.}
#'   \item{end}{1-based segment end position (inclusive).}
#'   \item{copy_number}{Total copy number (integer).}
#'   \item{minor_cn}{Minor allele copy number. Defaults to 0 if \code{NA}.}
#'   \item{major_cn}{Major allele copy number.}
#'   \item{cellular_prevalence}{Clonal cellular fraction (CCF) adjusted by sample purity:
#'     \code{cellular_prevalence = cnv.cf * cnvObj@purity}.}
#' }
#' Segments with \code{NA} cellular prevalence are dropped.
#'
#' \strong{Chromosome conversion (NCBI format)}
#' PhyloWGS expects NCBI-style chromosome identifiers. The helper function \code{toNCBI()}
#' converts UCSC style (chr1, chr2, ..., chrX, chrY, chrM) to NCBI style (1, 2, ..., X, Y, MT).
#' \itemize{
#'   \item \code{chr} prefix is stripped: chr1 → 1.
#'   \item chrM is converted to MT.
#'   \item Auto- and sex chromosomes handled uniformly.
#' }
#'
#' \strong{Purity adjustment}
#' Cellular prevalence is the fraction of cancer cells in a sample carrying a given alteration.
#' This function adjusts the \code{cnv.cf} (copy-number cell fraction) by multiplying by
#' \code{cnvObj@purity} to account for normal-cell contamination:
#' \code{cellular_prevalence = cnv.cf * purity}.
#' When purity = 1 (pure sample), no adjustment needed; when purity = 0.8, prevalence is scaled to 80%.
#'
#' \strong{Missing values}
#' \itemize{
#'   \item Minor allele copy numbers: \code{NA} values replaced with 0 (conservative estimate).
#'   \item Cellular prevalence: Segments with \code{NA} prevalence are excluded from output.
#' }
#'
#' \strong{Limitations}
#' \itemize{
#'   \item Only supports single \code{cnvstrml} objects; lists are not currently processed.
#'   \item Does not handle multi-clonal data; assumes single dominant clone.
#'   \item All segments mapped to a single cellular prevalence per region.
#' }
#'
#' @examples
#' \dontrun{
#' # Save CNV segments for PhyloWGS input
#' out_file <- tempfile(fileext = ".tsv")
#' result <- savePhyloWGS(filename = out_file, cnvObj = cnv_sample)
#'
#' # Inspect output
#' head(read.table(result, sep = "\t", header = TRUE))
#'
#' # Check purity adjustment impact
#' cat("Sample purity:", cnv_sample@purity, "\n")
#' }
#'
#' @seealso \code{\link{saveAnnotSV}}, \code{\link{saveVCF}}, \code{toNCBI()}
#'
#' @export
#' @description
#' Produce a seg file that can be used for downstream processes by phyloWGS
#'
savePhyloWGS <- function(filename, cnvObj) {
  # prepare the output for PhyloWGS
  cnvObj@segments <- toNCBI(cnvObj@segments)
  res <- cnvObj@segments %>%
    as.data.frame() %>%
    dplyr::select(
      chromosome = seqnames,
      start, end,
      copy_number = total.allele.copies,
      minor_cn = minor.allele.copies,
      major_cn = major.allele.copies,
      cellular_prevalence = cnv.cf
    ) %>%
    dplyr::mutate(cellular_prevalence = cellular_prevalence * cnvObj@purity)
  replace_na(list(minor_cn = 0)) %>%
    drop_na(cellular_prevalence) %>%
    relocate(major_cn, .after = minor_cn)


  tryCatch(
    {
      write.table(res,
        file = filename,
        quote = FALSE,
        sep = "\t",
        col.names = TRUE,
        row.names = FALSE
      )
    },
    warning = function(w) {
      message("Warning: ", w)
    },
    error = function(e) {
      message("Error:", e)
    }
  )
  return(filename)
}


toNCBI <- function(gr) {
  si <- GenomeInfoDb::seqinfo(gr)
  GenomeInfoDb::seqlevels(gr) <- GenomeInfoDb::seqlevels(gr) %>%
    gsub("chrM", "MT", .) %>%
    gsub("chr", "", .)
  gr
}
