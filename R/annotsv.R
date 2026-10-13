#' Save CNV segments to an AnnotSV-compatible BED-like file
#'
#' @description
#' Writes non-neutral CNV segments from a \code{cnvstrml} object to a tab-delimited
#' file suitable for downstream annotation with AnnotSV.
#'
#' @param cnvObj A \code{cnvstrml} object containing CNV segments for a single sample.
#' @param filename Character path to the output file to write.
#'
#' @return Character path to the written file on success; \code{NA_character_} on failure.
#'
#' @details
#' The output contains one header line (commented) followed by records with the columns:
#' \describe{
#'   \item{chrom}{Chromosome name.}
#'   \item{start}{1-based segment start.}
#'   \item{end}{Segment end (inclusive).}
#'   \item{SV_type}{Variant type derived from \code{cnv.annotation}: \code{DEL} for deletions,
#'     \code{DUP} for amplifications, otherwise \code{CNV}. Neutral segments are excluded.}
#'   \item{sample}{Sample identifier taken from \code{cnvObj@sample}.}
#'   \item{loh}{Logical LOH flag from \code{cnv.loh}.}
#'   \item{ccf}{Clonal cancer cell fraction from \code{cnv.ccf}.}
#' }
#'
#' This format mirrors the BED-like input expected by AnnotSV; see the AnnotSV documentation
#' for full specifications: https://lbgi.fr/AnnotSV/Documentation/README.AnnotSV_latest.pdf
#'
#' @examples
#' \dontrun{
#' out <- tempfile(fileext = ".tsv")
#' saveAnnotSV(cnvObj = cnv_sample, filename = out)
#' readLines(out, n = 3)
#' }
#'
#' @seealso \code{\link{annotateCNV}}, \code{\link{parseManta}}
#'
#' @export
saveAnnotSV <- function(cnvObj, filename) {
  # we store a line with the column names
  cat("#chrom\tstart\tend\tSV_type\tsample\tloh\tccf", file = filename, append = FALSE)

  # store the data from the cnvObj
  df <- cnvObj@segments %>%
    as.data.frame() %>%
    dplyr::select(
      chrom = seqnames,
      start,
      end,
      sv_type = cnv.annotation,
      loh = cnv.loh,
      ccf = cnv.ccf
    ) %>%
    dplyr::filter(sv_type != "NEUT") %>%
    dplyr::mutate(
      sv_type = dplyr::case_when(
        grepl(pattern = "DEL", sv_type) ~ "DEL",
        grepl(pattern = "AMP", sv_type) ~ "DUP",
        TRUE ~ "CNV"
      ),
      sample = cnvObj@sample
    ) %>%
    dplyr::relocate(
      chrom, start, end, sv_type, sample, loh, ccf
    )

  filename <- tryCatch(
    {
      write.table(df, file = filename, append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
      filename
    },
    error = function(e) {
      logger::log_error("Unable to store the file {filename}")
      return(NA_character_)
    }
  )
  return(filename)
}
