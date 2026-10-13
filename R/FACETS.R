#' create a cnvstrml object from a dataframe with facets data
#'
#' @param df (data.frame)  dataframe produced by FACETS
#' @param purity (numeric)  the purity produces by FACETS (a value between 0 - 1)
#' @param ploidy (numeric) the ploidy produced by FACETS
#' @param genome (character) the genome version assumed for this dataset
#' @param sample (character)  the name of the sample
#' @return an incomplete cnvstrml object
#'
parseFACETS_df <- function(
  df,
  purity = NA_real_,
  ploidy = NA_real_,
  sample = NA_character_,
  genome = "hg38"
) {
  logger::log_debug("Parsing FACETS dataframe with {nrow(df)} segments", namespace = "parser")
  # convert the input to a dataframe
  if (!is.data.frame(df)) {
    df <- as.data.frame(df)
  }
  # FACETS tend to call Chr X as chr23
  df[which(df$chrom == 23), "chrom"] <- "X"
  logger::log_debug("Converted chromosome 23 to X notation", namespace = "parser")

  # convert to GRanges
  # and normalize the chromosome name to UCSC
  seqinfo <- getSeqInfo(genome)
  logger::log_debug("Creating GRanges object for genome {genome}", namespace = "parser")
  facets.gr <- GenomicRanges::makeGRangesFromDataFrame(df,
    keep.extra.columns = TRUE,
    ignore.strand = TRUE,
    seqinfo = seqinfo
  ) %>%
    GenomeInfoDb::sortSeqlevels() %>%
    BiocGenerics::sort()
  GenomeInfoDb::seqlevelsStyle(facets.gr) <- "UCSC"


  # Add the flags
  logger::log_debug("Adding annotation flags and copy number information", namespace = "parser")

  mcols_new <- S4Vectors::mcols(facets.gr) %>%
    tibble::as_tibble() %>%
    dplyr::mutate(
      total.allele.copies = tcn.em,
      major.allele.copies = ifelse(is.na(lcn.em), tcn.em, tcn.em - lcn.em),
      minor.allele.copies = lcn.em,
      cnv.cf = cf.em
    )

  S4Vectors::mcols(facets.gr) <- mcols_new

  logger::log_debug("Creating cnvstrml object for sample {sample}", namespace = "parser")
  obj <- cnvstrml(
    segments = facets.gr,
    purity = as.numeric(purity),
    ploidy = as.numeric(ploidy),
    effective_ploidy = calculateEffectivePloidy(facets.gr),
    method = "FACETS",
    sample = sample,
    genome = genome
  )
  obj
}


#' parse a VCF file in BND format with FACETS data
#'
#' Such a file is created from the cnv_facets tool
#'
#' @param vcf (character) the input VCF file
#' @param genome (character)  the genome version assumed for this dataset
#' @param sample (character)  the name of the sample
#' @return an incomplete cnvstrml object
#'
parseFACETS_vcf <- function(
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

  dipLogR <- VariantAnnotation::meta(VariantAnnotation::header(ff.vcf))$dipLogR[1, 1] %>% as.numeric()
  emflags <- VariantAnnotation::meta(VariantAnnotation::header(ff.vcf))$emflags[, 1]
  estinsertsize <- VariantAnnotation::meta(VariantAnnotation::header(ff.vcf))$estinsertsize[1, 1] %>% as.numeric()
  purity_cval <- tryCatch(
    {
      VariantAnnotation::meta(VariantAnnotation::header(ff.vcf))$prepuritycval[1, 1] %>% as.numeric()
    },
    warning = function(w) {
      return(NA)
    },
    error = function(e) {
      return(NA)
    }
  )
  pre_cval <- VariantAnnotation::meta(VariantAnnotation::header(ff.vcf))$precval[1, 1] %>% as.numeric()
  proc_cval <- VariantAnnotation::meta(VariantAnnotation::header(ff.vcf))$proccval[1, 1] %>% as.numeric()
  nbhd <- VariantAnnotation::meta(VariantAnnotation::header(ff.vcf))$snpnbhd[1, 1] %>% as.numeric()
  matchednormal <- VariantAnnotation::meta(VariantAnnotation::header(ff.vcf))$matchednormal[1, 1] %>% as.logical()
  tumor <- VariantAnnotation::meta(VariantAnnotation::header(ff.vcf))$tumorsample[1, 1] %>% as.character()
  normal <- VariantAnnotation::meta(VariantAnnotation::header(ff.vcf))$normalsample[1, 1] %>% as.character()

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
    tidyr::replace_na(list(LCN_EM = 0, CF_EM = 0)) %>%
    dplyr::mutate(
      total.allele.copies = TCN_EM,
      minor.allele.copies = LCN_EM,
      major.allele.copies = TCN_EM - LCN_EM,
      cnv.totalmarkers = NUM_MARK,
      cnv.markers = NHET,
      cnv.lrr = CNLR_MEDIAN,
      cnv.cf = CF_EM
    )

  S4Vectors::mcols(ff.gr) <- mcols_new
  GenomicRanges::end(ff.gr) <- ff.gr$END

  idx <- which((ff.gr$SVLEN == GenomicRanges::width(ff.gr)) == FALSE)
  if (length(idx) > 0) {
    GenomicRanges::start(ff.gr)[idx] <- ff.gr$END[idx] - ff.gr$SVLEN[idx] + 1
  }

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

  qual.df <- list(
    flags = emflags,
    pre_cval = pre_cval,
    proc_cval = proc_cval,
    estinsertsize = estinsertsize,
    matchednormal = matchednormal,
    nbhd = nbhd,
    purity_cval = purity_cval,
    dipLogR = dipLogR,
    tumor = tumor,
    normal = normal
  )
  logger::log_debug("Creating cnvstrml object for sample {sample} with purity {purity} and ploidy {ploidy}", namespace = "parser")
  obj <- cnvstrml(
    segments = ff.gr,
    purity = as.numeric(purity),
    ploidy = as.numeric(ploidy),
    effective_ploidy = calculateEffectivePloidy(ff.gr),
    method = "FACETS",
    sample = sample,
    quality = qual.df,
    genome = genome
  )
  obj
}

#' Parse FACETS CNV and allelic imbalance output
#'
#' @description
#' High-level parser for FACETS (Fraction and Allele-specific Copy number Estimates from Tumor Sequencing)
#' output. FACETS is a tool for detecting somatic CNVs and allelic imbalances (including LOH) from
#' sequencing data (WGS, WES, targeted panels). It also estimates tumour purity and ploidy.
#' This function dispatches to VCF or data frame parsers and constructs a \code{cnvstrml} object
#' with standardized fields and CNV annotations.
#'
#' @param facets.data Character path to a FACETS VCF file (\code{.vcf}, \code{.vcf.gz}, \code{.vcf.bgz}),
#'   or a data frame of FACETS segments (e.g., from \code{read.table()} of FACETS output).
#' @param sample Character sample identifier. If a VCF file is supplied, the sample name is extracted
#'   from the VCF header; if not found, the filename (minus extension) is used. For data frames, must be
#'   supplied explicitly or defaults to \code{NA_character_}.
#' @param genome Character genome build (e.g., \code{"hg38"}) used to set \code{Seqinfo} and UCSC seqlevels.
#' @param methodName Character method description to store (default: \code{"FACETS"}).
#'
#' @return A \code{cnvstrml} object with:
#' \describe{
#'   \item{@segments}{A \code{GRanges} with metadata:
#'     \describe{
#'       \item{total.allele.copies}{From FACETS \code{TCN_EM} (VCF) or \code{tcn.em} (data frame).}
#'       \item{major.allele.copies}{Total copies minus minor allele copies.}
#'       \item{minor.allele.copies}{From FACETS \code{LCN_EM} (VCF) or \code{lcn.em}.}
#'       \item{cnv.cf}{Cancer cell fraction from \code{CF_EM}.}
#'       \item{cnv.lrr}{Log-ratio median (VCF only, from \code{CNLR_MEDIAN}).}
#'       \item{cnv.totalmarkers}{Total SNP markers (VCF only, from \code{NUM_MARK}).}
#'       \item{cnv.markers}{Heterozygous SNPs used (VCF only, from \code{NHET}).}
#'       \item{cnv.annotation}{CNV status (AMP, DEL, NEUT, etc.) assigned via \code{annotateCNV()}.}
#'     }}
#'   \item{@purity}{Purity estimate from VCF header or \code{NA} if from data frame.}
#'   \item{@ploidy}{Ploidy estimate from VCF header or \code{NA} if from data frame.}
#'   \item{@effective_ploidy}{Rounded ploidy, calculated via \code{calculateEffectivePloidy()}.}
#'   \item{@sample}{Sample name (from header, filename, or user-supplied).}
#'   \item{@method}{Method identifier, default \code{"FACETS"}.}
#'   \item{@quality}{List of metadata, including processing parameters from VCF header (e.g.,
#'     \code{emflags}, \code{dipLogR}, \code{estinsertsize}, \code{matchednormal}, \code{nbhd},
#'     \code{pre_cval}, \code{proc_cval}, \code{tumor}, \code{normal} for VCF inputs).}
#' }
#'
#' @details
#' \strong{Input dispatch logic}
#' \itemize{
#'   \item \emph{VCF file}: If \code{facets.data} is a character path ending in \code{.vcf},
#'     \code{.vcf.gz}, or \code{.vcf.bgz}, the function calls \code{parseFACETS_vcf()}.
#'     Segments with \code{FILTER} not in \{PASS, neutral\} are discarded.
#'   \item \emph{Data frame}: If \code{facets.data} is a data frame, \code{parseFACETS_df()} is used.
#'     Purity and ploidy are \code{NA} unless set separately; chromosomes are normalized
#'     (FACETS chr 23 → chr X).
#' }
#'
#' \strong{VCF processing}
#' \itemize{
#'   \item Purity and ploidy extracted from VCF header metadata (\code{PURITY}, \code{PLOIDY}; case-insensitive).
#'   \item Sample name extracted from VCF header or filename if not provided.
#'   \item Quality metrics stored: \code{emflags}, \code{dipLogR}, SNP neighbourhood size, calibration values, etc.
#'   \item Segment coordinates adjusted if \code{SVLEN} (from INFO) disagrees with interval width.
#' }
#'
#' \strong{Data frame processing}
#' \itemize{
#'   \item Expected columns: \code{chrom}, \code{start}, \code{end}, \code{tcn.em} (total copies),
#'     \code{lcn.em} (minor alleles, optional), \code{cf.em} (cancer fraction).
#'   \item Chr 23 automatically converted to X notation.
#'   \item Effective ploidy calculated from \code{tcn.em} distribution.
#' }
#'
#' \strong{Annotation and CCF}
#' \itemize{
#'   \item After parsing, CNV status is assigned via \code{annotateCNV()} using effective ploidy.
#'   \item CCF (clonal cancer cell fraction) is not computed; use \code{addCCF()} post-parsing if needed.
#'   \item Both VCF and data frame outputs include all detected segments (including NEUT).
#' }
#'
#' @examples
#' \dontrun{
#' # Parse FACETS VCF output (includes purity/ploidy in header)
#' cnv_vcf <- parseFACETS("sample_facets.vcf.gz")
#' cat("Purity:", cnv_vcf@purity, "Ploidy:", cnv_vcf@ploidy, "\n")
#'
#' # Parse FACETS data frame (purity/ploidy must be set separately)
#' df <- read.table("facets_segments.tsv", header = TRUE)
#' cnv_df <- parseFACETS(df, sample = "S1", genome = "hg38")
#'
#' # Inspect segments and annotations
#' head(cnv_vcf@segments)
#' table(S4Vectors::mcols(cnv_vcf@segments)$cnv.annotation)
#' }
#'
#' @seealso \code{\link{parseFACETS_vcf}}, \code{\link{parseFACETS_df}}, \code{\link{parseCNV}},
#'   \code{\link{annotateCNV}}, \code{\link{addCCF}}
#'
#' @export
#'
parseFACETS <- function(
  facets.data,
  sample = NA_character_,
  genome = "hg38",
  methodName = "FACETS"
) {
  logger::log_debug("FACETS root is {facets.data}", namespace = "parser")
  # dispatch to the right function

  if (is.data.frame(facets.data)) {
    # message("Preparing the segments from a data frame. Purity and Ploidy have to be set independently.")
    logger::log_debug("Parsing a data.frame", namespace = "parser")
    obj <- parseFACETS_df(facets.data, sample = sample, genome = genome)
  } else if (is.character(facets.data) && file.exists(facets.data) &&
    (endsWith(facets.data, "vcf") || endsWith(facets.data, "vcf.gz") || endsWith(facets.data, "vcf.bgz"))
  ) {
    logger::log_debug("Loading FACETS from {facets.data}", namespace = "parser")
    # message("Preparing the segments from the vcf file ", basename( facets.data),". Purity and Ploidy are expected to be found in the header of the file.")
    obj <- parseFACETS_vcf(facets.data, sample = sample, genome = genome)
  } else {
    logger::log_error("Unable to find the data. The input {facets.data} does not seem to be a data.frame or a file", namespace = "parser")
    stop("Exiting parseFACETS")
  }
  # obj=cnvstrmliner::calculateEffectivePloidy(obj)
  # logger::log_debug("Adding CNV annotations and calculating CCF", namespace = "parser")
  # mcols_new <- S4Vectors::mcols(obj@segments) %>% as.data.frame()
  # flags <- apply(mcols_new, 1, function(K) {
  #   cnvstrmliner::annotateCNV(
  #     major.cn = K["major.allele.copies"], minor.cn = K["minor.allele.copies"],
  #     ploidy = obj@effective_ploidy
  #   )
  # })

  obj@method <- methodName
  # S4Vectors::mcols(obj@segments) <- mcols_new %>%
  # dplyr::mutate(cnv.annotation = flags)
  # obj=cnvstrmliner::addCCF( obj )

  obj@quality$file <- c("segments" = facets.data)
  obj
}


# functions to save mFACETS objects
# these functions mirror the saveVCF functions


#' Create a VCF header to store the output from the mFACETS
#'
#' @param cnvobj (cnvstrmliner object) the object containing the data to store
#' @returns a string with the header of the vcf file.

make_header_mFACETS <- function(cnvobj) {
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
  header <- c(header, '##INFO=<ID=NHET,Number=1,Type=Integer,Description="Number of SNPs that are deemed heterozygous">')
  header <- c(header, '##INFO=<ID=CNLR_MEDIAN,Number=1,Type=Float,Description="Median log-ratio (logR) of the segment. logR is defined by the log-ratio of total read depth in the tumor versus that in the normal">')
  header <- c(header, '##INFO=<ID=CF_EM,Number=1,Type=Float,Description="Cellular fraction, fraction of DNA associated with the aberrant genotype. Set to 1 for normal diploid">')
  header <- c(header, '##INFO=<ID=TCN_EM,Number=1,Type=Integer,Description="Total copy number. 2 for normal diploid">')
  header <- c(header, '##INFO=<ID=LCN_EM,Number=1,Type=Integer,Description="Lesser (minor) copy number. 1 for normal diploid">')
  header <- c(header, '##INFO=<ID=CNV_ANN,Number=.,Type=String,Description="Annotation features assigned to this CNV">')
  header <- c(header, '##INFO=<ID=LOH,Number=.,Type=Boolean,Description="Is this segment predicted to be LOH">')
  header <- c(header, paste0("##cnvstrmliner=", packageVersion("cnvstrmliner")))

  extra <- c(
    "purity" = cnvobj@purity,
    "ploidy" = cnvobj@ploidy,
    "sex" = cnvobj@sex,
    "cnv_method" = cnvobj@method
  )

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

      names(X) <- gsub(pattern = " ", replacement = "", names(X))
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


recordToVcf_mFACETS <- function(segments) {
  stopifnot(class(segments) == "GRanges")
  stopifnot(length(segments) > 0)
  # Convert the annotated facets record to a VCF record.
  info <- data.frame(
    "SVTYPE" = paste("SVTYPE", S4Vectors::mcols(segments)[, "cnv.annotation"], sep = "="),
    "SVLEN" = paste("SVLEN", BiocGenerics::end(segments) - BiocGenerics::start(segments), sep = "="),
    "END" = paste("END", BiocGenerics::end(segments), sep = "="),
    "NUM_MARK" = paste("NUM_MARK", S4Vectors::mcols(segments)[, "cnv.totalmarkers"], sep = "="),
    "NHET" = paste("NHET", S4Vectors::mcols(segments)[, "cnv.markers"], sep = "="),
    "CNLR_MEDIAN" = paste("CNLR_MEDIAN", round(S4Vectors::mcols(segments)[, "cnv.lrr"], 3), sep = "="),
    "CF_EM" = paste("CF_EM", S4Vectors::mcols(segments)[, "cnv.cf"], sep = "="),
    "TCN_EM" = paste("TCN_EM", S4Vectors::mcols(segments)[, "total.allele.copies"], sep = "="),
    "LCN_EM" = paste("LCN_EM", S4Vectors::mcols(segments)[, "minor.allele.copies"], sep = "="),
    "CNV_ANN" = paste("CNV_ANN", S4Vectors::mcols(segments)[, "cnv.annotation"], sep = "="),
    "LOH" = paste("LOH", as.integer(S4Vectors::mcols(segments)[, "cnv.loh"]), sep = "=")
  ) %>%
    dplyr::mutate(SVTYPE = ifelse(grepl("DEL", SVTYPE), "DEL", ifelse(grepl("AMP", SVTYPE), "DUP", "CNV")))

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
#' in the FACETS/mFACETS format
#'
#' @param filename (character) filename to store the results
#' @param cnvObj (cnvstrmliner) object with the CNV segments
#' @returns the filename of the output or NULL if there was a problem
#' @export

writemFACETSVCF <- function(filename, cnvObj) {
  filename <- filename %>%
    gsub(".gz$", "", .) %>%
    gsub(".bgz$", "", .) %>%
    gsub(".vcf$", "", .) %>%
    gsub(".bcf$", "", .)

  vcf_tmp <- tempfile(pattern = paste0(basename(filename), "."), tmpdir = dirname(filename), fileext = ".vcf")

  header <- make_header_mFACETS(cnvObj)
  body <- recordToVcf_mFACETS(cnvObj@segments)
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
