#' Apply CNAqc Quality Control Analysis to CNV Data
#'
#' @description
#' This function applies the CNAqc (Copy Number Alteration Quality Control) pipeline
#' to assess the quality of copy number variant (CNV) calls by comparing them against
#' somatic mutation data. The analysis validates CNV segments by checking if the
#' variant allele frequencies (VAFs) of somatic mutations match the expected VAF
#' distributions for the inferred karyotype states.
#'
#' The function performs comprehensive quality control including:
#' \itemize{
#'   \item Peak analysis to validate karyotype calls against mutation VAFs
#'   \item Per-chromosome quality scoring
#'   \item Cancer cell fraction (CCF) estimation for each segment
#'   \item Calculation of the genome-wide QC pass coefficient
#' }
#'
#' @param cnvObj A \code{cnvstrmliner} object containing CNV segments with purity
#'   and ploidy information. The object must have valid purity value; samples
#'   without purity estimates will be skipped with a warning.
#' @param somatic.mutations Either a character string path to a VCF file containing
#'   somatic mutations, or a data.frame with mutation data. If a data.frame, it
#'   must contain the columns required by CNAqc: \code{chr}, \code{from}, \code{to},
#'   \code{ref}, \code{alt}, \code{DP} (depth), \code{NV} (variant reads), and
#'   \code{VAF} (variant allele frequency). Use \code{\link{loadMutationsFromVCF}}
#'   to properly format VCF data.
#' @param clonal Logical. If \code{TRUE} (default), all CNV segments are treated
#'   as clonal (CCF = 1). If \code{FALSE}, the function attempts to model
#'   subclonal copy number alterations and estimate clonal composition.
#' @param n_bootstrap Integer. Number of bootstrap iterations used for the Bmix
#'   algorithm in peak analysis. Higher values increase precision but require
#'   more computation time. Default is 10.
#' @param min_absolute_karyotype_mutations Integer. Minimum number of mutations
#'   required to test a specific karyotype state. Karyotypes with fewer mutations
#'   are excluded from quality assessment. Default is 100. See
#'   \code{CNAqc::analyze_peaks} for details.
#' @param min_segments_for_overfragmentation Integer. Minimum number of segments
#'   required on a chromosome arm to perform overfragmentation testing. Arms with
#'   fewer segments are not tested. Default is 20.
#' @param calculate_chromosome_scores Logical apply the QC algorithm on each
#'   chromosome independently.
#' @return A \code{cnvstrmliner} object with updated quality metrics. The returned
#'   object contains the following additions:
#'   \describe{
#'     \item{\code{@segments$cnv.cnaqc.pass}}{Logical vector indicating which
#'       segments pass QC (one value per segment)}
#'     \item{\code{@quality$cnaqcobj}}{Complete CNAqc result object}
#'     \item{\code{@quality$cnaqcobj$score}}{Overall quality score measuring
#'       deviation from expected VAF distributions based on purity}
#'     \item{\code{@quality$cnaqcobj$peaks_analysis$matches}}{Data frame showing
#'       which karyotypes PASS or FAIL the VAF tolerance threshold}
#'     \item{\code{@quality$cnaqcobj$peaks_analysis$QC}}{Overall QC flag (PASS/FAIL)}
#'     \item{\code{@quality$cnaqcobj$qc_coeff}}{Fraction of the genome on segments
#'       that pass quality control (0 to 1)}
#'     \item{\code{@quality$cnaqcobj$scores_chr}}{Named vector of quality scores
#'       per chromosome}
#'     \item{\code{@quality$cnaqcobj$mutations_chr}}{Named vector of mutation
#'       counts used per chromosome}
#'     \item{\code{@quality$cnaqcobj$ccf_estimate}}{Data frame with cancer cell
#'       fraction estimates per karyotype}
#'     \item{\code{@quality$cnaqcobj$seed}}{Random seed used for reproducibility}
#'   }
#'
#' @details
#' The CNAqc algorithm validates CNV calls by leveraging the relationship between
#' copy number states and variant allele frequencies of heterozygous somatic mutations.
#' For a given karyotype (Major:minor alleles), the expected VAF can be predicted
#' based on sample purity and ploidy. The function tests whether observed VAF
#' distributions match theoretical expectations.
#'
#' The analysis is performed both genome-wide and per-chromosome using parallel
#' processing via the \code{future} framework. The current future plan is logged
#' during execution.
#'
#' Quality control failures may indicate:
#' \itemize{
#'   \item Incorrect purity or ploidy estimates
#'   \item Poor quality CNV calls
#'   \item Incorrect karyotype assignments
#'   \item Contamination or sample quality issues
#'   \item Presence of significant subclonal populations (when \code{clonal = TRUE})
#' }
#'
#' @note
#' This function requires the CNAqc package to be installed and loaded. Verbose
#' output from CNAqc is suppressed during execution.
#'
#' @seealso
#' \code{\link{loadMutationsFromVCF}} for loading and formatting VCF mutation data
#'
#' \code{\link{qcTable}} for generating a summary table of QC metrics
#'
#' \code{\link{detectOverfragmentation}} for additional fragmentation analysis
#'
#' \code{CNAqc::analyze_peaks} for details on peak analysis parameters
#'
#' @examples
#' \dontrun{
#' # Load CNV data
#' cnv <- parseFACETS("sample.vcf.gz")
#'
#' # Apply QC with VCF file
#' cnv <- applyQCMetrics(cnv,
#'   somatic.mutations = "somatic_mutations.vcf",
#'   clonal = TRUE,
#'   n_bootstrap = 10
#' )
#'
#' # Check overall QC status
#' cnv@quality$cnaqcobj$peaks_analysis$QC
#'
#' # View QC pass rate
#' cnv@quality$cnaqcobj$qc_coeff
#'
#' # Apply QC with pre-loaded mutation data
#' mutations <- loadMutationsFromVCF("mutations.vcf", "sample_01")
#' cnv <- applyQCMetrics(cnv, somatic.mutations = mutations)
#'
#' # Use with subclonal analysis
#' cnv <- applyQCMetrics(cnv,
#'   somatic.mutations = "mutations.vcf",
#'   clonal = FALSE,
#'   n_bootstrap = 20
#' )
#' }
#'
#' @export
#'
applyQCMetrics <- function(
  cnvObj,
  somatic.mutations,
  clonal = TRUE,
  n_bootstrap = 10,
  min_absolute_karyotype_mutations = 100,
  min_segments_for_overfragmentation = 20,
  calculate_chromosome_scores = TRUE
) {
  logger::log_threshold(logger::WARN, namespace = "applyqcmetrics")
  rlang::local_options(rlib_message_verbosity = "quiet")
  fp <- attr(future::plan(), "class")
  logger::log_tictoc("Future plan is set to {fp} ", level = logger::DEBUG, namespace = "applyqcmetrics")
  purity <- cnvObj@purity
  if (is.na(purity)) {
    logger::log_warn("Sample {cnvObj@sample} does not have purity. QC metrics cannot be calculated for it")
    return(cnvObj)
  }
  if (is(somatic.mutations, "character")) {
    mut <- cnvstrmliner::loadMutationsFromVCF(somatic.mutations, cnvObj@sample, genome = GenomeInfoDb::genome(cnvObj@segments) %>% unique())
  } else {
    mut <- somatic.mutations
  }

  cnvs <- segmentsCNAqc(cnvObj@segments, cnvObj@purity, cnvObj@effective_ploidy)
  if (clonal == TRUE) {
    # let's make everything clonal
    cnvs <- cnvs %>% dplyr::select(-c(CCF, Major_2, minor_2))
  }
  # Peak analysis
  # VAFs peaks are known for mutations mapped to a segment with given minor/ major allele copies. CNAqc therefore computes expected peaks, and compares them to peaks detected from data.
  logger::log_tictoc("Initializing the object", level = logger::INFO, namespace = "applyqcmetrics")
  x <- suppressMessages(
    CNAqc::init(
      mutations = mut,
      cna = cnvs,
      purity = purity,
      sample = cnvObj@sample,
      ref = "hg38"
    )
  )
  seed <- 12323
  set.seed(seed)

  if (calculate_chromosome_scores == TRUE) {
    logger::log_tictoc("Starting the processing of the object by chromosome", level = logger::INFO, namespace = "applyqcmetrics")
    sink(nullfile()) # do this to avoid excessive empty lines
    x_chr <- suppressMessages(suppressWarnings(
      CNAqc::split_by_chromosome(x)
    ))
    sink()
    f2 <- listenv::listenv()
    for (cc in 1:length(x_chr)) {
      f2[[cc]] <- future::future(
        {
          o <- tryCatch(
            {
              CNAqc::analyze_peaks(x_chr[[cc]],
                n_bootstrap = n_bootstrap,
                min_absolute_karyotype_mutations = min_absolute_karyotype_mutations
              )
            },
            error = function(e) {
              logger::log_warn("Error in QC chromosome {names(x_chr)[cc]}. Chromosome will be ignored !")
              return(NULL)
            }
          )
          if (is.null(o)) {
            return(list("score" = NA, "mutations" = 0))
          }
          muts <- o$mutations %>%
            tidyr::drop_na(QC_PASS) %>%
            nrow()
          score <- o$peaks_analysis$score
          score <- ifelse(is.null(score), NA, score)
          list("score" = score, "mutations" = muts)
        },
        seed = TRUE
      )
    }
    f2 <- as.list(f2)
  }
  logger::log_tictoc("Starting the processing of the full object", level = logger::INFO, namespace = "applyqcmetrics")

  f <- future::future(
    {
      o <- tryCatch(
        {
          x %>%
            CNAqc::analyze_peaks(., n_bootstrap = n_bootstrap, min_absolute_karyotype_mutations = min_absolute_karyotype_mutations) %>%
            CNAqc::compute_CCF(.)
        },
        error = function(e) {
          logger::log_warn("Error in QC. No peaks and CCF will be calculated !")
          return(x)
        }
      )
      o
    },
    seed = TRUE
  )
  sink(nullfile())
  suppressMessages(suppressWarnings({
    x <- future::value(f)
  }))
  if (calculate_chromosome_scores == TRUE) {
    suppressMessages(suppressWarnings({
      per_chr <- lapply(f2, function(X) {
        future::value(X)
      }) %>% dplyr::bind_rows()
    }))
  }
  sink()


  segs <- dplyr::bind_rows(x$cna, dplyr::tibble(QC_PASS = logical())) %>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE, start.field = "from", end.field = "to", ignore.strand = TRUE)
  if (!"QC_PASS" %in% colnames(x$cna)) {
    x$cna["QC_PASS"] <- "UNKNOWN"
  }
  genome.stats <- x$cna %>%
    dplyr::group_by(Major, minor, QC_PASS) %>%
    dplyr::summarise(width = sum(to - from)) %>%
    tidyr::drop_na(QC_PASS)
  pass <- genome.stats %>%
    dplyr::filter(QC_PASS == TRUE) %>%
    dplyr::pull(width) %>%
    sum(., na.rm = TRUE)
  qc_coeff <- pass / genome.stats %>%
    dplyr::pull(width) %>%
    sum(na.rm = TRUE)

  logger::log_tictoc("Processing is completed", level = logger::INFO, namespace = "applyqcmetrics")
  # qc = sapply(xlist2, function(x){x$peaks_analysis$QC})
  # qc=x$peaks_analysis$QC
  # qc.fraction=sum(qc=='PASS')/length( qc )
  cnvObj@segments$cnv.cnaqc.pass <- segs$QC_PASS
  ccf_estimates <- lapply(x$CCF_estimates, function(X) {
    X$QC_table
  }) %>% dplyr::bind_rows()
  # get per chromosome stats


  cnvObj@quality[["cnaqcobj"]] <- x # the initial object
  # cnvObj@quality[['cnaqc']][['qc.samples']]=qc
  # cnvObj@quality[['cnaqc']][['qc.coeff']]=qc.fraction
  cnvObj@quality[["cnaqcobj"]][["qc_coeff"]] <- qc_coeff
  cnvObj@quality[["cnaqcobj"]][["ccf_estimate"]] <- ccf_estimates
  cnvObj@quality[["cnaqcobj"]][["seed"]] <- seed
  if (calculate_chromosome_scores == TRUE) {
    cnvObj@quality[["cnaqcobj"]][["mutations_chr"]] <- per_chr$mutations %>% magrittr::set_names(names(x_chr))
    cnvObj@quality[["cnaqcobj"]][["scores_chr"]] <- per_chr$score %>% magrittr::set_names(names(x_chr))
  } else {
    cnvObj@quality[["cnaqcobj"]][["mutations_chr"]] <- NA
    cnvObj@quality[["cnaqcobj"]][["scores_chr"]] <- NA
  }
  rlang::local_options(rlib_message_verbosity = "verbose")
  cnvObj
}

#' Load mutations from a VCF file and prepare a CNAqc compatible dataframe.
#'
#' @param somatic.fn (character) VCF file with _somatic_ mutations
#' @param sample (character) the name of the sample to load information from the VCF file.
#' @param genome (character) the version of the genome
#' @param filter (logical) filter mutations, and keep only the mutations with FILTER == PASS (default: TRUE)
#' @returns a dataframe compatible with CNAqc
#' @description
#' Loads the mutations from a vcf file and returns a dataframe needed for CNAqc
#' @export
#'
loadMutationsFromVCF <- function(somatic.fn, sample, genome = "hg38", filter = TRUE) {
  # germline=normalizePath(file.path('mFACETS','vcf',manifest[idx,'germline']))
  # germline.vcf=readVcf( germline, genome='hg38')%>% keepStandardChromosomes( pruning.mode='coarse')
  vcf <- VariantAnnotation::readVcf(somatic.fn, genome = genome) %>% GenomeInfoDb::keepStandardChromosomes(pruning.mode = "coarse")
  if (isTRUE(filter)) {
    keep <- which(SummarizedExperiment::rowRanges(vcf)$FILTER == "PASS")
    logger::log_debug(
      "Keeping {length(keep)} of {nrow(vcf)} mutations with FILTER == PASS",
      namespace = "cnvstrmliner"
    )
    # subset the VCF itself so that rowRanges() and geno() stay aligned
    vcf <- vcf[keep, ]
    if (length(keep) == 0) {
      logger::log_warn(
        "No mutation in {basename(somatic.fn)} has FILTER == PASS. Use filter = FALSE to keep all of them.",
        namespace = "cnvstrmliner"
      )
      return(data.frame(
        chr = character(), from = integer(), to = integer(),
        ref = character(), alt = character(),
        DP = integer(), NV = integer(), VAF = numeric()
      ))
    }
  }
  rr <- SummarizedExperiment::rowRanges(vcf) %>%
    as.data.frame() %>%
    dplyr::select(chr = seqnames, from = start, to = end, ref = REF, alt = ALT) %>%
    tibble::rownames_to_column("id")
  rr$alt <- rr$alt %>% sapply(., function(x) {
    as.character(x[[1]])
  })
  ad <- VariantAnnotation::geno(vcf)$AD[, sample] %>%
    magrittr::set_names(NULL)
  ad <- lapply(ad, function(x) {
    x[1:2]
  }) %>%
    do.call(rbind, .) %>%
    magrittr::set_colnames(c("ref", "alt")) %>%
    magrittr::set_rownames(rr$id) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("id") %>%
    dplyr::mutate(DP = alt + ref, NV = alt, VAF = NV / DP) %>%
    dplyr::select(id, DP, NV, VAF)

  mutations.somatic <- dplyr::inner_join(rr, ad, dplyr::join_by("id")) %>%
    dplyr::select(-id) %>%
    dplyr::mutate(chr = as.character(chr)) %>%
    dplyr::filter(ref %in% c("A", "C", "T", "G"), alt %in% c("A", "C", "T", "G")) %>%
    dplyr::filter(VAF > 0.05)
  mutations.somatic
}

#' Prepare CNV object
#' @param segments A GRanges object with the segments of an object.
#' @param purity Numeric purity of the sample.
#' @param ploidy Numeric ploidy of the sample. Default is 2.
#' @return A data frame formatted for CNAqc.
#' @description
#' Create a data frame with the CNV information needed for CNAqc
#'
#'
segmentsCNAqc <- function(segments, purity, ploidy = 2) {
  res <- segments %>%
    GenomeInfoDb::keepStandardChromosomes(pruning.mode = "coarse") %>%
    as.data.frame() %>%
    dplyr::select(
      "chr" = seqnames,
      "from" = start,
      "to" = end,
      "Major" = major.allele.copies,
      "minor" = minor.allele.copies,
      "CCF" = cnv.ccf
    ) %>%
    # dplyr::filter( ! (Major==minor & CCF ==1) ) %>%
    # dplyr::filter( CCF > 0 )%>%
    dplyr::mutate(
      chr = as.character(chr),
      "Major_2" = ceiling(ploidy / 2),
      "minor_2" = floor(ploidy / 2)
    )

  # sclust and BB
  if ("Clone1_A" %in% colnames(S4Vectors::mcols(segments))) {
    if ("Clone1_Fraction" %in% colnames(S4Vectors::mcols(segments))) {
      S4Vectors::mcols(segments)["Clone1_ccf"] <- S4Vectors::mcols(segments)["Clone1_Fraction"]
    }
    res <- segments %>%
      GenomeInfoDb::keepStandardChromosomes(pruning.mode = "coarse") %>%
      as.data.frame() %>%
      dplyr::mutate(
        "chr" = seqnames,
        "from" = start,
        "to" = end,
        "Major" = ifelse(!is.na(Clone1_A), Clone1_A, A),
        "minor" = ifelse(!is.na(Clone1_B), Clone1_B, B),
        "CCF" = ifelse(!is.na(Clone1_ccf), Clone1_ccf, 1),
        "Major_2" = ifelse(!is.na(Clone2_A), Clone2_A, NA),
        "minor_2" = ifelse(!is.na(Clone2_B), Clone2_B, NA)
      ) %>%
      dplyr::select(chr, from, to, Major, minor, CCF, Major_2, minor_2)
  }
  res
}

#' Generate a CNAqc Quality Control Summary Table
#'
#' Creates a comprehensive summary table containing CNAqc quality control metrics,
#' including karyotype quality checks and chromosome arm fragmentation analysis.
#' This function extracts and combines key QC information from a cnvstrmliner object
#' that has been processed with CNAqc analysis.
#'
#' @param cnvObj A cnvstrmliner object containing CNV segments and CNAqc results.
#'   The object must have been processed with \code{\link{applyQCMetrics}} or
#'   similar function that populates the CNAqc quality metrics.
#'
#' @return A data frame (tibble) containing:
#'   \itemize{
#'     \item \code{sample} - Sample identifier
#'     \item \code{method} - CNV calling method used (e.g., FACETS, Sequenza)
#'     \item \code{CNAqc_flag} - Overall CNAqc quality control flag (PASS/FAIL)
#'     \item \code{CNAqc_score} - Numerical quality score from CNAqc analysis
#'     \item Karyotype-specific QC results for standard karyotypes
#'     \item Chromosome arm fragmentation significance flags (one column per arm)
#'   }
#'
#' @details
#' The function requires that CNAqc analysis has been previously run on the object.
#' It combines two main quality control components:
#' \enumerate{
#'   \item Karyotype matching results - Shows which copy number states pass QC
#'   \item Chromosome arm fragmentation analysis - Indicates significant fragmentation per arm
#' }
#'
#' The resulting table provides a comprehensive view of the quality of CNV calls,
#' helping to identify samples with unreliable copy number profiles.
#'
#' @seealso \code{\link{applyQCMetrics}} for running CNAqc analysis
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Load CNV data and apply QC metrics
#' cnv <- parseFACETS("sample.vcf.gz")
#' cnv <- applyQCMetrics(cnv, somatic.mutations = "mutations.vcf")
#' cnv <- detectOverfragmentation(cnv)
#'
#' # Generate QC summary table
#' qc_summary <- qcTable(cnv)
#' print(qc_summary)
#' }
qcTable <- function(cnvObj) {
  cnaqc_metrics <- NA
  cnaqc_chr_metrics <- NA
  frag_metrics <- NA
  if ("cnaqcobj" %in% names(cnvObj@quality)) {
    if (!is.null(cnvObj@quality$cnaqcobj$peaks_analysis$matches)) {
      cnaqc_metrics <- cnvObj@quality$cnaqcobj$peaks_analysis$matches %>%
        dplyr::select(karyotype, QC) %>%
        dplyr::group_by(karyotype) %>%
        dplyr::summarise(significant = paste(unique(QC), collapse = "/"), .groups = "drop") %>%
        dplyr::mutate(karyotype = paste0("karyotype_score_", karyotype)) %>%
        tibble::column_to_rownames("karyotype") %>%
        t()
    }
    if (any(!is.na(cnvObj@quality$cnaqcobj$scores_chr))) {
      cnaqc_chr_metrics <- cnvObj@quality$cnaqcobj$scores_chr %>%
        as.data.frame() %>%
        magrittr::set_rownames(paste0("chromosome_score_", rownames(.))) %>%
        t()
    }
  }

  if ("fragmentation" %in% names(cnvObj@quality)) {
    frag_metrics <- cnvObj@quality$fragmentation$table %>%
      as.data.frame() %>%
      dplyr::select("seqnames", "arm", "significant") %>%
      dplyr::mutate(seqnames = paste(arm, seqnames, sep = ".")) %>%
      dplyr::select(-arm) %>%
      dplyr::mutate(seqnames = paste0("fragmentation_", seqnames)) %>%
      tibble::column_to_rownames("seqnames") %>%
      t()
  }

  dplyr::bind_cols(
    cnaqc_metrics,
    frag_metrics,
    cnaqc_chr_metrics
  ) %>%
    purrr::discard(~ is.na(.x)) %>%
    dplyr::mutate(
      sample = cnvObj@sample,
      method = cnvObj@method,
      segments = cnvObj@segments %>% length(),
      purity = cnvObj@purity,
      ploidy = cnvObj@ploidy,
      effective_ploidy = cnvObj@effective_ploidy,
      CNAqc_flag = cnvObj@quality$cnaqcobj$peaks_analysis$QC,
      CNAqc_score = cnvObj@quality$cnaqcobj$peaks_analysis$score,
      entropy = ifelse(!is.null(cnvObj@quality$entropy), cnvObj@quality$entropy, NA),
      gini = ifelse(!is.null(cnvObj@quality$gini), cnvObj@quality$gini, NA),
      mapd = ifelse(!is.null(cnvObj@quality$mapd), cnvObj@quality$mapd, NA),
      .before = 1
    )
}
