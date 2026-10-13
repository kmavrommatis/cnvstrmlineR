#' check the validity of the inputs for the cnvstrml object
#'
#' @param object the cnvstrml object to check
#' @return TRUE or errors
check_cnvstrml <- function(object) {
  errors <- character()
  object@sample <- as.character(object@sample)
  if (is.na(object@purity)) {
    object@purity <- 0
  }
  if (is.null(object@purity)) {
    object@purity <- 0
  }
  if (is.na(object@ploidy)) {
    object@ploidy <- 2
  }
  if (is.null(object@ploidy)) {
    object@ploidy <- 2
  }
  if (is.null(object@effective_ploidy)) {
    object@effective_ploidy <- round(object@ploidy)
  }
  if (object@purity < 0) {
    msg <- "Purity cannot be a negative number. It has to be between 0 and 1"
    errors <- c(errors, msg)
  }
  if (object@purity > 1) {
    msg <- "Purity cannot be a greater than 1. It has to be between 0 and 1"
    errors <- c(errors, msg)
  }
  if (object@ploidy < 0) {
    msg <- "Ploidy cannot be negative or 0. It has to be a positive number"
    errors <- c(errors, msg)
  }
  if (!object@genome %in% c("hg38", "hg19")) {
    msg <- "The genome should be either hg19 or hg38"
    errors <- c(errors, msg)
  }

  if (!is.na(object@sex)) {
    object@sex <- toupper(as.character(object@sex))

    if (!object@sex %in% c("XY", "XX", "M", "F")) {
      msg <- paste("Not a valid sex [", object@sex, "]. Sex can be either XY,M for male, or XX, F for female. Setting it to U") # errors=c(errors,msg)
      object@sex <- "U"
      logger::log_warn(msg)
    }
    object@sex_inferred <- FALSE
  }


  # if( abs( obj@ploidy - round(obj@ploidy)) > 0.35){
  #   message("The ploidy ", obj@ploidy," is ambiguous")
  #   newPloidy=cnvstrmliner::predictPloidy( obj@segments$cnv.lrr,
  #                                          obj@segments$total.allele.copies)
  #   message("it is more likely to be ", newPloidy)
  #
  #   flags=apply(mcols(obj@segments), 1, function(K){ cnvstrmliner::annotateCNV(major.cn = K['major.allele.copies'], minor.cn = K['minor.allele.copies'] ,
  #                                                                              ploidy=newPloidy)})
  #   loh=cnvstrmliner::annotateLOH( flags)
  #   mcols(obj@segments)[,'cnv.annotation']=flags
  #   mcols(obj@segments)[,'cnv.loh']=loh
  # }

  if (is.na(object@label)) {
    object@label <- paste0(object@sample, " (", object@method, ")")
  }

  if (length(errors) == 0) TRUE else errors
}


# setOldClass("GenomicRanges::GRanges")
#' create the a CNVstrml object
#'
#' cnvstrml is a structure that holds the information for the CNV predictions
#' and some additional information provided by the methods
#'
#' @slot segments set of CNV segments (GRanges) as predicted by a CNV calling method
#' @slot purity (real) the purity of the sample (default 1)
#' @slot ploidy (real) the ploidy of the sample (default 2)
#' @slot effective_ploidy (integer) the ploidy used to calculate the type of CNV (default 2)
#' @slot genome (character) the name of the genome to use (hg38, hg19)
#' @slot sample (character) name of sample
#' @slot method (character) the name of the method used to produce the data
#' @slot sex (character) the sex of the sample ('M', 'F')
#' @slot sex_inferred (logical) is the sex inferred or provided by the method
#' @slot segtype (character) the type of data stored in this object
#' @slot label (character) the label for this sample (default is sample + method)
#' @slot quality (list) metric(s) of the quality of segmentation (if provided by the methods)
#' @docType class
#' @export
#' @section Notes:
#' This is the main object provided by this package.
#'
#' Contains information about the CNV segments
#' the purity and ploidy, and effective ploidy of the samples,
#' and metadata about the method name, sample name and the
#' reference genome used
#'
#' The object contains segments which are produced by the CNV calling method and potentially
#' transformed using the functions within this package.
#'
#' The segment information of the object is ammended with the following  columns:
#'
#' `total.allele.copies`: integer of total allele copies for a segment
#' `major.allele.copies`: integer of major allele copies for a segment
#' `minor.allele.copies`: integer of minor allele copies for a segment
#' `cnv.annotation`     : A controlled list of annotation terms (e.g. NEUT, AMP, DEL). See [annotateCNV] for more information
#' `cnv.lrr`            : logR of this segment (0 == NEUTRAL)
#' `cnv.markers`        : Humber of heterozygous markers used in this segment
#' `cnv.totalmarkers`   : Number of markers (SNPs) included in this segment
#' `cnv.cf`            : Cancer cell fraction carrying this CNV. If not provided it is set to the purity of the sample
#' `cnv.loh`            : Is this a segment with Loss Of Heterozygosity ?
#' `cnv.markers_density`: Average density of cnv.markers (cnv.markers/width of segment)
#' `cnv.totalmarkers_density` : Average density of total markers
#' `cnv.fill`           : Was this segment included in the original prediction or added by the fillGaps function?
#'
#' When a new object is created the segments should contain the columns listed above
#' The cnv.annotation is added by the function  `cnvstrmliner::annotateCNV`
#' The effective ploidy is added by the function `cnvstrmliner::calculateEffectivPloidy`
#' cnv.loh is calculated automatically
#'
cnvstrml <- setClass("cnvstrml",
  # contains=c("GenomicRanges::GRanges"),
  slots = c(
    segments = "GRanges",
    purity = "numeric",
    ploidy = "numeric",
    effective_ploidy = "numeric",
    method = "character",
    sample = "character",
    genome = "character",
    sex = "character",
    sex_inferred = "logical",
    segtype = "character",
    label = "character",
    quality = "list"
  ),
  prototype = list(
    purity = 1,
    ploidy = 2,
    effective_ploidy = NA_real_,
    method = NA_character_,
    sample = NA_character_,
    genome = "hg38",
    sex = NA_character_,
    sex_inferred = TRUE,
    segtype = "single sample",
    label = NA_character_,
    quality = list()
  ),
  validity = check_cnvstrml
)


#' initialize the object
#'
#' check that purity is within limits
#' updates the sex based on the presence of chrY events if the user has not prodived that
#' adds the cnv.loh flag based on the number of minor alleles if it is not already present
#' adds the cnv.markers_density and cnv.totalmarkers_density
#' adds the effective_ploidy if it is not provided.
#'
#' @param .Object an .Object of class cnvstrml
#' @param ... additional arguments
#' @return an initialized .Objectect
#'
setMethod(
  "initialize", "cnvstrml",
  function(.Object, ...) {
    .Object <- callNextMethod()
    # callNextMethod()
    if (startsWith(.Object@segtype, "consensus")) {
      logger::log_info("This is a Consensus intermediate object")
      return(.Object)
    }

    # check if the proper columns are in teh segments
    check_for_columns <- c(
      "total.allele.copies", "major.allele.copies", "minor.allele.copies",
      "cnv.lrr", "cnv.markers", "cnv.totalmarkers", "cnv.cf"
    )
    if (!setequal(intersect(check_for_columns, colnames(S4Vectors::mcols(.Object@segments))), check_for_columns)) {
      logger::log_error("Object: ", .Object@sample, " ", .Object@segtype)
      logger::log_error("Missing columns from segment.")
      logger::log_error("The segments must contain the following columns: ", paste(check_for_columns, collapse = ","))
      logger::log_error("The current columns are: ", paste(colnames(S4Vectors::mcols(.Object@segments)), collapse = ","))
      stop("The segments must contain the following columns: ", paste(check_for_columns, collapse = ", "), ". Please check your input.")
    }

    # check and infer the sex
    logger::log_debug("checking (and if necessary inferring) sex")
    if (is.na(.Object@sex)) {
      logger::log_info("Predicting the sex of the sample")
      sex <- inferSex(.Object)
      .Object@sex <- toupper(sex)
    } else {
      logger::log_debug("Sex was provided with the parsed data")
    }

    # the effective ploidy has not been calculated
    # calculate it
    if (is.na(.Object@effective_ploidy)) {
      logger::log_info("Calculating the effective ploidy.")
      .Object@effective_ploidy <- calculateEffectivePloidy(.Object@segments)
    } else {
      logger::log_debug("User provided the effective ploidy")
    }


    if (!"cnv.annotation" %in% colnames(S4Vectors::mcols(.Object@segments)) | any(is.na(S4Vectors::mcols(.Object@segments)$cnv.annotation))) {
      logger::log_info("Adding annotation of cnv segments using {.Object@effective_ploidy} as the ploidy.")
      flags <- mapply(
        function(total, minor) {
          cnvstrmliner::annotateCNV(
            total.cn = total,
            minor.cn = minor,
            ploidy = .Object@effective_ploidy
          )
        },
        S4Vectors::mcols(.Object@segments)[, "total.allele.copies"],
        S4Vectors::mcols(.Object@segments)[, "minor.allele.copies"]
      )
      S4Vectors::mcols(.Object@segments)["cnv.annotation"] <- flags
    } else {
      logger::log_debug("User provided CNV annotation information in the cnv.annotation column")
    }


    # adjust sex segments. Remove any chrY segments if the sex is set to Female
    if (getSex(.Object) == "F") {
      .Object@segments[as.character(GenomeInfoDb::seqnames(.Object@segments)) %in% c("Y", "chrY")] <- NULL
    }

    # adjust the annotation (cnv.annotation) of the sex chromosomes

    phyletic <- which(as.character(GenomeInfoDb::seqnames(.Object@segments)) %in% c("chrX", "X", "chrY", "Y"))
    if (length(phyletic) > 0) {
      phyletic_ploidy <- ifelse(getSex(.Object) == "M", 1, 2) # ploidy for chrY is 1, chrX is 1 or 2 depending on sex
      logger::log_info("Adjusting annotation of segments located on phyletic chromosomes with ploidy {phyletic_ploidy}")
      flags <- mapply(
        function(total, minor) {
          cnvstrmliner::annotateCNV(
            total.cn = total,
            minor.cn = minor,
            ploidy = phyletic_ploidy
          )
        },
        S4Vectors::mcols(.Object@segments)[phyletic, "total.allele.copies"],
        S4Vectors::mcols(.Object@segments)[phyletic, "minor.allele.copies"]
      )
      S4Vectors::mcols(.Object@segments)$cnv.annotation[phyletic] <- flags
    } else {
      logger::log_info("There are no phyletic chromosomes in this sample")
    }

    # annotate the loh if the user has not already added it
    if (!"cnv.loh" %in% colnames(S4Vectors::mcols(.Object@segments)) | any(is.na(S4Vectors::mcols(.Object@segments)$cnv.loh))) {
      logger::log_info("Adding annotation of the LOH")
      loh <- mapply(
        function(chr, sex, minor, major) {
          cnvstrmliner::addLOH(chromosome = chr, sex = sex, minor.cn = minor, major.cn = major)
        },
        as.character(GenomeInfoDb::seqnames(.Object@segments)),
        .Object@sex,
        S4Vectors::mcols(.Object@segments)[, "minor.allele.copies"],
        S4Vectors::mcols(.Object@segments)[, "major.allele.copies"]
      )
      S4Vectors::mcols(.Object@segments)["cnv.loh"] <- loh
    } else {
      logger::log_debug("User provided LOH information in the cnv.loh column")
    }


    # annotate cancer cell fraction
    if (!"cnv.ccf" %in% colnames(S4Vectors::mcols(.Object@segments)) | any(is.na(S4Vectors::mcols(.Object@segments)$cnv.ccf))) {
      logger::log_info("Adding annotation for Cancer cell fraction (CCF).")
      ccf <- sapply(
        S4Vectors::mcols(.Object@segments)[, "cnv.cf"],
        function(cf, purity) {
          addCCF(cf, .Object@purity)
        }
      )
      S4Vectors::mcols(.Object@segments)[, "cnv.ccf"] <- ccf
    } else {
      logger::log_debug("User provided Cancer cell fraction(CCF) information in the cnv.ccf column")
    }

    # add the density of markers

    if ("cnv.markers" %in% colnames(S4Vectors::mcols(.Object@segments)) | any(is.na(S4Vectors::mcols(.Object@segments)$cnv.markers))) {
      logger::log_info("Adding the density of het markers")
      S4Vectors::mcols(.Object@segments)["cnv.markers_density"] <-
        S4Vectors::mcols(.Object@segments)[, "cnv.markers"] / BiocGenerics::width(.Object@segments)
    }
    if ("cnv.totalmarkers" %in% colnames(S4Vectors::mcols(.Object@segments)) | any(is.na(S4Vectors::mcols(.Object@segments)$cnv.totalmarkers))) {
      logger::log_info("Adding the density of all markers")
      S4Vectors::mcols(.Object@segments)["cnv.totalmarkers_density"] <-
        S4Vectors::mcols(.Object@segments)[, "cnv.totalmarkers"] / BiocGenerics::width(.Object@segments)
    }

    if (!"arm" %in% colnames(S4Vectors::mcols(.Object@segments)) | any(is.na(S4Vectors::mcols(.Object@segments)$arm))) {
      logger::log_info("Adding chromosome arm information")
      arms <- cnvstrmliner::addChromosomeAnnotation(.Object@segments, genome = .Object@genome)
      .Object@segments <- arms
    }

    invisible(.Object)
  }
)


#' show the contents of the object
#'
#' @param object Object of type cnvstrml
#'
#'
setMethod(
  "show", "cnvstrml",
  function(object) {
    cat("Contents of sample ", object@sample, "\n")
    cat(paste("Method used to produce data ", object@method, "\n"))
    cat(paste("Purity is set to ", object@purity, "\n"))
    cat(paste("Ploidy is set to ", object@ploidy, " \n"))
    cat(paste("Effective Ploidy is set to ", object@effective_ploidy, "\n"))
    cat(paste("Sex is set to ", object@sex, "(inferred: ", object@sex_inferred, ")\n"))
    if (length(object@quality) > 0) {
      nn <- names(object@quality)
      nn <- setdiff(nn, c("cnaqcobj", "fragmentation", "file"))
      cat("Quality of segmentation (parsed from the relevant files):\n")
      print(
        # object@quality %>% as.data.frame() %>% t()
        object@quality[nn]
      )
      if ("file" %in% names(object@quality)) {
        cat("Parsed the following file(s) \n")
        for (i in seq_along(object@quality$file)) {
          cat(names(object@quality$file)[i], ": ", object@quality$file[i], "\n", sep = "")
        }
      }
      if ("cnaqcobj" %in% names(object@quality)) {
        cat(paste(
          "CNAqc has been run on this sample with\n",
          "\tscore          : ", object@quality$cnaqcobj$peaks_analysis$score, "\n",
          "\tqc             : ", object@quality$cnaqcobj$peaks_analysis$QC, "\n"
        ))
      }
      if ("fragmentation" %in% names(object@quality)) {
        cat(paste(
          "Chromosome arm fragmentation has been run on this sample with\n",
          "\t", ifelse(is.na(object@quality$fragmentation$small_segment_size),
            paste("arm percentage cutoff for short segments ", object@quality$fragmentation$arm_size_percentage_cutoff),
            paste("small segment size shorter than ", object@quality$fragmentation$small_segment_size)
          ), "\n",
          "\t with minimum number of segments to test ", object@quality$fragmentation$minimum_segments_for_testing, "\n",
          "\tfragmented arms: ", object@quality$fragmentation$table %>% as.data.frame() %>% dplyr::filter(significant == TRUE) %>% nrow(), "\n"
        ))
      }
    }
    cat(paste(length(object@segments), " segments\n"))
  }
)


#' Get Supported CNV Calling Methods
#'
#' Returns a character vector of CNV calling methods that are supported by the
#' cnvstrmliner package for parsing and analysis.
#'
#' @return A character vector containing the names of supported CNV calling methods.
#'
#' @export
#'
#' @examples
#' # Get list of supported methods
#' methods <- getSupportedMethods()
#' print(methods)
#'
#' # Check if a method is supported
#' "FACETS" %in% getSupportedMethods()
getSupportedMethods <- function() {
  supported_methods <- c(
    "FACETS", "Sequenza", "SCLUST", "Battenberg",
    "FREEC", "CNVkit", "PureCN", "GATK",
    "Conseca",
    "CloneCNA",
    "TSV", "VCF", "SEG"
  )
  return(supported_methods)
}


#' Generate comprehensive summary statistics of copy number aberrations
#'
#' @description
#' Summarizes copy number aberrations by aggregating segments across three genomic levels
#' (genome, chromosome, and chromosomal arm). For each level, segments are grouped by
#' karyotype (major:minor allele copy numbers) and CNV annotation, and reported with the
#' number of segments, cumulative base-pair width, and normalized percentage of the
#' reference region affected.
#'
#' @param cnvObj A \code{cnvstrml} object containing CNV segmentation data with
#'   annotated segments including allelic composition (major and minor allele copies)
#'   and CNV classification. The object must have a valid genome build specified.
#'
#' @param level Character string specifying the summarization level. Must be one of:
#'   \describe{
#'     \item{genome}{Summarize across the entire genome (default)}
#'     \item{chromosome}{Summarize separately for each chromosome}
#'     \item{arm}{Summarize separately for each chromosomal arm (p and q)}
#'   }
#'
#' @return A data frame with columns depending on the specified \code{level}:
#'
#'   \strong{For level = "arm":}
#'   \describe{
#'     \item{seqnames}{Chromosome name (e.g., chr1, chr2)}
#'     \item{arm}{Chromosomal arm (p or q)}
#'     \item{karyotype}{Allelic composition as \code{major:minor} (e.g., \code{2:1})}
#'     \item{cnv.annotation}{CNV classification (AMP, DEL, NEUT, etc.)}
#'     \item{count}{Number of segments with this karyotype/annotation}
#'     \item{total_width}{Total genomic span in base pairs}
#'     \item{arm_length}{Total length of the chromosomal arm}
#'     \item{arm_percentage}{Percentage of arm covered (total_width/arm_length × 100)}
#'   }
#'
#'   \strong{For level = "chromosome":}
#'   \describe{
#'     \item{seqnames}{Chromosome name}
#'     \item{karyotype}{Allelic composition}
#'     \item{cnv.annotation}{CNV classification}
#'     \item{count}{Number of segments}
#'     \item{total_width}{Total genomic span in base pairs}
#'     \item{chr_length}{Total length of the chromosome}
#'     \item{chr_percentage}{Percentage of chromosome covered}
#'   }
#'
#'   \strong{For level = "genome":}
#'   \describe{
#'     \item{karyotype}{Allelic composition}
#'     \item{cnv.annotation}{CNV classification}
#'     \item{count}{Number of segments}
#'     \item{total_width}{Total genomic span in base pairs}
#'     \item{genome_length}{Total length of the genome}
#'     \item{genome_percentage}{Percentage of genome covered}
#'   }
#'
#'   Rows are sorted by descending total_width within each grouping, showing the most
#'   prevalent aberrations first.
#'
#' @details
#' \strong{Purpose and Applications}
#' \itemize{
#'   \item Quality assessment: identify samples with unusual CNV patterns
#'   \item Comparative analysis: compare CNV burden across samples or conditions
#'   \item Biological interpretation: determine which chromosomes/arms are affected
#'   \item Reporting: generate standardized metrics for publications or clinical reports
#'   \item Downstream analysis: provide input for burden analysis and correlations
#' }
#'
#' \strong{Key Metrics Explained}
#' \describe{
#'   \item{Karyotype}{Allelic composition in the format \code{major:minor}.
#'     Examples: \code{1:1} (diploid), \code{2:0} (copy-neutral LOH), \code{3:1}
#'     (gain with heterozygosity), \code{2:1} (single-copy gain), \code{1:0}
#'     (heterozygous deletion), \code{0:0} (homozygous deletion). Karyotype captures
#'     both total copy number and allelic imbalance.}
#'
#'   \item{CNV Annotation}{Biological classification of the aberration, such as
#'     AMP, AMP-LOH, NEUT, NEUT-LOH, HETDEL, HOMDEL, DEL, and INFER.}
#'
#'   \item{Count}{Number of distinct segments with the specified karyotype and
#'     annotation. High count suggests fragmentation or many focal events; low count suggests
#'     broad events or stable genomes. Interpret with total_width.}
#'
#'   \item{Total Width}{Cumulative base-pair size of all segments with the specified
#'     karyotype and annotation. Large values indicate substantial genomic territory affected;
#'     small values indicate focal aberrations. Useful for burden comparison across samples.}
#'
#'   \item{Percentage}{Proportion of the reference region (genome/chromosome/arm)
#'     affected: total_width divided by the relevant length times 100. Enables normalized
#'     comparisons across regions of different sizes and serves as a QC indicator.}
#'
#'   \item{Reference Length}{Total size of the reference region summarized
#'     (genome_length, chr_length, or arm_length). Provides context for absolute sizes.}
#' }
#'
#' \strong{Choosing the Summarization Level}
#' \describe{
#'   \item{genome}{Global view of CNV burden across the entire genome.
#'     Use for overall quality assessment, burden comparison, and outlier detection.
#'     Example: This sample has 15\% of genome amplified and 8\% deleted.
#' }
#'   \item{chromosome}{Aberration patterns per chromosome. Use to identify
#'     chromosome-specific aneuploidy and complex rearrangements. Example: Chr7 shows
#'     95\% amplification, Chr10 shows 90\% deletion.
#' }
#'   \item{arm}{Finest resolution per chromosomal arm (p and q). Use to detect
#'     arm-level events and arm-specific LOH. Example: 17p shows 100\% deletion (TP53 loss),
#'     while 17q is neutral.
#' }
#' }
#'
#' \strong{Interpreting Results in CNV Context}
#' \itemize{
#'   \item High percentage neutral (>80\%) genome-wide: normal/low-aberration sample
#'   \item Low percentage neutral (<30\%) genome-wide: highly rearranged genome or QC issue
#'   \item High count with low total_width: many focal events (e.g., ERBB2, MYC amplifications)
#'   \item Low count with high total_width: broad arm/chromosome events (aneuploidy)
#'   \item High NEUT-LOH percentage: copy-neutral LOH suggesting selective pressure
#'   \item Uneven arm patterns: centromeric breakpoints
#'   \item Consistent karyotype over large regions: early, stable clonal events
#'   \item Multiple karyotypes in same region: subclonality or over-segmentation
#' }
#'
#' @seealso
#' \code{\link{cnvMetrics}} for comparing CNV profiles across samples
#' \code{\link{cnvEntropy}} for assessing CNV complexity and fragmentation
#'
#' @examples
#' \dontrun{
#' # Genome-wide summary: overall burden assessment
#' genome_summary <- summaryTable(cnv_object, level = "genome")
#' # Check what % of genome is amplified
#' amp_burden <- sum(genome_summary$genome_percentage[genome_summary$cnv.annotation == "AMP"])
#'
#' # Chromosome-level summary: identify aneuploidies
#' chr_summary <- summaryTable(cnv_object, level = "chromosome")
#' # Find chromosomes with >80% aberration
#' aneuploidy <- chr_summary %>%
#'   filter(chr_percentage > 80, cnv.annotation != "NEUT")
#'
#' # Arm-level summary: detect arm-level events
#' arm_summary <- summaryTable(cnv_object, level = "arm")
#' # Classic 17p deletion in TP53
#' tp53_loss <- arm_summary %>%
#'   filter(seqnames == "chr17", arm == "p", cnv.annotation %in% c("DEL", "HETDEL"))
#'
#' # Identify most prevalent aberration pattern
#' top_pattern <- arm_summary %>%
#'   filter(cnv.annotation != "NEUT") %>%
#'   arrange(desc(total_width)) %>%
#'   head(5)
#' }
#'
#' @export

summaryTable <- function(cnvObj, level = c("genome", "chromosome", "arm")) {
  level <- match.arg(level)
  arms_lengths <- cnvstrmliner::chrom_arms[[cnvObj@genome]] %>% as.data.frame()
  chr_lengths <- arms_lengths %>%
    dplyr::group_by(seqnames) %>%
    dplyr::summarise(width = sum(width), .groups = c("keep")) %>%
    dplyr::ungroup()
  genome_length <- sum(arms_lengths$width)


  # calculate the summary table per chromosome arm
  df_raw <- as.data.frame(cnvObj@segments) %>%
    dplyr::select(
      seqnames, arm, width,
      major.allele.copies, minor.allele.copies,
      cnv.annotation
    ) %>%
    dplyr::filter(arm %in% c("p", "q")) %>%
    dplyr::mutate(karyotype = paste(round(major.allele.copies), round(minor.allele.copies), sep = ":")) %>%
    dplyr::group_by(seqnames, arm, karyotype, cnv.annotation) %>%
    dplyr::summarise(
      count = dplyr::n(),
      total_width = sum(width),
      .groups = c("keep")
    ) %>%
    dplyr::arrange(seqnames, arm, desc(total_width)) %>%
    dplyr::ungroup()
  if (level == "arm") {
    df <- df_raw %>%
      dplyr::left_join(arms_lengths %>% dplyr::select(seqnames, arm, arm_length = width),
        by = c("seqnames", "arm")
      ) %>%
      dplyr::mutate(arm_percentage = total_width / arm_length * 100) %>%
      dplyr::select(seqnames, arm, karyotype, cnv.annotation, count, total_width, arm_length, arm_percentage)
  } else if (level == "chromosome") {
    df <- df_raw %>%
      dplyr::group_by(seqnames, karyotype, cnv.annotation) %>%
      dplyr::summarise(
        count = sum(count),
        total_width = sum(total_width),
        .groups = c("keep")
      ) %>%
      dplyr::left_join(chr_lengths %>% dplyr::select(seqnames, chr_length = width),
        by = c("seqnames")
      ) %>%
      dplyr::mutate(chr_percentage = total_width / chr_length * 100) %>%
      dplyr::select(seqnames, karyotype, cnv.annotation, count, total_width, chr_length, chr_percentage) %>%
      dplyr::arrange(seqnames, desc(total_width)) %>%
      dplyr::ungroup()
  } else if (level == "genome") {
    df <- df_raw %>%
      dplyr::group_by(karyotype, cnv.annotation) %>%
      dplyr::summarise(
        count = sum(count),
        total_width = sum(total_width),
        .groups = c("keep")
      ) %>%
      dplyr::mutate(genome_length = genome_length) %>%
      dplyr::mutate(genome_percentage = total_width / genome_length * 100) %>%
      dplyr::select(karyotype, cnv.annotation, count, total_width, genome_length, genome_percentage) %>%
      dplyr::arrange(desc(total_width)) %>%
      dplyr::ungroup()
  }
  return(as.data.frame(df))
}
