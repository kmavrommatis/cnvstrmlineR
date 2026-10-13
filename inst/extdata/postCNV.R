#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(logger)
  library(dplyr)
  library(cnvstrmliner)
  library(VariantAnnotation)
})

option_list <- list(
  make_option(c("-i", "--input"),
    type = "character", default = NULL,
    help = "Path to the input CNV file",
    metavar = "character"
  ),
  make_option(c("-m", "--method"),
    type = "character", default = NULL,
    help = paste(
      "Method used to produce CNV results",
      paste(cnvstrmliner::getSupportedMethods(), collapse = "|")
    ),
    metavar = "character"
  ),
  make_option(c("-v", "--somatic_mutations"),
    type = "character",
    default = NULL,
    help = "somatic mutations VCF",
    metavar = "character"
  ),
  make_option(c("-c", "--logr"),
    type = "character", default = NULL,
    help = "File with logR values per target for Gini and MAPD calculation",
    metavar = "character",
    dest = "cnr_file"
  ),
  make_option(c("-d", "--bed"),
    type = "character", default = NULL,
    help = "BED file with covered regions for CNV entropy calculation and Gini and MAPD",
    metavar = "character",
    dest = "bed_file"
  ),
  make_option(c("-o", "--output"),
    type = "character", default = NULL,
    help = "Path to the output VCF file",
    metavar = "character"
  ),
  make_option(c("-s", "--sample"),
    type = "character",
    default = NA,
    help = "Sample name (optional)",
    metavar = "character",
    dest = "sample"
  )
)

parser <- OptionParser(
  option_list = option_list,
  description = "This is the workhorse script for post processing of CNV calls.
It can :
- Format CNV calls from various callers to VCF format.
  This script loads CNV calls from different methods and converts them to standardized VCF format
  using the cnvstrmliner package.
- Add CNAqc metrics.
  If a somatic mutations VCF is provided, it will add qc metrics with CNAqc.
- Compute CNV entropy.
  If a bed file with covered regions is provided, it will compute CNV entropy.
- Compute Gini and MAPD metrics.
  If a file with normalized coverage per target, it will compute Gini and MAPD metrics for CNV.
- Save output files.
  The processed CNV calls are saved in VCF format, along with optional QC statistics and
  RDS files containing the cnvstrml object.",
  epilogue = sprintf(
    "Developer: Konstantinos Mavrommatis\nDate: 8 January 2026\ncnvstrmliner version: %s",
    packageVersion("cnvstrmliner")
  )
)

args <- parse_args(parser)

# Validate required arguments
if (is.null(args$input) || !file.exists(args$input)) {
  stop("Error: --input is required. Please provide the path to the CNV file.")
} else {
  logger::log_info("CNV file                 : {args$input}")
}

if (is.null(args$method)) {
  stop("Error: --method is required. Please specify the method used (FACETS|Sequenza|Battenberg|CFREEC|CNVkit|GATK|PureCN|SCLUST|SEG).")
} else {
  logger::log_info("CNV method               : {args$method}")
}

if (is.null(args$output)) {
  stop("Error: --output is required. Please provide the path for the output VCF file.")
} else {
  logger::log_info("Output file              : {args$output}")
}

if (!is.na(args$sample)) {
  logger::log_info("Sample used              : {args$sample}")
}

# Validate input file exists


if (!is.null(args$somatic_mutations) && !file.exists(args$somatic_mutations)) {
  logger::log_error("Where is the somatic mutations VCF file? {args$somatic_mutations}")
  stop()
}
if (!is.null(args$somatic_mutations)) {
  logger::log_info("Using somatic mutations  : {args$somatic_mutations}")
}

if (!is.null(args$cnr_file) && !file.exists(args$cnr_file)) {
  logger::log_error("Where is the file with logR values per target (.cnr, .cns etc)? {args$cnr_file}")
  stop()
}
if (!is.null(args$cnr_file)) {
  logger::log_info("Using logR file          : {args$cnr_file}")
}

if (!is.null(args$bed_file) && !file.exists(args$bed_file)) {
  logger::log_error("Where is the BED file with the target coordinates? {args$bed_file}")
  stop()
}
if (!is.null(args$bed_file)) {
  logger::log_info("Using BED file           : {args$bed_file}")
}
# Validate method
valid_methods <- cnvstrmliner::getSupportedMethods()
if (!tolower(args$method) %in% tolower(valid_methods)) {
  logger::log_error("Invalid method: {args$method}")
  logger::log_error("Valid methods are: {paste(valid_methods, collapse = ', ')}")
  stop()
}


# Load CNV calls based on method
logger::log_info("Loading CNV calls from {args$method} in file {args$input} ")

cnv <- tryCatch(
  {
    parseCNV(
      cnv.data = args$input,
      cnv.method = args$method,
      sample = args$sample
    )
  },
  error = function(e) {
    logger::log_error("Failed to parse CNV file: {e$message}")
    stop()
  }
)

if (!is.na(args$sample) && cnv@sample != args$sample) {
  logger::log_warn("The sample name in the CNV file ({cnv@sample}) does not match the provided sample name ({args$sample}). Using the provided sample name.")
  cnv@sample <- args$sample
}

logger::log_info("Sample used              : {cnv@sample}")


# load the mutations if they are provided
entropy <- NA_real_

# add tje CNAqc metrics if we have somatic mutations
if (!is.null(args$somatic_mutations)) {
  logger::log_info("Loading somatic mutations from  {args$somatic_mutations}")
  vcfHeader <- scanVcfHeader(args$somatic_mutations)
  vcfSamples <- samples(vcfHeader)
  sample <- cnv@sample # this holds the final sample name used in the cnvstrml object
  if (!cnv@sample %in% vcfSamples) {
    logger::log_warn("The sample name {cnv@sample} is not found in the VCF file samples: {paste(vcfSamples, collapse=', ')}")
    sample <- NA_character_
  }
  if (is.na(sample)) {
    # try some heuristics to find the sample name as tumor or TUMOR or Tumor
    if ("tumor" %in% vcfSamples) {
      sample <- "tumor"
    } else if ("TUMOR" %in% vcfSamples) {
      sample <- "TUMOR"
    } else if ("Tumor" %in% vcfSamples) {
      sample <- "Tumor"
    }
    if (!is.na(sample)) {
      logger::log_warn("Using '{sample}' as sample name from VCF samples")
    }
  }

  # check if the samples in the VCF contain the sample name as substring
  # if there is a unique match use it
  if (is.na(sample)) {
    matched.samples <- vcfSamples[grepl(cnv@sample, vcfSamples)]
    if (length(matched.samples) == 1 && !is.na(matched.samples)) {
      sample <- matched.samples
      logger::log_warn("Using '{sample}' as sample name from VCF samples matching the input sample name")
    }
  }


  if (is.na(sample)) {
    logger::log_error("Cannot find a suitable sample name in the VCF file")
  } else {
    muts <- cnvstrmliner::loadMutationsFromVCF(
      args$somatic_mutations,
      sample = sample,
      genome = "hg38",
      filter = TRUE
    )

    cnv <- applyQCMetrics(cnv, muts, n_bootstrap = 10)
    cnv <- detectOverfragmentation(cnv)
  }
}


if (!args$bed_file %>% is.null()) {
  logger::log_info("Calculating CNV entropy using BED file: {args$bed_file}")
  entropy <- cnvEntropy(cnv, args$bed_file)
  cnv@quality$entropy <- entropy
  cnv@quality$bed_file <- args$bed_file
  logger::log_info("Calculated CNV entropy: {entropy}")
}

if (!args$cnr_file %>% is.null()) {
  logger::log_info("Calculating Gini and MAPD using BAM file: {args$cnr_file} and BED file: {args$bed_file}")

  coverage_tbl <- data.table::fread(args$cnr_file) # needs to have the column 'depth

  gini <- computeGini(coverage_tbl)
  mapd <- computeMAPD(coverage_tbl)
  cnv@quality$gini <- gini
  cnv@quality$mapd <- mapd
  cnv@quality$logRfile <- args$cnr_file
  logger::log_info("Calculated Gini: {gini} and MAPD: {mapd}")
}

# Create output directory if it doesn't exist
output_dir <- dirname(args$output)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  logger::log_info("Created output directory: {output_dir}")
}

# Save as VCF
logger::log_info("Saving CNV calls to VCF format: {args$output}")

tryCatch(
  {
    writeVCF(cnv, file = args$output)
    logger::log_info("Successfully saved CNV calls to             : {args$output}")
  },
  error = function(e) {
    logger::log_error("Failed to write VCF file: {e$message}")
    stop()
  }
)


output_stats <- args$output %>%
  gsub(pattern = "\\.vcf.gz$", replacement = "_cnvqc_stats.tsv") %>%
  gsub(pattern = "\\.vcf$", replacement = "_cnvqc_stats.tsv")

tryCatch(
  {
    tbl <- qcTable(cnv)
    write.table(tbl, file = output_stats, sep = "\t", quote = FALSE, row.names = FALSE)
    logger::log_info("Successfully saved statistics to            : {output_stats}")
  },
  error = function(e) {
    logger::log_error("Failed to write CNV QC stats file: {e$message}")
  }
)


output_summary <- args$output %>%
  gsub(pattern = "\\.vcf.gz$", replacement = "_cnvqc_summary.txt") %>%
  gsub(pattern = "\\.vcf$", replacement = "_cnvqc_summary.txt")

tryCatch(
  {
    summary_lines <- bind_rows(
      summaryTable(cnv),
      summaryTable(cnv, level = "chromosome"),
      summaryTable(cnv, level = "arm")
    )
    write.table(summary_lines, file = output_summary, sep = "\t", quote = FALSE, row.names = FALSE)
    logger::log_info("Successfully saved summary to               : {output_summary}")
  },
  error = function(e) {
    logger::log_error("Failed to write CNV QC summary file: {e$message}")
  }
)

if (!args$somatic_mutations %>% is.null()) {
  output_rds <- args$output %>%
    gsub(pattern = "\\.vcf.gz$", replacement = "_cnvqc.rds") %>%
    gsub(pattern = "\\.vcf$", replacement = "_cnvqc.rds")
  tryCatch(
    {
      saveRDS(cnv, file = output_rds)
      logger::log_info("Successfully saved CNV calls and QC to      : {output_rds}")
    },
    error = function(e) {
      logger::log_error("Failed to write CNV QC RDS file: {e$message}")
    }
  )
}

logger::log_info("Done!")
