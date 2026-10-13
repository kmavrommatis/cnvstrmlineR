#!/usr/bin/env  Rscript

suppressPackageStartupMessages(
    suppressWarnings({
        library(cnvstrmliner)
        library(future)
        library(future.apply)
        library(optparse)
        library(logger)
        library(GenomicRanges)
        library(VariantAnnotation)
    })
)
option_list <- list(
    make_option(c("-i", "--input_directory"),
        type = "character", default = NULL,
        help = "Path to the input directory with the CNV calls from the target method (FACETS or Sequenza). Expected to have multiple files. ",
        metavar = "character"
    ),
    make_option(c("-s", "--sample_name"),
        type = "character", default = NULL,
        help = "sample name. Will be used to select mutations from the VCF file and name the output",
        metavar = "character"
    ),
    make_option(c("-v", "--somatic_mutations"),
        type = "character", default = NULL,
        help = "somatic mutations VCF",
        metavar = "character"
    ),
    make_option(c("--method"),
        type = "character", default = "FACETS",
        help = "Method used to produce CNV results (FACETS|Sequenza) [default:FACETS]"
    ),
    make_option(c("--output_directory"),
        type = "character", default = ".",
        help = "output directory for the results [default: current working directory]"
    ),
    make_option(c("--threads"),
        type = "numeric", default = 1,
        help = "parallel threads to use [default: 1]"
    )
)

parser <- OptionParser(
    option_list = option_list,
    description = "Add CNAqc (Copy Number Aberration Quality Control) analysis to CNV calls.
==================================
It performs quality control analysis using the CNAqc framework.
When this script processes CNV calls from FACETS or Sequenza, it can check for alternative
solutions and/or focal solutions and process them independently.
Somatic mutations should be provided and will be integrated.
",
    epilogue = sprintf(
        "Developer: Konstantinos Mavrommatis\nDate: 8 January 2026\ncnvstrmliner version: %s",
        packageVersion("cnvstrmliner")
    )
) # Create the parser
args <- parse_args(parser) # Pass the parser to parse_args

logger::log_threshold(INFO, namespace = "script")
logger::log_threshold(INFO, namespace = "cna method")


# args <- list(
#     input_directory = "/data/",
#     sample_name = "MC3200-05_27_2008-B_L1.SSHT68_1",
#     somatic_mutations = "/data/MC3200-05_27_2008-B_L1.SSHT68_1-18348000780_L3.SSHT10_1.vcf.gz",
#     method = "Sequenza",
#     output_directory = "/data/res",
#     threads = 8
# )



method <- args$method
if (!tolower(method) %in% c("facets", "sequenza")) {
    stop("The methods that can be parsed are FACETS or Sequenza")
}
focal.size <- args$focal_size
focal.annotation <- args$focal_annotation

# Handle input_directory or input_file
if (!is.null(args$input_directory)) {
    cna_directory <- normalizePath(args$input_directory, mustWork = TRUE)
    cna_file <- NULL
    if (!dir.exists(cna_directory)) {
        logger::log_error("The input directory {args$input_directory} does not exist")
        stop()
    }
} else {
    # input_file is provided
    cna_file <- normalizePath(args$input_file, mustWork = TRUE)
    if (!file.exists(cna_file)) {
        logger::log_error("The input file {args$input_file} does not exist")
        stop()
    }
    # Set cna_directory to the directory containing the file for compatibility
    cna_directory <- dirname(cna_file)
}

if (!dir.exists(args$output_directory)) {
    dir.create(args$output_directory, recursive = TRUE, showWarnings = FALSE)
}
args$output_directory <- normalizePath(args$output_directory, mustWork = TRUE)


logger::log_tictoc("Processing files in directory {cna_directory}", level = INFO, namespace = "script")
logger::log_tictoc("Storing output in directory {args$output_directory}", level = INFO, namespace = "script")


if (is.null(args$sample_name)) {
    sample <- basename(cna_directory) %>%
        gsub(pattern = ".sorted.dedup.realigned.recal.hg38", replacement = "") %>%
        stringr::str_split_i(., pattern = "-", 1)
} else {
    sample <- args$sample_name
}

logger::log_tictoc("Working with sample {sample} processed by {method}", level = INFO, namespace = "cna method")

output.file <- file.path(
    args$output_directory,
    paste0(sample, "-", tolower(method), "-qc.Rdata")
)
if (!file.create(output.file)) {
    logger::log_warn("Unable to create the output file {output.file}")
    stop()
}


if (args$somatic_mutations %>% is.null() || !file.exists(args$somatic_mutations)) {
    logger::log_error("Where is the somatic mutations VCF file? {args$somatic_mutations}")
    stop()
}
logger::log_tictoc("Loading somatic mutations for {sample} from {args$somatic_mutations}", level = INFO, namespace = "cna method")

# check the sample names in the VCF file
vcfHeader <- VariantAnnotation::scanVcfHeader(args$somatic_mutations)
vcfSamples <- VariantAnnotation::samples(vcfHeader)
if (!sample %in% vcfSamples) {
    sample <- NA_character_
    logger::log_warn("The sample name {sample} is not found in the VCF file samples: {paste(vcfSamples, collapse=', ')}")
    # try some heuristics
    if ("tumor" %in% vcfSamples) {
        sample <- "tumor"
        logger::log_warn("Using 'tumor' as sample name")
    } else if ("TUMOR" %in% vcfSamples) {
        sample <- "TUMOR"
        logger::log_warn("Using 'TUMOR' as sample name")
    }

    # check if the samples in the VCF contain the sample name as substring
    # if there is a unique match use it
    if (is.na(sample)) {
        matched.samples <- vcfSamples[grepl(sample, vcfSamples)]
        if (length(matched.samples) == 1) {
            sample <- matched.samples
            logger::log_warn("Using '{sample}' as sample name from VCF samples matching the input sample name")
        }
    }

    if (is.na(sample)) {
        logger::log_error("Cannot find a suitable sample name in the VCF file")
        stop()
    }
}

muts <- cnvstrmliner::loadMutationsFromVCF(
    args$somatic_mutations,
    sample = sample,
    genome = "hg38",
    filter = TRUE
)
output.list <- c()
focal <- FALSE
alternatives <- FALSE
allfiles <- list.files(cna_directory, recursive = FALSE, full.names = TRUE)

# when we work with a directory
# i.e. expect to have either focal or alternative solutions or both

if (!is.null(cna_directory)) {
    if (any(grepl("focal.tgz$", allfiles) == TRUE)) {
        focal.file <- allfiles[grepl("focal.tgz$", allfiles)]
        focal <- TRUE
        logger::log_tictoc("Found focal solutions in {focal.file}", level = INFO, namespace = "cna method")
    } else {
        logger::log_tictoc("There aren't any focal events to process.", level = INFO, namespace = "cna method")
    }


    if (any(grepl("alternatives.tgz$", allfiles) == TRUE) && method == "Sequenza") {
        alternatives.file <- allfiles[grepl("alternatives.tgz$", allfiles)]
        alternatives <- TRUE
        logger::log_tictoc("Found alternative solutions in {alternatives.file}", level = INFO, namespace = "cna method")
    } else {
        logger::log_tictoc("There aren't any alternative solutions to process", level = INFO, namespace = "cna method")
    }
    cat(paste(allfiles, sep = ","))


    oplan <- plan(multisession, workers = args$threads)
    main.solution <- NA
    if (tolower(method) == "facets") {
        main.solution <- allfiles[grepl("vcf.gz$", allfiles)]
        if (length(main.solution) != 1) {
            logger::log_error("Cannot find the main solution seg file in {cna_directory}. ", namespace = "cna method")
            logger::log_error("Found {length(main.solution)} files matching the pattern. Expected exactly 1.", namespace = "cna method")
            stop()
        }
        logger::log_tictoc("Loading main solution {main.solution}", level = INFO, namespace = "cna method")
    }

    if (tolower(method) == "sequenza") {
        cat(paste(allfiles, sep = ","))
        logger::log_tictoc("Looking for the main solution in {cna_directory}", level = INFO, namespace = "cna method")
        main.solution <- allfiles[grepl("_segments.txt$", allfiles)]
        logger::log_tictoc("Loading main solution {main.solution}", level = INFO, namespace = "cna method")
        if (length(main.solution) != 1) {
            logger::log_error("Cannot find the main solution seg file in {cna_directory}. ", namespace = "cna method")
            logger::log_error("Found {length(main.solution)} files matching the pattern. Expected exactly 1.", namespace = "cna method")
            stop()
        }
    }
    main.cnv <- loadCNV(cnv.data = main.solution, cnv.method = method, sample = sample)
} else {
    stop("No input directory provided")
}
logger::log_tictoc("Applying QC metrics on main solution", level = INFO, namespace = "cna method")
main.cnv <- applyQCMetrics(main.cnv, muts, n_bootstrap = 1)

output.list <- ("main.cnv")

# Load alternatives if available
if (alternatives) {
    out_alt_dir <- file.path(tempdir(), sample, "alternatives")
    logger::log_tictoc("Loading alternative solutions from {alternatives.file}", level = INFO, namespace = "cna method")
    untar(file.path(cna_directory, alternatives.file), exdir = out_alt_dir)
    # parse the options under the focal directory
    dd <- out_alt_dir

    alt.solutions <- list.files(dd, pattern = "_segments.txt", recursive = TRUE, full.names = TRUE)

    alt.cnv <- lapply(alt.solutions, function(x) {
        a <- parseSequenza(x, sample = sample, methodName = basename(x))
        a <- applyQCMetrics(a, muts, n_bootstrap = 1)
        a
    })
    output.list <- c(output.list, "alt.cnv")
}

# load the focal solutions if available
if (focal) {
    out_focal_dir <- file.path(tempdir(), sample, "focal")
    logger::log_tictoc("Loading focal solutions from {focal.file}", level = INFO, namespace = "cna method")
    logger::log_tictoc("Unpacking focal solutions {file.path(cna_directory, focal.file)} to {out_focal_dir}", level = INFO, namespace = "cna method")
    untar(file.path(cna_directory, focal.file), exdir = out_focal_dir)

    # parse the options under the focal directory
    dd <- out_focal_dir
    if (method == "FACETS") {
        focal.solutions <- list.files(dd, pattern = "vcf.gz$", recursive = TRUE, full.names = TRUE)
    }
    if (method == "Sequenza") {
        focal.solutions <- list.files(dd, pattern = "_segments.txt$", recursive = TRUE, full.names = TRUE)
    }
    logger::log_tictoc("There are {length(focal.solutions)} focal solutions", level = INFO, namespace = "cna method")

    focal.cnv <- lapply(focal.solutions, function(x) {
        logger::log_tictoc("Loading focal solution {x}", level = INFO, namespace = "cna method")
        a <- if (method == "FACETS") {
            parseFACETS(x, sample = sample, methodName = basename(x))
        } else if (method == "Sequenza") {
            parseSequenza(x, sample = sample, methodName = basename(x))
        }
        if (is.na(a@purity)) {
            a@purity <- main.cnv@purity
        }
        a
    })


    mean_purity <- sapply(focal.cnv, function(X) {
        X@purity
    }) %>% mean(., na.rm = TRUE)
    logger::log_tictoc("Loaded focal solutions. Mean purity {mean_purity}", level = INFO, namespace = "cna method")
    focal.cnv <- future_lapply(focal.cnv, function(a) {
        if (is.na(a@purity)) {
            a@purity <- mean_purity
        }
        a <- applyQCMetrics(a, muts, n_bootstrap = 1)
        a <- detectOverfragmentation(a, arm_percentage_cutoff = 1e6, minimum_segments_for_testing = 20)
        a@segments <- cnvstrmliner:::mergeSegments(a@segments, minDistance = 1e5)
        a
    }, future.seed = TRUE)

    names(focal.cnv) <- sapply(focal.cnv, function(x) {
        x@sample
    })
    output.list <- c(output.list, "focal.cnv")
}
logger::log_tictoc("Storing results in {basename(output.file)}", level = INFO, namespace = "cna method")
save(list = output.list, file = file.path(args$output_directory, basename(output.file)))
logger::log_tictoc("Cleaning up temporary files", level = INFO, namespace = "cna method")
unlink(file.path(tempdir(), sample), recursive = TRUE)
logger::log_tictoc("Processing of sample {sample} from {cna_directory} finished successfully", level = INFO, namespace = "script")
logger::log_tictoc("Results are stored in {args$output_directory}", level = INFO, namespace = "script")
