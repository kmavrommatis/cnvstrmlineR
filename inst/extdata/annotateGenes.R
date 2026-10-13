#!/usr/bin/env Rscript
suppressPackageStartupMessages(
  suppressWarnings({
    library(optparse)
    library(dplyr)
    library(logger)
    library(cnvstrmliner)
    library(AnnotationHub)
    library(GenomeInfoDb)
    library(ensembldb)
    library(org.Hs.eg.db)
  })
)
# This is a utility script that
# is used to create annotation of genes and segments
# from CNV calls.
#
# the CNV calls should be in VCF format coming from cnvstrmliner
# the AnnotationHub should be cached in a directory.
#

option_list <- list(
  make_option(c("-i", "--input"),
    type = "character", default = NULL,
    help = "Path to the input file with the CNV calls as a vcf.gz file.",
    metavar = "character"
  ),
  make_option(c("-c", "--cache"),
    type = "character", default = Sys.getenv("ANNOTATION_HUB_CACHE"),
    help = "Path to the cache directory with AnnotationHub data.",
    metavar = "character"
  ),
  make_option(c("-s", "--sample"),
    type = "character", default = NULL,
    help = "The name of the sample",
    metavar = "character"
  ),
  make_option(c("-o", "--output_directory"),
    type = "character", default = ".",
    help = "Path to the directory with the annotations will be stored. ",
    metavar = "character"
  )
)

parser <- OptionParser(
  option_list = option_list,
  description = "Annotate CNV segments with overlapping genes using Ensembl annotations.
This script takes CNV calls in VCF format from cnvstrmliner and annotates them with
gene information from AnnotationHub (Ensembl database release 109 for hg38).",
  epilogue = sprintf(
    "Developer: Konstantinos Mavrommatis\nDate: 8 January 2026\ncnvstrmliner version: %s",
    packageVersion("cnvstrmliner")
  )
) # Create the parser
args <- parse_args(parser)

if (is.null(args$input)) {
  stop("Input file must be provided")
}
if (is.null(args$cache) |
  !dir.exists(args$cache)) {
  stop("A valid cache directory must be provided")
}

# load the CNV calls
logger::log_info("Loading CNV calls from {args$input}")
cnvObj <- cnvstrmliner::parseVCF(vcf = args$input)
if (!is.null(args$sample)) {
  cnvObj@sample <- args$sample
}


logger::log_info("Using cache directory {args$cache}")
setAnnotationHubOption("CACHE", args$cache)

ah <- AnnotationHub(localHub = TRUE)
# Query for the Ensembl database for Homo sapiens, release 109 (corresponds to hg38)
query_results <- query(ah, c("EnsDb", "Homo sapiens", "109"))

# Load the database
edb <- query_results[[1]]



# Get all genes with their symbols
# do a conversion of the chromosome names
# but not by using GenomeInfoDb::seqlevelsStyle
# which frequently fails with network issues.
genes_hg38 <- suppressWarnings(
  genes(edb,
    columns = c("gene_name", "gene_id", "symbol", "entrezid", "ensemblid"),
    return.type = "GRanges"
  )
) %>%
  keepStandardChromosomes(pruning.mode = "coarse") %>%
  as.data.frame() %>%
  dplyr::mutate(seqnames = dplyr::case_when(
    seqnames %in% c(seq(1, 22), "X", "Y") ~ paste0("chr", seqnames),
    seqnames == "M" ~ "MT",
    TRUE ~ seqnames
  )) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(entrezid = paste(entrezid, collapse = ", ")) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    gene = ifelse(gene_name == "", gene_id, gene_name)
  ) %>%
  dplyr::arrange(seqnames, gene, start) %>%
  dplyr::group_by(seqnames, gene) %>%
  dplyr::mutate(
    # Identify overlapping groups: new group when current start > previous end
    overlap_group = cumsum(dplyr::lag(end, default = 0L) < start)
  ) %>%
  dplyr::group_by(seqnames, gene, overlap_group) %>%
  dplyr::arrange(desc(width)) %>%
  dplyr::slice(1) %>% # Keep the longest gene in each overlap group
  dplyr::ungroup() %>%
  dplyr::select(-overlap_group) %>%
  makeGRangesFromDataFrame(., keep.extra.columns = TRUE) %>%
  GenomicRanges::sort()

# since we want to include ENTREZid we also query the NCBI gene database
# with the ensembl gene id and get the corresponding entrez id from org.Hs.eg.db
# this is because the ensembl database does not always have the entrez id for all genes
# and the org.Hs.eg.db has a more comprehensive mapping of ensembl gene id to

all_human_genes_by_ens <- select(
  x = org.Hs.eg.db,
  keys = genes_hg38$gene_id,
  columns = c("ENSEMBL", "SYMBOL", "ENTREZID"),
  keytype = "ENSEMBL"
) %>%
  dplyr::group_by(ENSEMBL) %>%
  dplyr::mutate(entrezid = paste(ENTREZID, collapse = ",")) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(gene = ifelse(SYMBOL == "", ENSEMBL, SYMBOL))

genes_hg38$entrezid <- all_human_genes_by_ens$entrezid[match(genes_hg38$gene, all_human_genes_by_ens$gene)]

# overlap with the segments
ov <- suppressWarnings(
  findOverlaps(
    query = cnvObj@segments,
    subject = genes_hg38
  )
)
# find the overlaps and map to each other
overlaps <- suppressWarnings(
  pintersect(
    cnvObj@segments[queryHits(ov)],
    genes_hg38[subjectHits(ov)]
  )
)
# Add the original gene index to track which overlap belongs to which gene
mcols(overlaps)$gene_idx <- subjectHits(ov)
mcols(overlaps)$gene_id <- genes_hg38$gene_id[subjectHits(ov)] # ensembl id
mcols(overlaps)$gene_name <- genes_hg38$gene_name[subjectHits(ov)] # hugo symbol
mcols(overlaps)$gene_entrez <- genes_hg38$entrezid[subjectHits(ov)] # entrez id
mcols(overlaps)$gene_start <- start(genes_hg38)[subjectHits(ov)]
mcols(overlaps)$gene_end <- end(genes_hg38)[subjectHits(ov)]
mcols(overlaps)$gene_width <- width(genes_hg38)[subjectHits(ov)]

# Calculate the weight (width of the overlap)
# svpluscnv uses genomic length for similar weighted averages
mcols(overlaps)$overlap_width <- width(overlaps)
# Perform the weighted average
df_overlaps <- as.data.frame(mcols(overlaps))
stats_per_gene <- aggregate(
  cbind(
    weighted_val = overlap_width * total.allele.copies,
    covered_width = overlap_width
  ) ~ gene_idx,
  data = df_overlaps,
  FUN = sum
)
stats_per_gene$weighted_avg <- stats_per_gene$weighted_val / stats_per_gene$covered_width
genes_hg38$weighted_cn <- NA
genes_hg38$weighted_cn[stats_per_gene$gene_idx] <- stats_per_gene$weighted_avg
mcols(overlaps)$weighted_cn <- genes_hg38$weighted_cn[subjectHits(ov)]
# now show the genes with the annotation of the cnv

# subject_ranges <- ov$genes_hg38
# chrom <- seqnames(subject_ranges) %>% as.character()
# st <- start(subject_ranges)
# en <- end(subject_ranges)
mcols(overlaps)$seqnames <- seqnames(genes_hg38)[subjectHits(ov)]
ann <- mcols(overlaps)[, c(
  "gene_id", # ensembl id
  "gene_name", # hugo symbol
  "gene_entrez", # entrez id
  "cnv.annotation",
  "cnv.loh",
  "cnv.lrr",
  "cnv.cf",
  "weighted_cn",
  "seqnames",
  "gene_start",
  "gene_end"
)] %>%
  as.data.frame() %>%
  dplyr::mutate(
    gene = ifelse(gene_name == "", gene_id, gene_name)
  ) %>%
  unique()

logger::log_info("Computing gene level CNV information")

gene_level <- cbind(
  "sample_id" = cnvObj@sample,
  "Ploidy" = cnvObj@effective_ploidy,
  "B.Allele.Freq" = NA,
  ann
) %>%
  dplyr::rename(
    chr = seqnames,
    start = gene_start,
    end = gene_end,
    hugo_symbol = gene_name,
    ensembl_id = gene_id,
    entrez_id = gene_entrez,
    gene_name = gene,
    Status = cnv.annotation,
    Ratio = cnv.lrr,
    Absolute.Copy.Number = weighted_cn,
    Sample.Tumor.Purity = cnv.cf,
    LOH = cnv.loh
  ) %>%
  dplyr::select(chr, start, end, gene_name, Status, Ratio, Absolute.Copy.Number, Sample.Tumor.Purity, LOH, hugo_symbol, entrez_id, ensembl_id) %>%
  dplyr::filter(!is.na(gene_name) & gene_name != "") %>%
  dplyr::mutate(
    chr = stringr::str_remove(chr, "chr"),
    chr = ifelse(chr == "M", "MT", chr),
    Status = dplyr::case_when(
      grepl("NEUT", Status) ~ "Neutral",
      grepl("HOMDEL", Status) ~ "Deep deletion",
      grepl("HETDEL", Status) ~ "Shallow deletion",
      grepl("AMP", Status) ~ "Amplification",
      TRUE ~ Status
    ),
    LOH = dplyr::case_when(
      (LOH == TRUE & Status == "Neutral") ~ "copy-neutral LOH",
      (LOH == TRUE & Status == "Shallow deletion") ~ "LOH: loss of heterozygosity",
      TRUE ~ ""
    )
  )

logger::log_info("Computing segment level CNV information")

segment_level <- cnvObj@segments %>%
  as.data.frame(row.names = NULL) %>%
  dplyr::rename(
    chr = seqnames,
    Status = cnv.annotation,
    Ratio = cnv.lrr,
    Absolute.Copy.Number = total.allele.copies,
    Sample.Tumor.Purity = cnv.cf,
    LOH = cnv.loh
  ) %>%
  dplyr::select(chr, start, end, Status, Ratio, Absolute.Copy.Number, Sample.Tumor.Purity, LOH, major.allele.copies, minor.allele.copies) %>%
  dplyr::mutate(
    chr = stringr::str_remove(chr, "chr"),
    chr = ifelse(chr == "M", "MT", chr),
    Genotype = paste0(
      strrep("A", major.allele.copies),
      strrep("B", minor.allele.copies)
    ),
    Status = dplyr::case_when(
      grepl("NEUT", Status) ~ "Neutral",
      grepl("HOMDEL", Status) ~ "Deep deletion",
      grepl("HETDEL", Status) ~ "Shallow deletion",
      grepl("AMP", Status) ~ "Amplification",
      TRUE ~ Status
    ),
    LOH = dplyr::case_when(
      (LOH == TRUE & Status == "Neutral") ~ "copy-neutral LOH",
      (LOH == TRUE & Status == "Shallow deletion") ~ "LOH: loss of heterozygosity",
      TRUE ~ ""
    )
  ) %>%
  dplyr::select(-c(major.allele.copies, minor.allele.copies))

if (!dir.exists(args$output_directory)) {
  dir.create(args$output_directory, recursive = TRUE)
}
gene_file <- file.path(args$output_directory, paste0(cnvObj@sample, "-", cnvObj@method, ".gene_level.tsv"))
data.table::fwrite(
  gene_level,
  file = gene_file,
  append = FALSE,
  quote = FALSE,
  sep = "\t"
)

logger::log_info("Gene level annotation stored in file {gene_file}")

seg_file <- file.path(args$output_directory, paste0(cnvObj@sample, "-", cnvObj@method, ".segment_level.tsv"))
data.table::fwrite(
  segment_level,
  file = seg_file,
  append = FALSE,
  quote = FALSE,
  sep = "\t"
)

logger::log_info("Segment level annotation stored in file {seg_file}")
