#' Base and tidyverse functions used across the package
#'
#' Imports for functions that are used unqualified in the package code so that
#' \code{R CMD check} can resolve them.
#'
#' @importFrom methods as is
#' @importFrom stats aggregate as.dist dbinom median quantile sd var weighted.mean
#' @importFrom utils packageVersion read.csv read.table write.table globalVariables
#' @importFrom dplyr across bind_rows desc relocate rowwise starts_with
#' @importFrom tidyr drop_na
#' @name cnvstrmliner-imports
#' @keywords internal
NULL

# Column names referenced inside data-masking verbs (dplyr/tidyr/data.table) and
# the magrittr placeholder are not visible to R's code analysis. Declaring them
# here keeps `R CMD check` quiet without changing behaviour.
utils::globalVariables(c(
  ".", "..col_select", "A", "B", "C", "FALSE", "TRUE", "L", "M",
  "ALT", "REF", "DP", "NV", "VAF", "alt", "ref", "id", "chr", "Chr", "chrom",
  "chromosome", "Chromosome", "chrom_gaps", "chr_length", "chr_percentage",
  "arm_length", "arm_percentage", "genome", "genome_percentage",
  "CALL", "CCF", "CF", "CF_EM", "CNLR_MEDIAN", "CNV_ANN", "CN_LR_MEDIAN",
  "CN", "CNt", "CopyNr", "Copy_Nr_Raw", "Cellularity",
  "Clone1_A", "Clone1_B", "Clone1_ccf", "Clone2_A", "Clone2_B",
  "LCN_EM", "LOH", "LogR", "MAJOR_CN", "MINOR_CN", "TOTAL_CN", "TCN_EM",
  "MEAN_LOG2_COPY_RATIO", "MedianRatio", "Major", "Major_2",
  "NHET", "NUM_MARK", "NUM_HET_MARK", "NUM_POINTS_COPY_RATIO", "N.BAF",
  "QC", "QC_PASS", "SVTYPE", "Start", "START_T", "End", "END_T",
  "StartPos", "EndPos", "Values", "WilcoxonRankSumTestPvalue",
  "KolmogorovSmirnovPvalue", "bin_idx", "ccf", "cellular_prevalence",
  "cf.em", "tcn.em", "cn", "cn1", "cn2", "cnv.cf", "cnv.cnaqc.pass", "cnv.loh",
  "cnv.lrr", "cnv.markers", "cnv.totalmarkers", "consistency",
  "copy number", "copy.number", "count", "depth", "depth.ratio",
  "e", "end.pos", "endpos", "flags", "frac1_A", "frac2_A",
  "freec_cn", "freec_cn1", "freec_cn2", "from", "to",
  "gatk_call", "gatk_cn", "gatk_mean_log2_copy_ratio", "genotype",
  "karyotype", "lcn.em", "logR", "loh", "lr",
  "mCN", "major_cn", "minor_cn", "major.allele.copies", "minor.allele.copies",
  "total.allele.copies", "minor", "minor_2",
  "nMaj1_A", "nMaj2_A", "nMin1_A", "nMin2_A", "n_SNPs", "n_len",
  "n_long", "n_short", "num.mark", "num.snps", "p_value", "probes", "pval",
  "query", "query_CCF", "query_CNV", "subject", "subject_CCF", "subject_CNV",
  "run_id", "seg", "seg_idx", "seg.mean", "segment", "segments", "segvalue",
  "segment.chr.start", "segment.end.start", "segment.start.start",
  "segment.type", "segment.type.start", "segment.width.start",
  "significant", "size_cutoff", "smaller", "snps", "start.pos", "startpos",
  "status", "strand", "support", "sv_type", "sv.support", "sv.support.class",
  "SV3", "SV3.chr", "SV3.end", "SV3.len", "SV3.partner", "SV3.partner.chr",
  "SV3.partner.end", "SV3.partner.start", "SV3.start", "SV3.type",
  "SV5", "SV5.chr", "SV5.end", "SV5.len", "SV5.partner", "SV5.partner.chr",
  "SV5.partner.end", "SV5.partner.start", "SV5.start", "SV5.type",
  "svend", "svstart", "tmp_matched", "tmp_total", "tmp_value", "tmp_weight",
  "total_bp", "total_segments", "total_width", "uncertainty", "val"
))
