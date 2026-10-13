# cnvstrmliner

<p align="center">
  <img src="man/figures/logo.png" alt="cnvstrmliner logo" width="400">
</p>

R package for parsing CNV tools' outputs

# Introduction

The proliferation of Next-Generation Sequencing (NGS) tools for predicting 
Copy Number Aberrations (CNAs) has led to a significant data integration 
challenge in bioinformatics. Each application, while serving the same core purpose, 
generates results in a unique output format, often spread across multiple files. 
Furthermore, the specific metrics and fields reported—such as segment start/end positions, 
log-ratio values, and associated p-values—are frequently named differently or 
are not uniformly present across all tool outputs. 
This lack of standardization creates substantial overhead for researchers, 
who must write custom parsing scripts for every new tool, leading to complex 
and error-prone downstream analysis. It is exceptionally difficult to compare 
results across different methods, integrate them into common analysis pipelines, 
or use them with existing visualization and interpretation software, 
even when using semi-standard formats like SEG or VCF which haven't achieved universal adoption.

A dedicated package capable of reading the diverse outputs of CNA prediction tools 
and translating them into a single, common, standardized format is crucial for 
streamlining genomic data analysis. Such a package would significantly reduce 
the need for custom scripting, allowing researchers to focus on biological 
interpretation rather than data wrangling. By providing a consistent data 
structure, it would immediately enable method comparison, consensus calling, 
and seamless integration with downstream analysis and visualization tools, 
irrespective of the initial prediction method. This unified approach promotes 
reproducibility and interoperability across different studies and laboratories, 
ultimately accelerating the pace of discovery in cancer genomics and other fields 
relying on CNA analysis.

`cnvstrmliner` is a package that tries to assume the role of this middle layer
that is capable of parsing multiple output formats and provide a single uniform one 
using standard Bioconductor objects.
Furthermore, some additional functionality is added to facilitate the comparisons and 
QC.


# Installation 

The github repo is

https://github.com/kmavrommatis/cnvstrmlineR

Before installing `cnvstrmliner` itself, install its GitHub-only dependencies
(see [Dependencies](#dependencies) below), then run:

```r
# potentially unload a previous copy of the package
detach("package:cnvstrmliner", unload = TRUE)

devtools::install_github("kmavrommatis/cnvstrmlineR", force = TRUE)
```

To install with the vignettes available locally via
`browseVignettes("cnvstrmliner")`, ask for them explicitly:

```r
devtools::install_github("kmavrommatis/cnvstrmlineR",
                         force = TRUE, build_vignettes = TRUE)
```

Note that building the vignettes takes a while — they execute the full parsing,
consensus and QC workflow on the bundled example data — and additionally
requires `knitr`, `rmarkdown` and `BiocStyle`. If you only need the package,
use the first form; the same articles are available online at
<https://kmavrommatis.github.io/cnvstrmlineR/>.

## Dependencies

`cnvstrmliner` depends on a number of external packages. Most of them are
installed automatically when you install `cnvstrmliner`, but two GitHub-only
packages are not resolved automatically and should be installed beforehand:

- **CNAqc** is a required, GitHub-only dependency (used by `applyQCMetrics`
  and listed in `Remotes:`) and must be installed first:

  ```r
  devtools::install_github("caravagnalab/CNAqc")
  ```

- **svpluscnv** is optional and adds additional plotting capabilities. It is
  not always installed automatically, and its own dependencies must be
  fetched from the CRAN archive first:

  ```r
  install.packages("https://cran.r-project.org/src/contrib/Archive/D3GB/D3GB_1.1.tar.gz",
                   repos = NULL)
  install.packages("https://cran.r-project.org/src/contrib/Archive/taRifx/taRifx_1.0.6.2.tar.gz",
                   repos = NULL)
  devtools::install_github("gonzolgarcia/svpluscnv")
  ```

## Supporting CNV algorithms and formats

`cnvstrmliner` is able to parse outputs from the following methods

| method | function | comments |
| ------- | ------  | ----     |
| [SEG](https://cnvkit.readthedocs.io/en/v0.9.9/fileformats.html#:~:text=The%20SEG%20format%20is%20the,chrom%20%E2%80%93%20chromosome%20name%20or%20ID) | parseSEG | |
| [Sclust](https://pubmed.ncbi.nlm.nih.gov/29844525/) | parseSCLUST | Requires updateSCLUSTclone for clone specific information |
| [ControlFreec](http://bioinformatics.oxfordjournals.org/content/28/3/423.long) | parseCfreec  | subclonal information is buggy, thus not supported in cnvstrmlinr |
| [FACETS](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5027494/) | parseFACETS | ensemblCall can be used to merge multiple FACETS runs in a single object (suitable for unmatched tumor samples) |
| [PureCN](https://scfbm.biomedcentral.com/articles/10.1186/s13029-016-0060-z) | parsePureCN || 
| [GATK CopyNumber](https://gatk.broadinstitute.org/hc/en-us/articles/360035535892-Somatic-copy-number-variant-discovery-CNVs) | parseGATK | |
| [Battenberg](https://github.com/Wedge-lab/battenberg) | parseBattenberg | |
| [Sequenza](https://sequenzatools.bitbucket.io/#/home) | parseSequenza | |
| [CloneCNA](https://pmc.ncbi.nlm.nih.gov/articles/PMC4990858/) | parseCloneCNA | |
| [CNVkit](https://cnvkit.readthedocs.io/en/stable/) | parseCNVkit ||
| [Conseca](https://github.com/GerkeLab/conseca) | parseConseca | Aggregated CNV calls combining ControlFreec and GATK results |
| TSV format | parseTSV | Tab-delimited files with columns: chromosome, start, end, cn1, cn2 |
| VCF format | parseVCF | Generic VCF files with CNV information in header and INFO fields |
| [manta](https://pubmed.ncbi.nlm.nih.gov/26647377/) | parseManta | |

For examples of use see the vignette of the package and the help pages of each function.

## Object description

Each cnvojb contains a GRanges slot (accessible with obj@segments) with normalized set of columns: 

* total.allele.copies: the total copies found in a segment (major.allele.copies + minor.allele.copies)
* major.allele.copies: the copies of the major allele
* minor.allele.copies: the copies of the minor allele
* cnv.totalmarkers: total number of markers SNPs used to test this segment
* cnv.markers: number of heterozygous SNPs used to test this segment
* cnv.annotation: annotation of the segment (NEUT, AMP, DEL etc)
* cnv.lrr: the log2 ratio of that segment. 
* cnv.cf: This is the cell fraction that has the event. Typically the purity of the tumor sample equals the highest cnv.cf.
* cnv.cff: This is the Cancer cell fraction. Events appearing in all cancer cells (CCF = 1) are considered clonal. To find the prevalence (cn.cf) of these cells in the sample divide this by the purity.
* cnv.loh: This is set to TRUE if the segment has lost its heterozygosity. Depending on the total copy number it can be a deletion, neutral LOH or even LOH in an amplified region. 

Additional slots contain the information for 
* purity: is the predicted overall purity of the sample.
* ploidy: is the ploidy of the sample
* effective_ploidy: is the integer ploidy of the sample which is either provided by the tool itself or reflects the state with the majority of basepairs.
* sex: the sex of the sample. 

Finally if the function applyQCmetrics has been run an additional slot called cnaqcobj is provided with additional metrics about the quality of the CNA calls.




## Exported Functions

### Parsers
- **parseBattenberg**: Parse Battenberg copy number output files integrating information from `_subclones.txt` and `_cellularity_ploidy.txt` (or `_rho_and_psi.txt`) to create a `cnvstrml` object.
- **parseCfreec**: Parse Control-FREEC output, combining information from `_info.txt`, `_CNVs`, `_ratio.txt`, and `_BAF.txt` files to construct a `cnvstrml` object with gap filling.
- **parseCNVkit**: Parse CNVkit copy number output (typically `.cns` files from the `call` command), with automatic inference of allele-specific copy numbers when missing.
- **parseConseca**: Parse Conseca copy number output (aggregated tool combining ControlFreec and GATK results), accepting data frames or files.
- **parseFACETS**: Parse FACETS copy number output from VCF or data frame format, extracting somatic CNVs, purity, and ploidy with subclone support.
- **parseGATK**: Parse GATK CopyNumber `.seg` files, inferring allele-specific copy numbers from log2 ratio values using heuristic cutoffs.
- **parseManta**: Parse Manta structural variant VCF files, extracting breakpoint ranges to create a `cnvstrml` object (no purity/ploidy information).
- **parsePureCN**: Parse PureCN output, integrating information from `<filename>.csv` (purity, ploidy) and `<filename>_loh.csv` (segments).
- **parseSCLUST**: Parse SCLUST output by combining information from `_cn_summary.txt`, `_subclonal_cn.txt`, and `_allelic_states.txt` with subclonal support.
- **parseSEG**: Parse generic SEG format files (`.seg` or `.cbs`), interpreting the last column as either copy numbers or logR values.
- **parseSequenza**: Parse Sequenza output by combining information from `_seqzinfo.txt`, `_confints_CP.txt`, and `_segments.txt`.
- **parseTSV**: Parse simple tab-delimited files with columns: `chromosome`, `start`, `end`, `cn1`, and `cn2` for custom CNV data.
- **parseCNV**: Universal parser that automatically dispatches to appropriate method-specific parser based on file format and content.
- **loadCNV**: Alias for parseCNV.
- **readCNV**: Alias for parseCNV.
- **parseVCF**: Parse VCF files containing CNV information.

### Annotation

- **annotateCNV**: Annotate copy number segments with descriptive CNV type labels (HOMDEL, HETDEL, NEUT, AMP, etc.) based on allele copies and ploidy.
- **addChromosomeAnnotation**: Add chromosome, arm (p/q), centromere, and telomere annotations to segments, splitting segments that cross boundaries.
- **addLOH**: Annotate loss of heterozygosity (LOH) status based on minor allele copy number, with correct handling of sex chromosomes.
- **addCCF**: Calculate and add Cancer Cell Fraction (CCF) to segments based on cellular fraction and sample purity.

### Consensus & Ensemble

- **findConsensus**: Find consensus CNV segments across multiple samples or methods, creating regions of consistent annotation with caller-specific metadata.
- **applyConsensus**: Apply pre-computed consensus segmentation to individual CNV objects for standardized breakpoint alignment.
- **ensemblCall**: Create ensemble CNV calls from multiple callers using majority voting, with ploidy adjustment and consensus value calculation.
- **mergeSegments**: Merge adjacent segments with identical copy number, annotation, and LOH status, using median/weighted averaging for values.
- **updateBattenbergClone**: Extract and create a new `cnvstrml` object for a specific Battenberg subclone.
- **updateSCLUSTClone**: Extract and create a new `cnvstrml` object for a specific SCLUST subclone.

### Comparison & Metrics

- **compareSegmentStatus**: Compare CNV status between two samples with detailed per-segment TP/FP/FN metrics at annotation, copy number, and CCF levels.
- **compareSegmentOverlap**: Compare segment overlaps between two `cnvstrml` objects and annotate both segment sets with best-match overlap metrics: `jaccard.index` (maximum per-segment Jaccard overlap) and `jaccard.status` (CNV annotation state associated with that maximum; `NO_OVERLAP` if no overlap exists).
- **getJaccardIndex**: Compute altered-genome Jaccard overlap statistics from a segment comparison table, including total altered bp in truth/prediction and matched fractions.
- **cnvMetrics**: Calculate pairwise Sørensen-Dice similarity metrics between CNV profiles using the CNVMetrics package, with optional state simplification.
- **cnvEntropy**: Calculate conditional entropy for CNV patterns across samples to quantify profile complexity and fragmentation.
- **computeStateSizes**: Compute total genomic size (bp) per CNV state (NEUTRAL, AMP, HETDEL, HOMDEL) for a `cnvstrml` object.
- **computeStateWeights**: Compute state-specific weights from query/subject state sizes for weighted metric aggregation.

### Segmentation comparison metrics

A family of functions that score one segmentation against another (a known
truth, or simply another caller) and, where no truth exists, score a family of
segmentations against each other. All of them accept either a `cnvstrml` object
or a bare `GRanges` of segments.

Most of them take a `bins` argument, which defines **where the assay could see
anything**: pass the capture target regions for exome/panel data, and a
`GenomicRanges::tileGenome()` grid for WGS. Bins with no assigned segment are
never imputed to neutral, so an empty segmentation cannot earn a good score by
saying nothing.

See the `metrics` vignette for a worked example of every function below, with a
figure per metric and guidance on when each is appropriate.

**State agreement (needs a reference)**

- **computeStateKappa**: Cohen's kappa between the copy-number state two segmentations assign to each bin. Chance-corrected, so an all-neutral caller scores exactly 0 regardless of how much of the reference is genuinely neutral. The default single-number quality score.
- **scoreNullCaller**: Score a synthetic all-neutral caller against a reference. Diagnostic rather than a metric: `null_p_o` is the raw agreement that doing no work achieves on that sample, and is the number every raw-agreement figure must be read against.
- **computeStateConfusion**: The bin-level confusion matrix kappa is read off, plus per-state recall, precision, F1 and support. Says *which* states were confused for which.
- **computeEventRecovery**: Splits the two error directions kappa merges — recall over genuinely altered reference territory, and precision of the query's altered calls. Width-weighted; segments carrying no bin are excluded from both numerator and denominator rather than counted as misses.

**Boundary and partition structure (state-blind)**

- **getSegmentBreakpoints**: Every segment start/end in a segmentation, de-duplicated, as width-1 `GRanges`. Chromosome-end positions are dropped by default since they record where the segmentation starts and stops rather than a detected transition.
- **computeBoundaryConcordance**: Precision, recall and F1 for breakpoint detection, counting a query breakpoint as correct when a reference breakpoint lies within `tolerance_bp` (default 300 kb).
- **computeSegmentationARI**: Adjusted Rand index between two segmentations treated as partitions of genomic bins. `partition_by = "segment"` scores segment *identity* (sensitive to over/under-segmentation); `partition_by = "state"` scores copy-number state.
- **computeSegmentationVI**: Variation of information between the two partitions — a true metric, and one that decomposes into an over-segmentation term and an under-segmentation term. Use the normalised `nvi` for anything compared across grids.
- **computeFragmentation**: Whether the query splits the reference or merges it, and by how much (`split_ratio`, `merge_ratio`, and bin-weighted mean partner counts). The direction ARI and VI cannot give.

**Per-segment annotation**

- **assignBinState**: Label every bin with the state of the segment it overlaps most. The projection every bin-level metric is built on; exported so a fixed reference can be projected once and reused across a sweep.
- **annotateBinStateConcordance**: Attach `state.concordance`, `state.n_matched` and `state.n_total` to each segment — the fraction of its bins whose state matches a reference. Size-neutral and state-aware.
- **annotateBinJaccard**: Attach `bin.jaccard.index`, `bin.jaccard.status` and the three bin counters to the segments on *both* sides, giving each segment's overlap with its best-matching counterpart measured in bins rather than base pairs (so capture deserts do not inflate it).
- **annotateBreakpointConsistency**: Attach `breakpoint.consistency.start` / `.end` to each segment from a consistency table, plus a one-row summary. Says how much of a segmentation rests on well-supported boundaries.

**Truth-free, across a family of segmentations**

- **computeBreakpointConsistency**: Cluster the boundaries of several segmentations of the same sample and report, per cluster, the fraction of segmentations placing a boundary there. Each segmentation votes at most once per cluster.
- **computeCallerAgreement**: Pairwise overlap of the altered-bin sets across a family of callers, plus summary statistics including a per-bin vote entropy. No reference required.
- **computeConsensusAgreement**: Leave-one-out agreement of each caller with the majority verdict of the others, with separate `captured`/`supported` directions so over- and under-calling are distinguishable.

  Note: agreement measures typicality, not correctness. Where callers share a failure mode, or where purity is too low for any caller to see anything, the consensus is wrong and these metrics are actively misleading.

**Rollups and the one-call driver**

- **summariseSegmentJaccard**: Roll per-segment Jaccard columns (from `annotateBinJaccard` or `compareSegmentOverlap`) up to one row per group: median, mean, width-weighted mean, fraction above a threshold, and fraction with no counterpart at all.
- **summariseStateConcordance**: Roll per-segment state concordance up to one row per group, reporting both the pooled (bin-weighted) accuracy and the mean of the per-segment fractions — they differ sharply and neither substitutes for the other.
- **summariseBreakpointConsistency**: Reduce a breakpoint consistency table to scalars (`mean_consistency`, `frac_high`, `frac_unanimous`, ...), optionally within groups.
- **segmentationMetrics**: Run all of the above genome-wide metrics for one query/subject pair and flatten them into a single row. The returned column set never depends on the input — a failed or skipped metric contributes `NA` — so `dplyr::bind_rows()` over a whole sweep cannot produce a ragged frame. This is the function analysis code should normally call.

### Quality Control

- **applyQCMetrics**: Apply CNAqc quality control pipeline to assess CNV call quality using somatic mutation VAF concordance.
- **detectOverfragmentation**: Detect chromosome arms with excessive short segments using binomial tests with Bonferroni correction.
- **filterLowQualitySegments**: Remove or flag segments failing quality criteria (e.g., marker support, size/noise thresholds) before downstream analysis.
- **qcTable**: Generate comprehensive QC summary table with metrics from CNAqc analysis.

### Plotting

- **plotSegments**: Visualize CNV segments across the genome with color-coded annotations, copy number tracks, and optional region highlighting for single or multiple samples.
- **plotGRanges**: Low-level plotting function for GRanges objects with flexible value and color mapping, ploidy baseline, and chromosome faceting.
- **plotCNVFreq**: Plot CNV frequency heatmap across multiple samples showing recurrent amplification and deletion patterns.

### Gap Filling

- **fillGaps**: Fill genomic gaps between CNV segments with neutral-state segments to create continuous genome-wide profiles.

### Export & Save

- **saveSEG** / **writeSEG**: Export CNV data to SEG format (tab-delimited) for IGV. The `segmean` argument selects the value written in the sixth column — absolute copy number (default), the caller's own log2 ratio (optionally purity-corrected), or the ploidy-relative log2 ratio — and `style = "gistic"` switches to GISTIC's column names. For GISTIC input prefer `exportGISTIC`; note the default writes absolute copy numbers, which GISTIC would misread as amplifications.
- **exportGISTIC**: Export CNV data to a tab-delimited GISTIC2.0 segment file (`Sample`, `Chromosome`, `Start`, `End`, `Num_Markers`, `Seg.CN`), where `Seg.CN` is the ploidy-relative log2 copy ratio `log2(total.allele.copies / ploidy)`. A thin preset over `saveSEG(segmean = "ploidy", minSegCN = -5, dropNA = TRUE, style = "gistic")`. Accepts a list of objects to build a cohort file.
- **writeVCF** / **saveVCF**: Export CNV data to bgzipped and tabix-indexed VCF 4.2 format.
- **savePhyloWGS**: Export CNV segments for PhyloWGS subclonal reconstruction analysis with NCBI chromosome format and purity adjustment.
- **saveAnnotSV**: Export CNV segments to AnnotSV-compatible BED format for structural variant annotation.

### Accessors & Utilities

- **getPloidy**: Retrieve ploidy value from a `cnvstrml` object (raw fractional or rounded effective ploidy).
- **getPurity**: Retrieve tumor purity value from a `cnvstrml` object.
- **getSex**: Retrieve sex information (M/F) from a `cnvstrml` object.
- **getSupportedMethods**: Get list of all CNV calling methods supported by the package parsers.
- **summaryTable**: Generate comprehensive summary statistics at genome, chromosome, or arm levels including CNV burden and aberration frequencies.
- **matchCNVSV**: Match CNV segment boundaries with structural variant breakpoints to validate and annotate CNV calls.

### External Metrics (require non-cnvstrml inputs)

- **computeCoverage**: Compute per-target coverage and genome-normalized log2 ratios from a BAM over BED/GRanges targets. Requires external BAM and BED/GRanges inputs.
- **computeGini**: Compute Gini coefficient of coverage uniformity from a coverage table (e.g., CNVkit `.cnr`). Requires a data frame with a `depth` column.
- **computeMAPD**: Compute Median Absolute Pairwise Difference (MAPD) from log2 ratio data (e.g., CNVkit `.cnr`). Requires columns `chromosome` and `log2`.
- **computeBIC**: Compute Bayesian Information Criterion (BIC) score for a segmentation model.
- **computeAdaptiveBIC**: Compute adaptive BIC variant for segmentation quality/model selection.
- **computeCVIC**: Compute cross-validated information criterion (CVIC) for segmentation comparison.
- **computeMBICseg**: Compute modified BIC metric tailored for segmented copy-number profiles.
- **computeMSQBIC**: Compute mean-squared-error-based BIC variant for segmentation model assessment.



## Utility scripts

### postCNV.R

**Location:** `inst/extdata/postCNV.R`

**Purpose:** Comprehensive post-processing script for CNV calls from various tools. Acts as a workhorse for standardizing CNV outputs and computing quality control metrics.

**Key Features:**
- **Format Conversion:** Converts CNV calls from any supported method to standardized VCF 4.2 format (bgzipped and tabix-indexed)
- **Quality Control:** Computes CNAqc metrics when somatic mutations VCF is provided
- **CNV Entropy:** Calculates genomic fragmentation metrics using covered regions from BED files
- **Coverage Metrics:** Computes Gini coefficient and MAPD (Median Absolute Pairwise Difference) from BAM files
- **Multi-format Output:** Saves results as VCF, RDS (cnvstrml object), and QC statistics tables

**Usage:**
```bash
./postCNV.R -i <input_file> -m <method> -o <output.vcf.gz> [options]
```

**Options:**
- `-i, --input`: Path to the input CNV file (required)
- `-m, --method`: CNV calling method - FACETS|Sequenza|SCLUST|Battenberg|FREEC|CNVkit|PureCN|GATK|Conseca|CloneCNA|TSV|VCF|SEG (required)
- `-o, --output`: Path to the output VCF file (required)
- `-s, --sample`: Sample name (optional, inferred from input if not provided)
- `-v, --somatic_mutations`: Somatic mutations VCF for CNAqc quality control (optional)
- `-b, --bam`: BAM file for Gini and MAPD coverage metrics calculation (optional)
- `-d, --bed`: BED file with covered regions for CNV entropy and coverage metrics (optional)
- `-h, --help`: Show help message

**Output Files:**
1. **`<sample>-<method>.vcf.gz`** - Standardized CNV calls in VCF 4.2 format (bgzipped)
2. **`<sample>-<method>.vcf.gz.tbi`** - Tabix index for the VCF file
3. **`<sample>-<method>.rds`** - R data file containing the cnvstrml object
4. **`<sample>-<method>_qc_stats.tsv`** - QC metrics table (if somatic mutations provided)
5. **`<sample>-<method>_entropy.tsv`** - CNV entropy metrics (if BED file provided)
6. **`<sample>-<method>_coverage_metrics.tsv`** - Gini and MAPD metrics (if BAM file provided)

**Example:**
```bash
# Basic conversion
./postCNV.R -i sample1_facets.vcf.gz -m FACETS -o sample1_processed.vcf.gz

# With full QC metrics
./postCNV.R -i sample1.cns -m CNVkit -o sample1_cnv.vcf.gz \
  -v sample1_somatic.vcf.gz \
  -b sample1.bam \
  -d targets.bed \
  -s SAMPLE001
```

---

### annotateGenes.R

**Location:** `inst/extdata/annotateGenes.R`

**Purpose:** Annotate CNV segments with overlapping gene information using Ensembl annotations from AnnotationHub.

**Key Features:**
- **Gene-Level Annotation:** Maps CNV segments to individual genes with copy number status. For genes overlapping multiple segments the weighted copy number is provided according to:
- 
  $$\text{Gene CN} = \frac{\sum (\text{Segment CN}_i \times \text{Overlap Width}_i)}{\text{Total Overlap Width}}$$

- **Segment-Level Annotation:** Provides segment-wise CNV information with genotype representation
- **Ensembl Integration:** Uses AnnotationHub with Ensembl database release 109 (hg38/GRCh38)
- **LOH Detection:** Annotates loss of heterozygosity including copy-neutral LOH events
- **Standardized Output:** Generates tab-delimited files compatible with downstream analysis tools

**Usage:**
```bash
./annotateGenes.R -i <input.vcf.gz> -c <cache_dir> -o <output_dir> [options]
```

**Options:**
- `-i, --input`: Path to CNV calls in VCF format from cnvstrmliner (required)
- `-c, --cache`: Path to AnnotationHub cache directory (required, or set via `ANNOTATION_HUB_CACHE` environment variable)
- `-o, --output_directory`: Output directory for annotation files (default: current directory)
- `-s, --sample`: Sample name (optional, extracted from VCF if not provided)
- `-h, --help`: Show help message

**Output Files:**
1. **`<sample>-<method>.gene_level.tsv`** - Gene-level CNV annotations
   - Columns: `chr`, `start`, `end`, `gene_name`, `Status`, `Ratio`, `Absolute.Copy.Number`, `Sample.Tumor.Purity`, `LOH`
   - Only includes genes with non-neutral CNV status
   - Status categories: Deep deletion, Shallow deletion, Amplification, Neutral
   - LOH categories: "copy-neutral LOH", "LOH: loss of heterozygosity", or empty

2. **`<sample>-<method>.segment_level.tsv`** - Segment-level CNV annotations
   - Columns: `chr`, `start`, `end`, `Status`, `Ratio`, `Absolute.Copy.Number`, `Sample.Tumor.Purity`, `LOH`, `Genotype`
   - Genotype format: allelic representation (e.g., "AAB" = 2 copies of A allele, 1 copy of B allele)
   - Excludes neutral segments for focused analysis

**Prerequisites:**
- AnnotationHub cache must be pre-downloaded with Ensembl release 109 for Homo sapiens
- Input VCF must be generated by cnvstrmliner (contains required metadata)

**Example:**
```bash
# Set cache directory and annotate CNVs
export ANNOTATION_HUB_CACHE=/path/to/annotation_cache
./annotateGenes.R -i sample1_cnv.vcf.gz -c $ANNOTATION_HUB_CACHE -o ./annotations/

# With explicit sample name
./annotateGenes.R -i sample1_cnv.vcf.gz -c /data/annotation_hub -o ./output/ -s SAMPLE001
```

**Output Example (gene_level.tsv):**
```
chr    start      end        gene_name    Status              Ratio    Absolute.Copy.Number    Sample.Tumor.Purity    LOH
1      100000     200000     GENE1        Amplification       0.58     5                       0.85                   
3      500000     600000     GENE2        Shallow deletion    -0.32    1                       0.85                   LOH: loss of heterozygosity
7      1000000    1200000    GENE3        Deep deletion       -2.15    0                       0.85                   
```