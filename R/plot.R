if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    "arm", "chrom_arms", "chrom.length", "cnv.annotation", "cnv.ccf", "end",
    "freq.gain", "freq.loss", "name",
    "segment_mid", "segment_width_bp", "seqnames", "seqlengths", "start", "state",
    "transparency", "value", "yvalues"
  ))
}


# Truncate a label to at most `width` characters, marking truncation with an ellipsis.
# `width` may be Inf, which leaves the label untouched.
.abbrevLabel <- function(x, width) {
  if (!is.finite(width) || nchar(x) <= width) {
    return(x)
  }
  paste0(substr(x, 1, width - 1), "\u2026")
}


#' Plot CNV segments across the genome with genomic coordinates and copy number status
#'
#' @description
#' Visualizes copy-number variation (CNV) segments along chromosomes with segments colored by
#' annotation status (amplification, deletion, neutral, LOH, etc.). Supports single or multiple
#' samples with customizable y-axis copy-number range, optional highlight regions, and chromosome
#' subsetting. Automatically computes subtitle with sample metadata (method, purity, ploidy, QC).
#'
#' @param cnvList A \code{cnvstrml} object or a list of \code{cnvstrml} objects. If a list,
#'   plots are stacked vertically with unified axes and legend.
#' @param highlightgr Optional \code{GRanges} or \code{cnvstrml} object specifying regions to
#'   highlight with a semi-transparent overlay (default: \code{NULL}, no highlight).
#' @param max.copies Integer maximum copy number to display on y-axis (default: 6).
#'   Segments with copy numbers exceeding this value are capped at \code{max.copies}.
#'   Y-axis range: 0 to \code{max.copies}.
#' @param main Character plot title (default: \code{NULL}, uses sample name from \code{cnvList@sample}).
#' @param subtitle Character plot subtitle. If \code{NULL}, auto-generated from sample metadata
#'   (method, purity, ploidy, effective ploidy, QC confidence if available).
#' @param chromosomes Character or character vector of chromosomes to plot
#'   (default: \code{NULL}, all standard chromosomes). Useful for focusing on specific regions.
#' @param strip Logical; if \code{TRUE} (default), displays chromosome names at the top of each facet.
#'   Set to \code{FALSE} to hide strip labels.
#' @param tiny.segment.frac Numeric fraction of chromosome span used to detect tiny segments
#'   that may be visually lost as lines (default: \code{0.001}). Segments with width
#'   \code{<= tiny.segment.frac * chromosome_length} are additionally marked with a diamond
#'   at segment midpoint. Set to \code{0} to disable midpoint diamonds. If \code{NULL},
#'   the threshold is estimated automatically from the number of plotted chromosomes
#'   (\code{auto_tiny_frac = 1.5 / (3000 / n_chromosomes)}).
#'
#' @param sample.label.width Integer maximum number of characters of the sample name shown in
#'   the per-sample label drawn beside each panel when more than three samples are plotted
#'   (default: 24). Longer names are truncated with an ellipsis so the label gutter cannot
#'   grow without bound and squeeze the panels. Set to \code{Inf} to disable truncation.
#'
#' @return A \code{ggplot2} object (or \code{patchwork} composition for multi-sample input).
#'   For single samples, returns a faceted plot with one facet per chromosome showing
#'   segments as colored horizontal bars with heights proportional to copy number.
#'   For multiple samples, returns a vertically stacked collection of plots with shared legend.
#'
#' @details
#' \strong{Segment coloring}
#' Segments are colored by \code{cnv.annotation} status according to a fixed palette:
#' \itemize{
#'   \item AMP (brown4): amplifications.
#'   \item AMP-LOH (cyan4): amplifications with loss of heterozygosity.
#'   \item NEUT (chartreuse3): neutral (copy-number unchanged).
#'   \item NEUT-LOH (cyan): neutral with LOH.
#'   \item HETDEL (blue): heterozygous deletion (single-copy loss).
#'   \item HETDEL-LOH (cyan2): heterozygous deletion with LOH.
#'   \item HOMDEL (darkblue): homozygous deletion.
#'   \item INFER (yellow): inferred segment (e.g., from gap-filling).
#' }
#'
#' \strong{Multi-sample visualization}
#' When \code{cnvList} is a list, plots are generated recursively for each sample and
#' combined using \code{patchwork::wrap_plots()}. Only the first sample shows strip labels;
#' subsequent plots suppress them for clarity. Legend and axes are collected across all plots.
#' With three or fewer samples each panel keeps its own title and subtitle. With more than
#' three samples the titles are dropped and each panel is instead identified by a compact
#' horizontal label in the left gutter reading \code{sample} on the first line and
#' \code{method (purity)} on the second; sample names are truncated at
#' \code{sample.label.width} characters. Because \code{patchwork} sizes the axis-title
#' gutter across the whole composition, panels stay vertically aligned even when the result
#' is stacked with other plots.
#'
#' \strong{Copy-number capping}
#' Segments with \code{total.allele.copies > max.copies} are capped for visualization,
#' though their original values are preserved. The y-axis range is \[0, max.copies\].
#'
#' \strong{Consensus sample restriction}
#' If \code{cnvList@segtype == "consensus"}, the function returns \code{NULL} with a message,
#' as consensus objects lack single-sample segment attribution.
#'
#' \strong{Highlight regions}
#' If \code{highlightgr} is provided (\code{GRanges} or \code{cnvstrml}), matching regions
#' are shaded with a light overlay, useful for highlighting recurrent mutations, target regions, or
#' alternate caller results.
#'
#' @examples
#' \dontrun{
#' # Plot single sample
#' p <- plotSegments(cnv_sample)
#' print(p)
#'
#' # Plot with title and custom copy-number range
#' p <- plotSegments(
#'   cnv_sample,
#'   max.copies = 8,
#'   main = "Sample S1 - WGS",
#'   subtitle = "Custom analysis"
#' )
#'
#' # Plot multiple samples stacked
#' samples <- list(sample1, sample2, sample3)
#' p <- plotSegments(samples, max.copies = 6)
#'
#' # Highlight specific chromosomes
#' p <- plotSegments(cnv_sample, chromosomes = c("chr1", "chr2", "chr17"))
#'
#' # Highlight regions of interest
#' rois <- GRanges(seqnames = "chr17", IRanges(start = 37844393, end = 37884915))
#' p <- plotSegments(cnv_sample, highlightgr = rois)
#' }
#'
#' @seealso \code{\link{plotGRanges}} for low-level plotting, \code{\link{cnvEntropy}} for
#'   fragmentation assessment.
#'
#' @export
plotSegments <- function(cnvList, highlightgr = NULL, max.copies = 6, main = NULL, subtitle = NULL, chromosomes = NULL, strip = TRUE, tiny.segment.frac = NULL,
                         sample.label.width = 24) {
  if (is.list(cnvList)) {
    n_samples <- length(cnvList)

    pp <- lapply(seq_along(cnvList), function(idx) {
      X <- cnvList[[idx]]
      strip <- TRUE
      if (idx > 1) {
        strip <- FALSE
      }

      # Keep original behavior for one to three panels.
      if (n_samples <= 3) {
        return(cnvstrmliner::plotSegments(X,
          highlightgr = highlightgr,
          max.copies = max.copies,
          chromosomes = chromosomes,
          main = main,
          subtitle = subtitle,
          strip = strip,
          tiny.segment.frac = tiny.segment.frac
        ))
      }

      p_i <- cnvstrmliner::plotSegments(X,
        highlightgr = highlightgr,
        max.copies = max.copies,
        chromosomes = chromosomes,
        main = "",
        subtitle = "",
        strip = strip,
        tiny.segment.frac = tiny.segment.frac
      )

      meta_label <- paste0(
        .abbrevLabel(X@sample, sample.label.width), "\n",
        X@method,
        " (", round(X@purity, digits = 2), ")"
      )

      p_i +
        ggplot2::ylab(meta_label) +
        ggplot2::theme(
          axis.title.y = ggplot2::element_text(
            angle = 0,
            hjust = 1,
            vjust = 0.5,
            size = ggplot2::rel(0.8),
            lineheight = 0.9
          )
        )
    })

    return(patchwork::wrap_plots(pp,
      ncol = 1
    ) +
      patchwork::plot_annotation(main, theme = ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))) +
      patchwork::plot_layout(
        axis_titles = if (n_samples > 3) "keep" else "collect",
        guides = "collect"
      ))
  }

  if (!is(cnvList, "cnvstrml")) {
    stop("This function plots only cnvstrml objects. You provided a ", class(cnvList))
  }
  aa <- cnvList@segments
  if (cnvList@segtype == "consensus") {
    message("This function can plot only segments for a single sample")
    return(NULL)
  }
  if (is.null(main)) {
    main <- cnvList@sample
  }
  if (is.null(subtitle)) {
    subtitle <- paste("method:", cnvList@method, "Purity:", round(cnvList@purity, digits = 2), "Ploidy:", round(cnvList@ploidy, digits = 2), "Eff Ploidy:", cnvList@effective_ploidy)
    qc <- NULL
    qc <- tryCatch(cnvList@quality$cnaqc$qc.coeff, error = function(e) {
      return(e)
    })
    if (!is.null(qc)) {
      subtitle <- paste(subtitle, "QC confidence (0-1):", qc)
    }
  }
  aa <- aa %>%
    GenomeInfoDb::keepStandardChromosomes(pruning.mode = "coarse")

  aa$total.allele.copies <- ifelse(aa$total.allele.copies > max.copies, max.copies, aa$total.allele.copies)
  aa$cnv.annotation <- ifelse((aa$cnv.loh == TRUE & aa$cnv.annotation == "HETDEL" & aa$total.allele.copies > 1), "HETDEL-LOH", aa$cnv.annotation)

  if (!is.null(highlightgr) && is(highlightgr, "cnvstrml")) {
    hl <- highlightgr@segments
  } else {
    hl <- highlightgr
  }

  colors <- c(
    "AMP" = "brown4", "AMP-LOH" = "cyan4",
    "NEUT" = "chartreuse3", "NEUT-LOH" = "cyan", "INFER" = "yellow",
    "HETDEL" = "blue", "HOMDEL" = "darkblue", "HETDEL-LOH" = "cyan2"
  )
  p <- cnvstrmliner::plotGRanges(
    gr = aa,
    ploidy = cnvList@effective_ploidy,
    highlightgr = hl,
    chromosomes = chromosomes,
    valuecol = "total.allele.copies",
    colorcol = "cnv.annotation",
    max.copies = max.copies,
    colors = colors,
    main = main,
    subtitle = subtitle,
    strip = strip,
    tiny.segment.frac = tiny.segment.frac
  )
  p + ggplot2::labs(y = "Copy number")
}

#' Plot the segments using CanvasExpress
#'
#' @param cnvList A cnvstrml object or list of cnvstrml objects used to create the plot.
#' @param ggplotObj A ggplot object created previously.
#' @param seqinfo Seqinfo object describing the plotted genome.
#' @param arms GRanges object with chromosome arm coordinates.
#' @description
#' if CanvasExpress is installed we can plot the ggplot object created previously
#'
#'

plotSegmentsCE <- function(cnvList, ggplotObj, seqinfo, arms) {
  plot_list <- if (is.list(ggplotObj)) ggplotObj else list(ggplotObj)

  d0 <- lapply(seq_along(plot_list), function(x) {
    d <- plot_list[[x]]$data %>% dplyr::select("seqnames", "start", "end", "total.allele.copies", "state")
    d["method"] <- plot_list[[x]]$labels$title
    d
  }) %>% dplyr::bind_rows()

  seqs <- d0 %>%
    dplyr::pull("seqnames") %>%
    unique()
  d0 <- data.frame(lapply(d0, function(x) {
    gsub("chr", "", x)
  }))

  slengths <- GenomeInfoDb::seqlengths(seqinfo)[which(GenomeInfoDb::seqnames(seqinfo) %in% seqs)]
  names(slengths) <- sapply(names(slengths), function(x) {
    gsub("chr", "", x)
  })
  armslengths <- end(arms[which(GenomeInfoDb::seqnames(arms) %in% seqs)])[seq(1, length(arms), by = 2)]
  names(armslengths) <- GenomeInfoDb::seqnames(arms[which(GenomeInfoDb::seqnames(arms) %in% seqs)])[seq(1, length(arms), by = 2)]
  names(armslengths) <- sapply(names(armslengths), function(x) {
    gsub("chr", "", x)
  })
  canvasXpress::canvasXpress(
    d0,
    graphType = "Scatter2D",
    scatterType = "manhattanCNV",
    colorBy = "state",
    segregateVariablesBy = "method",
    layoutTopology = "7X1",
    manhattanMarkerChromosomeNumber = "seqnames",
    manhattanRangeStart = "start",
    manhattanRangeEnd = "end",
    manhattanCNV = "total.allele.copies",
    manhattanRangeMiddleShow = FALSE,
    manhattanRangeLineWidth = 10,
    chromosomeLengths = as.list(slengths),
    armLengths = as.list(armslengths),
    xAxisTextScaleFontFactor = 0.6,
    yAxisTextScaleFontFactor = 0.6,
    zoomYDisable = TRUE
  )
}


#' Plot a GRanges object as genome-wide segments with flexible coloring and value mapping
#'
#' @description
#' Low-level plotting function for genomic ranges (\code{GRanges}). Maps segment values
#' (e.g., copy number) to y-coordinates and segment metadata (e.g., CNV status) to colors.
#' Produces a faceted genome browser-style plot with optional ploidy baseline, centromeric shading,
#' and region highlights. Used internally by \code{plotSegments()} and suitable for custom
#' genomic data visualization.
#'
#' @param gr A \code{GRanges} object with standardized metadata columns: at minimum
#'   \code{seqnames}, \code{start}, \code{end}. Optional columns: \code{cnv.annotation},
#'   \code{total.allele.copies}, \code{arm}, etc.
#' @param ploidy Numeric baseline copy number to plot as a horizontal dotted green line
#'   (default: \code{NULL}, no line). Typically 2 for diploid; highlights deviations.
#' @param highlightgr Optional \code{GRanges} specifying regions to shade with a light
#'   overlay (default: \code{NULL}, no highlight).
#' @param chromosomes Character or character vector of chromosomes to include
#'   (default: \code{NULL}, all standard chromosomes). Values are matched as-is;
#'   must respect seqlevels style (e.g., chr1 vs 1).
#' @param valuecol Character name of the metadata column to map to y-axis values
#'   (default: \code{NULL}, plots y = 0 for all segments). Common: \code{"total.allele.copies"},
#'   \code{"cnv.lrr"}, or any numeric column.
#' @param colorcol Character name of the metadata column to map to segment colors
#'   (default: \code{NULL}, uses default color). Common: \code{"cnv.annotation"} for status,
#'   \code{"arm"} for p/q arm. Values mapped via \code{colors} parameter.
#' @param colors Named character vector mapping \code{colorcol} values to RGB/named colors
#'   (default: \code{NULL}, uses ggplot2 automatic palette). Example:
#'   \code{c(HOMDEL = "darkblue", HETDEL = "blue", NEUT = "green", AMP = "red")}.
#' @param max.copies Integer; if non-\code{NULL}, caps y-values and y-axis limit at this value
#'   (default: \code{NULL}, axis scales to data range).
#' @param main Character plot title (default: \code{NULL}, no title).
#' @param subtitle Character plot subtitle (default: \code{NULL}, no subtitle).
#' @param strip Logical; if \code{TRUE} (default), shows chromosome facet labels rotated 90°.
#'   Set to \code{FALSE} to hide labels.
#' @param tiny.segment.frac Numeric fraction of chromosome span used to detect tiny segments
#'   that may be visually lost as lines (default: \code{0.001}). Segments with width
#'   \code{<= tiny.segment.frac * chromosome_length} are additionally marked with a diamond
#'   at segment midpoint. Set to \code{0} to disable midpoint diamonds. If \code{NULL},
#'   the threshold is estimated automatically from the number of plotted chromosomes
#'   (\code{auto_tiny_frac = 1.5 / (3000 / n_chromosomes)}).
#'
#' @return A \code{ggplot2} object; a faceted plot with one facet per chromosome.
#'   Segments rendered as horizontal bars (\code{geom_segment}) with positions determined
#'   by genomic coordinates and heights by \code{valuecol}. Facets scaled to chromosome
#'   size (free x-scales). No x-axis tick labels or breaks (genomic positions readable via
#'   segment positions).
#'
#' @details
#' \strong{Coordinate system}
#' \itemize{
#'   \item X-axis: genomic position (bp), scaled per chromosome.
#'   \item Y-axis: values from \code{valuecol} (default 0 if not specified).
#'   \item Facets: one per chromosome, with \code{space = "free_x"} to avoid empty gaps.
#' }
#'
#' \strong{Segment rendering}
#' Uses \code{geom_segment} with line width 4 to draw segments from (start, 0) to (end, value).
#' Segments are colored by \code{state} (mapped from \code{colorcol}).
#'
#' \strong{Centromeric/arm shading}
#' If arms are available (loaded via \code{loadCentromere()}), chromosome arms are shaded
#' with low alpha (0.1) to distinguish p and q arms visually.
#'
#' \strong{Highlight overlay}
#' If \code{highlightgr} is provided, regions are shaded with \code{geom_rect} spanning
#' from -Inf to +Inf on the y-axis (full height), alpha = 0.1.
#'
#' \strong{Chromosome subsetting}
#' Standard chromosomes are kept; sex chromosomes and contigs are handled via
#' \code{GenomeInfoDb::keepStandardChromosomes()}. Optionally filter further with \code{chromosomes}.
#'
#' @examples
#' \dontrun{
#' # Basic plot: y = total copy number, color by annotation
#' p <- plotGRanges(
#'   gr = segments,
#'   ploidy = 2,
#'   valuecol = "total.allele.copies",
#'   colorcol = "cnv.annotation",
#'   colors = c(AMP = "red", NEUT = "green", DEL = "blue"),
#'   max.copies = 6
#' )
#' print(p)
#'
#' # Custom plot: y = logR, color by arm, no ploidy line
#' p <- plotGRanges(
#'   gr = segments,
#'   valuecol = "cnv.lrr",
#'   colorcol = "arm",
#'   colors = c(p = "lightblue", q = "darkblue"),
#'   main = "Copy-ratio by arm"
#' )
#'
#' # Plot subset of chromosomes with highlights
#' rois <- GRanges(seqnames = "chr17", IRanges(37844393, 37884915))
#' p <- plotGRanges(
#'   gr = segments,
#'   ploidy = 2,
#'   valuecol = "total.allele.copies",
#'   colorcol = "cnv.annotation",
#'   chromosomes = c("chr1", "chr17"),
#'   highlightgr = rois
#' )
#' }
#'
#' @seealso \code{\link{plotSegments}} for convenience wrapper handling single/multi-sample
#'   \code{cnvstrml} objects, \code{\link{loadCentromere}} for arm annotations.
#'
#' @export
#'
plotGRanges <- function(gr,
                        ploidy = NULL,
                        highlightgr = NULL,
                        chromosomes = NULL,
                        valuecol = NULL,
                        colorcol = NULL,
                        colors = NULL,
                        max.copies = NULL,
                        main = NULL,
                        subtitle = NULL,
                        strip = TRUE,
                        tiny.segment.frac = NULL) {
  if (!is.null(chromosomes)) {
    if (!is.list(chromosomes)) {
      chromosomes <- c(chromosomes)
    }
  }


  aa <- gr %>%
    GenomeInfoDb::keepStandardChromosomes(pruning.mode = "coarse") %>%
    GenomeInfoDb::sortSeqlevels() %>%
    BiocGenerics::sort()

  if (!is.null(chromosomes)) {
    aa <- aa[which(as.vector(GenomeInfoDb::seqnames(aa)) %in% as.vector(chromosomes))]
  }

  n_chromosomes <- length(unique(as.character(GenomeInfoDb::seqnames(aa))))
  if (is.null(tiny.segment.frac)) {
    # Approximate a 3000 px total plotting width and mark segments below ~1.5 px.
    tiny.segment.frac <- if (n_chromosomes > 0) 1.5 / (3000 / n_chromosomes) else 0
  }
  if (!is.numeric(tiny.segment.frac) || length(tiny.segment.frac) != 1 || is.na(tiny.segment.frac) || tiny.segment.frac < 0) {
    stop("tiny.segment.frac must be NULL or a single non-negative numeric value")
  }

  if (!is.null(highlightgr)) {
    highlightgr <- highlightgr %>%
      GenomeInfoDb::keepStandardChromosomes(pruning.mode = "coarse") %>%
      GenomeInfoDb::sortSeqlevels() %>%
      BiocGenerics::sort()
    if (!is.null(chromosomes)) {
      highlightgr <- highlightgr[which(as.vector(GenomeInfoDb::seqnames(highlightgr)) %in% as.vector(chromosomes))]
    }
  }


  if (!is.null(valuecol)) {
    yvalues <- S4Vectors::mcols(aa)[, valuecol]
  } else {
    yvalues <- rep(0, length(aa))
  }

  if (!is.null(colorcol)) {
    .state <- S4Vectors::mcols(aa)[, colorcol]
  } else {
    .state <- rep("NEUT", length(aa))
  }

  if (!is.null(max.copies)) {
    yupperlim <- max.copies
  } else {
    yupperlim <- max(yvalues)
  }
  ylowerlim <- min(c(0, yvalues))
  if (length(unique(yvalues)) == 1) {
    ylowerlim <- 0
  }
  arms <- loadCentromere(gr)
  if (!is.null(chromosomes)) {
    arms <- arms[which(as.vector(GenomeInfoDb::seqnames(arms)) %in% as.vector(chromosomes))]
  }
  # S4Vectors::mcols(aa)[,'alpha']=ifelse( S4Vectors::mcols(aa)[,'arm']=='p', 0.6, 0.8)
  # S4Vectors::mcols(aa)[,'alpha']=ifelse( is.na(S4Vectors::mcols(aa)[,'arm']),0.3,S4Vectors::mcols(aa)[,'alpha'] )
  # message("Processing ", length(aa), " segments for plotting...")
  plot_df <- aa %>%
    as.data.frame(row.names = NULL) %>%
    dplyr::mutate(
      seqlengths = GenomeInfoDb::seqlengths(aa)[GenomeInfoDb::seqnames(aa) %>% as.character()],
      chrom.length = dplyr::if_else(
        !is.na(seqlengths) & seqlengths > 0,
        as.numeric(seqlengths),
        as.numeric(stats::ave(end, seqnames, FUN = max))
      ),
      values = yvalues,
      state = factor(.state, levels = names(colors)),
      cnv.annotation = factor(cnv.annotation, levels = names(colors)),
      yvalues = ifelse(yvalues > yupperlim, yupperlim, yvalues),
      transparency = pmax(0.3, ifelse(is.na(as.numeric(cnv.ccf)), 1, as.numeric(cnv.ccf))),
      segment_width_bp = pmax(1, end - start + 1),
      tiny_segment = tiny.segment.frac > 0 & segment_width_bp <= pmax(1, chrom.length * tiny.segment.frac),
      segment_mid = (start + end) / 2
    ) %>%
    tidyr::complete(state)

  has_tiny_segment <- any(plot_df$tiny_segment %in% TRUE, na.rm = TRUE)
  tiny_df <- plot_df[plot_df$tiny_segment %in% TRUE, , drop = FALSE]

  p <- plot_df %>%
    ggplot2::ggplot(ggplot2::aes(y = yvalues, x = start)) +
    ggplot2::geom_blank(ggplot2::aes(x = 0), na.rm = TRUE) +
    ggplot2::geom_blank(ggplot2::aes(x = seqlengths), na.rm = TRUE) +
    ggplot2::geom_rect(
      data = if (!is.null(arms)) as.data.frame(arms) else NULL,
      ggplot2::aes(
        xmin = start,
        xmax = end,
        ymin = -Inf,
        ymax = Inf,
        color = NA,
        fill = arm
      ),
      linewidth = 0,
      alpha = 0.1,
      inherit.aes = FALSE, na.rm = TRUE
    ) +
    ggplot2::geom_rect(
      data = if (!is.null(highlightgr)) as.data.frame(highlightgr) else NULL,
      ggplot2::aes(
        xmin = start,
        xmax = end,
        ymin = -Inf,
        ymax = Inf
      ),
      alpha = 0.1,
      inherit.aes = FALSE, na.rm = TRUE
    ) +
    ggplot2::geom_segment(ggplot2::aes(yend = yvalues, xend = end, color = state, alpha = transparency), linewidth = 4, na.rm = TRUE) +
    ggplot2::ylim(ylowerlim, yupperlim) +
    ggplot2::facet_grid(~seqnames, space = "free_x", scales = "free_x") +
    ggplot2::labs(title = main, subtitle = subtitle, y = NULL, x = NULL) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      legend.position = "right"
    )
  if (has_tiny_segment) {
    p <- p + ggplot2::geom_point(
      data = tiny_df,
      ggplot2::aes(x = segment_mid, y = yvalues, color = state, alpha = transparency, shape = "Small segments"),
      size = 2.4,
      na.rm = TRUE,
      show.legend = c(shape = TRUE, color = FALSE, alpha = FALSE)
    )
  }
  if (!is.null(ploidy)) {
    p <- p + ggplot2::geom_hline(yintercept = ploidy, linetype = "dotted", color = "green")
  }
  if (strip == TRUE) {
    p <- p + ggplot2::theme(strip.text.x.top = ggplot2::element_text(angle = 90))
  } else {
    p <- p + ggplot2::theme(strip.text = ggplot2::element_blank())
  }

  if (!is.null(colors)) {
    p <- p + ggplot2::scale_color_manual(
      values = colors,
      drop = FALSE
    )
  }
  if (has_tiny_segment) {
    p <- p + ggplot2::scale_shape_manual(values = c("Small segments" = 18), name = NULL)
  }
  p <- p + ggplot2::scale_alpha_identity(guide = "none")
  p
}


#' Plot the CNV frequencies across a population of samples
#'
#' @description
#' Generates a frequency plot showing the proportion of samples with copy number
#' gains (amplifications) and losses (deletions) across the genome. This visualization
#' is commonly used to identify recurrent copy number alterations in a cohort of samples,
#' which may indicate driver events or regions under selective pressure.
#'
#' @details
#' This function creates a genome-wide visualization where:
#' \itemize{
#'   \item \strong{Gains} (amplifications) are plotted as positive values (red bars)
#'   \item \strong{Losses} (deletions) are plotted as negative values (blue bars)
#'   \item The y-axis shows both the absolute number of samples and the percentage
#'   \item Dotted horizontal lines indicate the 50% threshold
#'   \item Chromosome arms (p and q) are distinguished by color shade
#' }
#'
#' \subsection{Input Requirements}{
#' The input must be a consensus cnvstrmliner object created by \code{\link{findConsensus}}
#' with \code{ignore.samples = TRUE}. This aggregates copy number information across
#' multiple samples into a single object where each segment contains copy number values
#' for all samples.
#'
#' For best results, it is recommended to run \code{\link{fillGaps}} on each sample
#' before creating the consensus to ensure complete genome coverage:
#' \preformatted{
#' cnvList <- lapply(cnvList, fillGaps)
#' consensus <- findConsensus(cnvList, ignore.samples = TRUE)
#' plotCNVFreq(consensus)
#' }
#' }
#'
#' \subsection{Ploidy Adjustment}{
#' When \code{ploidy = TRUE} (default), the function adjusts for sample-specific ploidy
#' by subtracting each sample's median copy number from its values. This centers the
#' data around 0, making gains positive and losses negative relative to each sample's
#' baseline. This is useful when samples have different overall ploidy levels.
#'
#' When \code{ploidy = FALSE}, the raw copy numbers are used with 2 (diploid) as the
#' neutral state. Segments with copy number > 2 are counted as gains, and < 2 as losses.
#' }
#'
#' @param cnvObj A cnvstrmliner object of type "consensus" created by
#'   \code{\link{findConsensus}} with \code{ignore.samples = TRUE}. The object must
#'   contain aggregated copy number data across multiple samples.
#' @param ploidy Logical. If \code{TRUE} (default), adjusts for sample-specific ploidy
#'   by centering each sample's copy numbers around its median. If \code{FALSE}, uses
#'   diploid (CN = 2) as the neutral baseline for all samples.
#' @param main Character string for the plot title. If \code{NULL} (default), uses
#'   "CNV frequency" as the title.
#' @param subtitle Character string for the plot subtitle. If \code{NULL} (default),
#'   no subtitle is displayed.
#'
#' @return A ggplot2 object displaying the CNV frequency plot. The plot can be further
#'   customized using standard ggplot2 functions or saved using \code{ggsave()}.
#'
#' @note
#' \itemize{
#'   \item This function only accepts consensus cnvstrmliner objects. Attempting to
#'     use a single-sample object will result in an error.
#'   \item The plot uses a dual y-axis showing both absolute sample counts (left) and
#'     percentages (right).
#'   \item Chromosome arm annotations (p/q) are automatically loaded and used to
#'     shade the bars differently for each arm.
#'   \item Standard chromosomes only are displayed (sex chromosomes and non-standard
#'     contigs are excluded).
#' }
#'
#' @seealso
#' \code{\link{findConsensus}} for creating consensus objects from multiple samples,
#' \code{\link{fillGaps}} for filling gaps between segments,
#' \code{\link{plotSegments}} for plotting individual sample CNV profiles
#'
#' @examples
#' \dontrun{
#' # Load multiple CNV samples
#' cnv1 <- parseFACETS("sample1_facets.vcf.gz")
#' cnv2 <- parseFACETS("sample2_facets.vcf.gz")
#' cnv3 <- parseFACETS("sample3_facets.vcf.gz")
#' cnvList <- list(cnv1, cnv2, cnv3)
#'
#' # Fill gaps for complete coverage
#' cnvList <- lapply(cnvList, fillGaps)
#'
#' # Create consensus object
#' consensus <- findConsensus(cnvList, ignore.samples = TRUE)
#'
#' # Plot frequency with ploidy adjustment (default)
#' plotCNVFreq(consensus)
#'
#' # Plot frequency without ploidy adjustment
#' plotCNVFreq(consensus, ploidy = FALSE)
#'
#' # Customize title and subtitle
#' plotCNVFreq(consensus,
#'   main = "CNV Frequency in Tumor Cohort",
#'   subtitle = "N = 50 samples"
#' )
#'
#' # Save the plot
#' p <- plotCNVFreq(consensus)
#' ggplot2::ggsave("cnv_frequency.pdf", p, width = 12, height = 4)
#' }
#'
#' @export
#' @importFrom ggplot2 ggplot aes geom_rect geom_blank geom_hline facet_grid
#' @importFrom ggplot2 scale_color_manual scale_y_continuous labs theme element_blank element_text unit


plotCNVFreq <- function(cnvObj, ploidy = TRUE, main = NULL, subtitle = NULL) {
  neutral <- 2
  # fill in the gaps
  # cnvList=lapply(cnvList, function(X){cnvstrmliner::fillGaps(X)})
  # cnvObj=cnvstrmliner::findConsensus( cnvList, method='quantitative', ignore.samples = TRUE)

  if (!startsWith(cnvObj@segtype, prefix = "consensus")) {
    logger::log_error("Function plotCNVFreq accepts as input a cnvstrmliner object from findConsensus")
    stop("Please provide a consensus cnvstrml object")
  }

  aa <- cnvObj@segments %>%
    GenomeInfoDb::keepStandardChromosomes(pruning.mode = "coarse") %>%
    GenomeInfoDb::sortSeqlevels() %>%
    BiocGenerics::sort()

  arms <- loadCentromere(aa)


  # find the frequencies of loss and gain
  df <- S4Vectors::mcols(aa) %>%
    as.data.frame() %>%
    dplyr::select(-tidyselect::starts_with("stats"))
  if (ploidy == TRUE) {
    # find ploidy per sample
    medploidy <- apply(df, 2, function(x) {
      median(x, na.rm = TRUE)
    })
    df <- df - medploidy
    neutral <- 0
  } else {
    neutral <- 2
  }
  freq.loss <- apply(df, 1, function(X) {
    sum(X < neutral)
  })
  freq.gain <- apply(df, 1, function(X) {
    sum(X > neutral)
  })


  totalSamples <- ncol(df)
  df <- df %>%
    dplyr::mutate(
      freq.loss = -freq.loss,
      freq.gain = freq.gain
    ) %>%
    dplyr::select(freq.loss, freq.gain)
  S4Vectors::mcols(aa) <- df
  ov <- GenomicRanges::findOverlaps(query = aa, subject = arms)
  S4Vectors::mcols(aa)[S4Vectors::queryHits(ov), "arm"] <- S4Vectors::mcols(arms)[S4Vectors::subjectHits(ov), "arm"]


  # dplyr::mutate(name=paste( name, arm)) %>%

  colors <- c(
    "freq.loss p" = "blue1", "freq.loss q" = "blue4", "freq.loss u" = "blue2", "freq.loss" = "blue2",
    "freq.gain p" = "firebrick1", "freq.gain q" = "firebrick4", "freq.gain u" = "firebrick2", "freq.gain" = "firebrick2"
  )
  p <- aa %>%
    as.data.frame(row.names = NULL) %>%
    dplyr::mutate(seqlengths = GenomeInfoDb::seqlengths(aa)[GenomeInfoDb::seqnames(aa) %>% as.character()]) %>%
    tidyr::pivot_longer(cols = c("freq.loss", "freq.gain")) %>%
    dplyr::mutate(name = paste(name, arm)) %>%
    ggplot2::ggplot(ggplot2::aes(group = name), na.rm = TRUE) +
    ggplot2::geom_blank(ggplot2::aes(x = 0), na.rm = TRUE) +
    ggplot2::geom_blank(ggplot2::aes(x = seqlengths), na.rm = TRUE) +
    ggplot2::geom_rect(ggplot2::aes(ymin = 0, ymax = value, xmin = start, xmax = end, fill = name), na.rm = TRUE) +
    ggplot2::scale_color_manual(values = colors, aesthetics = c("colour", "fill")) +
    ggplot2::scale_y_continuous("# samples",
      limits = c(-totalSamples - 2, totalSamples + 2), expand = c(0, 0), breaks = c(-totalSamples, -floor(totalSamples / 2), 0, ceiling(totalSamples / 2), totalSamples),
      sec.axis = ggplot2::sec_axis(name = "% samples", trans = ~ . / totalSamples * 100)
    ) +
    ggplot2::facet_grid(~seqnames, space = "free_x", scales = "free_x") +
    ggplot2::labs(title = ifelse(!is.null(main), main, "CNV frequency"), subtitle = subtitle, y = "loss", x = "chromosome") +
    ggplot2::geom_hline(ggplot2::aes(yintercept = totalSamples / 2), linetype = "dotted") +
    ggplot2::geom_hline(ggplot2::aes(yintercept = -totalSamples / 2), linetype = "dotted") +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      panel.spacing = unit(2, "points"),
      legend.position = "none",
      strip.text.x.top = ggplot2::element_text(angle = 90),
      axis.title.y.right = ggplot2::element_text(angle = 90)
    )


  p +
    ggplot2::ylab("Total copies") +
    ggplot2::xlab("Position on genome")
}

.centromere_cache <- new.env(parent = emptyenv())

#' helper function to load the centromere locations
#'
#' @param gr A GRanges object used to infer the genome build and retrieve centromere locations.
#'
#' @return A \code{GRanges} object with the centromere and chromosome arm
#'   annotations for the genome build of \code{gr}, or \code{NULL} if the build is
#'   not supported.
#'
#' @keywords internal
loadCentromere <- function(gr) {
  cache_env <- .centromere_cache

  # get the genome from the seqinfo
  genome <- GenomeInfoDb::genome(GenomeInfoDb::seqinfo(gr)) %>%
    unique() %>%
    as.character()
  genome_key <- genome[1]

  if (exists(genome_key, envir = cache_env, inherits = FALSE)) {
    return(get(genome_key, envir = cache_env, inherits = FALSE))
  }
  # if(!is.null(genome)){
  #   # load( file.path(
  #   #   system.file('extdata', package = 'cnvstrmliner'),
  #   #   'seqinfo',
  #   #   paste0( genome, ".gaps.rds")
  #   # ))
  #   # we work with gaps.centromeres
  #   gaps.centromere=chrom_gaps[[genome]][['centromere']]
  #   arms=GenomicRanges::reduce(gaps.centromere, min.gapwidth=1000000)%>% GenomicRanges::gaps()
  #   arms=arms[ GenomicRanges::strand( arms) == "*"]
  #   arms$arm=c("p","q")
  # }

  # if(!is.null( arms )){
  #   ov=GenomicRanges::findOverlaps( query=gr, subject=arms)
  #   S4Vectors::mcols(gr)[ S4Vectors::queryHits(ov) , 'arm']= S4Vectors::mcols(arms)[S4Vectors::subjectHits(ov), 'arm']
  # }else{
  #   S4Vectors::mcols(gr)[ S4Vectors::queryHits(ov) , 'arm']=NA
  # }

  arms <- chrom_arms[[genome_key]]
  assign(genome_key, arms, envir = cache_env)
  arms
}
