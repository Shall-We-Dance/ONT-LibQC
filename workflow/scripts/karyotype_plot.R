#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (!(length(args) %in% c(8, 9))) {
  stop("Usage: karyotype_plot.R <reads.tsv> <chromosomes.tsv> <genes.bed> <output.svg> <title> <window_bp> <mapq> <assembly> [gene_window_bp]")
}

reads_path <- args[[1]]
chromosomes_path <- args[[2]]
genes_path <- args[[3]]
output_path <- args[[4]]
plot_title <- args[[5]]
window_size <- as.numeric(args[[6]])
mapq <- as.integer(args[[7]])
assembly <- args[[8]]
gene_window_size <- if (length(args) >= 9) as.numeric(args[[9]]) else 250000
if (!is.finite(window_size) || window_size <= 0 || !is.finite(gene_window_size) || gene_window_size <= 0) stop("Window sizes must be positive")

suppressPackageStartupMessages({
  library(karyoploteR)
  library(GenomicRanges)
})

chromosomes <- read.delim(chromosomes_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
chromosomes$size <- as.numeric(chromosomes$size)
chromosome_names <- chromosomes$chromosome

genome <- GRanges(
  seqnames = chromosome_names,
  ranges = IRanges(start = rep(1, nrow(chromosomes)), end = chromosomes$size)
)
seqlengths(genome)[chromosome_names] <- chromosomes$size

reads <- read.delim(reads_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
if (nrow(reads) > 0) {
  mapped <- GRanges(
    seqnames = reads$chromosome,
    ranges = IRanges(start = reads$position, width = 1)
  )
  seqlevels(mapped) <- chromosome_names
  seqlengths(mapped)[chromosome_names] <- chromosomes$size
} else {
  mapped <- GRanges()
}

genes <- read.delim(
  genes_path, header = FALSE, sep = "\t", stringsAsFactors = FALSE,
  quote = "", comment.char = ""
)
if (ncol(genes) < 3) {
  stop("Gene BED must contain at least chromosome, start and end columns")
}
genes <- genes[genes[[1]] %in% chromosome_names, , drop = FALSE]
if (nrow(genes) > 0) {
  all.genes <- GRanges(
    seqnames = genes[[1]],
    ranges = IRanges(start = as.numeric(genes[[2]]) + 1, end = as.numeric(genes[[3]]))
  )
  seqlevels(all.genes) <- chromosome_names
  seqlengths(all.genes)[chromosome_names] <- chromosomes$size
} else {
  all.genes <- GRanges()
}

plot_height <- min(10, max(5.2, 3.8 + length(chromosome_names) * 0.22))
svg(output_path, width = 14, height = plot_height, pointsize = 9, bg = "white", onefile = FALSE)
on.exit(dev.off(), add = TRUE)
# Reserve an outer footer so the bin-size caption cannot overlap chr1 density.
par(oma = c(2.2, 0, 0, 0))
# Plot type 1 keeps the mapped-read signal in its own panel. Plot type 6 maps
# every data panel onto the ideogram, which makes gene and read density overlap.
plot_params <- getDefaultPlotParams(plot.type = 1)
plot_params$leftmargin <- 0.11
plot_params$rightmargin <- 0.04
plot_params$topmargin <- 55
plot_params$bottommargin <- 42
plot_params$ideogramheight <- 32
plot_params$data1height <- 78
plot_params$data1inmargin <- 7
plot_params$data1outmargin <- 5
kp <- plotKaryotype(
  genome = genome,
  chromosomes = chromosome_names,
  plot.type = 1,
  plot.params = plot_params,
  main = NULL,
  labels.plotter = NULL
)
kpAddChromosomeNames(kp, cex = 0.68, col = "#34495E")
kpAddMainTitle(kp, main = plot_title, cex = 1.05)
kpDataBackground(kp, data.panel = "ideogram", color = "#FFFFFFAA")
if (length(all.genes) > 0) {
  kpPlotDensity(
    kp, all.genes, window.size = gene_window_size, data.panel = "ideogram",
    col = "#3388FF", border = "#3388FF", r0 = 0.5, r1 = 1
  )
  kpPlotDensity(
    kp, all.genes, window.size = gene_window_size, data.panel = "ideogram",
    col = "#3388FF", border = "#3388FF", r0 = 0.5, r1 = 0
  )
}
kpAddBaseNumbers(
  kp,
  tick.dist = max(1, round(max(chromosomes$size) / 5, -max(0, nchar(round(max(chromosomes$size) / 5)) - 2))),
  minor.ticks = FALSE,
  add.units = TRUE,
  cex = 0.34,
  tick.col = "#607080"
)
kpDataBackground(kp, data.panel = 1, r0 = 0.08, r1 = 0.92, color = "#F5F8FA")
if (length(mapped) > 0) {
  kpPlotDensity(
    kp,
    data = mapped,
    window.size = window_size,
    data.panel = 1,
    r0 = 0.08,
    r1 = 0.92,
    col = "#E4575666",
    border = "#A82030",
    lwd = 0.65
  )
}
window_label <- if (window_size >= 1e6) sprintf("%.2g Mb", window_size / 1e6) else sprintf("%.0f kb", window_size / 1e3)
mtext(
  sprintf("%s | blue ideogram: gene density (%s) | red panel: mapped-read density (%s) | %s reads | MAPQ >= %d", assembly, sprintf("%.0f kb", gene_window_size / 1000), window_label, format(length(mapped), big.mark = ","), mapq),
  side = 1,
  line = 0.7,
  outer = TRUE,
  cex = 0.72,
  col = "#536475"
)
if (length(mapped) == 0) {
  mtext("No MAPQ-pass alignments on configured primary chromosomes", side = 1, line = 1, cex = 0.75, col = "#687986")
}
