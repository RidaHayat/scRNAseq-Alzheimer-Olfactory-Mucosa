#!/usr/bin/env Rscript

# Descriptive screen of genes highlighted in a recent independent olfactory
# epithelium AD study. This does not test disease effects: the local data have
# one sample per condition and lack cell-type annotations.
suppressPackageStartupMessages(library(Matrix))

input_dir <- Sys.getenv("DATA_DIR", "data")
output_dir <- file.path("results")
figure_dir <- file.path(output_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

gene_panel <- data.frame(
  gene = c("AKR1B1", "CERT1", "LBH", "HLA-DQA2", "EGR3", "GZMA", "MIF"),
  reported_context = c(
    "Olfactory sensory neurons", "Olfactory sensory neurons",
    "Olfactory sensory neurons", "Myeloid cells", "Myeloid cells",
    "Neuron-immune interaction", "Neuron-immune interaction"
  ),
  stringsAsFactors = FALSE
)

read_filtered_counts <- function(sample_name) {
  folder <- file.path(input_dir, sample_name)
  feature_path <- file.path(folder, "genes.tsv")
  if (!file.exists(feature_path)) feature_path <- file.path(folder, "features.tsv")
  barcode_path <- file.path(folder, "barcodes.tsv")
  matrix_path <- file.path(folder, "matrix.mtx")
  missing <- c(feature_path, barcode_path, matrix_path)[
    !file.exists(c(feature_path, barcode_path, matrix_path))
  ]
  if (length(missing)) stop("Missing 10x input file(s) under ", folder)

  features <- read.delim(feature_path, header = FALSE, stringsAsFactors = FALSE,
                         check.names = FALSE)
  symbols <- if (ncol(features) >= 2L) as.character(features[[2]]) else
    as.character(features[[1]])
  barcodes <- readLines(barcode_path, warn = FALSE)
  counts <- Matrix::readMM(matrix_path)
  if (nrow(counts) != length(symbols) || ncol(counts) != length(barcodes)) {
    stop("Feature/barcode dimensions do not match matrix for ", sample_name)
  }
  counts <- methods::as(counts, "CsparseMatrix")

  # Match CreateSeuratObject(min.cells=3, min.features=100) and the project's
  # QC rule (200 < detected genes < 5,000; mitochondrial fraction < 10%).
  keep_genes <- Matrix::rowSums(counts > 0) >= 3
  counts <- counts[keep_genes, , drop = FALSE]
  symbols <- symbols[keep_genes]
  detected <- Matrix::colSums(counts > 0)
  library_size <- Matrix::colSums(counts)
  mito <- grepl("^MT-", symbols)
  mito_counts <- if (any(mito)) Matrix::colSums(counts[mito, , drop = FALSE]) else
    rep(0, ncol(counts))
  percent_mito <- 100 * mito_counts / pmax(library_size, 1)
  keep_cells <- detected > 200 & detected < 5000 & percent_mito < 10
  counts <- counts[, keep_cells, drop = FALSE]
  list(counts = counts, symbols = symbols, library_size = library_size[keep_cells])
}

summarize_sample <- function(sample_name) {
  x <- read_filtered_counts(sample_name)
  rows <- lapply(seq_len(nrow(gene_panel)), function(i) {
    gene <- gene_panel$gene[i]
    idx <- which(x$symbols == gene)
    if (!length(idx)) {
      detected <- rep(NA_real_, ncol(x$counts))
      pct_detected <- mean_expression <- NA_real_
    } else {
      detected <- Matrix::colSums(x$counts[idx, , drop = FALSE])
      normalized <- log1p(10000 * detected / pmax(x$library_size, 1))
      pct_detected <- 100 * mean(detected > 0)
      mean_expression <- mean(normalized)
    }
    data.frame(
      sample = sample_name,
      condition_label = if (sample_name == "AD1") "AD1" else "Control2",
      gene = gene,
      reported_context = gene_panel$reported_context[i],
      qc_pass_cells = ncol(x$counts),
      pct_cells_detected = pct_detected,
      mean_log1p_normalized_expression = mean_expression,
      gene_found = length(idx) > 0L,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

results <- do.call(rbind, lapply(c("AD1", "Control2"), summarize_sample))
write.csv(results, file.path(output_dir, "external_signature_screen.csv"),
          row.names = FALSE)

pct <- matrix(results$pct_cells_detected, nrow = 2, ncol = nrow(gene_panel),
              byrow = TRUE,
              dimnames = list(c("AD1", "Control2"), gene_panel$gene))
available <- colSums(is.na(pct)) == 0
png(file.path(figure_dir, "published_gene_screen.png"),
    width = 2200, height = 1400, res = 220)
par(mar = c(10, 5, 4, 1))
barplot(pct[, available, drop = FALSE], beside = TRUE, las = 2,
        col = c("#B44A4A", "#4277A8"),
        ylab = "QC-passed cells with detected transcript (%)",
        main = "Published olfactory AD genes in the two local samples",
        ylim = c(0, max(5, max(pct, na.rm = TRUE) * 1.15)))
legend("topright", legend = c("AD1", "Control2"),
       fill = c("#B44A4A", "#4277A8"), bty = "n")
mtext("Genes absent from the input matrix are omitted. Detection is descriptive; n=1 sample/group.",
      side = 1, line = 8, cex = 0.9)
dev.off()

cat("Wrote results/external_signature_screen.csv and",
    "results/figures/published_gene_screen.png\n")
