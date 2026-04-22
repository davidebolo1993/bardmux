#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

usage <- function() {
  cat(
    "Usage: plot_cb_compare.R --summary SUMMARY.tsv --out-prefix PREFIX [--title TITLE]\n",
    "\n",
    "Required:\n",
    "  --summary FILE     Summary from compare_cb_assignments\n",
    "\n",
    "Optional:\n",
    "  --out-prefix STR   Output prefix (default: cb_compare)\n",
    "  --title STR        Figure title prefix\n",
    sep = ""
  )
}

parse_args <- function(a) {
  cfg <- list(summary = NULL, out_prefix = "cb_compare", title = "CB Comparison")
  i <- 1
  while (i <= length(a)) {
    key <- a[[i]]
    if (key %in% c("-h", "--help")) return(NULL)
    if (key == "--summary" && i + 1 <= length(a)) {
      cfg$summary <- a[[i + 1]]
      i <- i + 2
      next
    }
    if (key == "--out-prefix" && i + 1 <= length(a)) {
      cfg$out_prefix <- a[[i + 1]]
      i <- i + 2
      next
    }
    if (key == "--title" && i + 1 <= length(a)) {
      cfg$title <- a[[i + 1]]
      i <- i + 2
      next
    }
    stop(sprintf("Unknown or incomplete argument: %s", key))
  }
  if (is.null(cfg$summary)) return(NULL)
  cfg
}

split_blocks <- function(lines) {
  blocks <- list()
  cur <- character(0)
  for (ln in lines) {
    if (nchar(trimws(ln)) == 0) {
      if (length(cur) > 0) {
        blocks[[length(blocks) + 1]] <- cur
        cur <- character(0)
      }
    } else {
      cur <- c(cur, ln)
    }
  }
  if (length(cur) > 0) blocks[[length(blocks) + 1]] <- cur
  blocks
}

read_tsv_block <- function(block) {
  con <- textConnection(paste(block, collapse = "\n"))
  on.exit(close(con), add = TRUE)
  read.delim(con, sep = "\t", header = TRUE, stringsAsFactors = FALSE,
             check.names = FALSE, quote = "")
}

find_block <- function(blocks, header_prefix) {
  for (b in blocks) {
    if (length(b) >= 1 && startsWith(b[[1]], header_prefix)) return(b)
  }
  NULL
}

get_metric <- function(metrics, key, default = NA_real_) {
  idx <- which(metrics$metric == key)
  if (length(idx) == 0) return(default)
  as.numeric(metrics$value[[idx[[1]]]])
}

format_count_pct <- function(counts) {
  total <- sum(counts)
  pct <- if (total > 0) 100 * counts / total else rep(0, length(counts))
  sprintf("%s\n(%.1f%%)", format(counts, big.mark = ",", scientific = FALSE), pct)
}

cfg <- tryCatch(parse_args(args), error = function(e) {
  message(e$message)
  NULL
})
if (is.null(cfg)) {
  usage()
  quit(status = 1)
}

lines <- readLines(cfg$summary, warn = FALSE)
blocks <- split_blocks(lines)

metrics_block <- find_block(blocks, "metric\tvalue")
if (is.null(metrics_block)) {
  stop("Could not find metric block in summary")
}
metrics <- read_tsv_block(metrics_block)

edit_block <- find_block(blocks, "strict_edit_distance\t")
wf_status_block <- find_block(blocks, "wf_only_bardmux_status\t")

edit_df <- NULL
if (!is.null(edit_block) && length(edit_block) > 1) {
  edit_df <- read_tsv_block(edit_block)
}
wf_status_df <- NULL
if (!is.null(wf_status_block) && length(wf_status_block) > 1) {
  wf_status_df <- read_tsv_block(wf_status_block)
}

intersection <- get_metric(metrics, "intersection_assigned_ids", 0)
bardmux_only <- get_metric(metrics, "bardmux_only_assigned_ids", 0)
wf_only <- get_metric(metrics, "wf_only_assigned_ids", 0)
strict_concordant <- get_metric(metrics, "strict_concordant_ids", 0)
strict_discordant <- get_metric(metrics, "strict_discordant_ids", 0)
strict_concordance <- get_metric(metrics, "strict_concordance_pct", 0)
jaccard <- get_metric(metrics, "assigned_id_jaccard_pct", 0)
bardmux_vs_wf <- get_metric(metrics, "bardmux_vs_wf_assigned_coverage_pct", 0)
wf_vs_bardmux <- get_metric(metrics, "wf_vs_bardmux_assigned_coverage_pct", 0)

out_png <- paste0(cfg$out_prefix, ".overview.png")
png(out_png, width = 1800, height = 1300, res = 150)
par(mfrow = c(2, 2), mar = c(6, 5, 4, 2) + 0.1)

# Panel 1: Assigned-set composition
counts1 <- c(intersection, bardmux_only, wf_only)
names1 <- c("Intersection", "Bardmux only", "wf only")
cols1 <- c("#2a9d8f", "#457b9d", "#e76f51")
bp1 <- barplot(counts1, col = cols1, ylim = c(0, max(counts1, 1) * 1.28), las = 2,
               main = "Assigned Read-ID Set Composition", ylab = "Distinct read IDs")
text(bp1, counts1, labels = format_count_pct(counts1), pos = 3, cex = 0.85)

# Panel 2: Strict concordance
counts2 <- c(strict_concordant, strict_discordant)
names2 <- c("Concordant", "Discordant")
cols2 <- c("#2a9d8f", "#d62828")
bp2 <- barplot(counts2, col = cols2, ylim = c(0, max(counts2, 1) * 1.3), las = 2,
               main = sprintf("Strict CB Concordance (%.2f%%)", strict_concordance),
               ylab = "Intersection IDs (strict-comparable)")
text(bp2, counts2, labels = format_count_pct(counts2), pos = 3, cex = 0.85)

# Panel 3: wf-only status breakdown
if (!is.null(wf_status_df) && nrow(wf_status_df) > 0) {
  ord <- order(wf_status_df$count, decreasing = TRUE)
  wf_status_df <- wf_status_df[ord, , drop = FALSE]
  counts3 <- as.numeric(wf_status_df$count)
  labels3 <- wf_status_df$wf_only_bardmux_status
  cols3 <- rep("#6c757d", length(counts3))
  bp3 <- barplot(counts3, col = cols3, ylim = c(0, max(counts3, 1) * 1.32), las = 2,
                 main = "wf-only Assigned IDs: bardmux Status", ylab = "Distinct read IDs")
  pct3 <- ifelse(is.na(wf_status_df$pct_of_wf_only), 0, wf_status_df$pct_of_wf_only)
  lab3 <- sprintf("%s\n(%.1f%%)", format(counts3, big.mark = ",", scientific = FALSE), pct3)
  text(bp3, counts3, labels = lab3, pos = 3, cex = 0.8)
} else {
  plot.new()
  title("wf-only Assigned IDs: bardmux Status")
  text(0.5, 0.5, "No wf_only_bardmux_status section found", cex = 1)
}

# Panel 4: strict edit-distance concordance
if (!is.null(edit_df) && nrow(edit_df) > 0) {
  edit_df <- edit_df[order(edit_df$strict_edit_distance), , drop = FALSE]
  mat <- rbind(as.numeric(edit_df$concordant_ids), as.numeric(edit_df$discordant_ids))
  colnames(mat) <- as.character(edit_df$strict_edit_distance)
  bp4 <- barplot(mat, beside = FALSE, col = c("#2a9d8f", "#d62828"),
                 main = "Strict Concordance by bardmux Edit Distance",
                 xlab = "Edit distance", ylab = "Read IDs")
  legend("topright", legend = c("Concordant", "Discordant"),
         fill = c("#2a9d8f", "#d62828"), bty = "n")
} else {
  plot.new()
  title("Strict Concordance by bardmux Edit Distance")
  text(0.5, 0.5, "No strict_edit_distance section found", cex = 1)
}

mtext(sprintf("%s | Jaccard: %.2f%% | bardmux_vs_wf: %.2f%% | wf_vs_bardmux: %.2f%%",
              cfg$title, jaccard, bardmux_vs_wf, wf_vs_bardmux),
      side = 1, outer = FALSE, line = -1.2, cex = 0.9)

dev.off()

cat(sprintf("Wrote: %s\n", out_png))
