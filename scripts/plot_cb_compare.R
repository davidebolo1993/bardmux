#!/usr/bin/env Rscript

options(scipen = 999)

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
    "  --title STR        Figure title\n",
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

label_count_pct <- function(counts, total = NULL) {
  if (is.null(total)) total <- sum(counts)
  pct <- if (total > 0) 100 * counts / total else rep(0, length(counts))
  sprintf("%s\n(%.1f%%)", format(counts, big.mark = ",", scientific = FALSE), pct)
}

make_empty_plot <- function(ggplot2, title, subtitle, text) {
  ggplot2::ggplot() +
    ggplot2::theme_void(base_size = 13) +
    ggplot2::labs(title = title, subtitle = subtitle) +
    ggplot2::annotate("text", x = 0, y = 0, label = text, size = 4)
}

cfg <- tryCatch(parse_args(args), error = function(e) {
  message(e$message)
  NULL
})
if (is.null(cfg)) {
  usage()
  quit(status = 1)
}

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package 'ggplot2' is required. Install with install.packages('ggplot2').")
}
if (!requireNamespace("scales", quietly = TRUE)) {
  stop("Package 'scales' is required. Install with install.packages('scales').")
}

lines <- readLines(cfg$summary, warn = FALSE)
blocks <- split_blocks(lines)

metrics_block <- find_block(blocks, "metric\tvalue")
if (is.null(metrics_block)) {
  stop("Could not find metric block in summary")
}
metrics <- read_tsv_block(metrics_block)
mode <- metrics$value[metrics$metric == "comparison_mode"]
if (length(mode) > 0 && mode[[1]] != "bardmux_vs_wf") {
  stop("Summary is not bardmux_vs_wf mode. Use plot_bardmux_pair_compare.R for bardmux_vs_bardmux.")
}

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
strict_comparable <- get_metric(metrics, "strict_comparable_ids", 0)
strict_concordant <- get_metric(metrics, "strict_concordant_ids", 0)
strict_discordant <- get_metric(metrics, "strict_discordant_ids", 0)
strict_concordance <- get_metric(metrics, "strict_concordance_pct", 0)
jaccard <- get_metric(metrics, "assigned_id_jaccard_pct", 0)
bardmux_vs_wf <- get_metric(metrics, "bardmux_vs_wf_assigned_coverage_pct", 0)
wf_vs_bardmux <- get_metric(metrics, "wf_vs_bardmux_assigned_coverage_pct", 0)

# Panel 1
p1_df <- data.frame(
  category = factor(c("Intersection", "Bardmux only", "wf only"),
                    levels = c("Intersection", "Bardmux only", "wf only")),
  count = c(intersection, bardmux_only, wf_only)
)
p1_df$label <- label_count_pct(p1_df$count)

p1 <- ggplot2::ggplot(p1_df, ggplot2::aes(x = category, y = count, fill = category)) +
  ggplot2::geom_col(width = 0.75, color = "#1f1f1f") +
  ggplot2::geom_text(ggplot2::aes(label = label), vjust = -0.25, size = 3.8, lineheight = 0.95) +
  ggplot2::scale_fill_manual(values = c("Intersection" = "#2a9d8f", "Bardmux only" = "#457b9d", "wf only" = "#e76f51")) +
  ggplot2::scale_y_continuous(labels = scales::comma, expand = ggplot2::expansion(mult = c(0, 0.16))) +
  ggplot2::labs(
    title = "Assigned Read-ID Set Composition",
    subtitle = sprintf("Jaccard %.2f%% | bardmux_vs_wf %.2f%% | wf_vs_bardmux %.2f%%", jaccard, bardmux_vs_wf, wf_vs_bardmux),
    x = NULL,
    y = "Distinct read IDs"
  ) +
  ggplot2::theme_bw(base_size = 13) +
  ggplot2::theme(
    legend.position = "none",
    axis.text.x = ggplot2::element_text(angle = 20, hjust = 1),
    plot.title = ggplot2::element_text(face = "bold")
  )

# Panel 2
p2_df <- data.frame(
  category = factor(c("Concordant", "Discordant"), levels = c("Concordant", "Discordant")),
  count = c(strict_concordant, strict_discordant)
)
p2_df$label <- label_count_pct(p2_df$count, total = strict_comparable)

p2 <- ggplot2::ggplot(p2_df, ggplot2::aes(x = category, y = count, fill = category)) +
  ggplot2::geom_col(width = 0.7, color = "#1f1f1f") +
  ggplot2::geom_text(ggplot2::aes(label = label), vjust = -0.25, size = 3.8, lineheight = 0.95) +
  ggplot2::scale_fill_manual(values = c("Concordant" = "#2a9d8f", "Discordant" = "#d62828")) +
  ggplot2::scale_y_continuous(labels = scales::comma, expand = ggplot2::expansion(mult = c(0, 0.16))) +
  ggplot2::labs(
    title = sprintf("Strict CB Concordance (%.2f%%)", strict_concordance),
    subtitle = sprintf("Strict-comparable IDs: %s", format(strict_comparable, big.mark = ",", scientific = FALSE)),
    x = NULL,
    y = "Strict-comparable read IDs"
  ) +
  ggplot2::theme_bw(base_size = 13) +
  ggplot2::theme(
    legend.position = "none",
    plot.title = ggplot2::element_text(face = "bold")
  )

# Panel 3
if (!is.null(wf_status_df) && nrow(wf_status_df) > 0) {
  ord <- order(wf_status_df$count, decreasing = TRUE)
  wf_status_df <- wf_status_df[ord, , drop = FALSE]
  wf_status_df$wf_only_bardmux_status <- factor(wf_status_df$wf_only_bardmux_status,
                                                 levels = wf_status_df$wf_only_bardmux_status)
  wf_status_df$label <- sprintf("%s\n(%.1f%%)",
                                format(wf_status_df$count, big.mark = ",", scientific = FALSE),
                                wf_status_df$pct_of_wf_only)

  p3 <- ggplot2::ggplot(wf_status_df,
                        ggplot2::aes(x = wf_only_bardmux_status, y = count)) +
    ggplot2::geom_col(width = 0.8, fill = "#6c757d", color = "#1f1f1f") +
    ggplot2::geom_text(ggplot2::aes(label = label), vjust = -0.25, size = 3.8, lineheight = 0.95) +
    ggplot2::scale_y_continuous(labels = scales::comma, expand = ggplot2::expansion(mult = c(0, 0.16))) +
    ggplot2::labs(
      title = "wf-only Assigned IDs: bardmux Status",
      subtitle = sprintf("wf-only IDs: %s", format(wf_only, big.mark = ",", scientific = FALSE)),
      x = NULL,
      y = "Distinct read IDs"
    ) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 25, hjust = 1),
      plot.title = ggplot2::element_text(face = "bold")
    )
} else {
  p3 <- make_empty_plot(ggplot2,
                        "wf-only Assigned IDs: bardmux Status",
                        "No status section found",
                        "No wf_only_bardmux_status block found in summary")
}

# Panel 4
if (!is.null(edit_df) && nrow(edit_df) > 0) {
  edit_df <- edit_df[order(edit_df$strict_edit_distance), , drop = FALSE]

  ed <- as.character(edit_df$strict_edit_distance)
  p4_df <- data.frame(
    edit_distance = factor(rep(ed, each = 2), levels = ed),
    concordance = factor(rep(c("Concordant", "Discordant"), times = nrow(edit_df)),
                         levels = c("Concordant", "Discordant")),
    count = as.numeric(c(rbind(edit_df$concordant_ids, edit_df$discordant_ids)))
  )

  p4 <- ggplot2::ggplot(p4_df, ggplot2::aes(x = edit_distance, y = count, fill = concordance)) +
    ggplot2::geom_col(width = 0.75, color = "#1f1f1f") +
    ggplot2::scale_fill_manual(values = c("Concordant" = "#2a9d8f", "Discordant" = "#d62828")) +
    ggplot2::scale_y_continuous(labels = scales::comma, expand = ggplot2::expansion(mult = c(0, 0.08))) +
    ggplot2::labs(
      title = "Strict Concordance by bardmux Edit Distance",
      subtitle = "Discordance concentrated at edit distance = 2",
      x = "Edit distance",
      y = "Read IDs",
      fill = NULL
    ) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(
      legend.position = "top",
      plot.title = ggplot2::element_text(face = "bold")
    )
} else {
  p4 <- make_empty_plot(ggplot2,
                        "Strict Concordance by bardmux Edit Distance",
                        "No edit-distance section found",
                        "No strict_edit_distance block found in summary")
}

out_png <- paste0(cfg$out_prefix, ".overview.png")

# Arrange 2x2 without extra dependencies beyond ggplot2 + base grid.
png(out_png, width = 2300, height = 1700, res = 160)
grid::grid.newpage()
layout <- grid::grid.layout(nrow = 3, ncol = 2,
                            heights = grid::unit(c(0.08, 0.46, 0.46), "npc"),
                            widths = grid::unit(c(0.5, 0.5), "npc"))
grid::pushViewport(grid::viewport(layout = layout))

grid::grid.text(cfg$title,
                vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1:2),
                gp = grid::gpar(fontsize = 18, fontface = "bold"))

print(p1, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 1))
print(p2, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 2))
print(p3, vp = grid::viewport(layout.pos.row = 3, layout.pos.col = 1))
print(p4, vp = grid::viewport(layout.pos.row = 3, layout.pos.col = 2))

dev.off()

cat(sprintf("Wrote: %s\n", out_png))
