#!/usr/bin/env Rscript

options(scipen = 999)

args <- commandArgs(trailingOnly = TRUE)

usage <- function() {
  cat(
    "Usage: plot_bardmux_pair_compare.R --summary SUMMARY.tsv --out-prefix PREFIX [--title TITLE]\n",
    "\n",
    "Required:\n",
    "  --summary FILE     Summary from compare_cb_assignments in bardmux_vs_bardmux mode\n",
    "\n",
    "Optional:\n",
    "  --out-prefix STR   Output prefix (default: bardmux_pair_compare)\n",
    "  --title STR        Figure title\n",
    sep = ""
  )
}

parse_args <- function(a) {
  cfg <- list(summary = NULL, out_prefix = "bardmux_pair_compare", title = "bardmux vs bardmux")
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
if (is.null(metrics_block)) stop("Could not find metric block in summary")
metrics <- read_tsv_block(metrics_block)

mode <- metrics$value[metrics$metric == "comparison_mode"]
if (length(mode) == 0 || mode[[1]] != "bardmux_vs_bardmux") {
  stop("Summary is not bardmux_vs_bardmux mode")
}

trans_block <- find_block(blocks, "status_a\tstatus_b\tcount\tpct_of_shared_status_ids")
trans_df <- NULL
if (!is.null(trans_block) && length(trans_block) > 1) {
  trans_df <- read_tsv_block(trans_block)
}

assigned_intersection <- get_metric(metrics, "assigned_intersection_ids", 0)
a_only <- get_metric(metrics, "bardmux_a_only_assigned_ids", 0)
b_only <- get_metric(metrics, "bardmux_b_only_assigned_ids", 0)
strict_comp <- get_metric(metrics, "assigned_strict_comparable_ids", 0)
cb_match <- get_metric(metrics, "assigned_cb_match_ids", 0)
cb_mismatch <- get_metric(metrics, "assigned_cb_mismatch_ids", 0)
cb_match_pct <- get_metric(metrics, "assigned_cb_match_pct", 0)
jaccard <- get_metric(metrics, "assigned_id_jaccard_pct", 0)
a_vs_b_cov <- get_metric(metrics, "bardmux_a_vs_b_assigned_coverage_pct", 0)
b_vs_a_cov <- get_metric(metrics, "bardmux_b_vs_a_assigned_coverage_pct", 0)
shared_status <- get_metric(metrics, "shared_status_ids", 0)
anchor_both <- get_metric(metrics, "anchor_found_both_shared_ids", 0)
anchor_a_only <- get_metric(metrics, "anchor_found_a_only_shared_ids", 0)
anchor_b_only <- get_metric(metrics, "anchor_found_b_only_shared_ids", 0)
anchor_neither <- get_metric(metrics, "anchor_found_neither_shared_ids", 0)
status_exact_pct <- get_metric(metrics, "status_exact_match_pct", 0)

# Panel 1: assigned overlap
p1_df <- data.frame(
  category = factor(c("Intersection", "A only", "B only"),
                    levels = c("Intersection", "A only", "B only")),
  count = c(assigned_intersection, a_only, b_only)
)
p1_df$label <- label_count_pct(p1_df$count)

p1 <- ggplot2::ggplot(p1_df, ggplot2::aes(x = category, y = count, fill = category)) +
  ggplot2::geom_col(width = 0.75, color = "#1f1f1f") +
  ggplot2::geom_text(ggplot2::aes(label = label), vjust = -0.25, size = 3.8, lineheight = 0.95) +
  ggplot2::scale_fill_manual(values = c("Intersection" = "#2a9d8f", "A only" = "#457b9d", "B only" = "#f4a261")) +
  ggplot2::scale_y_continuous(labels = scales::comma, expand = ggplot2::expansion(mult = c(0, 0.16))) +
  ggplot2::labs(
    title = "Assigned Read-ID Overlap",
    subtitle = sprintf("Jaccard %.2f%% | A_vs_B %.2f%% | B_vs_A %.2f%%", jaccard, a_vs_b_cov, b_vs_a_cov),
    x = NULL,
    y = "Distinct read IDs"
  ) +
  ggplot2::theme_bw(base_size = 13) +
  ggplot2::theme(legend.position = "none", plot.title = ggplot2::element_text(face = "bold"))

# Panel 2: CB agreement
p2_df <- data.frame(
  category = factor(c("CB match", "CB mismatch"), levels = c("CB match", "CB mismatch")),
  count = c(cb_match, cb_mismatch)
)
p2_df$label <- label_count_pct(p2_df$count, strict_comp)

p2 <- ggplot2::ggplot(p2_df, ggplot2::aes(x = category, y = count, fill = category)) +
  ggplot2::geom_col(width = 0.72, color = "#1f1f1f") +
  ggplot2::geom_text(ggplot2::aes(label = label), vjust = -0.25, size = 3.8, lineheight = 0.95) +
  ggplot2::scale_fill_manual(values = c("CB match" = "#2a9d8f", "CB mismatch" = "#d62828")) +
  ggplot2::scale_y_continuous(labels = scales::comma, expand = ggplot2::expansion(mult = c(0, 0.16))) +
  ggplot2::labs(
    title = sprintf("Assigned CB Concordance (%.2f%%)", cb_match_pct),
    subtitle = sprintf("Strict-comparable assigned IDs: %s", format(strict_comp, big.mark = ",", scientific = FALSE)),
    x = NULL,
    y = "Read IDs"
  ) +
  ggplot2::theme_bw(base_size = 13) +
  ggplot2::theme(legend.position = "none", plot.title = ggplot2::element_text(face = "bold"))

# Panel 3: anchor outcomes on shared status IDs
p3_df <- data.frame(
  category = factor(c("Both anchor", "A anchor only", "B anchor only", "Neither anchor"),
                    levels = c("Both anchor", "A anchor only", "B anchor only", "Neither anchor")),
  count = c(anchor_both, anchor_a_only, anchor_b_only, anchor_neither)
)
p3_df$label <- label_count_pct(p3_df$count, shared_status)

p3 <- ggplot2::ggplot(p3_df, ggplot2::aes(x = category, y = count, fill = category)) +
  ggplot2::geom_col(width = 0.8, color = "#1f1f1f") +
  ggplot2::geom_text(ggplot2::aes(label = label), vjust = -0.25, size = 3.7, lineheight = 0.95) +
  ggplot2::scale_fill_manual(values = c("Both anchor" = "#2a9d8f", "A anchor only" = "#457b9d",
                                        "B anchor only" = "#f4a261", "Neither anchor" = "#6c757d")) +
  ggplot2::scale_y_continuous(labels = scales::comma, expand = ggplot2::expansion(mult = c(0, 0.16))) +
  ggplot2::labs(
    title = "Anchor Presence on Shared Read IDs",
    subtitle = sprintf("Status exact-match: %.2f%% on %s shared IDs", status_exact_pct,
                       format(shared_status, big.mark = ",", scientific = FALSE)),
    x = NULL,
    y = "Shared read IDs"
  ) +
  ggplot2::theme_bw(base_size = 13) +
  ggplot2::theme(
    legend.position = "none",
    axis.text.x = ggplot2::element_text(angle = 18, hjust = 1),
    plot.title = ggplot2::element_text(face = "bold")
  )

# Panel 4: top status transitions
if (!is.null(trans_df) && nrow(trans_df) > 0) {
  trans_df$transition <- paste0(trans_df$status_a, " -> ", trans_df$status_b)
  trans_df <- trans_df[order(trans_df$count, decreasing = TRUE), , drop = FALSE]
  top_n <- min(12, nrow(trans_df))
  trans_top <- trans_df[seq_len(top_n), , drop = FALSE]
  trans_top$transition <- factor(trans_top$transition, levels = rev(trans_top$transition))
  trans_top$label <- sprintf("%s (%.1f%%)", format(trans_top$count, big.mark = ",", scientific = FALSE),
                             trans_top$pct_of_shared_status_ids)

  p4 <- ggplot2::ggplot(trans_top, ggplot2::aes(x = transition, y = count)) +
    ggplot2::geom_col(width = 0.75, fill = "#8d99ae", color = "#1f1f1f") +
    ggplot2::geom_text(ggplot2::aes(label = label), hjust = -0.02, size = 3.5) +
    ggplot2::coord_flip() +
    ggplot2::scale_y_continuous(labels = scales::comma, expand = ggplot2::expansion(mult = c(0, 0.22))) +
    ggplot2::labs(
      title = "Top Status Transitions (A -> B)",
      subtitle = sprintf("Showing top %d transitions by count", top_n),
      x = NULL,
      y = "Read IDs"
    ) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
} else {
  p4 <- make_empty_plot(ggplot2,
                        "Top Status Transitions (A -> B)",
                        "No transition table found",
                        "No status transition block in summary")
}

out_png <- paste0(cfg$out_prefix, ".overview.png")
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
