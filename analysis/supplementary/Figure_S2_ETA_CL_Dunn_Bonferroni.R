# Figure S2: patient-level ETA(CL) with Dunn-Bonferroni comparisons.
# This script creates new outputs; it does not replace the original figure.
# ETA estimates are exploratory diagnostics, not independent validation.

.script_file <- local({
  ofiles <- vapply(sys.frames(), function(frame) {
    value <- frame$ofile
    if (is.null(value) || length(value) != 1L) "" else as.character(value)
  }, character(1))
  ofiles <- ofiles[nzchar(ofiles)]
  if (length(ofiles)) {
    normalizePath(tail(ofiles, 1L), winslash = "/", mustWork = FALSE)
  } else {
    file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
    if (!length(file_arg)) stop("Run with Rscript or source() so the script location is known.")
    normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = FALSE)
  }
})
source(file.path(dirname(.script_file), "../../scripts/project_paths.R"))
.project_root <- find_project_root(.script_file)

required_packages <- c("ggplot2", "dunn.test", "ragg")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Please install these packages first: ", paste(missing_packages, collapse = ", "))
}
library(ggplot2)

configured_input <- Sys.getenv("ATX_FIGURE_S2_INPUT_FILE", unset = "")
input_file <- if (!nzchar(configured_input)) {
  file.path(.project_root, "data/analysis_inputs/supplementary/7.fit")
} else if (grepl("^(?:/|[A-Za-z]:[/\\\\]|\\\\\\\\)", configured_input, perl = TRUE)) {
  path.expand(configured_input)
} else {
  file.path(.project_root, configured_input)
}
output_dir <- project_output_dir(
  .project_root, "ATX_FIGURE_S2_OUTPUT_DIR", "results/supplementary"
)
output_stem <- file.path(output_dir, "Figure_S2_ETA_CL_Dunn_Bonferroni")

fit <- read.table(input_file, header = TRUE, skip = 1, fill = TRUE)
stopifnot(all(c("ID", "ETA1", "rs1065852") %in% names(fit)))

# Check that one ID has one ETA and one genotype before removing repeat rows.
for (variable in c("ETA1", "rs1065852")) {
  counts <- tapply(fit[[variable]], fit$ID, function(x) length(unique(x)))
  if (any(counts != 1L)) {
    stop("More than one ", variable, " value for at least one patient.")
  }
}
plot_data <- fit[!duplicated(fit$ID), c("ID", "ETA1", "rs1065852")]
plot_data$genotype <- factor(
  plot_data$rs1065852, levels = c(1, 2, 3), labels = c("CC", "CT", "TT")
)
stopifnot(!anyNA(plot_data$genotype), all(is.finite(plot_data$ETA1)))

kw <- kruskal.test(ETA1 ~ genotype, data = plot_data)
dunn <- dunn.test::dunn.test(
  plot_data$ETA1, plot_data$genotype,
  method = "bonferroni", altp = TRUE, kw = FALSE, table = FALSE, list = TRUE
)

# With altp = TRUE, the two-sided adjusted p values are in altP.adjusted.
pairwise <- data.frame(
  comparison = gsub(" ", "", dunn$comparisons, fixed = TRUE),
  z = dunn$Z,
  p_raw = dunn$altP,
  p_adjusted = dunn$altP.adjusted,
  stringsAsFactors = FALSE
)
stopifnot(nrow(pairwise) == 3L, all(is.finite(pairwise$p_adjusted)))
print(kw)
print(pairwise, digits = 10)

brackets <- data.frame(
  comparison = c("CC-CT", "CT-TT", "CC-TT"),
  xmin = c(1, 2, 1), xmax = c(2, 3, 3),
  y = c(0.80, 0.99, 1.18)
)
comparison_index <- match(brackets$comparison, pairwise$comparison)
stopifnot(!anyNA(comparison_index))
brackets$p_adjusted <- pairwise$p_adjusted[comparison_index]
brackets$label <- vapply(brackets$p_adjusted, function(p) {
  if (p < 0.0001) {
    "italic(p)[adj] < '0.0001'"
  } else {
    sprintf("italic(p)[adj] == '%0.3f'", p)
  }
}, character(1))
brackets$xmid <- (brackets$xmin + brackets$xmax) / 2
stopifnot(max(plot_data$ETA1) < min(brackets$y) - 0.06)

group_sizes <- table(plot_data$genotype)
group_labels <- setNames(
  paste0(names(group_sizes), "\n(n = ", as.integer(group_sizes), ")"),
  names(group_sizes)
)
kw_label <- if (kw$p.value < 0.001) {
  "CYP2D6*10 (c.100C>T)  |  Kruskal-Wallis, p < 0.001"
} else {
  paste0("CYP2D6*10 (c.100C>T)  |  Kruskal-Wallis, p = ", formatC(kw$p.value, digits = 3, format = "f"))
}

p <- ggplot(plot_data, aes(x = genotype, y = ETA1, fill = genotype)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.7, linewidth = 0.5) +
  geom_point(
    position = position_jitter(width = 0.12, height = 0, seed = 20261013),
    size = 1.8, alpha = 0.6, color = "black"
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.4) +
  geom_segment(
    data = brackets, aes(x = xmin, xend = xmax, y = y, yend = y),
    inherit.aes = FALSE, linewidth = 0.42, color = "black"
  ) +
  geom_segment(
    data = brackets, aes(x = xmin, xend = xmin, y = y - 0.035, yend = y),
    inherit.aes = FALSE, linewidth = 0.42, color = "black"
  ) +
  geom_segment(
    data = brackets, aes(x = xmax, xend = xmax, y = y - 0.035, yend = y),
    inherit.aes = FALSE, linewidth = 0.42, color = "black"
  ) +
  geom_text(
    data = brackets, aes(x = xmid, y = y + 0.045, label = label),
    inherit.aes = FALSE, parse = TRUE, size = 3.8, family = "sans"
  ) +
  scale_x_discrete(labels = group_labels) +
  scale_y_continuous(breaks = c(-1, -0.5, 0, 0.5), expand = expansion(mult = c(0, 0))) +
  coord_cartesian(ylim = c(-1.15, 1.32)) +
  scale_fill_manual(values = c("CC" = "#4E79A7", "CT" = "#F28E2B", "TT" = "#E15759")) +
  labs(
    title = "Distribution of ETA(CL) by rs1065852 genotype",
    subtitle = kw_label,
    x = "Genotype", y = "ETA(CL)",
    caption = "Pairwise comparisons: Dunn's test with Bonferroni correction."
  ) +
  theme_bw(base_size = 13, base_family = "sans") +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "plain", size = 13.5, margin = margin(b = 6)),
    plot.subtitle = element_text(size = 11.5, margin = margin(b = 11)),
    plot.caption = element_text(size = 10.5, color = "black", hjust = 0, margin = margin(t = 12)),
    plot.caption.position = "plot",
    axis.title = element_text(size = 13),
    axis.title.x = element_text(margin = margin(t = 8)),
    axis.title.y = element_text(margin = margin(r = 8)),
    axis.text = element_text(size = 11.5, color = "black"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "#EBEBEB", linewidth = 0.35),
    panel.border = element_rect(color = "#333333", linewidth = 0.55),
    plot.margin = margin(12, 14, 12, 12)
  )

ggsave(paste0(output_stem, ".png"), p, device = ragg::agg_png,
       width = 6.8, height = 5.8, units = "in", dpi = 300, bg = "white")
ggsave(paste0(output_stem, ".tiff"), p, device = ragg::agg_tiff,
       width = 6.8, height = 5.8, units = "in", dpi = 600, bg = "white", compression = "lzw")
ggsave(paste0(output_stem, ".pdf"), p, device = grDevices::cairo_pdf,
       width = 6.8, height = 5.8, units = "in", bg = "white")
cat("Saved new figure files to: ", output_stem, "\n", sep = "")
