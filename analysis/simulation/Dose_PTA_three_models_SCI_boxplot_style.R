## Combined dose-PTA curves for three PopPK models
## Revised layout: legend to the right of Model I and styling aligned with
## the combined simulated-C2h boxplot figure.

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

required_packages <- c("dplyr", "ggplot2", "patchwork", "cowplot")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop(
    "Please install: ",
    paste(missing_packages, collapse = ", ")
  )
}

library(dplyr)
library(ggplot2)
library(patchwork)
library(cowplot)

## ---------------------------------------------------------------------------
## 1. Paths
## ---------------------------------------------------------------------------

simulation_root <- project_input_dir(
  .project_root, "ATX_SIMULATION_ROOT", "data/analysis_inputs/simulation"
)
output_folder <- project_output_dir(
  .project_root, "ATX_DOSE_PTA_OUTPUT_DIR", "results/simulation/dose_pta"
)

file_model1 <- file.path(
  simulation_root, "BW", "attainment_ModelI_BW_3strata.csv"
)
file_model2 <- file.path(
  simulation_root, "RS1", "attainment_ModelII_3strata.csv"
)
file_model3 <- file.path(
  simulation_root, "AS", "attainment_ModelAS_3strata.csv"
)

## ---------------------------------------------------------------------------
## 2. Settings
## ---------------------------------------------------------------------------

DOSE_BREAKS <- c(0.5, 0.8, 1.0, 1.1, 1.2, 1.3, 1.4)
WT_LEVELS <- c("17-<30", "30-<40", "40-60")
WT_LABELS <- c(
  "17-<30" = "17\u2013<30 kg",
  "30-<40" = "30\u2013<40 kg",
  "40-60" = "40\u201360 kg"
)
WT_COLORS <- c(
  "17-<30" = "#0072B2",
  "30-<40" = "#E69F00",
  "40-60" = "#009E73"
)

PTA_MAIN <- 75
PTA_SENSITIVITY <- 80

FIGURE_WIDTH_IN <- 180 / 25.4
FIGURE_HEIGHT_IN <- 225 / 25.4

## ---------------------------------------------------------------------------
## 3. Read and validate data
## ---------------------------------------------------------------------------

read_pta <- function(file, required_columns, model_name) {
  if (!file.exists(file)) {
    stop(model_name, " input file does not exist: ", file)
  }

  dat <- read.csv(
    file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  missing_columns <- setdiff(required_columns, names(dat))
  if (length(missing_columns) > 0L) {
    stop(
      model_name, " is missing column(s): ",
      paste(missing_columns, collapse = ", ")
    )
  }

  dat <- dat %>%
    mutate(
      WT_BAND = factor(WT_BAND, levels = WT_LEVELS),
      DOSEKG = as.numeric(DOSEKG),
      PTA = as.numeric(PTA)
    ) %>%
    filter(!is.na(WT_BAND), !is.na(DOSEKG), !is.na(PTA))

  if (any(dat$PTA < 0 | dat$PTA > 100)) {
    stop(model_name, " contains PTA outside 0%-100%.")
  }

  if (!isTRUE(all.equal(sort(unique(dat$DOSEKG)), DOSE_BREAKS))) {
    stop(model_name, " does not contain the expected seven-dose grid.")
  }

  dat
}

model1 <- read_pta(
  file_model1,
  c("WT_BAND", "DOSEKG", "PTA"),
  "Model I"
)

model2 <- read_pta(
  file_model2,
  c("WT_BAND", "Group", "DOSEKG", "PTA"),
  "Model II"
) %>%
  mutate(Group = factor(Group, levels = c("CC+CT", "TT")))

model3 <- read_pta(
  file_model3,
  c("WT_BAND", "AS_group", "DOSEKG", "PTA"),
  "Model III"
) %>%
  mutate(
    AS_group = factor(
      AS_group,
      levels = c("AS = 0.5", "AS = 1.25", "AS = 1.5", "AS = 2")
    )
  )

## ---------------------------------------------------------------------------
## 4. Panel function
## ---------------------------------------------------------------------------

make_panel <- function(dat, panel_title) {
  dat <- dat %>%
    mutate(Panel_title = panel_title)

  ggplot(
    dat,
    aes(
      x = DOSEKG,
      y = PTA,
      color = WT_BAND,
      group = WT_BAND
    )
  ) +
    geom_hline(
      yintercept = PTA_MAIN,
      color = "#444444",
      linetype = "dashed",
      linewidth = 0.45
    ) +
    geom_hline(
      yintercept = PTA_SENSITIVITY,
      color = "#888888",
      linetype = "dotted",
      linewidth = 0.42
    ) +
    geom_line(linewidth = 0.68) +
    geom_point(size = 1.75) +
    facet_wrap(~Panel_title) +
    scale_x_continuous(
      breaks = DOSE_BREAKS,
      limits = range(DOSE_BREAKS),
      expand = expansion(mult = c(0.025, 0.025))
    ) +
    scale_y_continuous(
      breaks = seq(0, 100, 20),
      limits = c(0, 100),
      expand = expansion(mult = c(0, 0.015))
    ) +
    scale_color_manual(
      values = WT_COLORS,
      limits = WT_LEVELS,
      labels = WT_LABELS,
      drop = FALSE,
      name = NULL
    ) +
    guides(
      color = guide_legend(
        ncol = 1,
        byrow = TRUE,
        override.aes = list(linewidth = 0.85, size = 2.2)
      )
    ) +
    theme_bw(base_size = 9.5, base_family = "sans") +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(
        color = "#E5E5E5",
        linewidth = 0.32
      ),
      panel.border = element_rect(
        color = "#444444",
        linewidth = 0.45
      ),
      strip.background = element_rect(
        fill = "#F2F2F2",
        color = "#555555",
        linewidth = 0.4
      ),
      strip.text = element_text(size = 9.5, face = "bold"),
      axis.title = element_blank(),
      axis.text = element_text(size = 8.3, color = "black"),
      axis.ticks = element_line(linewidth = 0.35, color = "black"),
      legend.position = "right",
      legend.text = element_text(size = 8.7),
      legend.key.width = grid::unit(10, "mm"),
      legend.spacing.y = grid::unit(1.5, "mm"),
      plot.tag = element_text(size = 10, face = "bold"),
      plot.tag.position = c(0.01, 0.99),
      plot.margin = margin(4, 5, 4, 5)
    )
}

## ---------------------------------------------------------------------------
## 5. Construct panels
## ---------------------------------------------------------------------------

pA <- make_panel(model1, "Model I")
pB <- make_panel(filter(model2, Group == "CC+CT"), "Model II: CC+CT")
pC <- make_panel(filter(model2, Group == "TT"), "Model II: TT")
pD <- make_panel(
  filter(model3, AS_group == "AS = 0.5"),
  "Model III: AS = 0.5"
)
pE <- make_panel(
  filter(model3, AS_group == "AS = 1.25"),
  "Model III: AS = 1.25"
)
pF <- make_panel(
  filter(model3, AS_group == "AS = 1.5"),
  "Model III: AS = 1.5"
)
pG <- make_panel(
  filter(model3, AS_group == "AS = 2"),
  "Model III: AS = 2"
)

## Model I occupies the upper-left half; the collected legend occupies the
## upper-right half. The remaining panels preserve the previous two-column grid.
figure_design <- c(
  area(t = 1, l = 1, b = 1, r = 2), # Model I
  area(t = 2, l = 1, b = 2, r = 2), # Model II: CC+CT
  area(t = 2, l = 3, b = 2, r = 4), # Model II: TT
  area(t = 3, l = 1, b = 3, r = 2), # Model III: AS 0.5
  area(t = 3, l = 3, b = 3, r = 4), # Model III: AS 1.25
  area(t = 4, l = 1, b = 4, r = 2), # Model III: AS 1.5
  area(t = 4, l = 3, b = 4, r = 4), # Model III: AS 2
  area(t = 1, l = 3, b = 1, r = 4)  # collected legend
)

panel_grid <- (
  pA + pB + pC + pD + pE + pF + pG + guide_area()
) +
  plot_layout(
    design = figure_design,
    heights = c(1, 1, 1, 1),
    guides = "collect"
  ) +
  plot_annotation(tag_levels = "A")

combined_plot <- ggdraw() +
  draw_plot(
    panel_grid,
    x = 0.055,
    y = 0.045,
    width = 0.935,
    height = 0.95
  ) +
  draw_label(
    "Dose (mg/kg/day)",
    x = 0.53,
    y = 0.012,
    size = 9,
    fontfamily = "sans"
  ) +
  draw_label(
    "Probability of target attainment (%)",
    x = 0.014,
    y = 0.52,
    angle = 90,
    size = 9,
    fontfamily = "sans"
  )

## ---------------------------------------------------------------------------
## 6. Export
## ---------------------------------------------------------------------------

output_stem <- file.path(
  output_folder,
  "Combined_Dose_PTA_three_models_boxplot_style"
)

ggsave(
  paste0(output_stem, ".png"),
  combined_plot,
  width = FIGURE_WIDTH_IN,
  height = FIGURE_HEIGHT_IN,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  paste0(output_stem, ".tiff"),
  combined_plot,
  width = FIGURE_WIDTH_IN,
  height = FIGURE_HEIGHT_IN,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

ggsave(
  paste0(output_stem, ".pdf"),
  combined_plot,
  width = FIGURE_WIDTH_IN,
  height = FIGURE_HEIGHT_IN,
  units = "in",
  device = cairo_pdf,
  bg = "white"
)

message("Completed. Files saved to: ", output_folder)
