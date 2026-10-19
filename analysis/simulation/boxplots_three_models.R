## Simulated steady-state C2h distributions by dose
## Model I (body weight), Model II (rs1065852), Model III (activity score)
##
## Box: 25th-75th percentiles
## Center line: median
## Whiskers: 5th-95th percentiles

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

required_packages <- c("dplyr", "ggplot2", "readr")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop(
    "Please install the following packages: ",
    paste(missing_packages, collapse = ", ")
  )
}

library(dplyr)
library(ggplot2)
library(readr)

## ---------------------------------------------------------------------------
## 1. Paths and analysis settings
## ---------------------------------------------------------------------------

simulation_root <- project_input_dir(
  .project_root, "ATX_SIMULATION_ROOT", "data/analysis_inputs/simulation"
)
output_dir <- project_output_dir(
  .project_root, "ATX_BOXPLOT_OUTPUT_DIR", "results/simulation/C2h_boxplots"
)

files <- tibble::tribble(
  ~Model,      ~File_type, ~Path,
  "Model I",   "Main",     file.path(simulation_root, "BW", "sim_BW.tab"),
  "Model I",   "Dose 1.3", file.path(simulation_root, "BW", "1.3", "sim_BW1.3.tab"),
  "Model II",  "Main",     file.path(simulation_root, "RS1", "sim_RS1.tab"),
  "Model II",  "Dose 1.3", file.path(simulation_root, "RS1", "1.3", "sim_RS1.3.tab"),
  "Model III", "Main",     file.path(simulation_root, "AS", "sim_AS.tab"),
  "Model III", "Dose 1.3", file.path(simulation_root, "AS", "1.3", "sim_AS1.3.tab")
)

missing_files <- files$Path[!file.exists(files$Path)]
if (length(missing_files) > 0) {
  stop("Missing simulation file(s):\n", paste(missing_files, collapse = "\n"))
}

TARGET <- 268
C2H_TAD <- 2
EXPECTED_SUBPROBLEMS <- 500L
CHUNK_SIZE <- 200000L

dose_levels <- c("0.5", "0.8", "1", "1.1", "1.2", "1.3", "1.4")
weight_levels <- c("17-<30", "30-<40", "40-60")
model_levels <- c("Model I", "Model II", "Model III")
subgroup_levels <- c(
  "Weight-based", "CC+CT", "TT",
  "AS = 0.5", "AS = 1.25", "AS = 1.5", "AS = 2"
)

audit_records <- list()

## ---------------------------------------------------------------------------
## 2. Stream one NONMEM table and retain only steady-state C2h records
## ---------------------------------------------------------------------------

read_c2h_file <- function(path, model, file_type) {
  if (model == "Model II") {
    tab_cols <- c(
      "ID", "TIME", "TAD", "AMT", "DV", "IPRED",
      "CL", "V", "WT", "rs_grp", "DOSEKG", "MDV"
    )
    grouping_column <- "rs_grp"
  } else {
    tab_cols <- c(
      "ID", "TIME", "TAD", "AMT", "DV", "IPRED",
      "CL", "V", "WT", "AS", "DOSEKG", "MDV"
    )
    grouping_column <- "AS"
  }

  con <- file(path, open = "r")
  on.exit(close(con), add = TRUE)

  retained_chunks <- list()
  retained_id <- 0L
  total_lines <- 0L
  table_count <- 0L

  repeat {
    txt <- readLines(con, n = CHUNK_SIZE, warn = FALSE)
    if (length(txt) == 0L) break

    total_lines <- total_lines + length(txt)
    table_count <- table_count + sum(grepl("^\\s*TABLE\\s+NO\\.", txt))

    keep <- nzchar(trimws(txt)) &
      !grepl("^\\s*TABLE\\s+NO\\.", txt) &
      !grepl("^\\s*ID\\s+TIME\\s+TAD\\s+AMT", txt)
    txt <- txt[keep]

    if (length(txt) > 0L) {
      block <- read.table(
        text = txt,
        header = FALSE,
        col.names = tab_cols,
        colClasses = "numeric",
        fill = FALSE,
        comment.char = "",
        stringsAsFactors = FALSE
      )

      block <- block[
        block$MDV == 0 & abs(block$TAD - C2H_TAD) < 1e-6,
        c("DV", "IPRED", "WT", grouping_column, "DOSEKG")
      ]

      if (nrow(block) > 0L) {
        retained_id <- retained_id + 1L
        retained_chunks[[retained_id]] <- block
      }
    }
  }

  if (length(retained_chunks) == 0L) {
    stop("No TAD = 2 and MDV = 0 records found in ", path)
  }

  c2h <- bind_rows(retained_chunks)
  rm(retained_chunks)

  if (model == "Model I") {
    c2h$Subgroup <- "Weight-based"
  } else if (model == "Model II") {
    c2h$Subgroup <- dplyr::case_when(
      c2h$rs_grp == 0 ~ "CC+CT",
      c2h$rs_grp == 1 ~ "TT",
      TRUE ~ NA_character_
    )
  } else {
    c2h$AS <- round(c2h$AS, 2)
    c2h$Subgroup <- dplyr::case_when(
      c2h$AS == 0.50 ~ "AS = 0.5",
      c2h$AS == 1.25 ~ "AS = 1.25",
      c2h$AS == 1.50 ~ "AS = 1.5",
      c2h$AS == 2.00 ~ "AS = 2",
      TRUE ~ NA_character_
    )
  }

  c2h <- c2h %>%
    mutate(
      Model = model,
      File_type = file_type,
      WT_GROUP = case_when(
        WT >= 17 & WT < 30 ~ "17-<30",
        WT >= 30 & WT < 40 ~ "30-<40",
        WT >= 40 & WT <= 60 + 1e-8 ~ "40-60",
        TRUE ~ NA_character_
      ),
      Dose = sub("\\.0$", "", sprintf("%.1f", round(DOSEKG, 1)))
    )

  if (anyNA(c2h$Subgroup)) {
    stop("Unrecognized subgroup value(s) in ", path)
  }
  if (anyNA(c2h$WT_GROUP)) {
    stop("Body weight outside 17-60 kg was found in ", path)
  }
  if (any(!is.finite(c2h$DV))) {
    stop("Non-finite simulated DV value(s) found in ", path)
  }

  observed_doses <- sort(unique(round(c2h$DOSEKG, 1)))
  expected_doses <- if (file_type == "Dose 1.3") {
    1.3
  } else {
    c(0.5, 0.8, 1.0, 1.1, 1.2, 1.4)
  }

  if (!identical(observed_doses, expected_doses)) {
    stop(
      "Unexpected dose grid in ", path, ". Found: ",
      paste(observed_doses, collapse = ", ")
    )
  }

  audit_records[[length(audit_records) + 1L]] <<- tibble(
    Model = model,
    File_type = file_type,
    File = path,
    Table_count = table_count,
    Expected_tables = EXPECTED_SUBPROBLEMS,
    C2h_records = nrow(c2h),
    Minimum_C2h = min(c2h$DV),
    Maximum_C2h = max(c2h$DV),
    Doses = paste(observed_doses, collapse = ", ")
  )

  if (table_count != EXPECTED_SUBPROBLEMS) {
    warning(
      basename(path), " contains ", table_count,
      " TABLE blocks; expected ", EXPECTED_SUBPROBLEMS, "."
    )
  }

  summary <- c2h %>%
    group_by(Model, WT_GROUP, Subgroup, Dose) %>%
    summarise(
      n = n(),
      q05 = quantile(DV, 0.05, na.rm = TRUE, names = FALSE),
      q25 = quantile(DV, 0.25, na.rm = TRUE, names = FALSE),
      median = median(DV, na.rm = TRUE),
      q75 = quantile(DV, 0.75, na.rm = TRUE, names = FALSE),
      q95 = quantile(DV, 0.95, na.rm = TRUE, names = FALSE),
      .groups = "drop"
    )

  rm(c2h)
  gc(verbose = FALSE)
  summary
}

summary_file <- file.path(
  output_dir, "Simulated_C2h_boxplot_summary_three_models.csv"
)
audit_file <- file.path(
  output_dir, "Simulated_C2h_boxplot_file_audit.csv"
)
reuse_summary <- tolower(Sys.getenv(
  "ATX_REUSE_BOXPLOT_SUMMARY", unset = "false"
)) == "true"

if (reuse_summary && file.exists(summary_file) && file.exists(audit_file)) {
  message("Reusing existing summary and audit files.")
  box_summary <- read_csv(summary_file, show_col_types = FALSE) %>%
    mutate(
      ## Repair summaries produced by an earlier script version in which
      ## the character label for 1.0 mg/kg/day was read as NA.
      Dose = if_else(is.na(Dose), "1", as.character(Dose))
    )
  audit <- read_csv(audit_file, show_col_types = FALSE)
} else {
  ## Read and immediately summarise each large file to limit memory use.
  summary_list <- vector("list", nrow(files))

  for (i in seq_len(nrow(files))) {
    message(
      "[", i, "/", nrow(files), "] Reading ",
      files$Model[i], " — ", files$File_type[i]
    )
    summary_list[[i]] <- read_c2h_file(
      path = files$Path[i],
      model = files$Model[i],
      file_type = files$File_type[i]
    )
  }

  box_summary <- bind_rows(summary_list)
  audit <- bind_rows(audit_records)
}

box_summary <- box_summary %>%
  mutate(
    Model = factor(Model, levels = model_levels),
    WT_GROUP = factor(WT_GROUP, levels = weight_levels),
    Dose = factor(Dose, levels = dose_levels),
    Subgroup = factor(Subgroup, levels = subgroup_levels)
  ) %>%
  arrange(Model, WT_GROUP, Dose, Subgroup) %>%
  mutate(Box_id = row_number())

## Every model × weight stratum must contain all seven doses.
grid_check <- box_summary %>%
  distinct(Model, WT_GROUP, Dose) %>%
  count(Model, WT_GROUP, name = "Dose_count")

if (any(grid_check$Dose_count != length(dose_levels))) {
  stop("At least one model × weight stratum does not contain all seven doses.")
}

write_csv(
  box_summary,
  summary_file
)
write_csv(
  audit,
  audit_file
)

## ---------------------------------------------------------------------------
## 3. Combined manuscript figure
## ---------------------------------------------------------------------------

subgroup_colors <- c(
  "Weight-based" = "#4C78A8",
  "CC+CT" = "#56B4E9",
  "TT" = "#D55E00",
  "AS = 0.5" = "#CC79A7",
  "AS = 1.25" = "#009E73",
  "AS = 1.5" = "#E69F00",
  "AS = 2" = "#0072B2"
)

model_labels <- c(
  "Model I" = "Model I",
  "Model II" = "Model II",
  "Model III" = "Model III"
)

weight_labels <- c(
  "17-<30" = "17\u2013<30 kg",
  "30-<40" = "30\u2013<40 kg",
  "40-60" = "40\u201360 kg"
)

y_upper <- ceiling(
  max(c(800, box_summary$q95), na.rm = TRUE) / 100
) * 100

p_box <- ggplot() +
  geom_hline(
    yintercept = TARGET,
    linetype = "dashed",
    color = "#444444",
    linewidth = 0.45
  ) +
  geom_boxplot(
    data = box_summary,
    aes(
      x = Dose,
      ymin = q05,
      lower = q25,
      middle = median,
      upper = q75,
      ymax = q95,
      fill = Subgroup,
      group = Box_id
    ),
    stat = "identity",
    width = 0.78,
    position = position_dodge2(width = 0.82, preserve = "single", padding = 0.08),
    color = "#2F2F2F",
    linewidth = 0.34
  ) +
  facet_grid(
    rows = vars(Model),
    cols = vars(WT_GROUP),
    switch = "y",
    labeller = labeller(
      Model = model_labels,
      WT_GROUP = weight_labels
    )
  ) +
  scale_fill_manual(
    name = NULL,
    values = subgroup_colors,
    drop = TRUE
  ) +
  scale_y_continuous(
    name = "Simulated steady-state C2h (ng/mL)",
    breaks = seq(0, y_upper, by = 200),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_x_discrete(
    name = "Dose (mg/kg/day)",
    drop = FALSE
  ) +
  coord_cartesian(ylim = c(0, y_upper), clip = "on") +
  guides(
    fill = guide_legend(order = 1, nrow = 1, byrow = TRUE)
  ) +
  theme_bw(base_size = 9.5, base_family = "sans") +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "#E5E5E5", linewidth = 0.32),
    panel.border = element_rect(color = "#444444", linewidth = 0.45),
    strip.background = element_rect(fill = "#F2F2F2", color = "#555555", linewidth = 0.4),
    strip.text.x = element_text(size = 9.5, face = "bold"),
    strip.text.y.left = element_text(size = 9.5, face = "bold", angle = 90),
    strip.placement = "outside",
    axis.title = element_text(size = 10.5, color = "black"),
    axis.text = element_text(size = 8.3, color = "black"),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "vertical",
    legend.box.just = "center",
    legend.title = element_text(size = 9.3),
    legend.text = element_text(size = 8.7),
    legend.key.width = grid::unit(1.15, "lines"),
    panel.spacing = grid::unit(0.55, "lines"),
    plot.margin = margin(5.5, 8, 5.5, 5.5)
  )

## ---------------------------------------------------------------------------
## 4. Export
## ---------------------------------------------------------------------------

ggsave(
  file.path(output_dir, "Simulated_C2h_boxplots_three_models.png"),
  p_box,
  width = 220,
  height = 225,
  units = "mm",
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(output_dir, "Simulated_C2h_boxplots_three_models.tiff"),
  p_box,
  width = 220,
  height = 225,
  units = "mm",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

ggsave(
  file.path(output_dir, "Simulated_C2h_boxplots_three_models.pdf"),
  p_box,
  width = 220 / 25.4,
  height = 225 / 25.4,
  units = "in",
  device = cairo_pdf,
  bg = "white"
)

message("Completed. Files saved to: ", output_dir)
