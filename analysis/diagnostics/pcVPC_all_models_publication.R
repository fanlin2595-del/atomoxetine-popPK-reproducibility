# ============================================================================
# Figure 2: all atomoxetine pc-VPC panels in ONE reproducible R script
# Reviewer revision: readable text at final print size, high-contrast curves,
# an explicit shared key, and vector PDF plus 600-dpi PNG/TIFF exports.
# Style approved from Figure2_pcVPC_preview_v4.png (2026-10-10): single-line
# axis/strip titles, compact fonts, aligned key, and pale outer CI ribbons.
#
# Layout: a = Model I / Model II / Model III; b = CC+CT / TT / shared key.
# Five equal-size panels, left-aligned rows, grey regular-weight title strips,
# and a dashed outside border. No image cropping, rsvg, or image stitching.
#
# Run from the project root: Rscript --vanilla analysis/diagnostics/pcVPC_all_models_publication.R
# Packages: ggplot2, xpose4, gtable. ragg is optional (recommended for raster).
# Original scripts, source data, simulation summaries, and old images are NOT
# changed. Existing files bearing the NEW Figure2_pcVPC_revised name are updated.
# ============================================================================

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

main_pcVPC <- function() {
  required <- c("ggplot2", "xpose4", "gtable")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop("Missing R packages: ", paste(missing, collapse = ", "),
         ". Install these packages before running this script.")
  }

  # --------------------------- Editable settings -----------------------------
  root <- project_input_dir(.project_root, "ATX_PCVPC_INPUT_DIR", "data/analysis_inputs/pcVPC")
  output_dir <- project_output_dir(.project_root, "ATX_PCVPC_OUTPUT_DIR", "results/pcVPC")
  output_name <- "Figure2_pcVPC_revised"

  # Physical size matters: do not export a very wide figure and then shrink it.
  # Text sizes below are points AT the final 180-mm-wide figure size.
  figure_width_mm  <- 180
  figure_height_mm <- 138
  raster_dpi <- 600
  font_family <- "Helvetica"
  font <- list(axis_title_x = 7.5, axis_title_y = 6.5, tick = 7.5, strip = 8,
               row_label = 10, key = 8.5, key_title = 9)

  # Keep the original statistical settings. The supplied PsN files contain
  # 95% confidence intervals (NOT 90%). Do not relabel them as 90% intervals.
  # The 5th/95th percentile curves bound a central 90% prediction interval.
  # This differs from the 95% CI around EACH simulated percentile; for an
  # equal-tailed CI its bounds are the 2.5th/97.5th percentiles of the simulated
  # percentile estimates. Read those CI bounds from PsN, not concentration PI bounds.
  ci <- 95
  lo.pi <- 5
  up.pi <- 95
  x_limits <- c(0, 25)
  x_breaks <- c(0, 2, 6, 12, 18, 24)
  y_limits <- c(0, 1500)
  y_breaks <- c(0, 500, 1000, 1500)

  palette <- list(observed = "#B52D28", simulated = "#005A8D",
                  points = "#333333", ci_outer = "#E2EDF5",
                  ci_median = "#91B8D6", strip = "#D9D9D9",
                  text = "#111111", border = "#4A4A4A")
  line_median <- 0.60
  line_outer <- 0.55
  point_size <- 1.10
  point_stroke <- 0.22

  # Panels occupy identical physical cells. The bottom row uses columns 1/2.
  geometry <- list(left = 5.5, top = 6, col_gap = 3.0, right = 4.0,
                   row_gap = 8.0, bottom = 6.0, frame_inset = 1.2,
                   title_gap = 1.8)
  cell_w <- (figure_width_mm - geometry$left - geometry$right -
               2 * geometry$col_gap) / 3
  cell_h <- (figure_height_mm - geometry$top - geometry$bottom -
               geometry$row_gap) / 2
  col_left <- geometry$left + (0:2) * (cell_w + geometry$col_gap)
  row_top <- geometry$top + (0:1) * (cell_h + geometry$row_gap)

  input_folders <- c(
    BW = file.path(root, "BW", "pcVPC"),
    RS = file.path(root, "RS1065852", "pcVPC", "pcvpc-overall"),
    AS = file.path(root, "AS", "pcVPC"),
    STRAT = file.path(root, "RS1065852", "pcVPC", "pcvpc-strat")
  )

  # ---------------------- Read existing PsN summaries -------------------------
  read_pcVPC <- function(folder) {
    results_file <- file.path(folder, "vpc_results.csv")
    tab_file <- file.path(folder, "vpctab")
    if (!file.exists(tab_file)) {
      candidates <- list.files(folder, pattern = "^vpctab", full.names = TRUE)
      if (length(candidates) != 1L) {
        stop("Cannot identify one vpctab file in: ", folder)
      }
      tab_file <- candidates[[1]]
    }
    if (!file.exists(results_file)) stop("Missing: ", results_file)
    vpc <- xpose4::read.npc.vpc.results(results_file)
    obs <- xpose4::read.vpctab(tab_file)@Data

    # PsN's prediction-corrected vpctab can legitimately name its corrected
    # concentration column DV. Check run metadata, rather than silently using
    # uncorrected observations in a plot labelled pc-VPC.
    if (!isTRUE(vpc$pred.corr)) {
      stop("Prediction correction is not confirmed in: ", results_file)
    }
    dv <- if ("PREDCORR_DV" %in% names(obs)) "PREDCORR_DV" else "DV"
    if (!all(c("ID", "TAD", dv) %in% names(obs))) {
      stop("Required observation columns are missing in: ", tab_file)
    }
    obs <- obs[!is.na(obs$TAD) & obs$TAD >= x_limits[1] &
                 obs$TAD <= x_limits[2] & !is.na(obs[[dv]]), , drop = FALSE]
    obs$pc_concentration <- obs[[dv]]
    tables <- if (is.data.frame(vpc$result.tables)) {
      list(vpc$result.tables)
    } else {
      Filter(is.data.frame, vpc$result.tables)
    }
    if (!length(tables)) stop("No continuous VPC tables in: ", results_file)
    list(obs = obs, tables = tables, strata_names = vpc$strata.names,
         dv_column = dv, prediction_corrected = vpc$pred.corr,
         source_files = c(results_file, tab_file))
  }

  getcol <- function(df, name) {
    if (!name %in% names(df)) {
      stop("Required VPC column not found: ", name,
           ". Check the chosen CI and percentile levels.")
    }
    df[[name]]
  }

  make_statistics <- function(df) {
    # Preserve the original scripts' BIN MIDPOINTS and percentile values.
    # No re-binning, re-estimation, prediction correction, or simulation here.
    out <- data.frame(
      mid = (getcol(df, "upper") + getcol(df, "lower")) / 2,
      obs_med = getcol(df, "50.real"),
      obs_lo = getcol(df, paste0(lo.pi, ".real")),
      obs_up = getcol(df, paste0(up.pi, ".real")),
      sim_med = getcol(df, "50.sim"),
      sim_lo = getcol(df, paste0(lo.pi, ".sim")),
      sim_up = getcol(df, paste0(up.pi, ".sim")),
      ci_lo_from = getcol(df, paste0(ci, ".CI.for.", lo.pi, ".from")),
      ci_lo_to = getcol(df, paste0(ci, ".CI.for.", lo.pi, ".to")),
      ci_up_from = getcol(df, paste0(ci, ".CI.for.", up.pi, ".from")),
      ci_up_to = getcol(df, paste0(ci, ".CI.for.", up.pi, ".to")),
      ci_md_from = getcol(df, paste0(ci, ".CI.for.50.from")),
      ci_md_to = getcol(df, paste0(ci, ".CI.for.50.to"))
    )
    out <- out[!is.na(out$mid) & out$mid >= x_limits[1] &
                 out$mid <= x_limits[2], , drop = FALSE]
    out[order(out$mid), , drop = FALSE]
  }

  sources <- lapply(input_folders, read_pcVPC)
  for (name in c("BW", "RS", "AS")) {
    if (length(sources[[name]]$tables) != 1L) {
      stop("Expected one overall VPC table for model: ", name)
    }
  }
  strata <- sources$STRAT
  if (!all(c("rs_grp", "strata_no") %in% names(strata$obs))) {
    stop("The stratified vpctab must identify rs_grp and strata_no.")
  }
  get_stratum <- function(code) {
    pattern <- paste0("^rs_grp\\s*==\\s*", code, "\\s*$")
    index <- which(grepl(pattern, trimws(strata$strata_names)))
    if (length(index) != 1L) stop("Cannot uniquely identify rs_grp = ", code)
    obs <- strata$obs[strata$obs$rs_grp == code & !is.na(strata$obs$rs_grp), , drop = FALSE]
    if (!nrow(obs) || any(obs$strata_no != index)) {
      stop("Observation group and VPC stratum do not agree for rs_grp = ", code)
    }
    list(obs = obs, stat = make_statistics(strata$tables[[index]]),
         n = length(unique(obs$ID)))
  }
  ccct <- get_stratum(0)
  tt <- get_stratum(1)

  panels <- list(
    BW = list(obs = sources$BW$obs, stat = make_statistics(sources$BW$tables[[1]]),
              title = "MODEL I (weight based)"),
    RS = list(obs = sources$RS$obs, stat = make_statistics(sources$RS$tables[[1]]),
              title = "MODEL II (rs1065852 genotype)"),
    AS = list(obs = sources$AS$obs, stat = make_statistics(sources$AS$tables[[1]]),
              title = "MODEL III (activity score)"),
    CCCT = list(obs = ccct$obs, stat = ccct$stat,
                title = sprintf("rs1065852 CC+CT (n = %d)", ccct$n)),
    TT = list(obs = tt$obs, stat = tt$stat,
              title = sprintf("rs1065852 TT (n = %d)", tt$n))
  )

  # -------------------------- Native vector plots ----------------------------
  draw_panel <- function(panel) {
    s <- panel$stat
    o <- panel$obs
    s$panel_title <- panel$title
    o$panel_title <- panel$title
    p <- ggplot2::ggplot() +
      # Ribbons FIRST, to avoid washing out the points and lines.
      ggplot2::geom_ribbon(data = s, ggplot2::aes(x = mid, ymin = ci_lo_from, ymax = ci_lo_to),
                           fill = palette$ci_outer, colour = NA) +
      ggplot2::geom_ribbon(data = s, ggplot2::aes(x = mid, ymin = ci_up_from, ymax = ci_up_to),
                           fill = palette$ci_outer, colour = NA) +
      ggplot2::geom_ribbon(data = s, ggplot2::aes(x = mid, ymin = ci_md_from, ymax = ci_md_to),
                           fill = palette$ci_median, colour = NA) +
      ggplot2::geom_point(data = o, ggplot2::aes(x = TAD, y = pc_concentration),
                          shape = 1, size = point_size, stroke = point_stroke,
                          colour = palette$points) +
      ggplot2::geom_line(data = s, ggplot2::aes(x = mid, y = sim_med),
                         colour = palette$simulated, linewidth = line_median) +
      ggplot2::geom_line(data = s, ggplot2::aes(x = mid, y = sim_lo),
                         colour = palette$simulated, linewidth = line_outer, linetype = "dashed") +
      ggplot2::geom_line(data = s, ggplot2::aes(x = mid, y = sim_up),
                         colour = palette$simulated, linewidth = line_outer, linetype = "dashed") +
      ggplot2::geom_line(data = s, ggplot2::aes(x = mid, y = obs_med),
                         colour = palette$observed, linewidth = line_median) +
      ggplot2::geom_line(data = s, ggplot2::aes(x = mid, y = obs_lo),
                         colour = palette$observed, linewidth = line_outer, linetype = "dashed") +
      ggplot2::geom_line(data = s, ggplot2::aes(x = mid, y = obs_up),
                         colour = palette$observed, linewidth = line_outer, linetype = "dashed") +
      ggplot2::scale_x_continuous(breaks = x_breaks,
                                   expand = ggplot2::expansion(mult = 0.04)) +
      ggplot2::scale_y_continuous(breaks = y_breaks,
                                   expand = ggplot2::expansion(mult = 0.04)) +
      # Coordinate limits crop only the displayed window, not the source data.
      ggplot2::coord_cartesian(xlim = x_limits, ylim = y_limits) +
      ggplot2::facet_wrap(ggplot2::vars(panel_title)) +
      ggplot2::labs(x = "Time after dose (h)",
                     y = "Prediction-corrected concentration (ng/mL)") +
      ggplot2::theme_bw(base_size = font$tick, base_family = font_family) +
      ggplot2::theme(
        text = ggplot2::element_text(colour = palette$text, face = "plain"),
        axis.title = ggplot2::element_text(face = "plain"),
        axis.title.x = ggplot2::element_text(size = font$axis_title_x,
                                             margin = ggplot2::margin(t = 3)),
        axis.title.y = ggplot2::element_text(size = font$axis_title_y,
                                             margin = ggplot2::margin(r = 4), lineheight = 1.0),
        axis.text = ggplot2::element_text(size = font$tick, colour = palette$text),
        axis.text.x = ggplot2::element_text(margin = ggplot2::margin(t = 2)),
        axis.text.y = ggplot2::element_text(margin = ggplot2::margin(r = 2)),
        axis.ticks = ggplot2::element_line(linewidth = 0.30, colour = palette$border),
        axis.ticks.length = grid::unit(0.8, "mm"),
        panel.border = ggplot2::element_rect(fill = NA, colour = palette$border, linewidth = 0.35),
        panel.grid.major = ggplot2::element_line(colour = "#E1E1E1", linewidth = 0.23),
        panel.grid.minor = ggplot2::element_line(colour = "#F0F0F0", linewidth = 0.15),
        strip.text = ggplot2::element_text(size = font$strip, face = "plain",
                                            lineheight = 1.0, margin = ggplot2::margin(2.5, 1, 2.5, 1)),
        strip.background = ggplot2::element_rect(fill = palette$strip, colour = NA),
        legend.position = "none", # A single explicit key is drawn in cell 6.
        plot.margin = ggplot2::margin(0, 0, 0, 0)
      )
    p
  }
  plots <- lapply(panels, draw_panel)
  grobs <- lapply(plots, function(p) {
    g <- ggplot2::ggplotGrob(p)
    strip_row <- max(g$layout$b[grepl("^strip-t", g$layout$name)])
    gtable::gtable_add_rows(g, grid::unit(geometry$title_gap, "mm"), pos = strip_row)
  })

  # Equalise axes/strips as well as cell dimensions: all FIVE plotting frames
  # are the same size. Title rectangles share the plot frame's left/right edges.
  common_widths <- do.call(grid::unit.pmax, lapply(grobs, function(g) g$widths))
  common_heights <- do.call(grid::unit.pmax, lapply(grobs, function(g) g$heights))
  grobs <- lapply(grobs, function(g) {
    g$widths <- common_widths
    g$heights <- common_heights
    g
  })

  # ------------------------------ Shared key ---------------------------------
  key_labels <- c("Observed median", "Observed 5th/95th", "Simulated median",
                  "Simulated 5th/95th", sprintf("%d%% CIs: median", ci),
                  sprintf("%d%% CIs: 5th/95th", ci), "Observations")
  draw_key <- function() {
    grid::grid.text("", x = grid::unit(0, "mm"), y = grid::unit(cell_h - 2.1, "mm"),
                    just = "left", gp = grid::gpar(fontfamily = font_family,
                    fontsize = font$key_title, col = palette$text, fontface = "plain"))
    key_y <- cell_h - 11 - (0:6) * 6.0
    for (i in seq_along(key_labels)) {
      if (i <= 4L) {
        grid::grid.lines(x = grid::unit(c(0, 7), "mm"), y = grid::unit(rep(key_y[i], 2), "mm"),
                         gp = grid::gpar(col = if (i <= 2L) palette$observed else palette$simulated,
                         lwd = if (i %% 2L == 1L) line_median * 72.27 / 25.4 else line_outer * 72.27 / 25.4,
                         lty = if (i %% 2L == 1L) "solid" else "dashed"))
      } else if (i %in% c(5L, 6L)) {
        grid::grid.rect(x = grid::unit(3.5, "mm"), y = grid::unit(key_y[i], "mm"),
                        width = grid::unit(7, "mm"), height = grid::unit(2.7, "mm"),
                        gp = grid::gpar(col = NA, fill = if (i == 5L) palette$ci_median else palette$ci_outer))
      } else {
        grid::grid.points(x = grid::unit(3.5, "mm"), y = grid::unit(key_y[i], "mm"), pch = 1,
                          size = grid::unit(point_size, "mm"),
                          gp = grid::gpar(col = palette$points, lwd = 0.7))
      }
      grid::grid.text(key_labels[i], x = grid::unit(9, "mm"),
                      y = grid::unit(key_y[i], "mm"), just = "left",
                      gp = grid::gpar(fontfamily = font_family, fontsize = font$key,
                                     col = palette$text, fontface = "plain"))
    }
  }

  draw_figure <- function() {
    grid::grid.newpage()
    grid::grid.rect(gp = grid::gpar(fill = "white", col = NA))
    inset <- geometry$frame_inset
    grid::grid.rect(width = grid::unit(figure_width_mm - 2 * inset, "mm"),
                    height = grid::unit(figure_height_mm - 2 * inset, "mm"),
                    gp = grid::gpar(fill = NA, col = "#8C8C8C", lwd = 0.6, lty = "dashed"))
    positions <- list(c(1, 1), c(1, 2), c(1, 3), c(2, 1), c(2, 2))
    for (i in seq_along(grobs)) {
      row <- positions[[i]][1]
      col <- positions[[i]][2]
      grid::pushViewport(grid::viewport(x = grid::unit(col_left[col], "mm"),
                        y = grid::unit(figure_height_mm - row_top[row], "mm"),
                        width = grid::unit(cell_w, "mm"), height = grid::unit(cell_h, "mm"),
                        just = c("left", "top"), clip = "off"))
      grid::grid.draw(grobs[[i]])
      grid::popViewport()
    }
    for (row in 1:2) {
      grid::grid.text(c("a", "b")[row], x = grid::unit(2.4, "mm"),
                      y = grid::unit(figure_height_mm - row_top[row] - 2.1, "mm"),
                      just = c("left", "center"), gp = grid::gpar(fontfamily = font_family,
                      fontsize = font$row_label, fontface = "bold", col = palette$text))
    }
    # Align the shared key with Model III's plotting frame, not its full cell.
    # Measure on the active export device so PDF/PNG/TIFF remain aligned.
    panel <- grobs[[3]]$layout[grepl("^panel", grobs[[3]]$layout$name), ][1, ]
    plot_offset <- grid::convertWidth(sum(common_widths[seq_len(panel$l - 1)]),
                                     "mm", valueOnly = TRUE)
    grid::pushViewport(grid::viewport(x = grid::unit(col_left[3] + plot_offset, "mm"),
                      y = grid::unit(figure_height_mm - row_top[2], "mm"),
                      width = grid::unit(cell_w - plot_offset, "mm"), height = grid::unit(cell_h, "mm"),
                      just = c("left", "top"), clip = "off"))
    draw_key()
    grid::popViewport()
  }

  # ------------------------------- Export ------------------------------------
  if (!dir.exists(output_dir)) stop("Output folder does not exist: ", output_dir)
  base <- file.path(output_dir, output_name)
  use_ragg <- requireNamespace("ragg", quietly = TRUE)
  render <- function(path, kind) {
    if (kind == "pdf") {
      grDevices::pdf(path, width = figure_width_mm / 25.4,
                     height = figure_height_mm / 25.4, family = font_family,
                     useDingbats = FALSE, onefile = FALSE, bg = "white")
    } else if (use_ragg && kind == "png") {
      ragg::agg_png(path, width = figure_width_mm, height = figure_height_mm,
                    units = "mm", res = raster_dpi, background = "white")
    } else if (use_ragg && kind == "tiff") {
      ragg::agg_tiff(path, width = figure_width_mm, height = figure_height_mm,
                     units = "mm", res = raster_dpi, background = "white", compression = "lzw")
    } else if (kind == "png") {
      grDevices::png(path, width = figure_width_mm, height = figure_height_mm,
                     units = "mm", res = raster_dpi, bg = "white")
    } else {
      grDevices::tiff(path, width = figure_width_mm, height = figure_height_mm,
                      units = "mm", res = raster_dpi, bg = "white", compression = "lzw")
    }
    on.exit(grDevices::dev.off(), add = TRUE)
    draw_figure()
  }
  for (format in c("pdf", "png", "tiff")) {
    path <- paste0(base, ".", format)
    render(path, format)
    message("Saved: ", path)
  }
  all_sources <- unlist(lapply(sources, function(x) x$source_files), use.names = FALSE)
  saveRDS(list(plots = plots, panels = panels, grobs = grobs,
               settings = list(width_mm = figure_width_mm, height_mm = figure_height_mm,
                 dpi = raster_dpi, font = font, palette = palette, ci = ci,
                 percentiles = c(lo.pi, 50, up.pi), geometry = geometry),
               source_checksums = tools::md5sum(all_sources),
               session = utils::sessionInfo()), paste0(base, ".rds"))
  message("Finished. CI = ", ci, "%; n(CC+CT) = ", ccct$n, "; n(TT) = ", tt$n,
          ". Print at ", figure_width_mm, " mm width; do not shrink to a single column.")
  invisible(plots)
}

main_pcVPC()
