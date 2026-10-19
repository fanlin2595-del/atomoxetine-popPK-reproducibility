# ============================================================================
# Atomoxetine GOF: one standalone script, four models x four diagnostics.
# Final layout: no in-figure legend; original fit files/scripts/images unchanged.
# Run this single script to generate PNG/TIFF (600 dpi), PDF, and a plot bundle.
# Required packages: install.packages("ggplot2"); optionally install.packages("ragg").
# Run from the project root: Rscript --vanilla analysis/diagnostics/GOF_all_models_publication.R
# Uses the source scripts' MDV == 0 selection and LOESS span = 1, se = FALSE.
# 2026-10-10: match the approved pc-VPC typography and restrained grey strips.
# Required: ggplot2. Optional: ragg (recommended for raster output).
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

main_gof <- function() {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Install ggplot2 first.")

  # ---------------------------- User settings -------------------------------
  input_root <- project_input_dir(.project_root, "ATX_GOF_INPUT_DIR", "data/analysis_inputs/GOF")
  output_dir <- project_output_dir(.project_root, "ATX_GOF_OUTPUT_DIR", "results/GOF")
  output_name <- "Figure1_GOF_revised"
  output_formats <- c("png", "pdf", "tiff")
  W <- 180
  H <- 187 # Remove the old legend area while retaining all panel dimensions.
  dpi <- 600
  family <- "Helvetica"
  font <- list(axis = 7.5, tick = 7.0, panel = 8, tag = 8, model = 9)
  colors <- list(point = "#4F81BD", smooth = "#005A8D", reference = "#222222",
                 threshold = "#777777", grid = "#E1E1E1", minor = "#F0F0F0",
                 strip = "#D9D9D9", border = "#4A4A4A", text = "#111111")
  point_size <- 0.70
  point_stroke <- 0.16
  point_alpha <- 0.75
  smooth_width <- 0.40
  reference_width <- 0.28
  threshold_width <- 0.23
  loess_span <- 1

  # Shared axes are display changes only. coord_cartesian does not censor
  # observations before smoothing. Never alter/remove outliers to fit an axis.
  concentration_limits <- c(0, 1800)
  concentration_breaks <- c(0, 500, 1000, 1500)
  residual_limits <- c(-3.5, 3.5)
  residual_breaks <- -3:3
  time_limits <- c(0, 25)
  time_breaks <- c(0, 2, 6, 12, 18, 24)
  pred_residual_limits <- c(0, 750)
  pred_residual_breaks <- c(0, 200, 400, 600)

  models <- data.frame(
    id = c("BASE", "BW", "RS", "AS"),
    title = c("Base Model", "MODEL I (weight based)",
              "MODEL II (rs1065852 genotype)", "MODEL III (activity score)"),
    file = c("BASE/5-M1.fit", "BW/7.fit", "RS1065852/15.fit", "AS/12-M1.fit"),
    stringsAsFactors = FALSE
  )
  models$path <- file.path(input_root, models$file)
  if (!all(file.exists(models$path))) stop("Missing fit file: ",
                                          paste(models$path[!file.exists(models$path)], collapse = ", "))
  checksums <- tools::md5sum(models$path)

  read_fit <- function(path) {
    first <- readLines(path, n = 1, warn = FALSE)
    skip <- if (grepl("^\\s*TABLE NO", first)) 1L else 0L
    fit <- utils::read.table(path, header = TRUE, skip = skip, check.names = TRUE)
    required <- c("ID", "TIME", "MDV", "DV", "PRED", "IPRED", "CWRES", "TAD")
    if (!all(required %in% names(fit))) stop("Required columns missing in ", path)
    obs <- fit[!is.na(fit$MDV) & fit$MDV == 0, , drop = FALSE]
    if (!nrow(obs)) stop("No MDV == 0 rows in ", path)
    if (!all(vapply(obs[required], function(x) all(is.finite(x)), logical(1)))) {
      stop("Non-finite diagnostic values in ", path, "; inspect the source data.")
    }
    # NONMEM commonly writes DV twice; detect inconsistent duplicate columns.
    if ("DV.1" %in% names(obs) && !isTRUE(all.equal(obs$DV, obs$DV.1))) {
      stop("Duplicated DV columns disagree in ", path)
    }
    obs
  }
  observations <- setNames(lapply(models$path, read_fit), models$id)
  reference <- observations[[1]][c("ID", "TIME", "TAD", "DV")]
  for (id in names(observations)) {
    if (!identical(observations[[id]][names(reference)], reference)) {
      stop("Observation identities/order differ for model ", id, "; check inputs.")
    }
  }
  check_range <- function(x, limits, label) {
    if (any(x < limits[1] | x > limits[2])) stop(label, " lies outside its display limits.")
  }
  for (o in observations) {
    check_range(unlist(o[c("DV", "PRED", "IPRED")]), concentration_limits, "Concentration")
    check_range(o$CWRES, residual_limits, "CWRES")
    check_range(o$TAD, time_limits, "TAD")
    check_range(o$PRED, pred_residual_limits, "PRED in residual plot")
  }

  # ------------------------------- Plots ------------------------------------
  theme_gof <- ggplot2::theme_bw(base_size = font$tick, base_family = family) +
    ggplot2::theme(
      text = ggplot2::element_text(colour = colors$text, face = "plain"),
      axis.title = ggplot2::element_text(size = font$axis, face = "plain"),
      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 3)),
      axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 3)),
      axis.text = ggplot2::element_text(size = font$tick, colour = colors$text),
      axis.text.x = ggplot2::element_text(margin = ggplot2::margin(t = 2)),
      axis.text.y = ggplot2::element_text(margin = ggplot2::margin(r = 2)),
      axis.ticks = ggplot2::element_line(linewidth = 0.25, colour = colors$border),
      axis.ticks.length = grid::unit(0.7, "mm"),
      panel.border = ggplot2::element_rect(fill = NA, colour = colors$border, linewidth = 0.3),
      panel.grid.major = ggplot2::element_line(colour = colors$grid, linewidth = 0.20),
      panel.grid.minor = ggplot2::element_line(colour = colors$minor, linewidth = 0.13),
      plot.title = ggplot2::element_text(size = font$panel, hjust = 0.5,
                                         face = "plain", margin = ggplot2::margin(b = 3)),
      legend.position = "none",
      plot.margin = ggplot2::margin(0, 0, 0, 0)
    )
  titles <- c("DV vs PRED", "DV vs IPRED", "CWRES vs TAD", "CWRES vs PRED")
  panel_specs <- list(
    list(x = "PRED", y = "DV", xlab = "PRED (ng/mL)", ylab = "DV (ng/mL)"),
    list(x = "IPRED", y = "DV", xlab = "IPRED (ng/mL)", ylab = "DV (ng/mL)"),
    list(x = "TAD", y = "CWRES", xlab = "Time after dose (h)", ylab = "CWRES"),
    list(x = "PRED", y = "CWRES", xlab = "PRED (ng/mL)", ylab = "CWRES")
  )

  draw_panel <- function(o, index) {
    spec <- panel_specs[[index]]
    p <- ggplot2::ggplot(o, ggplot2::aes(x = .data[[spec$x]], y = .data[[spec$y]]))
    if (index <= 2) {
      p <- p + ggplot2::geom_abline(intercept = 0, slope = 1,
                                     colour = colors$reference, linewidth = reference_width)
      xlim <- ylim <- concentration_limits
      xb <- yb <- concentration_breaks
      # Retain the original scale limits: all observations fall within them,
      # but negative smoothed DV estimates are suppressed just as before.
      scale_x_limits <- scale_y_limits <- c(0, max(c(o$DV, o[[spec$x]])))
    } else {
      p <- p + ggplot2::geom_hline(yintercept = 0, colour = colors$reference,
                                     linewidth = reference_width) +
        ggplot2::geom_hline(yintercept = c(-2, 2), colour = colors$threshold,
                             linewidth = threshold_width, linetype = "dashed")
      xlim <- if (index == 3) time_limits else pred_residual_limits
      xb <- if (index == 3) time_breaks else pred_residual_breaks
      ylim <- residual_limits
      yb <- residual_breaks
      scale_x_limits <- if (index == 3) time_limits else NULL
      scale_y_limits <- NULL
    }
    p + ggplot2::geom_point(shape = 1, size = point_size, stroke = point_stroke,
                             alpha = point_alpha, colour = colors$point) +
      ggplot2::geom_smooth(method = "loess", formula = y ~ x, se = FALSE, na.rm = TRUE,
                             span = loess_span, colour = colors$smooth,
                             linewidth = smooth_width) +
      ggplot2::scale_x_continuous(breaks = xb, limits = scale_x_limits,
                                  expand = ggplot2::expansion(mult = 0.04)) +
      ggplot2::scale_y_continuous(breaks = yb, limits = scale_y_limits,
                                  expand = ggplot2::expansion(mult = 0.04)) +
      ggplot2::coord_cartesian(xlim = xlim, ylim = ylim) +
      ggplot2::labs(x = spec$xlab, y = spec$ylab, title = titles[index]) + theme_gof
  }
  plots <- unlist(lapply(observations, function(o) lapply(1:4, function(i) draw_panel(o, i))),
                  recursive = FALSE)

  # Use an in-memory graphics device to prevent accidental Rplots.pdf output.
  grDevices::pdf(NULL, family = family)
  grobs <- lapply(plots, ggplot2::ggplotGrob)
  widths <- do.call(grid::unit.pmax, lapply(grobs, function(g) g$widths))
  heights <- do.call(grid::unit.pmax, lapply(grobs, function(g) g$heights))
  grobs <- lapply(grobs, function(g) { g$widths <- widths; g$heights <- heights; g })
  grDevices::dev.off()

  # ------------------------------- Layout -----------------------------------
  # Preserve original 2 x 2 model blocks, each containing a 2 x 2 GOF quartet.
  left <- 4
  right <- 4
  top <- 5
  bottom <- 5
  block_col_gap <- 7
  block_row_gap <- 7
  strip_h <- 4.5
  strip_gap <- 3
  panel_col_gap <- 2.5
  panel_row_gap <- 4
  block_w <- (W - left - right - block_col_gap) / 2
  block_h <- (H - top - bottom - block_row_gap) / 2
  cell_w <- (block_w - panel_col_gap) / 2
  cell_h <- (block_h - strip_h - strip_gap - panel_row_gap) / 2
  model_left <- left + c(0, 1, 0, 1) * (block_w + block_col_gap)
  model_top <- top + c(0, 0, 1, 1) * (block_h + block_row_gap)

  draw_figure <- function() {
    grid::grid.newpage()
    grid::grid.rect(gp = grid::gpar(fill = "white", col = NA))
    grid::grid.rect(width = grid::unit(W - 2.4, "mm"), height = grid::unit(H - 2.4, "mm"),
                    gp = grid::gpar(fill = NA, col = "#8C8C8C", lwd = 0.6, lty = "dashed"))
    panel_idx <- which(grepl("^panel", grobs[[1]]$layout$name))[1]
    panel_layout <- grobs[[1]]$layout[panel_idx, ]
    plot_left <- grid::convertWidth(sum(widths[seq_len(panel_layout$l - 1)]), "mm", valueOnly = TRUE)
    plot_right <- grid::convertWidth(sum(widths[seq.int(panel_layout$r + 1, length(widths))]),
                                     "mm", valueOnly = TRUE)
    strip_w <- block_w - plot_left - plot_right
    for (m in 1:4) {
      # Model strips align with the outermost plot frames of their own quartet.
      sx <- model_left[m] + plot_left
      sy <- H - model_top[m] - strip_h / 2
      grid::grid.rect(x = grid::unit(sx + strip_w / 2, "mm"), y = grid::unit(sy, "mm"),
                       width = grid::unit(strip_w, "mm"), height = grid::unit(strip_h, "mm"),
                       gp = grid::gpar(fill = colors$strip, col = NA))
      grid::grid.text(models$title[m], x = grid::unit(sx + strip_w / 2, "mm"),
                       y = grid::unit(sy, "mm"),
                       gp = grid::gpar(fontfamily = family, fontsize = font$model,
                                       fontface = "plain", col = colors$text))
      for (j in 1:4) {
        col <- (j - 1) %% 2
        row <- (j - 1) %/% 2
        x <- model_left[m] + col * (cell_w + panel_col_gap)
        y <- model_top[m] + strip_h + strip_gap + row * (cell_h + panel_row_gap)
        grid::pushViewport(grid::viewport(x = grid::unit(x, "mm"), y = grid::unit(H-y, "mm"),
                           width = grid::unit(cell_w, "mm"), height = grid::unit(cell_h, "mm"),
                           just = c("left", "top"), clip = "off"))
        grid::grid.draw(grobs[[(m - 1) * 4 + j]])
        grid::grid.text(LETTERS[j], x = grid::unit(0, "mm"), y = grid::unit(cell_h - 1.4, "mm"),
                         just = c("left", "center"),
                         gp = grid::gpar(fontfamily = family, fontsize = font$tag,
                                         fontface = "bold", col = colors$text))
        grid::popViewport()
      }
    }
    # No in-figure legend: define points and line styles in the figure caption.
  }

  # ------------------------------- Export -----------------------------------
  if (!dir.exists(output_dir)) stop("Output directory does not exist: ", output_dir)
  output_base <- file.path(output_dir, output_name)
  use_ragg <- requireNamespace("ragg", quietly = TRUE)
  render <- function(format) {
    path <- paste0(output_base, ".", format)
    if (format == "pdf") {
      grDevices::pdf(path, width=W/25.4, height=H/25.4, family=family,
                      useDingbats=FALSE, onefile=FALSE)
    } else if (format == "png" && use_ragg) {
      ragg::agg_png(path, width=W, height=H, units="mm", res=dpi, background="white")
    } else if (format == "tiff" && use_ragg) {
      ragg::agg_tiff(path, width=W, height=H, units="mm", res=dpi,
                      background="white", compression="lzw")
    } else if (format == "png") {
      grDevices::png(path, width=W, height=H, units="mm", res=dpi, bg="white")
    } else if (format == "tiff") {
      grDevices::tiff(path, width=W, height=H, units="mm", res=dpi, bg="white", compression="lzw")
    } else stop("Unsupported export format: ", format)
    on.exit(grDevices::dev.off(), add=TRUE)
    draw_figure()
    message("Saved: ", path)
  }
  for (format in output_formats) render(format)
  saveRDS(list(models=models, observations=observations, plots=plots, grobs=grobs,
               settings=list(width_mm=W, height_mm=H, dpi=dpi, font=font, colors=colors,
                             loess_span=loess_span, cell_w=cell_w, cell_h=cell_h,
                             legend=FALSE),
               source_checksums=checksums, session=utils::sessionInfo()), paste0(output_base,".rds"))
  stopifnot(identical(checksums, tools::md5sum(models$path)))
  message("Source files unchanged. Models: ", paste(models$id, collapse=", "),
          "; observations/model: ", paste(vapply(observations, nrow, integer(1)), collapse=", "),
          "; subjects/model: ", paste(vapply(observations, function(o)length(unique(o$ID)), integer(1)), collapse=", "))
  invisible(plots)
}

main_gof()
