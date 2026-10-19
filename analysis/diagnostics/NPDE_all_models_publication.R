# ============================================================================
# Atomoxetine: three final-model NPDE quartets, stacked vertically; no legend.
# Standalone script. Source tables and old scripts/figures are never modified.
# Dependencies: npde, ggplot2; optional ragg for raster export.
# Run: Rscript NPDE_all_models_publication.R
# Optional output directory as first argument (default: input_root).
#
# IMPORTANT: all four diagnostics explicitly use NPDE, not NPD.
# In npde 3.5, plot() defaults to NPD; changing an axis label does not change
# the diagnostic. Never extract $npd: it can partially match a column "npde".
#
# Reference envelopes retain the original script's settings:
#   Q-Q: pointwise 90% normal-reference envelope, 1,000 normal samples.
#   Histogram/scatter: 95% reference envelopes (npde pi.size = 0.95).
#   Scatter curves: 2.5th, 50th and 97.5th percentiles (vpc.interval = 0.95).
# These are NOT the 5th/50th/95th percentiles of the separate pc-VPC figure.
# No fitted model parameters, simulation files, or manuscript P values change.
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

main_npde <- function() {
  for (pkg in c("npde", "ggplot2")) {
    if (!requireNamespace(pkg, quietly = TRUE)) stop("Install package: ", pkg)
  }
  input_root <- project_input_dir(.project_root, "ATX_NPDE_INPUT_DIR", "data/analysis_inputs/NPDE")
  args <- commandArgs(trailingOnly = TRUE)
  output_dir <- if (length(args)) {
    requested <- path.expand(args[[1]])
    absolute <- grepl("^(?:/|[A-Za-z]:[/\\\\]|\\\\\\\\)", requested, perl = TRUE)
    if (absolute) requested else file.path(.project_root, requested)
  } else {
    project_output_dir(.project_root, "ATX_NPDE_OUTPUT_DIR", "results/NPDE")
  }
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_name <- "Figure3_NPDE_revised"
  output_formats <- c("png", "pdf", "tiff")
  W <- 140                         # final page width, mm
  H <- 280                         # final page height, mm
  dpi <- 600
  family <- "Helvetica"
  font <- list(axis = 7.5, tick = 7, panel = 8, model = 9, tag = 8)
  colors <- list(line = "#005A8D", point = "#303030", reference = "#737373",
                 qq = "#D5EAF6", band = "#DFECF5", median = "#FBE4E9",
                 outside = "#F2BDC7", strip = "#D9D9D9", text = "#111111")
  qq_envelope <- 0.90
  reference_envelope <- 0.95
  scatter_interval <- 0.95
  seed <- 123L
  nsim <- 1000L
  models <- data.frame(
    id = c("BW", "RS1065852", "AS"),
    title = c("MODEL I (weight based)", "MODEL II (rs1065852 genotype)",
              "MODEL III (activity score)"), stringsAsFactors = FALSE)
  if (!dir.exists(output_dir)) stop("Output directory does not exist: ", output_dir)

  # ---------------------- Read and validate source data -----------------------
  source_files <- unlist(lapply(models$id, function(id) {
    file.path(input_root, id, c("vpc_original.npctab.dta", "vpc_simulation.1.npctab.dta"))
  }), use.names = FALSE)
  if (!all(file.exists(source_files))) stop("Missing source: ",
                        paste(source_files[!file.exists(source_files)], collapse = ", "))
  source_checksums <- tools::md5sum(source_files)
  objects <- observations <- diagnostic_data <- list()
  for (id in models$id) {
    message("Computing NPDE: ", id)
    obs <- utils::read.table(file.path(input_root, id, "vpc_original.npctab.dta"), header = TRUE)
    sim <- utils::read.table(file.path(input_root, id, "vpc_simulation.1.npctab.dta"), header = TRUE)
    required <- c("ID", "MDV", "TAD", "DV")
    if (!all(required %in% names(obs)) || !all(required %in% names(sim))) {
      stop("Required ID, MDV, TAD or DV column missing for ", id)
    }
    obs <- obs[required]
    if (nrow(sim) != nsim * nrow(obs)) stop("Expected 1,000 simulation replicates for ", id)
    for (field in c("ID", "MDV", "TAD")) {
      if (!isTRUE(all.equal(sim[[field]], rep(obs[[field]], nsim)))) {
        stop("Simulation order differs from observed order: ", id, ", ", field)
      }
    }
    set.seed(seed)
    x <- npde::autonpde(
      namobs = obs, namsim = sim[c("ID", "TAD", "DV")],
      iid = 1, ix = 3, iy = 4, imdv = 2, boolsave = FALSE, verbose = FALSE,
      units = list(x = "h", y = "ng/mL"), decorr.method = "cholesky", ties = TRUE)
    # Keep the package's non-missing observation mask, not dosing/MDV rows.
    keep <- x@data@not.miss
    result <- x@results@res[keep, , drop = FALSE]
    values <- result[["npde"]]       # exact column matching is intentional
    if (is.null(values) || any(!is.finite(values))) stop("Invalid NPDE for ", id)
    observed <- obs[!is.na(obs$MDV) & obs$MDV == 0, , drop = FALSE]
    if (length(values) != nrow(observed)) stop("Observation count mismatch: ", id)
    objects[[id]] <- x
    observations[[id]] <- observed
    diagnostic_data[[id]] <- data.frame(ID = observed$ID, TAD = observed$TAD,
       prediction = result[["ypred"]], NPDE = values)
    rm(sim)
  }

  # -------------------------- Native R plots ---------------------------------
  plot_theme <- ggplot2::theme_classic(base_size = font$tick, base_family = family) +
    ggplot2::theme(
      text = ggplot2::element_text(face = "plain", colour = colors$text),
      axis.title = ggplot2::element_text(size = font$axis, face = "plain"),
      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 3)),
      axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 3)),
      axis.text = ggplot2::element_text(size = font$tick, colour = colors$text),
      axis.text.x = ggplot2::element_text(margin = ggplot2::margin(t = 2)),
      axis.text.y = ggplot2::element_text(margin = ggplot2::margin(r = 2)),
      axis.line = ggplot2::element_line(linewidth = 0.25, colour = "#333333"),
      axis.ticks = ggplot2::element_line(linewidth = 0.25, colour = "#333333"),
      axis.ticks.length = grid::unit(0.7, "mm"),
      panel.grid = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(size = font$panel, face = "plain", hjust = 0.5,
                                        margin = ggplot2::margin(b = 3)),
      legend.position = "none", plot.margin = ggplot2::margin(0, 0, 0, 0))

  # Use do.call: npde 3.5's plotting method parses match.call and needs evaluated
  # strings, not a literal variable name, for plot.type and which.
  package_plot <- function(x, type, extra = list()) {
    common <- list(x = x, plot.type = type, which = "npde", verbose = FALSE,
      vpc.interval = scatter_interval, pi.size = reference_envelope,
      approx.pi = TRUE, bin.method = "equal", bin.number = 10,
      fill.bands = colors$band, alpha.bands = 1,
      fill.med = colors$median, alpha.med = 1,
      fill.outliers.med = colors$outside, fill.outliers.bands = colors$outside,
      col.lobs = colors$line, lwd.lobs = 0.40, lty.lobs = 1,
      col.pobs = colors$point, pch.pobs = 1, size.pobs = 0.55,
      col.ther = colors$reference, lwd.ther = 0.30, lty.ther = 5,
      lty.med = 1, lty.bands = 1, grid = FALSE)
    p <- do.call(graphics::plot, utils::modifyList(common, extra))
    if (!inherits(p, "ggplot")) stop("npde plotting did not return a ggplot object.")
    # The package adds a single-level facet ("all"/"none"). Remove it so every
    # plot has the same gtable layout and no superfluous facet title.
    p <- p + ggplot2::facet_null()
    # npde does not pass point stroke/alpha through all plotting paths.
    # Style the point layer only; data/quantiles/reference envelopes are intact.
    for (i in seq_along(p$layers)) {
      if (inherits(p$layers[[i]]$geom, "GeomPoint")) {
        p$layers[[i]]$aes_params$stroke <- 0.18
        p$layers[[i]]$aes_params$alpha <- 0.80
      }
    }
    p
  }

  # Use common display limits without censoring observations or changing bins.
  # Stop if new data would fall outside these limits.
  residual_limits <- c(-3.5, 3.5)
  prediction_limits <- c(0, 550)
  time_limits <- c(0, 20)
  for (d in diagnostic_data) {
    if (any(abs(d$NPDE) > 3.5) || any(d$TAD < 0 | d$TAD > 20) ||
        any(d$prediction < 0 | d$prediction > 550)) {
      stop("An observation lies outside display limits; expand limits before drawing.")
    }
  }
  plots <- bands <- list()
  for (id in models$id) {
    x <- objects[[id]]
    values <- diagnostic_data[[id]]$NPDE
    qq <- data.frame(theo = stats::qnorm(stats::ppoints(length(values))), samp = sort(values))
    set.seed(seed)
    sim_q <- replicate(nsim, sort(stats::rnorm(length(values))))
    tails <- (1 - qq_envelope) / 2
    qq_band <- data.frame(theo = qq$theo,
      lo = apply(sim_q, 1, stats::quantile, probs = tails),
      med = apply(sim_q, 1, stats::quantile, probs = 0.5),
      hi = apply(sim_q, 1, stats::quantile, probs = 1 - tails))
    bands[[id]] <- qq_band
    p1 <- ggplot2::ggplot() +
      ggplot2::geom_ribbon(data = qq_band, ggplot2::aes(x = theo, ymin = lo, ymax = hi),
                            fill = colors$qq) +
      ggplot2::geom_line(data = qq_band, ggplot2::aes(x = theo, y = med),
                          colour = "#6F9DB6", linewidth = 0.25) +
      ggplot2::geom_abline(intercept = 0, slope = 1, colour = colors$reference, linewidth = 0.28) +
      ggplot2::geom_line(data = qq, ggplot2::aes(x = theo, y = samp),
                          colour = colors$line, linewidth = 0.42) +
      ggplot2::scale_x_continuous(breaks = -3:3) +
      ggplot2::scale_y_continuous(breaks = -3:3) +
      ggplot2::coord_cartesian(xlim = c(-3.2, 3.2), ylim = residual_limits, expand = FALSE) +
      ggplot2::labs(x = "Theoretical quantiles", y = "NPDE quantiles", title = "Normal Q-Q plot") +
      plot_theme
    p2 <- package_plot(x, "hist", list(breaks.x = 7, breaks.y = 4, col.ther = "#404040", lty.ther = 2))
    # Exact tests guard against accidental reversion to NPD/default plotting.
    if (!isTRUE(all.equal(sort(p2$data$x), sort(values)))) stop("Histogram is not NPDE: ", id)
    p2 <- p2 + ggplot2::scale_x_continuous(breaks = -3:3) +
      ggplot2::scale_y_continuous(breaks = c(0, 20, 40, 60)) +
      ggplot2::coord_cartesian(xlim = residual_limits, expand = FALSE) +
      ggplot2::labs(x = "NPDE", y = "Frequency", title = "NPDE histogram") + plot_theme
    scatter <- lapply(c("x.scatter", "pred.scatter"), function(type) {
      p <- package_plot(x, type, list(breaks.x = if (type == "x.scatter") 3 else 6, breaks.y = 6))
      point_idx <- which(vapply(p$layers, function(l) inherits(l$geom, "GeomPoint"), logical(1)))
      if (length(point_idx) != 1L) stop("Unexpected scatter layers in npde package.")
      point_data <- p$layers[[point_idx]]$data
      if (!isTRUE(all.equal(as.numeric(point_data$y1), as.numeric(values)))) {
        stop("Scatter values are not NPDE: ", id, ", ", type)
      }
      p + ggplot2::scale_y_continuous(breaks = -3:3) +
        ggplot2::scale_x_continuous(breaks = if (type == "x.scatter") c(0, 2, 6, 12, 18)
                                                        else seq(0, 500, 100)) +
        ggplot2::coord_cartesian(xlim = if (type == "x.scatter") time_limits else prediction_limits,
                                   ylim = residual_limits, expand = FALSE) +
        ggplot2::labs(x = if (type == "x.scatter") "Time after dose (h)"
                         else "Predicted concentration (ng/mL)", y = "NPDE",
                     title = if (type == "x.scatter") "NPDE vs time" else "NPDE vs predicted concentration") +
        plot_theme
    })
    plots[[id]] <- c(list(p1, p2), scatter)
  }

  # ------------------------------- Layout -----------------------------------
  # Three model blocks in a single column; each has a 2 x 2 quartet.
  # Grey model strips use regular type and align with the outer plotting axes.
  flat_plots <- unlist(plots, recursive = FALSE)
  grDevices::pdf(NULL, family = family)
  grobs <- lapply(flat_plots, ggplot2::ggplotGrob)
  if (!all(vapply(grobs, function(g) identical(g$layout$name, grobs[[1]]$layout$name), logical(1)))) {
    stop("Plot table layouts differ; cannot safely align widths/heights.")
  }
  widths <- do.call(grid::unit.pmax, lapply(grobs, function(g) g$widths))
  heights <- do.call(grid::unit.pmax, lapply(grobs, function(g) g$heights))
  grobs <- lapply(grobs, function(g) { g$widths <- widths; g$heights <- heights; g })
  grDevices::dev.off()
  left <- right <- 4
  top <- bottom <- 5
  block_gap <- 7
  strip_h <- 4.5
  strip_gap <- 3
  col_gap <- 5
  row_gap <- 4
  block_w <- W - left - right
  block_h <- (H - top - bottom - 2 * block_gap) / 3
  cell_w <- (block_w - col_gap) / 2
  cell_h <- (block_h - strip_h - strip_gap - row_gap) / 2

  draw_figure <- function() {
    grid::grid.newpage()
    grid::grid.rect(gp = grid::gpar(fill = "white", col = NA))
    grid::grid.rect(width = grid::unit(W - 2.4, "mm"), height = grid::unit(H - 2.4, "mm"),
                   gp = grid::gpar(fill = NA, col = "#8C8C8C", lwd = 0.6, lty = "dashed"))
    panel <- grobs[[1]]$layout[which(grobs[[1]]$layout$name == "panel")[1], ]
    plot_left <- grid::convertWidth(sum(widths[seq_len(panel$l - 1)]), "mm", valueOnly = TRUE)
    plot_right <- grid::convertWidth(sum(widths[seq.int(panel$r + 1, length(widths))]), "mm", valueOnly = TRUE)
    strip_w <- block_w - plot_left - plot_right
    for (m in seq_len(nrow(models))) {
      model_top <- top + (m - 1) * (block_h + block_gap)
      sx <- left + plot_left + strip_w / 2
      sy <- H - model_top - strip_h / 2
      grid::grid.rect(x = grid::unit(sx, "mm"), y = grid::unit(sy, "mm"),
                     width = grid::unit(strip_w, "mm"), height = grid::unit(strip_h, "mm"),
                     gp = grid::gpar(fill = colors$strip, col = NA))
      grid::grid.text(models$title[m], x = grid::unit(sx, "mm"), y = grid::unit(sy, "mm"),
                     gp = grid::gpar(fontfamily = family, fontsize = font$model,
                                    fontface = "plain", col = colors$text))
      for (j in 1:4) {
        cx <- left + ((j - 1) %% 2) * (cell_w + col_gap)
        cy <- model_top + strip_h + strip_gap + ((j - 1) %/% 2) * (cell_h + row_gap)
        grid::pushViewport(grid::viewport(x = grid::unit(cx, "mm"), y = grid::unit(H - cy, "mm"),
                          width = grid::unit(cell_w, "mm"), height = grid::unit(cell_h, "mm"),
                          just = c("left", "top"), clip = "off"))
        grid::grid.draw(grobs[[(m - 1) * 4 + j]])
        grid::grid.text(LETTERS[j], x = grid::unit(0, "mm"), y = grid::unit(cell_h - 1.4, "mm"),
                       just = c("left", "center"),
                       gp = grid::gpar(fontfamily = family, fontsize = font$tag, fontface = "bold"))
        grid::popViewport()
      }
    }
  }

  # ------------------------------- Export -----------------------------------
  output_base <- file.path(output_dir, output_name)
  render <- function(format) {
    path <- paste0(output_base, ".", format)
    if (format == "pdf") {
      grDevices::pdf(path, width = W/25.4, height = H/25.4, family = family,
                     useDingbats = FALSE, onefile = FALSE)
    } else if (format == "png" && requireNamespace("ragg", quietly = TRUE)) {
      ragg::agg_png(path, width = W, height = H, units = "mm", res = dpi, background = "white")
    } else if (format == "tiff" && requireNamespace("ragg", quietly = TRUE)) {
      ragg::agg_tiff(path, width = W, height = H, units = "mm", res = dpi,
                    background = "white", compression = "lzw")
    } else if (format == "png") {
      grDevices::png(path, width = W, height = H, units = "mm", res = dpi, bg = "white")
    } else if (format == "tiff") {
      grDevices::tiff(path, width = W, height = H, units = "mm", res = dpi,
                     bg = "white", compression = "lzw")
    } else stop("Unsupported output format: ", format)
    on.exit(grDevices::dev.off(), add = TRUE)
    draw_figure()
    message("Saved: ", path)
  }
  for (format in output_formats) render(format)
  if (!identical(tools::md5sum(source_files), source_checksums)) stop("Source data changed during this run.")
  saveRDS(list(plots = plots, npde_objects = objects, diagnostic_data = diagnostic_data,
               qq_bands = bands, metric = "NPDE", nsim = nsim,
               qq_reference_level = qq_envelope, other_reference_level = reference_envelope,
               scatter_quantile_probabilities = c((1-scatter_interval)/2, 0.5, (1+scatter_interval)/2),
               dimensions_mm = c(width = W, height = H), source_checksums = source_checksums,
               package_versions = c(npde = as.character(utils::packageVersion("npde")),
                                    ggplot2 = as.character(utils::packageVersion("ggplot2"))),
               session_info = utils::sessionInfo()), paste0(output_base, ".rds"))
  message("Finished. All 12 panels use NPDE. Original data and figures were not modified.")
}

main_npde()
