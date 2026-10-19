# Newly written, deterministic fictional event-grid generator. NOT study code.
# No clinical files or historical simulation inputs are read.
# This prepares files only; NONMEM must be run separately in a licensed environment.
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Run using Rscript.")
root <- dirname(dirname(normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)))
run_root <- file.path(root, "simulation_demo_runs")
if (dir.exists(run_root)) {
  stop("simulation_demo_runs already exists. Use a fresh package or move old runs aside; no files changed.")
}

configs <- list(
  list(model = "sim_model_I_BW", file = "synthetic_sim_BW.dat", field = "AS", values = 1.25),
  list(model = "sim_model_II_rs1065852", file = "synthetic_sim_RS.dat", field = "rs_grp", values = c(0, 1)),
  list(model = "sim_model_III_AS", file = "synthetic_sim_AS.dat", field = "AS", values = c(0.5, 1.25, 2))
)
stopifnot(all(file.exists(vapply(configs, function(x) {
  file.path(root, "simulation", paste0(x$model, ".mod"))
}, character(1)))))

# Arbitrary small grid, not the study's weight distribution or dose range.
# DOSEKG is mg/kg/day here because every demo regimen is once daily.
# BW-only AS is a dummy placeholder; that model does not use AS in its equations.
weights <- c(20, 30.7, 50)
doses <- c(0.8, 1.3)
sample_times <- c(2, 26, 674, 698)
# Dose at time 0 plus 29 additional doses, every 24 h (last dose: 696 h).
# Includes early and late samples: DO NOT pool all MDV=0 rows as steady-state C2h.
interval <- 24
addl <- 29L
data_dir <- file.path(root, "data", "simulation")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(run_root)
manifest <- c(
  "FICTIONAL SIMULATION INPUTS - NOT CLINICAL DATA OR THE HISTORICAL SIMULATION DESIGN.",
  "Generator is deterministic and reads no clinical data; NONMEM seed remains 20260717.",
  "Three BW values (20, 30.7, 50 kg), two illustrative doses (0.8, 1.3 mg/kg/day).",
  "These doses are test inputs, not clinical dosing recommendations.",
  "Dose at t=0; II=24 h; ADDL=29. Observation times: 2, 26, 674, 698 h; TAD=2 h.",
  "Early and late sampling are both present. Steady state is not imposed or verified.",
  "AS values and genotype codes are explicit demo conventions, not a genetic translation algorithm.",
  "DV=0 is an input placeholder on all rows, not a measured or simulated zero concentration.",
  "The unchanged source controls use 500 simulation subproblems, not estimation."
)
for (j in seq_along(configs)) {
  conf <- configs[[j]]
  grid <- expand.grid(WT = weights, DOSEKG = doses, covariate = conf$values,
                      KEEP.OUT.ATTRS = FALSE)
  schema <- c("ID", "TIME", "AMT", "ADDL", "II", "DV", "MDV", "EVID",
              "WT", conf$field, "DOSEKG", "TAD")
  rows <- lapply(seq_len(nrow(grid)), function(i) {
    g <- grid[i, ]
    id <- 910000L + 1000L * j + i
    dose <- c(id, 0, g$WT * g$DOSEKG, addl, interval, 0, 1, 1,
              g$WT, g$covariate, g$DOSEKG, 0)
    observations <- lapply(sample_times, function(time) {
      c(id, time, 0, 0, 0, 0, 0, 0, g$WT, g$covariate, g$DOSEKG, 2)
    })
    as.data.frame(do.call(rbind, c(list(dose), observations)))
  })
  dat <- do.call(rbind, rows)
  names(dat) <- schema
  rownames(dat) <- NULL
  stopifnot(all(is.finite(as.matrix(dat))), nrow(dat) == 5L * nrow(grid))
  path <- file.path(data_dir, conf$file)
  con <- file(path, open = "wt")
  writeLines(paste0("@", paste(schema, collapse = " ")), con)
  write.table(dat, con, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ")
  close(con)
  target <- file.path(run_root, conf$model)
  dir.create(target)
  stopifnot(file.copy(path, file.path(target, conf$file), overwrite = FALSE),
            file.copy(file.path(root, "simulation", paste0(conf$model, ".mod")),
                      file.path(target, "model.mod"), overwrite = FALSE))
  manifest <- c(manifest, paste(conf$model, ":", nrow(grid), "fictional subjects;",
                               nrow(dat), "records;", 4L * nrow(grid), "observation placeholders; MD5",
                               unname(tools::md5sum(path))))
}
manifest <- c(manifest, "NONMEM was not run. No PTA or study figures were calculated.")
writeLines(manifest, file.path(data_dir, "SYNTHETIC_SIMULATION_MANIFEST.txt"))
writeLines(capture.output(sessionInfo()), file.path(run_root, "R_SESSION_INFO.txt"))
cat(paste(manifest, collapse = "\n"), "\n")
