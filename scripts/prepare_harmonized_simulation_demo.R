# Create compact, wholly fictional NONMEM dose-simulation input grids.
# Reads no clinical files and does not create simulated concentration outputs.
# Run from any directory: Rscript --vanilla scripts/prepare_simulation_demo.R

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Run this script using Rscript.")
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
root <- dirname(dirname(script_path))
out_dir <- file.path(root, "data", "simulation_demo_combined")
if (dir.exists(out_dir)) {
  stop("data/simulation_demo_combined already exists. No files were changed.")
}

# A compact whole-kg grid demonstrates the model-input layout and all three
# BW strata. It intentionally does not reproduce the study's finer grid or
# any patient-level covariate distribution.
weights <- 17:60
dose_rates <- c(0.5, 0.8, 1.0, 1.1, 1.2, 1.3, 1.4) # mg/kg/day, demo only
observation_tad <- c(2, 12, 24)
n_prior_doses <- 9L
dosing_interval <- 24

make_scenario <- function(name, levels, level_field) {
  combinations <- expand.grid(
    WT = weights,
    level = levels,
    DOSEKG = dose_rates,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  combinations <- combinations[order(combinations$WT, combinations$level,
                                     combinations$DOSEKG), , drop = FALSE]
  id_start <- switch(name, BW = 970000L, rs1065852 = 980000L,
                     AS = 990000L, CYP2D6 = 960000L)
  combinations$ID <- seq_len(nrow(combinations)) + id_start

  records <- vector("list", nrow(combinations) * (length(observation_tad) + 1L))
  k <- 0L
  for (i in seq_len(nrow(combinations))) {
    x <- combinations[i, , drop = FALSE]
    id <- x$ID
    wt <- x$WT
    rate <- x$DOSEKG
    level <- x$level
    # Ten daily doses total; the final implied dose is at 216 h.
    k <- k + 1L
    records[[k]] <- data.frame(
      ID = id, TIME = 0, AMT = wt * rate,
      ADDL = n_prior_doses, II = dosing_interval,
      DV = 0, MDV = 1, WT = wt,
      scenario_level = level, DOSEKG = rate, TAD = 0,
      check.names = FALSE
    )
    for (tad in observation_tad) {
      k <- k + 1L
      records[[k]] <- data.frame(
        ID = id, TIME = n_prior_doses * dosing_interval + tad,
        AMT = 0, ADDL = 0, II = 0,
        DV = 0, MDV = 0, WT = wt,
        scenario_level = level, DOSEKG = rate, TAD = tad,
        check.names = FALSE
      )
    }
  }
  dat <- do.call(rbind, records[seq_len(k)])
  names(dat)[names(dat) == "scenario_level"] <- level_field
  dat <- dat[order(dat$ID, dat$TIME, -dat$MDV), , drop = FALSE]
  rownames(dat) <- NULL
  dat
}

grids <- list(
  BW = make_scenario("BW", 1.25, "AS"),
  rs1065852 = make_scenario("rs1065852", 0:1, "rs_grp"),
  AS = make_scenario("AS", c(0.5, 1.25, 1.5, 2), "AS"),
  CYP2D6 = make_scenario("CYP2D6", 1:2, "CYP2D6")
)

covariate_fields <- c(BW = "AS", rs1065852 = "rs_grp", AS = "AS",
                      CYP2D6 = "CYP2D6")
combined <- do.call(rbind, lapply(names(grids), function(scenario) {
  dat <- grids[[scenario]]
  covariate <- covariate_fields[[scenario]]
  value <- if (scenario == "BW") rep(NA_real_, nrow(dat)) else dat[[covariate]]
  data.frame(
    SCENARIO = scenario,
    ID = dat$ID, TIME = dat$TIME, AMT = dat$AMT,
    ADDL = dat$ADDL, II = dat$II, DV = dat$DV, MDV = dat$MDV,
    WT = dat$WT,
    COVARIATE = if (scenario == "BW") "none" else covariate,
    COVARIATE_VALUE = value,
    DOSEKG = dat$DOSEKG, TAD = dat$TAD,
    check.names = FALSE
  )
}))
rownames(combined) <- NULL

dir.create(out_dir, recursive = TRUE)
combined_path <- file.path(out_dir, "simulation_demo_all_scenarios.csv")
con <- file(combined_path, open = "wt")
writeLines(paste0("@", paste(names(combined), collapse = ",")), con)
write.table(combined, con, sep = ",", quote = FALSE, na = "",
            row.names = FALSE, col.names = FALSE)
close(con)

summary_lines <- c(
  "SYNTHETIC DEMONSTRATION ONLY — NOT OBSERVED PATIENT DATA.",
  "No clinical, genetic, or other source dataset was read by this generator.",
  "The combined file is a harmonized input grid; DV=0 is a simulation placeholder, not a concentration.",
  "SCENARIO and COVARIATE/COVARIATE_VALUE identify which model scenario and covariate level each row belongs to.",
  "For the BW-only scenario, COVARIATE_VALUE is blank because no stratifying covariate is used.",
  "Once-daily oral dose; AMT = WT × DOSEKG; 9 additional doses at II=24 h (10 total doses).",
  paste0("Post-dose observation times (TAD, h): ", paste(observation_tad, collapse = ", ")),
  paste0("Body weights (kg): ", min(weights), "–", max(weights), " in 1-kg steps."),
  paste0("Candidate dose rates (mg/kg/day): ", paste(dose_rates, collapse = ", ")),
  "rs_grp codes 0/1 and CYP2D6 codes 1/2 are arbitrary demonstration categories.",
  "BW uses AS=1.25 as an unused placeholder; it is not a phenotype assignment.",
  "AS scenario values are artificial positive values, not a reconstructed clinical distribution.",
  "The files do not reproduce the original 0.1-kg grid or study simulation results.",
  "No random concentration values or empirical estimates are generated by this script.",
  "The harmonized combined file is for inspection/filtering, not a direct single-control-stream NONMEM input.",
  "Filter one scenario and map its covariate field to the matching control stream before use.",
  "NONMEM is not run here."
)
for (scenario in names(grids)) {
  dat <- grids[[scenario]]
  summary_lines <- c(
    summary_lines, "",
    paste0(scenario, ": fictional IDs = ", length(unique(dat$ID)),
           "; rows = ", nrow(dat),
           "; dosing records = ", sum(dat$MDV == 1),
           "; simulation observation records = ", sum(dat$MDV == 0))
  )
}
writeLines(summary_lines, file.path(out_dir, "SYNTHETIC_SIMULATION_MANIFEST.txt"))
cat(paste(summary_lines, collapse = "\n"), "\n", sep = "")
