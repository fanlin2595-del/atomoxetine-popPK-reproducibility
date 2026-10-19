# Static checks and fictional event-history checks. NONMEM is NOT invoked.
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Run using Rscript.")
root <- dirname(dirname(normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)))
models <- c("sim_model_I_BW", "sim_model_II_rs1065852", "sim_model_III_AS")
files <- c("synthetic_sim_BW.dat", "synthetic_sim_RS.dat", "synthetic_sim_AS.dat")
fields <- c("AS", "rs_grp", "AS")
ns <- c(6L, 12L, 18L)
expected_theta <- list(c(14.1, 69.4, 1), c(17.1, 67.6, 1, -0.645), c(15.6, 67.8, 1, 0.585))
expected_omega <- list(c(0.19, 0.359), c(0.156, 0.329), c(0.168, 0.31))
expected_sigma <- c(0.0718, 0.0712, 0.0712)
active_lines <- function(lines) {
  out <- gsub("[[:space:]]", "", toupper(sub(";.*$", "", lines)))
  out[nzchar(out)]
}
block <- function(lines, tag) {
  start <- which(lines == paste0("$", tag))
  stopifnot(length(start) == 1L)
  nexts <- which(seq_along(lines) > start & startsWith(lines, "$"))
  end <- if (length(nexts)) nexts[1] - 1L else length(lines)
  lines[seq.int(start + 1L, end)]
}
fixed_values <- function(lines, tag) {
  x <- block(lines, tag)
  stopifnot(all(grepl("FIX", x, fixed = TRUE)))
  as.numeric(gsub("[()]|FIX", "", x))
}
for (j in seq_along(models)) {
  ctl <- readLines(file.path(root, "simulation", paste0(models[j], ".mod")), warn = FALSE)
  stopifnot(all(utf8ToInt(paste(ctl, collapse = "\n")) < 128L))
  a <- active_lines(ctl)
  schema <- c("ID", "TIME", "AMT", "ADDL", "II", "DV", "MDV", "EVID",
              "WT", fields[j], "DOSEKG", "TAD")
  expected_input <- paste0("$INPUT", paste(toupper(schema), collapse = ""))
  stopifnot(expected_input %in% a,
            paste0("$DATA", toupper(files[j]), "IGNORE=@") %in% a,
            "$SIMULATION(20260717)ONLYSIMNSUB=500" %in% a,
            "$SUBROUTINEADVAN2TRANS2" %in% a,
            "IPRED=F" %in% a, "Y=F*(1+EPS(1))" %in% a,
            !any(grepl("^\\$(EST|COV)", a)),
            identical(fixed_values(a, "THETA"), expected_theta[[j]]),
            identical(fixed_values(a, "OMEGA"), expected_omega[[j]]),
            identical(fixed_values(a, "SIGMA"), expected_sigma[j]))
  expected_table <- paste0("$TABLEIDTIMETADAMTDVIPREDCLVWT", toupper(fields[j]), "DOSEKGMDV")
  stopifnot(expected_table %in% a,
            paste0("NOPRINTONEHEADERNOAPPENDFILE=", toupper(models[j]), ".TAB") %in% a)
  if (j == 1L) stopifnot(!any(grepl("AS", block(a, "PK"), fixed = TRUE)))
  if (j == 2L) stopifnot("IF(RS_GRP.EQ.1)THEN" %in% a)
  if (j == 3L) stopifnot("TVCL=TVCL*(AS/1.25)**THETA(4)" %in% a)

  path <- file.path(root, "data", "simulation", files[j])
  if (!file.exists(path)) stop("Missing synthetic grid. Run prepare_simulation_demo.R first.")
  stopifnot(identical(readLines(path, n = 1L), paste0("@", paste(schema, collapse = " "))))
  d <- read.table(path, skip = 1L, col.names = schema, check.names = FALSE)
  stopifnot(ncol(d) == 12L, nrow(d) == 5L * ns[j], length(unique(d$ID)) == ns[j],
            all(is.finite(as.matrix(d))), all(d$DV == 0),
            all(d$ID %in% (910000L + 1000L * j + seq_len(ns[j]))))
  if (fields[j] == "rs_grp") stopifnot(all(d$rs_grp %in% 0:1))
  if (fields[j] == "AS") stopifnot(all(d$AS > 0))
  for (subject in split(d, d$ID)) {
    dose <- subject[subject$EVID == 1, , drop = FALSE]
    obs <- subject[subject$MDV == 0, , drop = FALSE]
    stopifnot(nrow(dose) == 1L, nrow(obs) == 4L, all(diff(subject$TIME) > 0),
              dose$TIME == 0, dose$ADDL == 29, dose$II == 24, dose$MDV == 1,
              abs(dose$AMT - dose$WT * dose$DOSEKG) < 1e-8,
              identical(obs$TIME, c(2L, 26L, 674L, 698L)) ||
                isTRUE(all.equal(obs$TIME, c(2, 26, 674, 698))),
              all(obs$EVID == 0), all(obs$AMT == 0), all(obs$ADDL == 0), all(obs$II == 0))
    history <- dose$TIME + seq.int(0, dose$ADDL) * dose$II
    reconstructed <- vapply(obs$TIME, function(t) t - max(history[history <= t]), numeric(1))
    stopifnot(all(reconstructed == obs$TAD), all(obs$TAD == 2))
  }
  run <- file.path(root, "simulation_demo_runs", models[j])
  if (dir.exists(file.path(root, "simulation_demo_runs"))) {
    stopifnot(identical(readLines(file.path(run, "model.mod")), ctl),
              identical(unname(tools::md5sum(path)),
                        unname(tools::md5sum(file.path(run, files[j])))))
  }
}
cat("PASS: 3 simulation controls, fixed parameters, unchanged NSUB/seed, 12-column inputs,\n",
    "fictional event chronology, dose amounts, TAD and isolated run copies.\n",
    "NOT TESTED: NONMEM execution, original simulation grid, steady-state selection, PTA or figures.\n",
    sep = "")
