# Synthetic NONMEM demonstration, package v2; not historical study code.
# Reads only the supplied control streams. NEVER reads clinical data.
# Run: Rscript scripts/prepare_demo.R, from any working directory.

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Run this script using Rscript.")
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
root <- dirname(dirname(script_path))
run_root <- file.path(root, "demo_runs")
if (dir.exists(run_root)) {
  stop("demo_runs already exists. Use a fresh package or move old runs aside; no files were changed.")
}

model_names <- c("model_I_BW", "model_II_rs1065852", "model_III_AS", "model_IV_CYP2D6")
model_paths <- file.path(root, "models", paste0(model_names, ".mod"))
stopifnot(all(file.exists(model_paths)))

get_block <- function(lines, tag) {
  lines <- trimws(sub(";.*$", "", lines))
  start <- which(grepl(paste0("^\\$", tag, "($|[[:space:]])"), lines))
  stopifnot(length(start) == 1L)
  following <- which(seq_along(lines) > start & grepl("^\\$", lines))
  finish <- if (length(following)) following[1] - 1L else length(lines)
  body <- lines[seq.int(start + 1L, finish)]
  body[nzchar(body)]
}
get_theta <- function(lines) {
  vapply(get_block(lines, "THETA"), function(x) {
    fixed <- grepl("FIX", x, ignore.case = TRUE)
    parts <- strsplit(gsub("[()]", "", x), ",", fixed = TRUE)[[1]]
    numeric_parts <- suppressWarnings(as.numeric(trimws(parts)))
    numeric_parts <- numeric_parts[is.finite(numeric_parts)]
    if (fixed || length(numeric_parts) == 1L) numeric_parts[1] else numeric_parts[2]
  }, numeric(1))
}

source_II <- readLines(model_paths[2], warn = FALSE)
theta <- get_theta(source_II)
omega <- as.numeric(get_block(source_II, "OMEGA"))
sigma <- as.numeric(get_block(source_II, "SIGMA"))
stopifnot(length(theta) == 4L, length(omega) == 2L, length(sigma) == 1L,
          all(omega > 0), sigma > 0, theta[3] == 1)

# Oral one-compartment concentration by superposition, in ng/mL.
# Equal KA and KE are handled by the limiting expression.
concentration <- function(t, dose_times, amounts, cl, v, ka) {
  lag <- t - dose_times
  active <- lag >= 0
  lag <- lag[active]
  amounts <- amounts[active]
  ke <- cl / v
  delta <- ka - ke
  kernel <- if (abs(delta) < 1e-8) {
    lag * exp(-ka * lag)
  } else {
    (exp(-ke * lag) - exp(-ka * lag)) / delta
  }
  1000 * sum(amounts / v * ka * kernel)
}

RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
set.seed(4587316)
n_subjects <- 96L # Arbitrary demonstration size, not the study sample size.
all_records <- vector("list", n_subjects)
n_excluded <- 0L
expected_observations <- 9L * n_subjects
schema <- c("ID", "DAT2", "TIME", "AMT", "ADDL", "II", "DV", "MDV",
            "WT", "SEX", "AGE", "rs1065852", "rs_grp", "rs5030865",
            "rs16947", "AS", "CYP2D6", "TAD")

for (i in seq_len(n_subjects)) {
  id <- 900000L + i
  genotype <- (i - 1L) %% 3L
  rs_grp <- as.integer(genotype == 2L)
  phenotype <- if (rs_grp == 1L) 2L else 1L
  activity <- c(2, 1.25, 0.5)[genotype + 1L]
  regimen <- ((i - 1L) %/% 3L) %% 3L # 0 morning, 1 evening, 2 BID, demo only.
  interval <- if (regimen == 2L) 12 else 24
  offset <- if (regimen == 1L) 12 else 0
  weight <- round(runif(1, 17, 60), 1)
  age <- sample(6:18, 1)
  sex <- sample(1:2, 1)
  daily_dose <- round(weight * runif(1, 0.6, 1.2), 1)
  dose1 <- daily_dose * interval / 24
  dose2 <- 1.1 * dose1 # Arbitrary day-7 change, not a dosing recommendation.
  n_doses_epoch <- 168 / interval
  epoch1 <- offset + seq(0, 168 - interval, by = interval)
  epoch2 <- offset + seq(168, 336 - interval, by = interval)
  dose_times <- c(epoch1, epoch2)
  dose_amounts <- c(rep(dose1, length(epoch1)), rep(dose2, length(epoch2)))
  eta <- rnorm(2, sd = sqrt(omega))
  cl <- theta[1] * (weight / 30.7)^0.75 * exp(theta[4] * rs_grp) * exp(eta[1])
  v <- theta[2] * (weight / 30.7) * exp(eta[2])

  make_record <- function(time, amount, addl, ii, dv, mdv, tad) {
    setNames(c(id, 0, time, amount, addl, ii, dv, mdv, weight, sex, age,
               genotype, rs_grp, 0, 0, activity, phenotype, tad), schema)
  }
  rows <- list(
    make_record(offset, dose1, n_doses_epoch - 1, interval, 0, 1, 0),
    make_record(offset + 168, dose2, n_doses_epoch - 1, interval, 0, 1, 0)
  )
  for (time in offset + c(2, 8, 12.5, 26, 74, 170, 176, 314, 320)) {
    pred <- concentration(time, dose_times, dose_amounts, cl, v, theta[3])
    dv <- pred * (1 + rnorm(1, sd = sqrt(sigma)))
    # Demo M1 preprocessing. Do not clamp artificial BLQ values to the LLOQ.
    if (dv < 5) {
      n_excluded <- n_excluded + 1L
      next
    }
    last_dose <- max(dose_times[dose_times <= time])
    rows[[length(rows) + 1L]] <- make_record(time, 0, 0, 0, round(dv, 6), 0, time - last_dose)
  }
  subject <- as.data.frame(do.call(rbind, rows))
  all_records[[i]] <- subject[order(subject$TIME, -subject$AMT), , drop = FALSE]
}

dat <- do.call(rbind, all_records)
rownames(dat) <- NULL
stopifnot(identical(names(dat), schema), ncol(dat) == 18L,
          length(unique(dat$ID)) == n_subjects, all(is.finite(as.matrix(dat))),
          all(dat$DV[dat$MDV == 0] >= 5),
          all((dat$rs_grp == 1) == (dat$CYP2D6 == 2)),
          sum(dat$MDV == 0) + n_excluded == expected_observations)

dir.create(file.path(root, "data"), showWarnings = FALSE)
data_path <- file.path(root, "data", "synthetic_demo.dat")
data_con <- file(data_path, open = "wt")
writeLines(paste0("@", paste(schema, collapse = " ")), data_con)
write.table(dat, data_con, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ")
close(data_con)

manifest <- c(
  "SYNTHETIC DEMONSTRATION ONLY - not observed patient data.",
  "Generator reads control streams only; no clinical dataset is accessed.",
  "Seed: 4587316", paste("Fictional subjects:", n_subjects),
  paste("Explicit dose records:", sum(dat$MDV == 1)),
  paste("Artificial observations retained:", sum(dat$MDV == 0)),
  paste("Artificial observations excluded at DV < 5:", n_excluded),
  paste("Generating THETA starting values:", paste(theta, collapse = ", ")),
  paste("Generating OMEGA variances:", paste(omega, collapse = ", ")),
  paste("Generating SIGMA variance:", sigma),
  paste("Input MD5:", unname(tools::md5sum(data_path))),
  "All four demo copies preserve the active code of the prepared v2 models.",
  "NONMEM estimation was NOT performed by this script."
)
writeLines(manifest, file.path(root, "data", "SYNTHETIC_MANIFEST.txt"))

dir.create(run_root)
for (j in seq_along(model_names)) {
  target <- file.path(run_root, model_names[j])
  dir.create(target)
  ctl <- readLines(model_paths[j], warn = FALSE)
  header <- "; SYNTHETIC DEMONSTRATION RUN - NOT A REPRODUCTION OF STUDY RESULTS."
  writeLines(c(header, ctl), file.path(target, "model.mod"))
  stopifnot(file.copy(data_path, file.path(target, "synthetic_demo.dat"), overwrite = FALSE))
}
writeLines(capture.output(sessionInfo()), file.path(run_root, "R_SESSION_INFO.txt"))
cat(paste(manifest, collapse = "\n"), "\n")
cat("Prepared four isolated run directories. NONMEM execution is a separate step.\n")
