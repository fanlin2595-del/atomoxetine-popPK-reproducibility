# Package v2: static and synthetic-event checks. Does NOT invoke NONMEM.
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Run using Rscript.")
root <- dirname(dirname(normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)))
paths <- file.path(root, "models", c("model_I_BW.mod", "model_II_rs1065852.mod",
                                    "model_III_AS.mod", "model_IV_CYP2D6.mod"))
expected_input <- "$INPUT ID DAT2=DROP TIME AMT ADDL II DV MDV WT SEX AGE rs1065852 rs_grp rs5030865 rs16947 AS CYP2D6 TAD"
controls <- lapply(paths, readLines, warn = FALSE)
stopifnot(length(controls) == 4L)
for (ctl in controls) {
  stopifnot(all(utf8ToInt(paste(ctl, collapse = "\n")) < 128L))
  active <- trimws(gsub("[[:space:]]+", " ", sub(";.*$", "", ctl)))
  stopifnot(expected_input %in% active,
            "$DATA synthetic_demo.dat IGNORE=@" %in% active,
            "$SUBROUTINE ADVAN2 TRANS2" %in% active,
            any(grepl("^\\(1,[[:space:]]*FIX\\)$", active)),
            "Y=F*(1+EPS(1))" %in% active,
            "IWRES=IRES/(IPRED+DEL)" %in% active,
            "METHOD=1 INTERACTION" %in% active,
            !any(grepl("/Users/|/Volumes/|ATX0623", active)))
}
table_lines <- vapply(controls, function(x) grep("^\\$TABLE", x, value = TRUE), character(1))
outputs <- sub(".*FILE=", "", table_lines)
stopifnot(length(unique(outputs)) == 4L)

data_path <- file.path(root, "data", "synthetic_demo.dat")
if (!file.exists(data_path)) stop("Missing synthetic input. Run prepare_demo.R first.")
header <- readLines(data_path, n = 1L)
stopifnot(startsWith(header, "@ID "))
columns <- strsplit(sub("^@", "", header), " ", fixed = TRUE)[[1]]
dat <- read.table(data_path, skip = 1L, col.names = columns, check.names = FALSE)
stopifnot(ncol(dat) == 18L, all(is.finite(as.matrix(dat))),
          length(unique(dat$ID)) == 96L, all(dat$ID %in% 900001:900096),
          all(dat$DAT2 == 0), all(dat$AS > 0),
          all(dat$DV[dat$MDV == 0] >= 5), all(dat$TAD >= 0),
          all((dat$rs_grp == 1) == (dat$CYP2D6 == 2)))

for (subject in split(dat, dat$ID)) {
  stopifnot(all(diff(subject$TIME) >= 0))
  dosing <- subject[subject$AMT > 0, , drop = FALSE]
  stopifnot(nrow(dosing) == 2L, all(dosing$MDV == 1),
            all(dosing$II %in% c(12, 24)), all(dosing$ADDL >= 0),
            all(dosing$ADDL == floor(dosing$ADDL)))
  expanded <- unlist(lapply(seq_len(nrow(dosing)), function(i) {
    dosing$TIME[i] + seq.int(0, dosing$ADDL[i]) * dosing$II[i]
  }))
  stopifnot(!anyDuplicated(expanded), all(diff(expanded) > 0))
  obs <- subject[subject$MDV == 0, , drop = FALSE]
  stopifnot(all(obs$AMT == 0))
  reconstructed <- vapply(obs$TIME, function(t) t - max(expanded[expanded <= t]), numeric(1))
  stopifnot(isTRUE(all.equal(reconstructed, obs$TAD, check.attributes = FALSE)))
}

# Conditional structural equivalence, using the SAME explicit demonstration parameters.
# Does not compare independent fitted estimates.
cl_II <- 17.1 * (dat$WT / 30.7)^0.75 * exp(-0.645 * (dat$rs_grp == 1))
cl_IV <- 17.1 * (dat$WT / 30.7)^0.75 * exp(-0.645 * (dat$CYP2D6 == 2))
stopifnot(identical(cl_II, cl_IV))

normalize_active <- function(x) {
  x <- gsub("[[:space:]]", "", toupper(sub(";.*$", "", x)))
  x[nzchar(x)]
}

# Matching active blocks in latest Models II/IV, except their group selector.
record_block <- function(lines, tag) {
  active <- trimws(sub(";.*$", "", lines))
  start <- which(grepl(paste0("^\\$", tag, "($|[[:space:]])"), active))
  stopifnot(length(start) == 1L)
  next_record <- which(seq_along(active) > start & grepl("^\\$", active))
  end <- if (length(next_record)) next_record[1] - 1L else length(active)
  normalize_active(active[seq.int(start, end)])
}
for (tag in c("THETA", "OMEGA", "SIGMA", "ERROR", "ESTIMATION", "COV")) {
  stopifnot(identical(record_block(controls[[2]], tag), record_block(controls[[4]], tag)))
}
pk_II <- record_block(controls[[2]], "PK")
pk_IV <- record_block(controls[[4]], "PK")
pk_II <- gsub("RS_GRP.EQ.1", "GROUP.EQ.1", pk_II, fixed = TRUE)
pk_IV <- gsub("CYP2D6.EQ.2", "GROUP.EQ.1", pk_IV, fixed = TRUE)
stopifnot(identical(pk_II, pk_IV))
stopifnot(identical(record_block(controls[[3]], "THETA")[5], "(-3,0.585,3)"),
          identical(record_block(controls[[2]], "SIGMA")[2], "0.0712"))

run_root <- file.path(root, "demo_runs")
if (dir.exists(run_root)) {
  for (j in seq_along(paths)) {
    model_name <- sub("\\.mod$", "", basename(paths[j]))
    run_path <- file.path(run_root, model_name)
    copy <- readLines(file.path(run_path, "model.mod"), warn = FALSE)
    expected <- controls[[j]]
    stopifnot(identical(normalize_active(copy), normalize_active(expected)),
              identical(unname(tools::md5sum(file.path(run_path, "synthetic_demo.dat"))),
                        unname(tools::md5sum(data_path))))
  }
}
cat("PASS: four control streams, 18-column synthetic schema, dose histories, TAD,\n",
    "conditional II/IV structural equivalence and prepared demo copies.\n",
    "NOT TESTED: NM-TRAN compilation, NONMEM estimation/covariance, or original-data results.\n",
    sep = "")
