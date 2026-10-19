# ============================================================
# Model II Ka敏感性分析：Ka = 0.5、1、2 h^-1
#
# 文件对应：
#   0.5_sim.tab + 0.5_sim1.3.tab -> Ka = 0.5 h^-1
#   sim_RS1.tab + 1.3/sim_RS1.3.tab -> Ka = 1 h^-1（原模型/参照）
#   2_sim.tab + 2_sim1.3.tab     -> Ka = 2 h^-1
#
# 主要比较：
#   1. 每1 kg体重层的稳态C2h中位数及相对Ka=1的变化
#   2. PTA(C2h >= 268 ng/mL)及其百分点变化
#   3. P(268 <= C2h <= 500 ng/mL)及其百分点变化
#   4. PTA >= 75%时的最低候选剂量及其变化
#
# 三次模拟的随机实现不同，因此进行独立Monte Carlo分布比较，
# 不进行逐ID配对相减。
# ============================================================

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

library(dplyr)
library(tidyr)
library(ggplot2)

folder_ka <- project_input_dir(
  .project_root, "ATX_KA_INPUT_DIR", "data/analysis_inputs/ka_sensitivity"
)
folder_rs1 <- project_input_dir(
  .project_root, "ATX_RS1_INPUT_DIR", "data/analysis_inputs/simulation/RS1"
)
output_dir <- project_output_dir(
  .project_root, "ATX_KA_OUTPUT_DIR", "results/simulation/Ka_sensitivity"
)

file_map <- tibble(
  KA = c(0.5, 1, 2),
  KA_CODE = c("KA05", "KA1", "KA2"),
  FILES = list(
    c(
      file.path(folder_ka, "0.5_sim.tab"),
      file.path(folder_ka, "0.5_sim1.3.tab")
    ),
    c(
      file.path(folder_rs1, "sim_RS1.tab"),
      file.path(folder_rs1, "1.3", "sim_RS1.3.tab")
    ),
    c(
      file.path(folder_ka, "2_sim.tab"),
      file.path(folder_ka, "2_sim1.3.tab")
    )
  )
)

TARGET <- 268
UPPER_TARGET <- 500
PTA_THRESHOLD <- 75
C2H_TAD <- 2
WT_MIN <- 17
WT_MAX <- 60
EXPECTED_SUBPROBLEMS <- 500L
EXPECTED_DOSES <- c(0.5, 0.8, 1.0, 1.1, 1.2, 1.3, 1.4)
CHUNK_SIZE <- 200000L

tab_cols <- c(
  "ID", "TIME", "TAD", "AMT", "DV", "IPRED",
  "CL", "V", "WT", "rs_grp", "DOSEKG", "MDV"
)

# ============================================================
# 1. 分块读取单个NONMEM输出，只保留C2h记录
# ============================================================
read_c2h <- function(tab_file, ka_value, ka_code) {
  cat("\n============================================================\n")
  cat("读取 Ka =", ka_value, ":", tab_file, "\n")

  if (!file.exists(tab_file)) {
    stop("文件不存在：", tab_file)
  }

  con <- file(tab_file, open = "r")
  chunks <- list()
  chunk_id <- 0L
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
        c("DV", "IPRED", "WT", "rs_grp", "DOSEKG")
      ]

      if (nrow(block) > 0L) {
        chunk_id <- chunk_id + 1L
        chunks[[chunk_id]] <- block
      }
    }

    if (total_lines %% 2000000L < CHUNK_SIZE) {
      cat("已扫描", format(total_lines, big.mark = ","), "行\n")
    }
  }

  close(con)

  if (length(chunks) == 0L) {
    stop("未在文件中找到TAD=2且MDV=0的记录：", tab_file)
  }

  dat <- bind_rows(chunks) %>%
    rename(C2h = DV, C2h_IPRED = IPRED) %>%
    mutate(
      KA = ka_value,
      KA_CODE = ka_code,
      Group = factor(
        rs_grp,
        levels = c(0, 1),
        labels = c("CC+CT", "TT")
      ),
      WT_BIN_LOWER = pmin(floor(WT + 1e-8), WT_MAX - 1),
      WT_BIN_UPPER = WT_BIN_LOWER + 1,
      WT_GROUP = if_else(
        WT_BIN_LOWER == WT_MAX - 1,
        paste0(WT_MAX - 1, "-", WT_MAX),
        paste0(WT_BIN_LOWER, "-<", WT_BIN_UPPER)
      ),
      WT_GROUP = factor(
        WT_GROUP,
        levels = c(
          paste0(WT_MIN:(WT_MAX - 2), "-<", (WT_MIN + 1):(WT_MAX - 1)),
          paste0(WT_MAX - 1, "-", WT_MAX)
        )
      )
    )

  rm(chunks)
  gc()

  cat("SUBPROBLEMS:", table_count,
      "（预期", EXPECTED_SUBPROBLEMS, "）\n")
  cat("C2h记录数:", format(nrow(dat), big.mark = ","), "\n")
  cat("体重范围:",
      round(min(dat$WT, na.rm = TRUE), 1), "-",
      round(max(dat$WT, na.rm = TRUE), 1), "kg\n")

  if (min(dat$WT, na.rm = TRUE) < WT_MIN - 1e-8 ||
      max(dat$WT, na.rm = TRUE) > WT_MAX + 1e-8) {
    warning("Ka=", ka_value, "包含", WT_MIN, "-", WT_MAX,
            " kg范围外的体重记录。")
  }

  if (max(dat$WT, na.rm = TRUE) < WT_MAX - 1e-8) {
    warning("Ka=", ka_value, "的最大体重小于", WT_MAX,
            " kg；请确认是否仍使用旧的17-55 kg模拟文件。")
  }

  if (table_count != EXPECTED_SUBPROBLEMS) {
    warning("Ka=", ka_value, "检测到", table_count,
            "个SUBPROBLEMS，预期", EXPECTED_SUBPROBLEMS, "。")
  }

  dat
}

# ============================================================
# 2. 每个Ka场景分别读取、汇总，然后释放原始大数据
# ============================================================
summary_1kg_list <- vector("list", nrow(file_map))
summary_01kg_list <- vector("list", nrow(file_map))
summary_3strata_list <- vector("list", nrow(file_map))

for (i in seq_len(nrow(file_map))) {
  scenario_files <- file_map$FILES[[i]]

  c2h_parts <- lapply(scenario_files, function(current_file) {
    read_c2h(
      tab_file = current_file,
      ka_value = file_map$KA[i],
      ka_code = file_map$KA_CODE[i]
    )
  })

  c2h <- bind_rows(c2h_parts)
  rm(c2h_parts)
  gc()

  observed_doses <- sort(unique(round(c2h$DOSEKG, 6)))
  if (!isTRUE(all.equal(observed_doses, EXPECTED_DOSES))) {
    stop(
      "Ka=", file_map$KA[i],
      "的候选剂量网格不完整。检测到：",
      paste(observed_doses, collapse = ", "),
      "；预期：", paste(EXPECTED_DOSES, collapse = ", ")
    )
  }

  invalid_cells <- c2h %>%
    count(WT, rs_grp, DOSEKG, name = "n_cell") %>%
    filter(n_cell != EXPECTED_SUBPROBLEMS)

  if (nrow(invalid_cells) > 0L) {
    stop(
      "Ka=", file_map$KA[i],
      "存在模拟次数不等于", EXPECTED_SUBPROBLEMS,
      "的体重×基因型×剂量单元；",
      "请检查补充文件是否缺失记录或与原文件存在剂量重叠。"
    )
  }

  summary_01kg_list[[i]] <- c2h %>%
    group_by(KA, KA_CODE, WT, rs_grp, Group, DOSEKG) %>%
    summarise(
      n = n(),
      median_C2h = median(C2h, na.rm = TRUE),
      median_IPRED = median(C2h_IPRED, na.rm = TRUE),
      PTA_268 = mean(C2h >= TARGET, na.rm = TRUE) * 100,
      PTA_268_500 = mean(
        C2h >= TARGET & C2h <= UPPER_TARGET,
        na.rm = TRUE
      ) * 100,
      .groups = "drop"
    )

  summary_1kg_list[[i]] <- c2h %>%
    group_by(
      KA, KA_CODE,
      WT_BIN_LOWER, WT_BIN_UPPER, WT_GROUP,
      rs_grp, Group, DOSEKG
    ) %>%
    summarise(
      n = n(),
      WT_min = min(WT),
      WT_max = max(WT),
      median_C2h = median(C2h, na.rm = TRUE),
      median_IPRED = median(C2h_IPRED, na.rm = TRUE),
      p5 = quantile(C2h, 0.05, na.rm = TRUE),
      p95 = quantile(C2h, 0.95, na.rm = TRUE),
      PTA_268 = mean(C2h >= TARGET, na.rm = TRUE) * 100,
      PTA_268_500 = mean(
        C2h >= TARGET & C2h <= UPPER_TARGET,
        na.rm = TRUE
      ) * 100,
      .groups = "drop"
    )

  # 与正文剂量分析保持一致的三档体重分层。
  summary_3strata_list[[i]] <- c2h %>%
    mutate(
      WT_STRATUM = case_when(
        WT >= 17 & WT < 30 ~ "17-<30",
        WT >= 30 & WT < 40 ~ "30-<40",
        WT >= 40 & WT <= 60 ~ "40-60",
        TRUE ~ NA_character_
      ),
      WT_STRATUM = factor(
        WT_STRATUM,
        levels = c("17-<30", "30-<40", "40-60")
      )
    ) %>%
    filter(!is.na(WT_STRATUM)) %>%
    group_by(KA, KA_CODE, WT_STRATUM, rs_grp, Group, DOSEKG) %>%
    summarise(
      n = n(),
      WT_min = min(WT),
      WT_max = max(WT),
      median_C2h = median(C2h, na.rm = TRUE),
      median_IPRED = median(C2h_IPRED, na.rm = TRUE),
      p5 = quantile(C2h, 0.05, na.rm = TRUE),
      p95 = quantile(C2h, 0.95, na.rm = TRUE),
      PTA_268 = mean(C2h >= TARGET, na.rm = TRUE) * 100,
      PTA_268_500 = mean(
        C2h >= TARGET & C2h <= UPPER_TARGET,
        na.rm = TRUE
      ) * 100,
      .groups = "drop"
    )

  rm(c2h)
  gc()
}

summary_01kg <- bind_rows(summary_01kg_list) %>%
  mutate(across(c(median_C2h, median_IPRED, PTA_268, PTA_268_500),
                ~ round(.x, 2)))

summary_1kg <- bind_rows(summary_1kg_list) %>%
  mutate(across(c(median_C2h, median_IPRED, p5, p95,
                  PTA_268, PTA_268_500), ~ round(.x, 2)))

summary_3strata <- bind_rows(summary_3strata_list) %>%
  mutate(across(c(median_C2h, median_IPRED, p5, p95,
                  PTA_268, PTA_268_500), ~ round(.x, 2)))

rm(summary_01kg_list, summary_1kg_list, summary_3strata_list)
gc()

write.csv(summary_01kg, file.path(output_dir, "Ka_summary_C2h_0.1kg.csv"), row.names = FALSE)
write.csv(summary_1kg, file.path(output_dir, "Ka_summary_C2h_1kg.csv"), row.names = FALSE)
write.csv(summary_3strata,
          file.path(output_dir, "Ka_summary_C2h_3strata.csv"),
          row.names = FALSE)

# 检查Ka=0.5和Ka=2是否与参照Ka=1使用完全相同的
# 体重×基因型×剂量模拟网格；若不一致，后续比较会产生NA。
grid_keys <- c("WT_BIN_LOWER", "WT_BIN_UPPER", "WT_GROUP",
               "rs_grp", "Group", "DOSEKG")

reference_grid <- summary_1kg %>%
  filter(KA_CODE == "KA1") %>%
  select(all_of(grid_keys)) %>%
  distinct()

grid_audit <- bind_rows(lapply(c("KA05", "KA2"), function(code) {
  current_grid <- summary_1kg %>%
    filter(KA_CODE == code) %>%
    select(all_of(grid_keys)) %>%
    distinct()

  tibble(
    KA_CODE = code,
    missing_strata_vs_KA1 = nrow(
      anti_join(reference_grid, current_grid, by = grid_keys)
    ),
    extra_strata_vs_KA1 = nrow(
      anti_join(current_grid, reference_grid, by = grid_keys)
    )
  )
}))

write.csv(grid_audit, file.path(output_dir, "Ka_simulation_grid_audit.csv"), row.names = FALSE)
print(grid_audit)

if (any(grid_audit$missing_strata_vs_KA1 > 0 |
        grid_audit$extra_strata_vs_KA1 > 0)) {
  warning(
    "三种Ka场景的体重×基因型×剂量网格不一致；",
    "请确认三个Ka场景均包含17-60 kg及全部七个候选剂量后再解释比较结果。"
  )
}

# ============================================================
# 3. 以Ka=1为参照，比较C2h中位数和PTA
# ============================================================
comparison <- summary_1kg %>%
  select(
    WT_BIN_LOWER, WT_BIN_UPPER, WT_GROUP,
    rs_grp, Group, DOSEKG, KA_CODE,
    median_C2h, median_IPRED, PTA_268, PTA_268_500
  ) %>%
  pivot_wider(
    names_from = KA_CODE,
    values_from = c(median_C2h, median_IPRED, PTA_268, PTA_268_500),
    names_sep = "_"
  ) %>%
  mutate(
    delta_median_KA05 = median_C2h_KA05 - median_C2h_KA1,
    pct_median_KA05 = 100 * delta_median_KA05 / median_C2h_KA1,
    delta_median_KA2 = median_C2h_KA2 - median_C2h_KA1,
    pct_median_KA2 = 100 * delta_median_KA2 / median_C2h_KA1,

    delta_IPRED_KA05 = median_IPRED_KA05 - median_IPRED_KA1,
    pct_IPRED_KA05 = 100 * delta_IPRED_KA05 / median_IPRED_KA1,
    delta_IPRED_KA2 = median_IPRED_KA2 - median_IPRED_KA1,
    pct_IPRED_KA2 = 100 * delta_IPRED_KA2 / median_IPRED_KA1,

    delta_PTA268_KA05_pp = PTA_268_KA05 - PTA_268_KA1,
    delta_PTA268_KA2_pp = PTA_268_KA2 - PTA_268_KA1,
    delta_PTAwindow_KA05_pp = PTA_268_500_KA05 - PTA_268_500_KA1,
    delta_PTAwindow_KA2_pp = PTA_268_500_KA2 - PTA_268_500_KA1
  ) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2))) %>%
  arrange(rs_grp, WT_BIN_LOWER, DOSEKG)

write.csv(comparison, file.path(output_dir, "Ka_comparison_vs_Ka1_1kg.csv"), row.names = FALSE)

# 同时生成用于正文三档体重分层图片的Ka比较结果。
comparison_3strata <- summary_3strata %>%
  select(
    WT_STRATUM, rs_grp, Group, DOSEKG, KA_CODE,
    median_C2h, median_IPRED, PTA_268, PTA_268_500
  ) %>%
  pivot_wider(
    names_from = KA_CODE,
    values_from = c(median_C2h, median_IPRED, PTA_268, PTA_268_500),
    names_sep = "_"
  ) %>%
  mutate(
    delta_median_KA05 = median_C2h_KA05 - median_C2h_KA1,
    pct_median_KA05 = 100 * delta_median_KA05 / median_C2h_KA1,
    delta_median_KA2 = median_C2h_KA2 - median_C2h_KA1,
    pct_median_KA2 = 100 * delta_median_KA2 / median_C2h_KA1,

    delta_IPRED_KA05 = median_IPRED_KA05 - median_IPRED_KA1,
    pct_IPRED_KA05 = 100 * delta_IPRED_KA05 / median_IPRED_KA1,
    delta_IPRED_KA2 = median_IPRED_KA2 - median_IPRED_KA1,
    pct_IPRED_KA2 = 100 * delta_IPRED_KA2 / median_IPRED_KA1,

    delta_PTA268_KA05_pp = PTA_268_KA05 - PTA_268_KA1,
    delta_PTA268_KA2_pp = PTA_268_KA2 - PTA_268_KA1,
    delta_PTAwindow_KA05_pp = PTA_268_500_KA05 - PTA_268_500_KA1,
    delta_PTAwindow_KA2_pp = PTA_268_500_KA2 - PTA_268_500_KA1
  ) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2))) %>%
  arrange(rs_grp, WT_STRATUM, DOSEKG)

write.csv(
  comparison_3strata,
  file.path(output_dir, "Ka_comparison_vs_Ka1_3strata.csv"),
  row.names = FALSE
)

# ============================================================
# 4. 推荐剂量及其相对Ka=1的变化
# ============================================================
# 推荐规则：在候选剂量中选择PTA(C2h>=268)首次达到75%的最低剂量。
# 若1.4 mg/kg/day仍不达标，则保留为NA，不外推超说明书剂量。

recommend_observed <- summary_1kg %>%
  filter(PTA_268 >= PTA_THRESHOLD) %>%
  group_by(
    KA, KA_CODE,
    WT_BIN_LOWER, WT_BIN_UPPER, WT_GROUP,
    rs_grp, Group
  ) %>%
  arrange(DOSEKG, .by_group = TRUE) %>%
  slice(1L) %>%
  ungroup() %>%
  transmute(
    KA, KA_CODE,
    WT_BIN_LOWER, WT_BIN_UPPER, WT_GROUP,
    rs_grp, Group,
    rec_dose_mgkg = DOSEKG,
    PTA_at_rec = PTA_268,
    window_PTA_at_rec = PTA_268_500
  )

recommend_grid <- summary_1kg %>%
  distinct(WT_BIN_LOWER, WT_BIN_UPPER, WT_GROUP, rs_grp, Group) %>%
  crossing(file_map %>% select(KA, KA_CODE))

recommend_all <- recommend_grid %>%
  left_join(
    recommend_observed,
    by = c(
      "KA", "KA_CODE",
      "WT_BIN_LOWER", "WT_BIN_UPPER", "WT_GROUP",
      "rs_grp", "Group"
    )
  )

dose_comparison <- recommend_all %>%
  select(
    WT_BIN_LOWER, WT_BIN_UPPER, WT_GROUP,
    rs_grp, Group, KA_CODE,
    rec_dose_mgkg, PTA_at_rec, window_PTA_at_rec
  ) %>%
  pivot_wider(
    names_from = KA_CODE,
    values_from = c(rec_dose_mgkg, PTA_at_rec, window_PTA_at_rec),
    names_sep = "_"
  ) %>%
  mutate(
    dose_change_KA05_vs_KA1 = rec_dose_mgkg_KA05 - rec_dose_mgkg_KA1,
    dose_change_KA2_vs_KA1 = rec_dose_mgkg_KA2 - rec_dose_mgkg_KA1,
    status_KA05_vs_KA1 = case_when(
      is.na(rec_dose_mgkg_KA1) & is.na(rec_dose_mgkg_KA05) ~
        "Neither attained",
      is.na(rec_dose_mgkg_KA1) & !is.na(rec_dose_mgkg_KA05) ~
        "Attained only at Ka=0.5",
      !is.na(rec_dose_mgkg_KA1) & is.na(rec_dose_mgkg_KA05) ~
        "Not attained at Ka=0.5",
      TRUE ~ "Both attained"
    ),
    status_KA2_vs_KA1 = case_when(
      is.na(rec_dose_mgkg_KA1) & is.na(rec_dose_mgkg_KA2) ~
        "Neither attained",
      is.na(rec_dose_mgkg_KA1) & !is.na(rec_dose_mgkg_KA2) ~
        "Attained only at Ka=2",
      !is.na(rec_dose_mgkg_KA1) & is.na(rec_dose_mgkg_KA2) ~
        "Not attained at Ka=2",
      TRUE ~ "Both attained"
    )
  ) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2))) %>%
  arrange(rs_grp, WT_BIN_LOWER)

write.csv(recommend_all, file.path(output_dir, "Ka_recommended_dose_long_1kg.csv"), row.names = FALSE)
write.csv(dose_comparison,
          file.path(output_dir, "Ka_recommended_dose_comparison_vs_Ka1.csv"),
          row.names = FALSE)

# ============================================================
# 4B. 与正文一致的三档体重分层剂量比较
# ============================================================
recommend_3strata_observed <- summary_3strata %>%
  filter(PTA_268 >= PTA_THRESHOLD) %>%
  group_by(KA, KA_CODE, WT_STRATUM, rs_grp, Group) %>%
  arrange(DOSEKG, .by_group = TRUE) %>%
  slice(1L) %>%
  ungroup() %>%
  transmute(
    KA, KA_CODE, WT_STRATUM, rs_grp, Group,
    rec_dose_mgkg = DOSEKG,
    PTA_at_rec = PTA_268,
    window_PTA_at_rec = PTA_268_500
  )

recommend_3strata_grid <- summary_3strata %>%
  distinct(WT_STRATUM, rs_grp, Group) %>%
  crossing(file_map %>% select(KA, KA_CODE))

recommend_3strata_all <- recommend_3strata_grid %>%
  left_join(
    recommend_3strata_observed,
    by = c("KA", "KA_CODE", "WT_STRATUM", "rs_grp", "Group")
  )

dose_comparison_3strata <- recommend_3strata_all %>%
  select(
    WT_STRATUM, rs_grp, Group, KA_CODE,
    rec_dose_mgkg, PTA_at_rec, window_PTA_at_rec
  ) %>%
  pivot_wider(
    names_from = KA_CODE,
    values_from = c(rec_dose_mgkg, PTA_at_rec, window_PTA_at_rec),
    names_sep = "_"
  ) %>%
  mutate(
    dose_change_KA05_vs_KA1 =
      rec_dose_mgkg_KA05 - rec_dose_mgkg_KA1,
    dose_change_KA2_vs_KA1 =
      rec_dose_mgkg_KA2 - rec_dose_mgkg_KA1,
    status_KA05_vs_KA1 = case_when(
      is.na(rec_dose_mgkg_KA1) & is.na(rec_dose_mgkg_KA05) ~
        "Neither attained",
      is.na(rec_dose_mgkg_KA1) & !is.na(rec_dose_mgkg_KA05) ~
        "Attained only at Ka=0.5",
      !is.na(rec_dose_mgkg_KA1) & is.na(rec_dose_mgkg_KA05) ~
        "Not attained at Ka=0.5",
      TRUE ~ "Both attained"
    ),
    status_KA2_vs_KA1 = case_when(
      is.na(rec_dose_mgkg_KA1) & is.na(rec_dose_mgkg_KA2) ~
        "Neither attained",
      is.na(rec_dose_mgkg_KA1) & !is.na(rec_dose_mgkg_KA2) ~
        "Attained only at Ka=2",
      !is.na(rec_dose_mgkg_KA1) & is.na(rec_dose_mgkg_KA2) ~
        "Not attained at Ka=2",
      TRUE ~ "Both attained"
    )
  ) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2))) %>%
  arrange(rs_grp, WT_STRATUM)

write.csv(
  recommend_3strata_all,
  file.path(output_dir, "Ka_recommended_dose_long_3strata.csv"),
  row.names = FALSE
)
write.csv(
  dose_comparison_3strata,
  file.path(output_dir, "Ka_recommended_dose_comparison_3strata_vs_Ka1.csv"),
  row.names = FALSE
)

# ============================================================
# 5. C2h中位数百分比变化热图
# ============================================================
median_change_plot_data <- bind_rows(
  comparison_3strata %>%
    transmute(
      WT_STRATUM, Group, DOSEKG,
      Scenario = "Ka = 0.5 vs 1",
      Percent_change = pct_median_KA05
    ),
  comparison_3strata %>%
    transmute(
      WT_STRATUM, Group, DOSEKG,
      Scenario = "Ka = 2 vs 1",
      Percent_change = pct_median_KA2
    )
) %>%
  mutate(
    Scenario = factor(
      Scenario,
      levels = c("Ka = 0.5 vs 1", "Ka = 2 vs 1")
    )
  )

change_limit <- max(abs(median_change_plot_data$Percent_change), na.rm = TRUE)
change_limit <- max(5, ceiling(change_limit / 5) * 5)

p_median_change <- ggplot(
  median_change_plot_data,
  aes(x = factor(DOSEKG), y = WT_STRATUM, fill = Percent_change)
) +
  geom_tile(color = "white", linewidth = 0.15) +
  geom_text(aes(label = sprintf("%.1f", Percent_change)), size = 3.0) +
  scale_fill_gradient2(
    low = "#2166ac", mid = "white", high = "#b2182b",
    midpoint = 0,
    limits = c(-change_limit, change_limit),
    name = "Change (%)"
  ) +
  scale_y_discrete(
    limits = c("17-<30", "30-<40", "40-60"),
    labels = c(
      "17-<30" = "17\u2013<30",
      "30-<40" = "30\u2013<40",
      "40-60" = "40\u201360"
    ),
    drop = FALSE
  ) +
  facet_grid(Scenario ~ Group) +
  labs(
    x = "Dose (mg/kg/day)",
    y = "Body-weight stratum (kg)",
    title = "Percent Change in Median C2h Relative to Ka = 1"
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 7),
    axis.title = element_text(size = 11),
    plot.title = element_text(size = 13, hjust = 0.5),
    panel.grid = element_blank()
  )

ggsave(file.path(output_dir, "Ka_median_C2h_percent_change_heatmap.png"),
       plot = p_median_change,
       width = 12, height = 6.5, dpi = 300)

# ============================================================
# 6. 三个Ka场景的推荐剂量热图
# ============================================================
dose_plot_data <- recommend_3strata_all %>%
  mutate(
    KA_label = factor(
      KA_CODE,
      levels = c("KA05", "KA1", "KA2"),
      labels = c("Ka = 0.5", "Ka = 1", "Ka = 2")
    ),
    dose_label = if_else(
      is.na(rec_dose_mgkg),
      "NA",
      sprintf("%.1f", rec_dose_mgkg)
    )
  )

p_dose <- ggplot(
  dose_plot_data,
  aes(x = KA_label, y = WT_STRATUM, fill = rec_dose_mgkg)
) +
  geom_tile(color = "white", linewidth = 0.2) +
  geom_text(aes(label = dose_label), size = 2.4) +
  scale_fill_gradientn(
    colours = c("#f7fbff", "#9ecae1", "#4292c6", "#084594"),
    limits = c(0.5, 1.4),
    na.value = "#d9d9d9",
    name = "Dose\n(mg/kg/day)"
  ) +
  scale_y_discrete(
    limits = c("17-<30", "30-<40", "40-60"),
    labels = c(
      "17-<30" = "17\u2013<30",
      "30-<40" = "30\u2013<40",
      "40-60" = "40\u201360"
    ),
    drop = FALSE
  ) +
  facet_wrap(~ Group, ncol = 2) +
  labs(
    x = NULL,
    y = "Body-weight stratum (kg)",
    title = paste0(
      "Lowest Dose Achieving PTA ≥ ", PTA_THRESHOLD,
      "% Under Alternative Ka Values"
    )
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 11, face = "bold"),
    axis.text.x = element_text(size = 9),
    axis.text.y = element_text(size = 7),
    axis.title = element_text(size = 11),
    plot.title = element_text(size = 13, hjust = 0.5),
    panel.grid = element_blank()
  )

ggsave(file.path(output_dir, "Ka_recommended_dose_heatmap.png"),
       plot = p_dose,
       width = 9, height = 5.5, dpi = 300)

# ============================================================
# 7. 打印敏感性分析摘要
# ============================================================
summary_report <- tibble(
  Comparison = c("Ka=0.5 vs 1", "Ka=2 vs 1"),
  Median_abs_percent_change_DV = c(
    median(abs(comparison$pct_median_KA05), na.rm = TRUE),
    median(abs(comparison$pct_median_KA2), na.rm = TRUE)
  ),
  Maximum_abs_percent_change_DV = c(
    max(abs(comparison$pct_median_KA05), na.rm = TRUE),
    max(abs(comparison$pct_median_KA2), na.rm = TRUE)
  ),
  Median_abs_PTA_change_pp = c(
    median(abs(comparison$delta_PTA268_KA05_pp), na.rm = TRUE),
    median(abs(comparison$delta_PTA268_KA2_pp), na.rm = TRUE)
  ),
  Maximum_abs_PTA_change_pp = c(
    max(abs(comparison$delta_PTA268_KA05_pp), na.rm = TRUE),
    max(abs(comparison$delta_PTA268_KA2_pp), na.rm = TRUE)
  ),
  Maximum_abs_selected_dose_change = c(
    max(abs(dose_comparison$dose_change_KA05_vs_KA1), na.rm = TRUE),
    max(abs(dose_comparison$dose_change_KA2_vs_KA1), na.rm = TRUE)
  ),
  Weight_RS_strata_with_dose_change = c(
    sum(!is.na(dose_comparison$dose_change_KA05_vs_KA1) &
          dose_comparison$dose_change_KA05_vs_KA1 != 0),
    sum(!is.na(dose_comparison$dose_change_KA2_vs_KA1) &
          dose_comparison$dose_change_KA2_vs_KA1 != 0)
  ),
  Weight_RS_strata_with_attainability_change = c(
    sum(!dose_comparison$status_KA05_vs_KA1 %in%
          c("Both attained", "Neither attained")),
    sum(!dose_comparison$status_KA2_vs_KA1 %in%
          c("Both attained", "Neither attained"))
  ),
  Main_weight_strata_with_dose_change = c(
    sum(!is.na(dose_comparison_3strata$dose_change_KA05_vs_KA1) &
          dose_comparison_3strata$dose_change_KA05_vs_KA1 != 0),
    sum(!is.na(dose_comparison_3strata$dose_change_KA2_vs_KA1) &
          dose_comparison_3strata$dose_change_KA2_vs_KA1 != 0)
  ),
  Main_weight_strata_with_attainability_change = c(
    sum(!dose_comparison_3strata$status_KA05_vs_KA1 %in%
          c("Both attained", "Neither attained")),
    sum(!dose_comparison_3strata$status_KA2_vs_KA1 %in%
          c("Both attained", "Neither attained"))
  )
) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2)))

write.csv(summary_report, file.path(output_dir, "Ka_sensitivity_summary.csv"), row.names = FALSE)

cat("\n=== Ka敏感性分析摘要 ===\n")
print(as.data.frame(summary_report))

cat("\n=== 正文三档体重分层的最低达标剂量比较 ===\n")
print(as.data.frame(dose_comparison_3strata))

cat("\n=== 完成 ===\n")
cat("输出目录：", output_dir, "\n")
cat("输出文件：\n")
cat("  1. Ka_summary_C2h_0.1kg.csv\n")
cat("  2. Ka_summary_C2h_1kg.csv\n")
cat("  3. Ka_comparison_vs_Ka1_1kg.csv\n")
cat("  4. Ka_recommended_dose_long_1kg.csv\n")
cat("  5. Ka_recommended_dose_comparison_vs_Ka1.csv\n")
cat("  6. Ka_sensitivity_summary.csv\n")
cat("  7. Ka_median_C2h_percent_change_heatmap.png\n")
cat("  8. Ka_recommended_dose_heatmap.png\n")
cat("  9. Ka_simulation_grid_audit.csv\n")
cat(" 10. Ka_summary_C2h_3strata.csv\n")
cat(" 11. Ka_recommended_dose_long_3strata.csv\n")
cat(" 12. Ka_recommended_dose_comparison_3strata_vs_Ka1.csv\n")
cat(" 13. Ka_comparison_vs_Ka1_3strata.csv\n")
