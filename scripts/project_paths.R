# Shared path helpers for portable, project-relative analysis scripts.

find_project_root <- function(script_file) {
  if (missing(script_file) || is.null(script_file) || !nzchar(script_file)) {
    stop("Provide the path to the running R script.")
  }

  candidate <- normalizePath(dirname(script_file), winslash = "/", mustWork = FALSE)
  repeat {
    is_root <- file.exists(file.path(candidate, "README.md")) &&
      dir.exists(file.path(candidate, "models")) &&
      dir.exists(file.path(candidate, "analysis"))
    if (is_root) return(candidate)

    parent <- dirname(candidate)
    if (identical(parent, candidate)) break
    candidate <- parent
  }

  stop("Could not locate the project root from script: ", script_file)
}

project_input_dir <- function(project_root, env_var, relative_path) {
  configured <- Sys.getenv(env_var, unset = "")
  if (!nzchar(configured)) return(file.path(project_root, relative_path))

  configured <- path.expand(configured)
  is_absolute <- grepl("^(?:/|[A-Za-z]:[/\\\\]|\\\\\\\\)", configured, perl = TRUE)
  if (is_absolute) configured else file.path(project_root, configured)
}

project_output_dir <- function(project_root, env_var, relative_path) {
  configured <- Sys.getenv(env_var, unset = "")
  output_dir <- if (!nzchar(configured)) {
    file.path(project_root, relative_path)
  } else {
    configured <- path.expand(configured)
    is_absolute <- grepl("^(?:/|[A-Za-z]:[/\\\\]|\\\\\\\\)", configured, perl = TRUE)
    if (is_absolute) configured else file.path(project_root, configured)
  }

  if (!dir.exists(output_dir) && !dir.create(output_dir, recursive = TRUE)) {
    stop("Could not create output directory: ", output_dir)
  }
  normalizePath(output_dir, winslash = "/", mustWork = TRUE)
}
