#' Copy the bundled Quarto analysis report template
#'
#' Copies `analysis-report`, a Quarto document with lab-style defaults for
#' HTML and PDF output, into a directory of the current project.
#'
#' @param path Directory into which the template files are copied.
#' @param overwrite Replace existing files.
#'
#' @return The paths of the copied files, invisibly.
#'
#' @export
copy_report_template <- function(path = ".", overwrite = FALSE) {
  source <- system.file("quarto", "analysis-report", package = "clinstats")
  if (!nzchar(source)) {
    cli::cli_abort("The bundled template could not be located.")
  }
  files <- list.files(source, recursive = TRUE, full.names = TRUE)
  target <- file.path(path, "analysis-report", sub(
    paste0("^", source, "/?"), "", files
  ))
  dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
  copied <- file.copy(files, target, overwrite = overwrite)
  if (!all(copied) && !overwrite) {
    cli::cli_abort(c(
      "Some template files already exist.",
      i = "Use {.code overwrite = TRUE} to replace them."
    ))
  }
  invisible(target)
}
