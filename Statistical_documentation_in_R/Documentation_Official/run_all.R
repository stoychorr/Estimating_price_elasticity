# run_all.R
# Execute the entire pipeline from a clean R session and render paper.Rmd.

runner_directory <- function() {
  frame_files <- vapply(
    sys.frames(),
    function(frame) {
      if (is.null(frame$ofile)) NA_character_ else as.character(frame$ofile)
    },
    character(1)
  )
  frame_files <- frame_files[!is.na(frame_files)]
  if (length(frame_files) > 0L) {
    return(dirname(normalizePath(tail(frame_files, 1L), mustWork = TRUE)))
  }

  command_file <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(command_file) > 0L) {
    return(dirname(normalizePath(sub("^--file=", "", command_file[1L]), mustWork = TRUE)))
  }
  normalizePath(getwd(), mustWork = TRUE)
}

runner_dir <- runner_directory()
scripts <- sprintf("%02d_%s.R", 1:6, c(
  "import", "cleaning", "descriptives", "baseline", "robustness", "figures"
))

for (script in scripts) {
  message("Running ", script)
  source(file.path(runner_dir, script), chdir = TRUE)
}

if (!requireNamespace("rmarkdown", quietly = TRUE)) {
  stop("Install rmarkdown to render paper.Rmd.", call. = FALSE)
}

rmarkdown::render(
  input = file.path(runner_dir, "paper.Rmd"),
  output_file = "paper.html",
  output_dir = runner_dir,
  params = list(run_pipeline = FALSE),
  envir = new.env(parent = globalenv()),
  quiet = FALSE
)

message("Pipeline complete")
