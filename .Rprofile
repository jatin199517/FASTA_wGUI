# Project-local .Rprofile.
# Loaded automatically once at R session startup (when this project is opened),
# so the error handler below stays active for the WHOLE session -- including
# scripts run partially / line-by-line. No sourcing or extra script needed.

local({
  # Capture the project directory at startup so the log path is fixed even
  # after scripts call setwd() into a subfolder.
  proj_dir <- tryCatch(normalizePath(getwd()), error = function(e) getwd())
  log_file <- file.path(proj_dir, "R_console_errors.log")

  # Interactive sessions only: leave Rscript/headless/knit runs on R's default
  # error behaviour. Logs IN ADDITION to the normal console error (does not
  # suppress it).
  if (interactive()) {
    options(error = function() {
      msg   <- geterrmessage()
      calls <- sys.calls()
      tb <- if (length(calls)) {
        paste(sprintf("%2d: %s", seq_along(calls),
                      vapply(calls, function(cc) {
                        txt <- paste(deparse(cc), collapse = " ")
                        if (nchar(txt) > 500) paste0(substr(txt, 1, 500), " ...") else txt
                      }, character(1))),
              collapse = "\n")
      } else "(no call stack captured)"

      entry <- paste0(
        "\n==== ERROR  ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "  ====\n",
        "Working dir: ", getwd(), "\n",
        "Message:   ", msg,
        "Traceback:\n", tb, "\n",
        "=================================================\n"
      )
      # Append; never let logging raise a new error.
      try(cat(entry, file = log_file, append = TRUE), silent = TRUE)
    })

    try(cat(sprintf("\n##### SESSION START %s #####\n",
                    format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
            file = log_file, append = TRUE), silent = TRUE)

    message("R error logging active -> ", log_file)
  }
})
