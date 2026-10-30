# Run with Rscript run_all.R; no Python process is called.
# See README.md for numerical conventions and scope.

args = commandArgs(FALSE)

f = sub("^--file=", "", args[grepl("^--file=", args)])

if (length(f)) setwd(dirname(normalizePath(f)))

if (!requireNamespace("jsonlite", quietly = TRUE)) stop("Install jsonlite: install.packages(\"jsonlite\")")

options(digits = 16, warn = 1)

dir.create("results", showWarnings = FALSE)

dir.create("validation", showWarnings = FALSE)

manifest <- read.csv('validation/input_manifest.csv', stringsAsFactors=FALSE)
stopifnot(all(file.exists(manifest$path)))
stopifnot(identical(unname(tools::md5sum(manifest$path)), manifest$md5))
source("R/core.R")

source("R/reconstruct.R")

source("R/vertical.R")

source("R/selection.R")

source("R/movement.R")

source("R/tables.R")

source("R/validate.R")

source("R/plots.R")

writeLines(capture.output(sessionInfo()), "validation/sessionInfo.txt")

