packages <- c(
  "data.table",
  "ggplot2",
  "ggpubr",
  "optparse",
  "qqman",
  "R.utils",
  "rjson",
  "stringi",
  "openxlsx"
)

to_install <- packages[!sapply(packages, requireNamespace, quietly = TRUE)]
if (length(to_install) > 0) {
  install.packages(to_install, dependencies = TRUE, repos = "http://cran.rstudio.com/")
}

failed <- packages[!sapply(packages, requireNamespace, quietly = TRUE)]
if (length(failed) > 0) stop(paste("Failed to install:", paste(failed, collapse = ", ")))
