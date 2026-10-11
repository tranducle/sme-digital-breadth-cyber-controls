# Install the R packages used by the analysis scripts.
pkgs <- c("survey", "digest", "jsonlite")
missing <- setdiff(pkgs, rownames(installed.packages()))
if (length(missing) > 0) install.packages(missing, repos = "https://cloud.r-project.org")
cat("R packages available:", paste(pkgs, collapse = ", "), "\n")
