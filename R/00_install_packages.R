# Install the R packages used by 01_build_dgg.R and 02_validate_dgg.R.
# Reproducible alternative: renv::init() in the repository root, then renv::snapshot()
# to write renv.lock (commit it), and renv::restore() on another machine.
pkgs <- c("dggridR", "sf", "data.table")
missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing)
cat("R packages ready:", paste(pkgs, collapse = ", "), "\n")
# Tested with: GEOS 3.13.0, GDAL 3.8.5, PROJ 9.5.1 (as reported by sf on load).
