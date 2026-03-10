# Install package dependencies from DESCRIPTION if any are missing.
# This ensures tests run without renv when deps are not yet installed.
.find_pkg_root <- function() {
  wd <- getwd()
  for (i in 1:10) {
    if (file.exists(file.path(wd, "DESCRIPTION"))) return(wd)
    wd <- dirname(wd)
  }
  getwd()
}

pkg_root <- .find_pkg_root()
desc_file <- file.path(pkg_root, "DESCRIPTION")
if (file.exists(desc_file)) {
  desc <- read.dcf(desc_file, fields = c("Imports", "Suggests"))
  fields <- c(desc[, "Imports"], desc[, "Suggests"])
  fields <- fields[!is.na(fields) & nzchar(fields)]
  if (length(fields) > 0) {
    raw <- trimws(strsplit(paste(fields, collapse = ","), ",")[[1]])
    pkgs <- unique(sub("^([^ (]+).*", "\\1", raw))  # strip version constraints
    pkgs <- pkgs[!pkgs %in% c("R", "")]  # drop R and empty
    missing <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
    if (length(missing) > 0) {
      message("Installing missing test dependencies: ", paste(missing, collapse = ", "))
      install.packages(missing, repos = "https://cloud.r-project.org", quiet = TRUE)
    }
  }
}
