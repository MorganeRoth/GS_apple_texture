cat("Loading packages...\n")

# Get the names of all installed packages
inst_pkgs <- installed.packages()[, "Package"]

# Install required packages from CRAN
req_pkgs_cran <- setdiff(pkgs_cran, inst_pkgs)
if (length(req_pkgs_cran) > 0) {
  install.packages(req_pkgs_cran, clean = TRUE)
}

# Install required packages from Bioconductor
req_pkgs_bioc <- setdiff(pkgs_bioc, inst_pkgs)
if (length(req_pkgs_bioc) > 0) {
  if (getRversion() < "3.5") {
    source("https://bioconductor.org/biocLite.R")
    BiocInstaller::biocLite(req_pkgs_bioc)
  } else {
    install.packages("BiocManager")
    BiocManager::install(req_pkgs_bioc, update = FALSE)
  }
}

# Load necessary packages
shh <- lapply(c(setdiff(pkgs_cran, "config"), pkgs_bioc), library,
              character.only = TRUE)

# Purge obsolete variables
rm(inst_pkgs, req_pkgs_cran, req_pkgs_bioc, shh)

cat("Packages loaded!\n")
