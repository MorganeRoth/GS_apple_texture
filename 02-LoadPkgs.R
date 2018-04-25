cat("Loading packages...\n")

# Install required packages
inst_pkgs <- installed.packages()[, "Package"]
req_pkgs <- setdiff(nec_pkgs, inst_pkgs)
if (length(req_pkgs) > 0) {
  install.packages(req_pkgs, clean = TRUE, quiet = TRUE)
}

# Load necessary packages
shh <- lapply(nec_pkgs, library, character.only = TRUE)

# Purge obsolete variables
rm(inst_pkgs, req_pkgs, shh)

cat("Packages loaded!\n")
