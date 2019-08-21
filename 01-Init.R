cat("Initialising...\n")

# Delete all objects
rm(list = setdiff(ls(), "linked_path"))

# Copy sample configuration file if it doesn't exist
if (!file.exists("config.yml")) {
  file.copy("config.yml.sample", "config.yml")
}

# Get user
if (.Platform$OS.type == "windows") {
  user <- Sys.getenv("USERNAME")
} else {
  user <- Sys.getenv("USER")
}

# Define paths to use
idir <- file.path(linked_path, "R_input")
odir <- file.path(linked_path, "R_output/20190816")

# Create odir if not already existent
dir.create(odir, showWarnings = FALSE, recursive = TRUE)

# Delete all files and folders in odir
rmdirs <- dir(odir)
# if (length(rmdirs) > 0) {
#   unlink(file.path(odir, rmdirs), recursive = TRUE)
# }

# set seed for reproducible results

set.seed(1234)

# Purge obsolete variables
rm(user, rmdirs)

cat("Initialised!\n")
