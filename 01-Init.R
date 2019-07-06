cat("Initialising...\n")

# Delete all objects
rm(list = setdiff(ls(), "linked_path"))

# Get user
user <- Sys.getenv("USER")

# Define paths to use
idir <- file.path(linked_path, "input")
odir <- file.path(linked_path, "output", user)

# Create odir if not already existent
dir.create(odir, showWarnings = FALSE, recursive = TRUE)

# Delete all files and folders in odir
rmdirs <- dir(odir)
# if (length(rmdirs) > 0) {
#   unlink(file.path(odir, rmdirs), recursive = TRUE)
# }

# Purge obsolete variables
rm(user, rmdirs)

cat("Initialised!\n")
