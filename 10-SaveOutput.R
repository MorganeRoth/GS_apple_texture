cat("Saving output...\n")

# Create subfolder with reformatted date and time.
path <- file.path(odir, format(Sys.time(), "%Y-%m-%d_%H%M%S"))
dir.create(path, showWarnings = FALSE)

# Your code goes here

# Purge obsolete variables
rm(path)

cat("Output saved!\n")
