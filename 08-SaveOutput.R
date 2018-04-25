cat("Save output...\n")

# Create subfolder with reformatted date and time.
path <- paste0(odir, ts_done, "/")
dir.create(path, showWarnings = FALSE)

# Your code goes here

# Purge obsolete variables
rm(path)

cat("Output saved!\n")
