cat("Write to DB...\n")

# Create subfolder with reformatted date and time.
path <- paste0(odir, ts_done, "/")
dir.create(path, showWarnings = FALSE)

# Your code goes here

# Purge obsolete variables
rm(path)

cat("Written to DB!\n")
