# Linked path on file share
linked_path <- "../../mnt/path/on/file/share/"

# Initialise
source("01-Init.R")

# Declare necessary packages
pkgs_cran <- c("config")
pkgs_bioc <- c()

# Load packages
source("02-LoadPkgs.R")
# Load functions
source("03-LoadFunc.R")

# Configuration
config <- config::get()

# Wrangle
## Import
source("04-ImportData.R")
## Tidy
source("05-TidyData.R")

# Wrangle/Explore
## Transform
source("06-TransformData.R")

# Explore
## Vizualise
source("07-VisualiseData.R")
## Model
source("08-ModelData.R")

# Communicate
## Save Output
source("09-SaveOutput.R")
