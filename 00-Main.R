# Linked path on file share
linked_path <- "~/mnt/agroscope_os/2/2/6/1/1/4/6021/GenSel_Costa_collab"
# For local analyses use
# linked_path <- "."

# Initialise
source("01-Init.R")

# Declare necessary packages
pkgs_cran <- c("config", "magrittr", "lme4","doBy", "snpStats", "snpReady", "readxl","rrBLUP","reshape2",
               "ggplot2","parallel","dplyr","gplots", "data.table", "corrplot","FactoMineR","factoextra","missMDA", "purr",
               "MM4LMM", "nnet", "impute", "snpStats")
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
## Model phenos and genos for downstream analyses
source("05-ModelPhenos.R")
source("06-ModelGenos.R.R")
# source("07-VisualiseData.R")
## Do analysis of variance to assess genotypic contribution to trait
source("08-VarianceAnalysis.R")

# Explore
## Vizualise
source("07-VisualiseData.R")
## Model
source("08-ModelData.R")

# Communicate
## Save Output
source("09-SaveOutput.R")
