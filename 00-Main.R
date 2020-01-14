# Linked path on file share
linked_path <- "~/mnt/agroscope_os/2/2/6/1/1/4/6021/GenSel_Costa_collab"
# For local analyses use
# linked_path <- "."

# Initialise
source("01-Init.R")

# Declare necessary packages
pkgs_cran <- c("config", "magrittr", "lme4","doBy", "snpStats", "snpReady", "readxl","rrBLUP","reshape2",
               "ggplot2","parallel","dplyr","gplots", "data.table", "corrplot","FactoMineR","factoextra","missMDA", "purrr",
               "MM4LMM", "nnet", "impute", "snpStats", "robustHD", "adegenet", "agridat", "STPGA", "plyr",
               "EMMREML", "viridis")
pkgs_bioc <- c("chopsticks")

# Load packages
source("02-LoadPkgs.R")
# Load functions
source("03-LoadFunc.R")

# Configuration
config <- config::get()

## if you want to start directly with predictions go to each specific script
## and load data from there
# Wrangle
## Import
source("04-ImportData.R")
## Model phenos and genos for downstream analyses
source("05-ModelPhenos.R")
source("06-ModelGenos.R")
# source("07-VisualiseData.R")
## Do analysis of variance to assess genotypic contribution to trait
source("08-VarianceAnalysis.R")

## Use different prediction scenarios
source("09-PredictCOLtoFAM.R")
source("10-Predict_COLLtoCOLL.R")
source("11-Predict_COLL_propFAM_toFAM.R")
source("12-Predict_COLL_relatedFAM_toFAM.R")
source("13-Predict_COLL_to_FAM_relatedness.R")
source("14-Predict_STPGA_package.R")

# Communicate
## Save Output
source("09-SaveOutput.R")
