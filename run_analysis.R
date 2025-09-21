# Activate renv
renv::activate()


# Attach packages
library(dplyr)
library(data.table)
library(testthat)


# Load functions
fun_files <- list.files("src", pattern = "\\.R$", full.names = TRUE)
invisible(lapply(fun_files, source))


# Load data
load("data/census_full.rda")


# Split hierarchical key and key variable columns and set B
hkey <- colnames(census_full)[1:3]
key <- colnames(census_full)[4:8]
mask_threshold <- 5


# Make full table (rds format)
save_full_tb(data = census_full,
             hkey = hkey,
             key = key,
             mask_threshold = mask_threshold,
             output_path = "artifacts/full_table.rds")


# Save aggregated table (csv format)
save_agg_tb(hkey_level = 2,
          key = key[1:4],
          input_path = "artifacts/full_table.rds",
          output_table_path = "results/agg_table.csv",
          output_infoloss_path = "results/infoloss.csv")


# Run tests
test_dir("tests/testthat")
