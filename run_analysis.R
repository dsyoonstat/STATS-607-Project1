# Activate renv
renv::activate()


# Attach packages
library(dplyr)
library(data.table)


# Load functions
fun_files <- list.files("src", pattern = "\\.R$", full.names = TRUE)
invisible(lapply(fun_files, source))


# Load data
load("data/census_full.rda")

# Split hierarchical key and key variable columns
hkey = colnames(census_full)[1:3]
key = colnames(census_full)[4:8]


# Make full table (rds format)
savefulltb(census_full, hkey = hkey, key = key, B = 5, 
           output.path = "artifacts/fulltable.rds")


# Save aggregated table (csv format)
saveaggtb(hkey.level = 2, key = key[1:4], input.path = "artifacts/fulltable.rds", 
          output.table.path = "results/aggtable.csv", output.infoloss.path = "results/infoloss.csv")
