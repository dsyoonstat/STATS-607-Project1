# Activate renv
renv::activate()


# Attach packages
library(dplyr)
library(data.table)


# Load functions
fun_files <- list.files("src", pattern = "\\.R$", full.names = TRUE)
invisible(lapply(fun_files, source))


# Load data
load("census_full.rda")
hkey = colnames(census_full)[1:3]
key = colnames(census_full)[4:8]

# Make fulltable
savefulltb(census_full, hkey = hkey, key = key, B = 5)
saveaggtb(hkey.level = 2, key = key[1:4])
