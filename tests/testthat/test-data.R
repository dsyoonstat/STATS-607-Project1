test_that("columns 1-3 satisfy the hierarchical code structure (2-5-7 digits)", {
  
  # Load data (only if run separately)
  # data_path <- test_path("..", "..", "data", "census_full.rda")
  # load(data_path)
  
  
  # Convert integer hierarchical values into strings
  hkey_1 <- sprintf("%02d", as.integer(census_full[[1]]))
  hkey_2 <- sprintf("%05d", as.integer(census_full[[2]]))
  hkey_3 <- sprintf("%07d", as.integer(census_full[[3]]))
  
  
  # Compute booleans for hierarchical structure
  ok12 <- substr(hkey_2, 1, 2) == hkey_1
  ok23 <- substr(hkey_3, 1, 5) == hkey_2
  
  
  # Find rows with unvalid hierarchical structure
  bad_idx <- which(!(ok12 & ok23) | is.na(ok12) | is.na(ok23))
  
  expect_true(
    length(bad_idx) == 0,
    paste0("Hierarchy violated in rows: ", paste(head(bad_idx, 10), collapse = ", "))
  )
})
