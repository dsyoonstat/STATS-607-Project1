test_that("infoloss are within the theoretical range [-[3B/2],[3B/2]]", {
  
  # Load infoloss.csv
  infoloss_path <- test_path("..", "..", "results", "infoloss.csv")
  infoloss <- read.csv(infoloss_path)
  
  
  # Get infoloss values from infoloss.csv
  suppressWarnings({
    loss_num <- as.numeric(infoloss$Loss)
  })
  is_num <- !is.na(loss_num)
  tested_vals <- loss_num[is_num]
  
  
  # Compute theoretical lower and upper bounds
  # B should be same with the input B used in savefulltb
  B <- 5
  lower <- -floor(3 * B / 2)
  upper <-  floor(3 * B / 2)
  
  
  # Check whether the infoloss values are within the theoretical range
  in_range <- (tested_vals >= lower) & (tested_vals <= upper)
  bad_vals <- tested_vals[!in_range]
  
  expect_true(
    all(in_range),
    paste0(
      "Loss out of range [", lower, ", ", upper, "]: ",
      paste(head(bad_vals, 10), collapse = ", ")
    )
  )
})
