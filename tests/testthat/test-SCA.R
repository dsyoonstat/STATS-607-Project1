# Load apply_SCA.R (only if run separately)
# SCA_path <- testthat::test_path("..", "..", "src", "apply_SCA.R")
# source(SCA_path)


test_that("SCA replaces values < B only with {0, B} and leaves values >=B unchanged", {
  
  # Define input and apply SCA
  SCA_input <- sample(0:10, size = 20, replace=T)
  SCA_output <- apply_SCA(SCA_input, mask_threshold = 5)
  
  
  # Check whether values < B are changed into {0,B}
  expect_true(all(SCA_output[SCA_input < 5] %in% c(0, 5)))
  
  
  # Check whether values >= B are unchanged
  expect_equal(SCA_output[SCA_input >= 5], SCA_input[SCA_input >= 5])
})