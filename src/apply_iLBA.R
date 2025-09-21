#' Information Loss Bounded Aggregation (iLBA)
#'
#' Applies iLBA algorithm when aggregating cell frequencies.
#'
#' @param input_summary Numeric vector of length 3: c(n_small, n_SCA_B, sum_small)
#'   - n_small: number of original cells ≤ B, K in the paper
#'   - n_SCA_B: number of cells equal to B after SCA, k in the paper
#'   - sum_small: true frequency from upper table, f_i in the paper
#' @param mask_threshold Threshold for masking (default = 3).
#' @return Numeric vector c(sum_masked, type1, type2) is returned.
#' @references Park, M.J., Kim, H.J., & Kwon, S. (2024).
#' *Journal of the Korean Statistical Society*. Springer Nature.
#'
apply_iLBA <- function(input_summary, mask_threshold = 3) {
  
  # Divide input summary
  n_small <- input_summary[1]
  n_SCA_B <- input_summary[2]
  sum_small <- input_summary[3]
  
  
  # Initialize masked sum as n_SCA_B * mask_threshold
  sum_masked <- n_SCA_B * mask_threshold
  
  
  # Create type 1,2 indicators boundary cases of the iLBA algorithm
  type1 <- 0
  type2 <- 0
  
  
  # Apply the iLBA algorithm
  if (n_small > 1) {
    
    # Compute default masked sum for n_small > 1
    quotient <- (sum_small - 1) %/% mask_threshold
    sum_masked <- quotient * mask_threshold + trunc(mask_threshold / 2) + 1
    
    
    # Compute lower and upper bounds inferred by SCA and system
    lower_SCA <- n_SCA_B
    upper_SCA <- n_SCA_B + n_small * (mask_threshold - 1)
    lower_sys <- quotient * mask_threshold + 1
    upper_sys <- (quotient + 1) * mask_threshold


    # Compute masked sum for boundary cases
    if (quotient < 0) {
      sum_masked <- 0
    } else if (lower_sys < lower_SCA) {
      sum_masked <- sum_masked + mask_threshold
      type1 <- 1
    } else if (upper_sys > upper_SCA) {
      sum_masked <- sum_masked - mask_threshold
      type2 <- 1
    }


    # Deal with special case
    if (sum_masked > 0 && sum_masked < mask_threshold) {
      sum_masked <- mask_threshold      
    }
  }
  
  
  # Return masked sum and type 1,2 indicators
  return(c(sum_masked, type1, type2))
}
