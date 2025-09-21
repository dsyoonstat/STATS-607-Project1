#' Small Cell Adjustment (SCA)
#'
#' Applies SCA algorithm to cell frequencies below a threshold B.
#'
#' @param freq Numeric vector of cell frequencies.
#' @param mask_threshold Threshold for masking (default = 3).
#' @return Numeric vector of masked cell frequencies is returned.
#' @references Hundepool, A., Domingo-Ferrer, J., Franconi, L., Giessing, S.,
#' Nordholt, E. S., Spicer, K., & de Wolf, P.-P. (2012).
#' *Statistical Disclosure Control*. Wiley.
#'
apply_SCA <- function(freq, mask_threshold = 3) {
  
  # Convert input frequency into vector
  if (!is.vector(freq)) freq <- as.vector(freq)
  
  
  # Return SCA masked frequency
  changer <- ifelse(runif(length(freq)) < freq / mask_threshold, mask_threshold, 0)
  return(ifelse(freq < mask_threshold, changer, freq))
}
