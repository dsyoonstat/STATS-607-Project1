SCA <- function(freq, B = 3) {
  if (!is.vector(freq)) freq <- as.vector(freq)
  if (!sum(freq < B)) return(freq)
  n <- length(freq)
  changer <- ifelse(runif(n) < freq / B, B, 0)
  return(ifelse(freq < B, changer, freq))
}