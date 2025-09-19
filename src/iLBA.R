iLBA <- function(x, B = 3) {
  nLessEqOrigin <- x[1] # K
  nEqSCA <- x[2]        # k
  SumSCA <- x[3]        # f_i
  Masked <- nEqSCA * B
  type1 <- type2 <- 0
  if (nLessEqOrigin > 1) {
    a <- (SumSCA - 1) %/% B
    Masked <- a * B + trunc(B / 2) + 1
    lbdSCA <- nEqSCA
    ubdSCA <- nEqSCA + nLessEqOrigin * (B - 1)
    lbdSystem <- a * B + 1
    ubdSystem <- (a + 1) * B
    if (a < 0) {
      Masked <- 0
    } else if (lbdSystem < lbdSCA) {
      Masked <- Masked + B
      type1 <- 1
    } else if (ubdSystem > ubdSCA) {
      Masked <- Masked - B
      type2 <- 1
    }
    if (Masked > 0 && Masked < B) Masked <- B
  }
  return(c(Masked, type1, type2))
}