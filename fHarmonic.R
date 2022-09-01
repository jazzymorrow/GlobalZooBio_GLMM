# Harmonic is to fit Hour and DOY

fHarmonic <- function (theta, k = 4) {
  X <- matrix(0, length(theta), 2 * k)
  nam <- as.vector(outer(c("c", "s"), 1:k, paste, sep = ""))
  dimnames(X) <- list(names(theta), nam)
  m <- 0
  for (j in 1:k) {
    X[, (m <- m + 1)] <- cos(j * theta)
    X[, (m <- m + 1)] <- sin(j * theta)
  }
  X
}

#### Notes ####
# creates empty matrix (length(theta) x 2k)
# columns are c1 and s1 and contain:
#c1 = cos(theta), s1 = sin(theta)

# why would k be greater than 1? extends the wave, multiple cycles 


