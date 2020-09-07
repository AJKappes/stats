# Data generation
# Bayesian modeling of latent ordinal data parameters

set.seed(7)

get_x <- function(k, nsamp) {
  
  X <- matrix(rnorm(nsamp*k), nrow = nsamp, ncol = k)
  return(X)
  
}

get_B <- function(p) {
  
  betas <- matrix(rnorm(p, sd = 3))
  return(betas)
  
}

# linear latent construction
get_y <- function(X, beta) {
  
  y <- rnorm(1, sd = 10) + X%*%B_p + rnorm(nsamp)
  
}

get_z <- function(latent) {
  
  y0 <- min(latent)
  y33 <- quantile(latent, .33)
  y67 <- quantile(latent, .67)
  y100 <- max(latent)
  
  z <- sapply(latent,
              function(z) (z >= y0) + (z > y33) + (z > y67)) %>%
    unname()
  
  return(z)
  
}

