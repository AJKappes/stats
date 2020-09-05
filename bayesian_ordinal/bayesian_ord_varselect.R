# Variable selection
# Bayesian modeling of latent ordinal data parameters

remove(list = objects())
setwd('~/research/stats/bayesian_ordinal/')

library(tidyverse)

source('bayesian_ord_datagen.R')
source('bayesian_ord_gibbs_samp.R')


# Generate data -----------------------------------------------------------
# X: choose number of aggregate covariates, k
#    choose sample size nsamp
#    the number of actual covariates that generate y is p
#
# Betas: generate parameters based on the number of covariates p
#        sample from N(0, 3)
#
# y: X and B compute y in the form y = c + XB + e
#    c sampled from N(0, 10)
#    e sampled from N(0, 1)
#
# z: observed ordinal data, generated from inital y

nsamp <- 500
k <- 25
# initial p only used to generate y
p <- 5
# aggregate X data matrix
X_k <- get_x(k, nsamp)

# random sample covariates to get initial p X covariates
p_idx <- sample(1:k, p)

# create p covariate X data matrix, generate B and y
X <- X_k[, p_idx]
B_p <- get_B(p)
y <- get_y(X, B_p)

# create observed ordinal z
z <- get_z(y)


# Gibbs and var samp ------------------------------------------------------
# Gibbs function returns:
#   prediction accuracy
#   estimated beta
#   delta parameters

ngibb <- 100
bench <- gibbs_fn(ngibb, X, p)
bench_acc <- bench$accuracy

cat('True X covariates indice:', p_idx,
    '\nBenchmark ordinal value prediction accuracy:', bench_acc,
    '\n')

# variable sampling
nvsamp <- 10
acc_vec <- c()
p_idx_list <- vector(mode = 'list', length = nvsamp)

vsamp_fn <- function(k) {
  
  p_samp <- sample(1:k, 1)
  p_idx_samp <- sample(1:k, p_samp)
  X_samp <- X_k[, p_idx_samp]
  gibbs_samp <- gibbs_fn(ngibb, X_samp, p_samp)
  
  out <- list(gibbs = gibbs_samp,
              j_samp = p_idx_samp)
  return(out)
  
}

for (iter in 1:nvsamp) {
  
  cat('\nGibbs function iteration', iter,
      '\n')
  
  s <- vsamp_fn(k)
  acc_vec[iter] <- s$gibbs$accuracy
  p_idx_list[[iter]] <- s$j_samp
  
}


# Testing -----------------------------------------------------------------



