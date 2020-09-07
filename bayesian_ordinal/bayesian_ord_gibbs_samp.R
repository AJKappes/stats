# Gibbs sampling routine
# Bayesian modeling of latent ordinal data prameters

library(tidyverse)
library(MCMCpack)
library(TruncatedNormal)


# Gibbs function ----------------------------------------------------------

gibbs_fn <- function(ngibb, X, p) {
  
  # sampling progress
  samp_prog <- seq(0, ngibb, 10)
  
  # initialize parameters
  betas <- rep(1, p)
  mu0 <- rep(0, p)
  Sigma0 <- diag(length(betas))
  y_std = 1
  a <- rep(2, max(z))
  u <- c(0, rdirichlet(1, a))
  
  # beta marginal hyperparameter fixed
  Sigma1 <- solve(1/100*Sigma0 + t(X)%*%X)
  
  # importance sampling A matrix construction
  A <- diag(length(unique(z)))
  A[lower.tri(A)] <- 1
  
  y_sample <- matrix(0, nrow = ngibb, ncol = nsamp)
  U_sample <- matrix(0, nrow = ngibb, ncol = length(u))
  betas_sample <- matrix(0, nrow = ngibb, ncol = length(betas))
  
  for (iter in 1:ngibb) {
    
    # latent variable draw
    y_g <- sapply(1:nsamp,
                  function(i) rtnorm(1,
                                     lb = qnorm(sum(u[1:z[i]])),
                                     ub = qnorm(sum(u[1:(z[i] + 1)])),
                                     # linear
                                     mu = X[i, ]%*%betas,
                                     # nonlinear
                                     # mu = fn_(X[i, ], betas),
                                     sd = y_std))
    
    # linear weight draw
    ## conditional parameters - update
    
    ## linear
    mu1 <- Sigma1%*%(t(X)%*%y_g + 1/100*Sigma0%*%mu0)
    
    ## nonlinear
    # mu1 <- Sigma1%*%(t(fn_(X, 'cond'))%*%y_g + 1/100*Sigma0%*%mu0)
    
    betas <- mvrnorm(mu = mu1, Sigma = Sigma1)
    
    # importance sampling for u draw
    ## uniform sampling bounds
    ## distribution function(y) for latent y being within bounds
    bounds <- sort(unique(z))
    Cmin <- sapply(bounds, function(c) min(pnorm(y_g[z == c])))
    Cmax <- sapply(bounds, function(c) max(pnorm(y_g[z == c])))
    
    # sample u from uniform for close boundary value approximation
    ## q1 ~ U(Cmax1, Cmin2), q2 ~ U(Cmax2, Cmin3), q3 ~ U(Cmax3, 1)
    nimport <- 100
    q1 <- runif(nimport, min = Cmax[1], max = Cmin[2])
    q2 <- runif(nimport, min = Cmax[2], max = Cmin[3])
    # q3 <- runif(nimport, min = Cmax[3], max = 1)
    q3 <- 1
    Q <- rbind(q1, q2, q3)
    U_import <- solve(A)%*%Q
    dens <- apply(U_import, 2, function(j) ddirichlet(j, alpha = a))
    p_draw <- dens/sum(dens)
    U_draw <- sample(1:ncol(U_import), 1, prob = p_draw)
    u <- c(0, U_import[, U_draw])
    
    # store draws
    y_sample[iter, ] <- y_g
    U_sample[iter, ] <- u
    betas_sample[iter, ] <- betas
    
    if (iter %in% samp_prog) {
      
      cat('Sampling iteration', iter, 'complete',
          '\n')
      
    }
    
  }
  
  beta_gibb <- colMeans(betas_sample)
  y_hat <- X%*%matrix(beta_gibb)
  
  u_hat <- colMeans(U_sample)[2:4]
  u_bounds <- sapply(1:(length(u_hat) - 1), function(i) sum(u_hat[1:i]))
  delta <- c(min(y_hat),
             sapply(1:length(u_bounds), function(i) qnorm(u_bounds[i])))
  obs_bounds <- c(min(y), quantile(y, .33), quantile(y, .67))
  z_hat <- sapply(y_hat,
                  function(z) (z >= delta[1]) + (z > delta[2]) + (z > delta[3])) %>% 
    unname()
  
  acc_vec <- map2_dbl(z, z_hat, function(i, ip) ifelse(i == ip, 1, 0))
  acc <- mean(acc_vec)
  
  out <- list(accuracy = acc,
              zpred = z_hat,
              beta = beta_gibb,
              delta = delta)
  return(out)
  
}