# Variable Selection Methods
# Bayesian Inference for Ordinal Categorical Data
# Reference: Dr. Yuan Wang, WSU

remove(list = objects())
library(tidyverse)
library(MCMCpack)
library(TruncatedNormal)


### Generate data -----------------------------------------------------
set.seed(7)
nsamp <- 400

# initialize 20 covariates
# randomly pick 5 to generate latent response
k <- 20
p <- 5
X_k <- matrix(rnorm(nsamp*k), nrow = nsamp, ncol = k)
X <- X_k[, sample(1:dim(X_k)[2], p)]
B_p <- matrix(runif(p, 1, 20))

# linear latent construction
y <- runif(1, 1, 20) + X%*%B_p + rnorm(nsamp)

dt <- cbind(y, X) %>% 
  data.frame()
colnames(dt)[1] <- 'y'

# nonlinear latent construction
L <- 5

# normalize y to (0,1)
fn_ <- function(X, B) {
  
  if (length(B) > 1) {
    
    fx <- L/(1 + exp(-X%*%B))
    
  } else if (is.numeric(B) & length(B) == 1) {
    
    fx <- L/(1 + exp(-B*X))  
    
  } else if (B == 'cond') {
    
    fx <- L/(1 + exp(-X))
    
  }
  
  return(fx)
  
}

y_nl <- fn_(X[, 1], 3) + rnorm(nsamp)
# plot(X[, 1], y_nl)

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

# observed ordinal data
z <- get_z(y)


### parameters ---------------------------------------------------------
## betas: linear weights
##        assume mvn(mu_0, Sigma_0) prior, (0, I)          
## v: latent variable boundaries F(delta)
##    assume dir(a) prior, pdf v_j^{a_j - 1}, v in (0, 1)

# linear betas
# betas <- rep(0, p)
# nonlinear betas
betas <- rep(1, p)
mu0 <- rep(0, p)
Sigma0 <- diag(length(betas))
y_std = 1
a <- rep(2, max(z))
u <- c(0, rdirichlet(1, a))

# conditional parameters - fixed
## Sigma0 inverse defined by 1/100*Sigma0

## linear
Sigma1 <- solve(1/100*Sigma0 + t(X)%*%X)

## nonlinear
# Sigma1 <- solve(1/100*Sigma0 + t(fn_(X, 'cond'))%*%fn_(X, 'cond'))

## importance sampling construction
A <- diag(length(unique(z)))
A[lower.tri(A)] <- 1


### Gibbs sampling y, betas, u ------------------------------------------
## y | betas, u ~ N(X'betas, 1)
## betas | y, u ~ mvn(mu_1, Sigma_1)
##                linear
##                  mu_1 = Sigma_1^{-1}(sum(y_iX_i) + Sigma_0^{-1}mu_0)
##                  Sigma_1^{-1} = (Sigma_0^{-1} + sum(X_iX_i'))^{-1}
##                nonlinear  
##                  mu_1 = 
##                  Sigma_1 = 
## u | y, betas ~ trunc dir sum_1^{z_i-1}p_k < F(y_i) <= sum_{1}^z_i

nsim <- 1000
y_sample <- matrix(0, nrow = nsim, ncol = nsamp)
U_sample <- matrix(0, nrow = nsim, ncol = length(u))
betas_sample <- matrix(0, nrow = nsim, ncol = length(betas))
for (iter in 1:nsim) {
  
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
  
  cat('Sampling iteration', iter, 'complete',
      '\n')
  
}


### Sampling results ------------------------------------------------
beta_gibb <- colMeans(betas_sample)

## linear least squares check
beta <- coef(lm(y ~ ., data = dt))[2:(p+1)] %>%
  unname()


## nonlinear least squares check
# y_g_est <- colMeans(y_sample)
# nlmod <- nls(y_g_est ~ L/(1 + exp(-b1*X[, 1])),
#              start = list(b1 = 1))
# beta <- coef(nlmod) %>% unname()

## linear approx for nonlinear fn
# beta <- coef(lm(y_g ~ X[, 1]))[2] %>% unname()

# plot sampling values
# plot(1:nsim, betas_sample[, 1], type = 'l', main = 'Beta_1 Sampling Values')

## prediction
# linear
y_hat <- X%*%matrix(beta_gibb)

# nonlinear
# y_hat <- fn_(X[, 1], beta_gibb)
# mspe <- 1/nsamp*sum((y_nl - y_hat)^2)
# plot(X[, 1], y_nl, main = 'Prediction Overlay (red)')
# lines(X[, 1], y_hat, type = 'p', col = 'red')

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

# print results
cat('Beta_1 Estimates:',
    '\n (Gibbs sampling)', beta_gibb,
    '\n (Least squares check)', beta,
    '\n',
    '\nDelta Estimates (J = 3):', delta,
    '\n',
    '\nOrdinal value prediction accuracy rate:', acc,
    '\n')


### Param conv plot ---------------------------------------------------------

## convergence
beta_conv <- sapply(1:nsim, function(i) mean(betas_sample[1:i, 1]))

u_conv <- matrix(0, nrow = nsim, ncol = length(u_hat))
for (j in 1:ncol(u_conv)) {

  conv <- sapply(1:nsim, function(i) mean(U_sample[1:i, j + 1]))
  u_conv[, j] <- conv

}

# par(mfrow = c(2, 2), oma = c(0, 0, 2, 0))
# plot(1:nsim, beta_conv, main = 'b1')
# plot(1:nsim, u_conv[, 1], main = 'u1')
# plot(1:nsim, u_conv[, 2], main = 'u2')
# plot(1:nsim, u_conv[, 3], main = 'u3')
# mtext('Sampling Estimator Convergence', outer = TRUE)

 
