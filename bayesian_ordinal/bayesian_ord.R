# Replication and nonlinear extension
# Bayesian Inference for Ordinal Categorical Data
# Reference: Dr. Yuan Wang, WSU

remove(list = objects())
library(tidyverse)
library(MCMCpack)
library(TruncatedNormal)


### Generate data -----------------------------------------------------
set.seed(7)
nsamp <- 400
p <- 4
X <- matrix(rnorm(nsamp*p), nrow = nsamp, ncol = p)

# linear latent construction
y <- 5 + 3*X[, 1] + 4*X[, 1]^2 + rnorm(nsamp)

# nonlinear latent construction
fn_ <- function(X, B = FALSE) {
  
  if (B) {
    
    if (length(B) > 1) {
      
      fx <- max(y)/(1 + exp(-X%*%B))
      
    } else {
      
      fx <- max(y)/(1 + exp(-B*X))
      
    }
  
  } else {
    
    fx <- max(y)/(1 + exp(-X))
    
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
z <- get_z(y_nl)


### parameters ---------------------------------------------------------
## betas: linear weights
##        assume mvn(mu_0, Sigma_0) prior, (0, I)          
## v: latent variable boundaries F(delta)
##    assume dir(a) prior, pdf v_j^{a_j - 1}, v in (0, 1)

# linear betas
# betas <- rep(0, p)
# nonlinear betas
betas <- rep(0, p)
mu0 <- rep(0, p)
Sigma0 <- diag(length(betas))
y_std = 1
a <- rep(2, max(z))
u <- c(0, rdirichlet(1, a))

# conditional parameters - fixed
## Sigma0 inverse defined by 1/100*Sigma0
## linear
# Sigma1 <- solve(1/100*Sigma0 + t(X)%*%X)
## nonlinear
Sigma1 <- solve(1/100*Sigma0 + t(fn_(X))%*%fn_(X))

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
U_sample <- matrix(0, nrow = nsim, ncol = length(u))
betas_sample <- matrix(0, nrow = nsim, ncol = length(betas))
for (iter in 1:nsim) {
  
  # latent variable draw
  y_g <- sapply(1:nsamp,
                function(i) rtnorm(1,
                                   lb = qnorm(sum(u[1:z[i]])),
                                   ub = qnorm(sum(u[1:z[i] + 1])),
                                   # linear
                                   # mu = X[i, ]%*%betas,
                                   # nonlinear
                                   mu = fn_(X[i, ], betas),
                                   sd = y_std))
  
  # linear weight draw
  ## conditional parameters - update
  ## linear
  # mu1 <- Sigma1%*%(t(X)%*%y_g + solve(Sigma0)%*%mu0)
  # nonlinear
  mu1 <- Sigma1%*%(t(fn_(X))%*%y_g + solve(Sigma0)%*%mu0)
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
  U_sample[iter, ] <- u
  betas_sample[iter, ] <- betas
  
}
 

### Sampling results ------------------------------------------------
beta_gibb <- colMeans(betas_sample)[1]

## linear least squares check
# beta <- coef(lm(y_g ~ X[, 1] + X[, 1]^2))[2] %>% unname()

## nonlinear least squares check
nlmod <- nls(y_g ~ max(y)/(1 + exp(-b1*X[, 1])),
             start = list(b1 = 1))
beta <- coef(nlmod) %>% unname()


## linear approx for nonlinear fn
# beta <- coef(lm(y_g ~ X[, 1]))[2] %>% unname()

# plot sampling values
plot(1:nsim, betas_sample[, 1], type = 'l', main = 'Beta_1 Sampling Values')

## prediction
y_nl_hat <- fn_(X[, 1], beta_gibb)
mspe <- 1/nsamp*sum((y_nl - y_nl_hat)^2)
plot(X[, 1], y_nl, main = 'Prediction Overlay (red)')
lines(X[, 1], y_nl_hat, type = 'p', col = 'red')

z_hat <- get_z(y_nl_hat)
acc_vec <- map2_dbl(z, z_hat,
                    function(i, ip) ifelse(i == ip, 1, 0))
acc <- mean(acc_vec)

# print results
cat('Beta_1 Estimates:',
    '\n (Gibbs sampling)', beta_gibb,
    '\n (Least squares check)', beta,
    '\n',
    '\nOrdinal value prediction accuracy rate:', acc,
    '\n')


### TODO -------------------------------------------------------------










