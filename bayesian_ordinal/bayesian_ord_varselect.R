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
nvsamp <- 25

vsamp_fn <- function(k) {
  
  p_samp <- sample(1:k, 1)
  p_idx_samp <- sample(1:k, p_samp)
  X_samp <- X_k[, p_idx_samp]
  gibbs_samp <- gibbs_fn(ngibb, X_samp, p_samp)
  
  out <- list(gibbs = gibbs_samp,
              j_samp = p_idx_samp)
  return(out)
  
}


acc_idx_fn <- function(k, nvsamp) {
  
  acc_vec <- c()
  z_idx_list <- vector(mode = 'list', length = nvsamp)
  p_idx_list <- vector(mode = 'list', length = nvsamp)
  
  iter <- 1
  while (iter <= nvsamp) {
    
    cat('\nGibbs function iteration', iter,
        '\n')
    
    s <- vsamp_fn(k)
    acc_vec[iter] <- s$gibbs$accuracy
    z_idx_list[[iter]] <- s$gibbs$zpred
    p_idx_list[[iter]] <- s$j_samp
    
    iter <- iter + 1
    
  }
  
  m_idx <- which(acc_vec == max(acc_vec))
  z_m <- z_idx_list[[m_idx]]
  p_m <- p_idx_list[[m_idx]]
  
  out <- list(maxpred = acc_vec[m_idx],
              z_maxpred = z_m,
              p_maxpred = p_m)
  return(out)
  
}

# get sampling results
samp_res <- acc_idx_fn(k, nvsamp)
samp_z <- samp_res$z_maxpred

# TP and FP rate construction}
j <- 1:max(z)
zmax_idx <- lapply(j, function(i) which(samp_z == i))
zpred_list <- lapply(j, function(i) z[zmax_idx[[i]]])
tp <- sapply(j, function(i) sum(zpred_list[[i]] == i))
tpr <- map2_dbl(tp, j, function(r, i) r/length(zpred_list[[i]])) %>% 
  round(3)

# confusion matrix construction
conf_mat <- diag(tp)
conf_mat_r <- diag(tpr)
for (c in j) {
  
  jfp <- j[j != c]
  
  fp <- sapply(jfp, function(i) sum(zpred_list[[c]] == i))
  fpr <- sapply(fp, function(r) r/length(zpred_list[[c]])) %>% 
    round(3)
  
  conf_mat[c, jfp] <- fp
  conf_mat_r[c, jfp] <- fpr
  
}

jnames <- paste('J', j, sep = '')
jprednames <- paste(jnames, '_hat', sep = '')
conf_tbl <- data.frame(cbind(jprednames, conf_mat))
colnames(conf_tbl) <- c('', jnames)
conf_r_tbl <- data.frame(cbind(jprednames, conf_mat_r))
colnames(conf_r_tbl) <- c('', jnames)



# Testing -----------------------------------------------------------------

