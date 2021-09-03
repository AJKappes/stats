# Variable selection
# Bayesian modeling of latent ordinal data parameters

remove(list = objects())
setwd('~/research/stats/bayesian_ordinal/')

library(tidyverse)
library(gridExtra)
library(ggplot2)

source('bayesian_ord_datagen.R')
source('bayesian_ord_gibbs_samp.R')


# Functions ---------------------------------------------------------------


# random variable sampling
vsamp_fn <- function(k) {
  
  p_samp <- sample(1:k, 1)
  p_idx_samp <- sample(1:k, p_samp)
  X_samp <- X_k[, p_idx_samp]
  gibbs_samp <- gibbs_fn(ngibb, X_samp, p_samp, z, nsamp)
  
  out <- list(gibbs = gibbs_samp,
              j_samp = p_idx_samp)
  return(out)
  
}

# fixed variable sampling
vsamp_fix_fn <- function(k, p) {
  
  p_samp <- sample(1:k, p)
  X_samp <- X_k[, p_samp]
  gibbs_samp <- gibbs_fn(ngibb, X_samp, length(p_samp), z, nsamp)
  
  out <- list(gibbs = gibbs_samp,
              j_samp = p_samp)
  return(out)
  
}

acc_idx_fn <- function(k, nvsamp, vs, p = FALSE) {
  
  acc_vec <- c()
  b_nonzero_count <- c()
  b_nonzero <- vector(mode = 'list', length = nvsamp)
  tp_rate <- c()
  fp_rate <- c()
  b_zero_idx <- vector(mode = 'list', length = nvsamp)
  beta_samp_list <- vector(mode = 'list', length = nvsamp)
  beta_list <- vector(mode = 'list', length = nvsamp)
  beta_se_list <- vector(mode = 'list', length = nvsamp)
  z_idx_list <- vector(mode = 'list', length = nvsamp)
  p_idx_list <- vector(mode = 'list', length = nvsamp)
  
  iter <- 1
  while (iter <= nvsamp) {
    
    cat('\nGibbs function iteration', iter,
        '\n')
    
    if (vs == 'random') {
      
      s <- vsamp_fn(k)
      
    } else if (vs == 'fixed') {
      
      s <- vsamp_fix_fn(k, p)
      
    } else {
      
      stop('must designate vs = (random, fixed)')
      
    }
    
    j_samp <- s$j_samp
    b_samp <- s$gibbs$beta_samp
    b <- s$gibbs$beta
    tcrit <- qt(1 - .05/2, nsamp - length(b))
    b_se <- sqrt(s$gibbs$beta_var)
    lb <- b - tcrit*b_se
    ub <- b + tcrit*b_se
    b_zero_idx <- map2_dbl(lb, ub,
                           function(l, u)
                             (l < 0) & (0 < u))
    
    j_nonzero_idx <- which(b_zero_idx == 0)
    j_nonzero <- j_samp[j_nonzero_idx]
    b_nonzero[[iter]] <- j_nonzero_idx
    tp_rate[iter] <- length(j_nonzero[j_nonzero %in% p_idx])/length(p_idx)
    fp_rate[iter] <- length(j_nonzero[!j_nonzero %in% p_idx])/length(p_idx)
    b_nonzero_count[iter] <- length(b) - sum(b_zero_idx)
    beta_samp_list[[iter]] <- b_samp
    
    acc_vec[iter] <- s$gibbs$accuracy
    
    beta_list[[iter]] <- b
    beta_se_list[[iter]] <- b_se
    z_idx_list[[iter]] <- s$gibbs$zpred
    p_idx_list[[iter]] <- s$j_samp
    
    iter <- iter + 1
    
  }
  
  
  m_idx <- which(acc_vec == max(acc_vec))[1]
  m_nonzero_idx <- which(b_nonzero_count == max(b_nonzero_count))
  b_m <- beta_list[[m_idx]]
  b_m_nonzero_data <- beta_samp_list[[m_nonzero_idx]]
  b_se_m <- beta_se_list[[m_idx]]
  z_m <- z_idx_list[[m_idx]]
  p_m <- p_idx_list[[m_idx]]
  p_m_nonzero <- p_idx_list[[m_nonzero_idx]]
  p_m_nonzero_signf <- b_nonzero[[m_nonzero_idx]]
  
  out <- list(beta_samp_m_nonzero = b_m_nonzero_data,
              accuracy = acc_vec,
              max_idx = m_idx,
              maxpred = acc_vec[m_idx],
              b_nonzero_count = b_nonzero_count,
              b_maxpred = b_m,
              b_se_maxpred = b_se_m,
              z_maxpred = z_m,
              p_maxpred = p_m,
              p_max_nonzero = p_m_nonzero,
              b_max_nonzero_signf = b_nonzero[[m_nonzero_idx]],
              p_list = p_idx_list,
              b_nonzero = b_nonzero,
              nonzero_truepos_rate = tp_rate,
              nonzero_falsepos_rate = fp_rate,
              tp_rate_avg = mean(tp_rate),
              fp_rate_avg = mean(fp_rate))
  
  return(out)
  
}

# # get sampling results
# samp_res <- acc_idx_fn(k, nvsamp)
# samp_z <- samp_res$z_maxpred
# 
# # TP and FP rate construction}
# j <- 1:max(z)
# zmax_idx <- lapply(j, function(i) which(samp_z == i))
# zpred_list <- lapply(j, function(i) z[zmax_idx[[i]]])
# tp <- sapply(j, function(i) sum(zpred_list[[i]] == i))
# tpr <- map2_dbl(tp, j, function(r, i) r/length(zpred_list[[i]])) %>% 
#   round(3)
# 
# # confusion matrix construction
# conf_mat <- diag(tp)
# conf_mat_r <- diag(tpr)
# for (c in j) {
#   
#   jfp <- j[j != c]
#   
#   fp <- sapply(jfp, function(i) sum(zpred_list[[c]] == i))
#   fpr <- sapply(fp, function(r) r/length(zpred_list[[c]])) %>% 
#     round(3)
#   
#   conf_mat[c, jfp] <- fp
#   conf_mat_r[c, jfp] <- fpr
#   
# }
# 
# jnames <- paste('J', j, sep = '')
# jprednames <- paste(jnames, '_hat', sep = '')
# conf_tbl <- data.frame(cbind(jprednames, conf_mat))
# colnames(conf_tbl) <- c('', jnames)
# conf_r_tbl <- data.frame(cbind(jprednames, conf_mat_r))
# colnames(conf_r_tbl) <- c('', jnames)

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


nsamp <- 100
k <- 20
# initial p only used to generate y
p <- 7
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


ngibb <- 1000
# bench <- gibbs_fn(ngibb, X, p, z, nsamp)
# bench_acc <- bench$accuracy

# cat('True X covariates indice:', p_idx,
#     '\nBenchmark ordinal value prediction accuracy:', bench_acc,
#     '\n')

# variable sampling iterations
nvsamp <- 30
p_samp <- 10

samp_res <- acc_idx_fn(k, nvsamp, 'fixed', p_samp)
maxpred_acc <- samp_res$maxpred
maxpred_vs <- samp_res$p_maxpred[samp_res$p_maxpred %in% p_idx]
nonzero_vs <- samp_res$p_max_nonzero[samp_res$p_max_nonzero %in% p_idx]
b_true <- which(samp_res$p_max_nonzero %in% p_idx)

cat('Fixed sampling', p_samp, 'covariates from', k, 'possible',
    '\nX observations:', nsamp,
    '\nGibbs sampling chain length:', ngibb,
    '\nNumber of variable selection simulations:', nvsamp,
    '\nTrue X_i covariates:', p_idx,
    '\nTrue variable selection from max prediction accuracy:', maxpred_vs,
    '\nTrue variable selection from nonzero Gibbs beta 95%CI:', nonzero_vs,
    '\nNumber of true covariates selected from max prediction:', length(maxpred_vs),
    '\nNumber of true covariates selected from nonzero CI:', length(nonzero_vs),
    '\nMax prediction accuracy:', maxpred_acc,
    '\nAverage variable selection true/false positive rates:',
    '\n     true positive:', samp_res$tp_rate_avg,
    '\n     false positive:', samp_res$fp_rate_avg)



# Param plots -------------------------------------------------------------


plot_fun <- function(idx) {
  
  p <- ggplot(data = beta_chain, aes(beta_chain[, 1], beta_chain[, idx])) +
    geom_line() +
    xlab(colnames(beta_chain)[1]) +
    ylab(colnames(beta_chain)[idx]) +
    theme_bw()
  
  return(p)
  
}

beta_chain <- data.frame(cbind(1:ngibb, samp_res$beta_samp_m_nonzero))
colnames(beta_chain) <- c('iter',
                          paste('B', 1:p_samp, sep = ''))
colnames(beta_chain)[!names(beta_chain) %in%
                       'iter'][b_true] <- paste('B', b_true, '_true', sep = '')

b_plots <- lapply(2:dim(beta_chain)[2],
                  function(i)
                    plot_fun(i))

grid.arrange(grobs = b_plots,
             top = 'B Convergence (_true indicates true X_i)',
             nrow = 5)
