# Logistic variable selection comparison

remove(list = objects())
library(tidyverse)
# RF VIM functions for party::cforest object

# Data --------------------------------------------------------------------

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

# var select sim ----------------------------------------------------------

nvsamp <- 100

tp_rate_max <- c()
fp_rate_max <- c()

tp_rate_min <- c()
fp_rate_min <- c()

tp_rate_qr <- c()
fp_rate_qr <- c()

tp_rate_rf <- c()
fp_rate_rf <- c()

iter <- 1
samp_prog <- seq(0, nvsamp, by = nvsamp/10)
while (iter <= nvsamp) {
  
  # number of aggregate covariates
  k <- 20
  k_vec <- 1:k
  # number of covariates to generate y
  p <- 10
  # X sample size
  nsamp <- 100
  # aggregate X data matrix
  X <- matrix(rnorm(nsamp*k), nrow = nsamp, ncol = k)
  # random sample covariates to get initial X covariates
  p_idx <- sample(k_vec, p)
  # create p covariate X data matrix, generate B and y
  X_p <- X[, p_idx]
  B_p <- matrix(rnorm(p, sd = 3))
  y <- rnorm(1, sd = 10) + X_p%*%B_p + rnorm(nsamp)
  # create observed ordinal z
  z <- get_z(y)
  
  # data construction for model function
  # convert z to binary
  X_mod <- X %>%
    as_tibble(.name_repair = ~ paste('x', seq(1, dim(X)[2], 1), sep = ''))
  z_mod_max <- ifelse(z < max(z), 0, 1)
  z_mod_min <- ifelse(z == min(z), 0, 1)
  df_rf <- data.frame(cbind(factor(z, ordered = TRUE), X_mod))%>% 
    rename(z = factor.z..ordered...TRUE.)
  
  
  lhs <- paste(names(X_mod), collapse = '+')
  fx_max <- paste('z_mod_max', lhs, sep = '~')
  fx_min <- paste('z_mod_min', lhs, sep = '~')
  fx_z <- paste('z', lhs, sep = '~')
  
  # no intercept lhs specification
  # fx_noint <- paste(fx_terms, '-1', sep = '')
  
  ### logistic binary spec - overdispersion ###
  mod_zmax <- glm(as.formula(fx_max),
                  control = list(epsilon = 1e-5, maxit = 100),
                  family = 'quasibinomial',
                  data = cbind(z_mod_max, X_mod))
  
  mod_zmin <- glm(as.formula(fx_min),
                  control = list(epsilon = 1e-5, maxit = 100),
                  family = 'quasibinomial',
                  data = cbind(z_mod_min, X_mod))
  
  modzmax_xvec <- which(summary(mod_zmax)$coefficients[2:dim(X_mod)[2], 4] < .05) %>% 
    unname()
  
  modzmin_xvec <- which(summary(mod_zmin)$coefficients[2:dim(X_mod)[2], 4] < .05) %>% 
    unname()
  
  # tp and fp rates for max and min specifications
  tp_rate_max[iter] <- sum(modzmax_xvec %in% p_idx)/p
  fp_rate_max[iter] <- sum(!modzmax_xvec %in% p_idx)/p
  
  tp_rate_min[iter] <- sum(modzmin_xvec %in% p_idx)/p
  fp_rate_min[iter] <- sum(!modzmin_xvec %in% p_idx)/p
  
  ### proportional odds model ###
  
  # Hessian singular...
  
  # MASS::polr(as.formula(fx_ord), data = df_ord, method = 'logistic')
  # ordinal::clm(as.formula(fx_ord), data = df_ord, link = 'logit')
  
  ### quantile regression evaluates threshold at 2 => median ###
  
  mod_qr <- quantreg::rq(as.formula(fx_z), data = cbind(z, X_mod), tau = 0.5)
  qr_bounds <- summary(mod_qr)$coefficients[, 2:3]
  # provides var intervals 1 contain 0, 0 not
  modqr_inter <- map2_dbl(qr_bounds[-1, 1], qr_bounds[-1, 2],
                          function(l, u) (l < 0) & (u > 0))
  modqr_xvec <- which(modqr_inter == 0) %>%
    unname()
  
  tp_rate_qr[iter] <- sum(modqr_xvec %in% p_idx)/p
  fp_rate_qr[iter] <- sum(!modqr_xvec %in% p_idx)/p
    
  ### random forest ###
  
  # creates random forest
  mod_rf <- party::cforest(as.formula(fx_z), data = df_rf)
                           # controls = cforest_control(ntree = 1000))
  # error metric, OOB diff bw pre and post permutation average
  # covariates are not highly correlated
  # conditional importance not crucial
  # want < 0 importance measures
  #  implies permuted effect on acc > non permuted effect
  #  response values depend on covariate
  rf_vim <- party::varimp(mod_rf)
  rf_xvec <- which(rf_vim < 0) %>%
    unname()
  
  tp_rate_rf[iter] <- sum(rf_xvec %in% p_idx)/p
  fp_rate_rf[iter] <- sum(!rf_xvec %in% p_idx)/p
  
  if (iter %in% samp_prog) {
    
    cat('Iteration', iter, 'complete',
        '\n')
    
  }
  
  iter <- iter + 1
  
}

cat('(true, false) positive rate means',
    '\n   logistic max(z) == 1: (', c(mean(tp_rate_max), mean(fp_rate_max)), ')',
    '\n   logistic z > 1 == 1: (', c(mean(tp_rate_min), mean(fp_rate_min)), ')',
    '\n   quant reg at median: (', c(mean(tp_rate_qr), mean(fp_rate_qr)), ')',
    '\n   RF var importance: (', c(mean(tp_rate_rf), mean(fp_rate_rf)), ')',
    '\n')

# rate summary stats results from varselectX_gibbs.R
tp_rate_gibbs <- c(0.925, 0.700, 1.000, 0.085)
fp_rate_gibbs <- c(0.445, 0.100, 0.900, 0.204)

tp_list <- list(tp_rate_max, tp_rate_min, tp_rate_qr, tp_rate_rf)
fp_list <- list(fp_rate_max, fp_rate_min, fp_rate_qr, fp_rate_rf)

tp_ss <- sapply(tp_list,
                function(i) c(mean(i), min(i), max(i), sd(i)))
tp_ss_full <- cbind(tp_rate_gibbs, tp_ss) %>% unname() %>% round(3)

fp_ss <- sapply(fp_list,
                function(i) c(mean(i), min(i), max(i), sd(i)))
fp_ss_full <- cbind(fp_rate_gibbs, fp_ss) %>% unname() %>% round(3)

ss_table <- rbind(tp_ss_full, fp_ss_full)
colnames(ss_table) <- c('Gibbs', 'Logistic (max(z))', 'Logistic (min(z))',
                        'Quantile Reg', 'Random Forest')
row.names(ss_table) <- rep(c('mean', 'min', 'max', 'sd'), 2)
