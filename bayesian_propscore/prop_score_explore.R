# Propensity scrore
remove(list = objects())
library(tidyverse)


# data --------------------------------------------------------------------

set.seed(7)
N <- 100
k <- 5

X <- matrix(rnorm(N*k), ncol = k)

t_B <- seq(.5, by = .5, length.out = k) %>% 
  matrix()
t <- X%*%t_B + rnorm(N, sd = 2)
tperc_vals <- c(.2, .4, .6, .8)
t_quin <- quantile(t, tperc_vals)

y_a <- 3 
y <- t*y_a + rnorm(N)


# helper funs -------------------------------------------------------------

# get dose effect values and treatment block values
dose_block_vals <- function(t_quin, y) {
  
  dval <- c()
  tidx <- list()
  for (i in 1:(length(t_quin) + 1)) {
    
    if (i == 1) {
      
      idx <- which(t <= t_quin[i])
      tidx[[i]] <- idx
      dval[i] <- mean(y[idx])
      
    }
    
    if (i > 1 & i < length(t_quin) + 1) {
      
      idx <- which(t <= t_quin[i] & t > t_quin[i - 1])
      tidx[[i]] <- idx
      dval[i] <- mean(y[idx])
      
    }
    
    if (i == length(t_quin) + 1) {
      
      idx <- which(t > t_quin[i - 1])
      tidx[[i]] <- idx
      dval[i] <- mean(y[idx])
      
    }
    
  }
  
  de <- sapply(2:(length(t_quin) + 1),
               function(i) dval[i] - dval[i - 1])

  Rblock <- sapply(1:(length(t_quin) + 1),
                   function(i) median(t[tidx[[i]]]))
  
  out <- list(de = de, Rblock = Rblock)
  return(out)

}

# propensity score function
# currently normal
Rvals <- function(x) {
  
  R <- dnorm(x, mean = t_center, sd = sqrt(t_spread))
  return(R)
  
}

# true dose effect values for comparison to estimated
actual_de <- dose_block_vals(t_quin, y)[['de']]

# estimation --------------------------------------------------------------

dft <- data.frame(t = t, X)

nboot <- 100
# nboot <- 1
d_te_vec <- array(NA, c(1, length(tperc_vals), nboot))
d_te_ipw_vec <- array(NA, c(1, length(tperc_vals), nboot))
samp_prog <- seq(0, nboot, by = nboot/10)
for (b in 1:nboot) {
  
  # bootstrap data
  idx_boot <- sample(1:N, size = N, replace = TRUE)
  # idx_boot <- 1:N
  dft_boot <- dft[idx_boot, ]
  yboot <- y[idx_boot]
  t <- dft_boot$t
  
  t_quin <- quantile(t, tperc_vals)
  Rblock <- dose_block_vals(t_quin, yboot)[['Rblock']]
  
  # T | X 
  # parameters for propensity score values
  # Pr(T | X) but in continuous sense, not binary for T = {0, 1}
  t_mod <- lm(t ~ ., data = dft_boot)
  res <- residuals(t_mod)
  t_spread <- 1/(N-k+1)*sum(res^2)
  t_center <- fitted(t_mod)
  
  # y | T, R dose response function
  # condt'l expectation for outcome given treatment and propensity score
  dftr <- data.frame(y = yboot, t = t, R = Rvals(t)) %>% 
    mutate(T2 = t^2, R2 = R^2, TR = t*R)
  
  Y_tr_params <- lm(y ~ ., data = dftr) %>% 
    coefficients() %>%
    matrix()
  
  Y_tr <- lapply(1:length(Rblock),
                 function(i)
                   matrix(cbind(1, Rblock[i], Rvals(Rblock[i]), Rblock[i]^2,
                                Rvals(Rblock[i])^2, Rblock[i]*Rvals(Rblock[i])),
                          nrow = N, ncol = ncol(dftr)))
  
  # average potential treatment effect for blocked levels
  # covariate adjustment GPS estimation
  # average potential treatment effect
  d_te <- sapply(2:length(Rblock),
                 function(i)
                   mean(Y_tr[[i]]%*%Y_tr_params) - mean(Y_tr[[i - 1]]%*%Y_tr_params))
  
  
  # inverse propensity weights estimation
  W_t <- dnorm(t, mean = mean(t), sd = sd(t))/Rvals(t)
  ipw_mod_param <- lm(y ~ t, data = dftr, weights = W_t) %>% 
    coefficients() %>% 
    matrix()
  
  # expected potential treatment effect
  d_te_ipw <- sapply(2:length(Rblock),
                     function(i)
                       c(1, Rblock[i])%*%ipw_mod_param - c(1, Rblock[i - 1])%*%ipw_mod_param)
  
  d_te_vec[, , b] <- d_te
  d_te_ipw_vec[, , b] <- d_te_ipw
  
  if (b %in% samp_prog) {
    
    cat('Bootstrap iteration', b, 'complete',
        '\n')
    
  }
  
}


# bootstrap tables --------------------------------------------------------

d_te_ss <- sapply(1:length(tperc_vals),
                  function(i)
                    c(mean = mean(d_te_vec[, i, ]),
                      quantile(d_te_vec[ , i, ], .05),
                      quantile(d_te_vec[ , i, ], .95)))

d_te_ipw_ss <- sapply(1:length(tperc_vals),
                      function(i) c(mean = mean(d_te_ipw_vec[, i, ]),
                                    quantile(d_te_ipw_vec[ , i, ], .05),
                                    quantile(d_te_ipw_vec[ , i, ], .95)))

boot_te_table <- rbind(actual_de, d_te_ss, d_te_ipw_ss) %>% 
  round(3)
colnames(boot_te_table) <- c('TE(21)', 'TE(32)', 'TE(43)', 'TE(54)')
rownames(boot_te_table)[1] <- c('actual effect')





