
# Simulation study for Bayesian wavelet shrinkage estimators
# - Priors: Raised Cosine and symmetric Beta (a,b)
# - Shrinkage via spike-and-slab with level-dependent alpha(j)
# - Numerical integration by deterministic grid-based quadrature
# - Outputs: MSEs, MedAEs, AMSEs and AMAEs across test functions, sample sizes and SNRs


library(wavethresh)
library(dplyr)

sample_sizes <- c(128, 512, 1024, 2048)
SNRs <- c(1, 3, 6, 9)
M <- 200

calc_MSE <- function(est, true) {
  mean((est - true)^2)
}

calc_MedAE <- function(est, true) {
  abs_error <- abs(true - est)
  median(abs_error)
}

alpha_level <- function(j, J0, gamma = 2) {
  if (j < J0) {
    return(0)
  } else {
    return(1 - 1 / ((j - J0 + 1)^gamma))
  }
}

#--------------------------------------
# Raised Cosine and Beta distributions
#--------------------------------------
rc_dist <- function(x, params) {
  ifelse(
    x >= (params$mu - params$tau) & x <= (params$mu + params$tau),
    1 / (2 * params$tau) * (1 + cos((x - params$mu) * pi / params$tau)),
    0
  )
}

beta_dist <- function(x, params) {
  ifelse(
    abs(x) <= params$m,
    (x + params$m)^(params$a - 1) * (params$m - x)^(params$b - 1) /
      ((2 * params$m)^(params$a + params$b - 1) * beta(params$a, params$b)),
    0
  )
}

#-------------------------------------------------------------------
# General shrinkage function using grid-based numerical quadrature
#-------------------------------------------------------------------
shrinkage_grid <- function(d, alpha, params_dist, func_dist, sigma, n_grid = 200) {
  n <- length(d)
  shrink <- numeric(n)
  
  if (identical(func_dist, rc_dist)) {
    lower <- -params_dist$tau
    upper <-  params_dist$tau
  } else if (identical(func_dist, beta_dist)) {
    lower <- -params_dist$m
    upper <-  params_dist$m
  } else {
    stop("unknown density function in shrinkage_grid")
  }
  
  theta_grid <- seq(lower, upper, length.out = n_grid)
  delta <- (upper - lower) / (n_grid - 1)
  
  prior_vals <- func_dist(theta_grid, params_dist)
  
  for (i in seq_len(n)) {
    lik_vals <- dnorm(d[i], mean = theta_grid, sd = sigma)
    
    num_int <- sum(theta_grid * prior_vals * lik_vals) * delta
    den_int <- sum(prior_vals * lik_vals) * delta
    
    num <- (1 - alpha) * num_int
    den <- alpha * dnorm(d[i], mean = 0, sd = sigma) + (1 - alpha) * den_int
    
    shrink[i] <- num / den
  }
  
  shrink
}

#-------------------------------------------------------
# Function to apply shrinkage to the DWT coefficients
#-------------------------------------------------------
wd_shrink <- function(y, params_dist, func_dist, J0 = 1, gamma = 2, n_grid = 200) {
  
  ywd <- wd(y, filter.number = 10, family = "DaubExPhase")
  
  # detail coefficients at each level
  getD_level <- function(i) {
    accessD(ywd, level = i)
  }
  
  # empirical coefficients from resolution level 1 onward
  coef_emp <- c(unlist(lapply(1:(nlevelsWT(ywd) - 1), getD_level)))
  
  # estimate of m (or tau, for the RC prior)
  m_hat <- max(abs(coef_emp))
  
  # estimate of sigma_hat from the finest-level detail coefficients
  FineCoefs <- accessD(ywd, level = (nlevelsWT(ywd) - 1))
  sigma_hat <- mad(FineCoefs)
  
  # update the distribution parameters
  params_dist$m <- m_hat
  params_dist$tau <- m_hat
  
  # apply shrinkage to the detail coefficients, with alpha(j)
  for (j in 1:(nlevelsWT(ywd) - 1)) {
    D_j <- getD_level(j)
    alpha_j <- alpha_level(j, J0 = J0, gamma = gamma)
    D_j_shrink <- shrinkage_grid(D_j, alpha = alpha_j,
                                 params_dist, func_dist,
                                 sigma = sigma_hat, n_grid = n_grid)
    ywd <- putD(ywd, level = j, v = D_j_shrink)
  }
  
  ywd
}

test_functions <- list(
  bumps = function(n) DJ.EX(n = n)$bumps,
  blocks = function(n) DJ.EX(n = n)$blocks,
  doppler = function(n) DJ.EX(n = n)$doppler,
  heavi = function(n) DJ.EX(n = n)$heavi
)

Results_df <- data.frame(
  Method = character(),
  Function = character(),
  Sample_size = numeric(),
  SNR = numeric(),
  MSE = numeric(),
  MedAE = numeric(),
  stringsAsFactors = FALSE
)

## Prior parameters (m/tau will be estimated inside wd_shrink)
params_rc <- list(mu = 0, tau = NA)
params_beta_1 <- list(a = 1, b = 1, m = NA)
params_beta_5 <- list(a = 5, b = 5, m = NA)


#---------------------
# Simulation loop
#---------------------
J0_use <- 1
gamma_use <- 2  

beginning <- Sys.time()

for (func in names(test_functions)) {
  for (n in sample_sizes) {
    for (SNR in SNRs) {
      for (i in 1:M) {
        set.seed(271079 + 2024*i)
        
        v <- test_functions[[func]](n)
        sd_sinal <- sd(v)
        sig <- sd_sinal / SNR
        e <- rnorm(n, mean = 0, sd = sig)
        y <- v + e
        
        ## RC
        ywd_shrink_rc <- wd_shrink(y,
                                   params_dist = params_rc,
                                   func_dist   = rc_dist,
                                   J0 = J0_use,
                                   gamma = gamma_use,
                                   n_grid = 200)
        ywr_rc <- wr(ywd_shrink_rc)
        MSE_rc <- calc_MSE(ywr_rc, v)
        MedAE_rc <- calc_MedAE(ywr_rc, v)
        Results_df <- rbind(
          Results_df,
          data.frame(Method = "RC", Function = func,
                     Sample_size = n, SNR = SNR,
                     MSE = MSE_rc, MedAE = MedAE_rc)
        )
        
        ## Beta (a=1,b=1)
        ywd_shrink_beta_1 <- wd_shrink(y,
                                       params_dist = params_beta_1,
                                       func_dist   = beta_dist,
                                       J0 = J0_use,
                                       gamma = gamma_use,
                                       n_grid = 200)
        ywr_beta_1 <- wr(ywd_shrink_beta_1)
        MSE_beta_1 <- calc_MSE(ywr_beta_1, v)
        MedAE_beta_1 <- calc_MedAE(ywr_beta_1, v)
        Results_df <- rbind(
          Results_df,
          data.frame(Method = "Beta_1", Function = func,
                     Sample_size = n, SNR = SNR,
                     MSE = MSE_beta_1, MedAE = MedAE_beta_1)
        )
        
        ## Beta (a=5,b=5)
        ywd_shrink_beta_5 <- wd_shrink(y,
                                       params_dist = params_beta_5,
                                       func_dist   = beta_dist,
                                       J0 = J0_use,
                                       gamma = gamma_use,
                                       n_grid = 200)
        ywr_beta_5 <- wr(ywd_shrink_beta_5)
        MSE_beta_5 <- calc_MSE(ywr_beta_5, v)
        MedAE_beta_5 <- calc_MedAE(ywr_beta_5, v)
        Results_df <- rbind(
          Results_df,
          data.frame(Method = "Beta_5", Function = func,
                     Sample_size = n, SNR = SNR,
                     MSE = MSE_beta_5, MedAE = MedAE_beta_5)
        )
      }
    }
  }
}

ending <- Sys.time()
total_time <- ending - beginning
print(total_time)

## AMSEs
AMSEs <- Results_df %>%
  group_by(Method, Function, Sample_size, SNR) %>%
  summarise(AMSE = mean(MSE), .groups = "drop")
AMSEs$AMSE <- round(AMSEs$AMSE, 3)

## AMAEs
AMAEs <- Results_df %>%
  group_by(Method, Function, Sample_size, SNR) %>%
  summarise(AMAE = mean(MedAE), .groups = "drop")
AMAEs$AMAE <- round(AMAEs$AMAE, 3)


get_sd_MSE <- function(method, func, n, snr) {
  subset <- Results_df[
    Results_df$Method == method &
      Results_df$Function == func &
      Results_df$Sample_size == n &
      Results_df$SNR == snr,
  ]
  sd(subset$MSE)
}

get_sd_MedAE <- function(method, func, n, snr) {
  subset <- Results_df[
    Results_df$Method == method &
      Results_df$Function == func &
      Results_df$Sample_size == n &
      Results_df$SNR == snr,
  ]
  sd(subset$MedAE)
}

