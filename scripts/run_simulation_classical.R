
# Simulation study for classical wavelet shrinkage estimators
# - Methods: soft-thresholding with Universal threshold, FDR, CV and SURE
# - Noise standard deviation estimated via MAD of the finest-scale wavelet detail coefficients
# - DWT: Daubechies (DaubExPhase), filter.number = 10
# - Outputs: MSE, MedAE, AMSE and AMAE across test functions, sample sizes and SNRs


library(wavethresh)
library(dplyr)

sample_sizes <- c(128, 512, 1024, 2048)
SNRs <- c(1, 3, 6, 9)
M <- 200

calc_MSE <- function(est, true) {
  MSE <- mean((est - true)^2)
  return(MSE)
}

calc_MedAE <- function(est, true) {
  abs_error <- abs(true - est)
  median_error <- median(abs_error)
  return(median_error)
}

test_functions <- list(
  bumps = function(n) DJ.EX(n = n)$bumps,
  blocks = function(n) DJ.EX(n = n)$blocks,
  doppler = function(n) DJ.EX(n = n)$doppler,
  heavi = function(n) DJ.EX(n = n)$heavi
)

Classical_results_df <- data.frame(
  Method = character(),
  Function = character(),
  Sample_size = numeric(),
  SNR = numeric(),
  MSE = numeric(),
  MedAE = numeric(),
  stringsAsFactors = FALSE
)


## Simulation loop
beginning <- Sys.time()

for (func in names(test_functions)) {
  for (n in sample_sizes) {
    for (SNR in SNRs) {
      
      for (j in 1:M) {
        
        set.seed(271079 + 2024*j)
        
        v <- test_functions[[func]](n)
        sd_sinal <- sd(v)
        sigma <- sd_sinal / SNR
        e <- rnorm(n, mean = 0, sd = sigma)
        y <- v + e
        
        ## DWT
        ywd <- wd(y, filter.number = 10, family = "DaubExPhase")
        FineCoefs <- accessD(ywd, lev = (nlevelsWT(ywd) - 1))
        s <- mad(FineCoefs)
        
        ## Soft Thresholding  
        # Universal
        utDJ <- s * sqrt(2 * log(n))
        univ_soft <- threshold(ywd, levels = 1:(nlevelsWT(ywd) - 1), policy = "manual",
                               type = "soft", value = utDJ)
        
        ywruniv_soft <- wr(univ_soft)
        MSE_univ <- calc_MSE(ywruniv_soft, v)
        MedAE_univ <- calc_MedAE(ywruniv_soft, v)
        Classical_results_df <- rbind(
          Classical_results_df, 
          data.frame(Method = "Universal", Function = func, 
                     Sample_size = n, SNR = SNR, 
                     MSE = MSE_univ, MedAE = MedAE_univ)
        )
        
        # FDR
        fdr_soft <- threshold(ywd, levels = 1:(nlevelsWT(ywd) - 1), policy = "fdr", 
                              dev = madmad, type = "soft")
        
        ywrfdr_soft <- wr(fdr_soft)
        MSE_fdr <- calc_MSE(ywrfdr_soft, v)
        MedAE_fdr <- calc_MedAE(ywrfdr_soft, v)
        Classical_results_df <- rbind(
          Classical_results_df, 
          data.frame(Method = "FDR", Function = func, 
                     Sample_size = n, SNR = SNR, 
                     MSE = MSE_fdr, MedAE = MedAE_fdr)
        )
        
        # Cross validation
        ywdcvT_soft <- threshold(ywd, levels = 1:(nlevelsWT(ywd) - 1), policy = "cv", 
                                 dev = madmad, type = "soft")
        
        ywrcv_soft <- wr(ywdcvT_soft)
        MSE_cv <- calc_MSE(ywrcv_soft, v)
        MedAE_cv <- calc_MedAE(ywrcv_soft, v)
        Classical_results_df <- rbind(
          Classical_results_df, 
          data.frame(Method = "CV", Function = func, 
                     Sample_size = n, SNR = SNR, 
                     MSE = MSE_cv, MedAE = MedAE_cv)
        )
        
        # SURE
        ywdsure <- threshold(ywd, levels = 1:(nlevelsWT(ywd) - 1), policy = "sure", 
                             type = "soft")
        
        ywr_sure <- wr(ywdsure)
        MSE_sure <- calc_MSE(ywr_sure, v)
        MedAE_sure <- calc_MedAE(ywr_sure, v)
        Classical_results_df <- rbind(
          Classical_results_df, 
          data.frame(Method = "SURE", Function = func, 
                     Sample_size = n, SNR = SNR, 
                     MSE = MSE_sure, MedAE = MedAE_sure)
        )
      }
    }
  }
}

ending <- Sys.time()
total_time <- ending - beginning
print(total_time)

## AMSEs
classical_AMSEs <- Classical_results_df %>%
  group_by(Method, Function, Sample_size, SNR) %>%
  summarise(AMSE = mean(MSE), .groups = "drop")
classical_AMSEs$AMSE <- round(classical_AMSEs$AMSE, 3)

## AMAEs
classical_AMAEs <- Classical_results_df %>%
  group_by(Method, Function, Sample_size, SNR) %>%
  summarise(AMAE = mean(MedAE), .groups = "drop")
classical_AMAEs$AMAE <- round(classical_AMAEs$AMAE, 3)


get_sd_classical_MSE <- function(method, func, n, snr) {
  subset <- Classical_results_df[
    Classical_results_df$Method == method &
      Classical_results_df$Function == func &
      Classical_results_df$Sample_size == n &
      Classical_results_df$SNR == snr,
  ]
  sd(subset$MSE)
}

get_sd_classical_MedAE <- function(method, func, n, snr) {
  subset <- Classical_results_df[
    Classical_results_df$Method == method &
      Classical_results_df$Function == func &
      Classical_results_df$Sample_size == n &
      Classical_results_df$SNR == snr,
  ]
  sd(subset$MedAE)
}
 # Example
round(get_sd_classical_MedAE("Universal", "heavi", 1024, 3), 3)


