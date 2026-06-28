###############################################
### FUNCTION TO RUN CCA FOR THE PARALLEL MODEL
###############################################
run_cca <- function(miss_list) {
  
  first_cc <- na.omit(miss_list[[1]])
  
  tmp_fit <- vglm(affairs_mod ~ children + religiousness + rating,
                  data = first_cc,
                  family = cumulative(parallel = TRUE))
  
  tmp_coeffs     <- coef(tmp_fit)           
  tmp_intercepts <- tmp_fit@misc$zeta        
  
  param_names <- c(names(tmp_coeffs), names(tmp_intercepts))
  K <- length(param_names)
  
  betas <- matrix(NA_real_, nrow = M, ncol = K)
  SEs   <- matrix(NA_real_, nrow = M, ncol = K)
  
  colnames(betas) <- param_names
  colnames(SEs)   <- param_names
  
  ###################################
  ### LOOP THROUGH M DATASETS
  ###################################
  for (i in 1:M) {
    
    dat_cc <- na.omit(miss_list[[i]])
    
    fit <- vglm(affairs_mod ~ children + religiousness + rating,
                data = dat_cc,
                family = cumulative(parallel = TRUE))
    
    coeffs     <- coef(fit)
    intercepts <- fit@misc$zeta
    all_params <- c(coeffs, intercepts)
    
    sum_fit    <- summary(fit)
    coef_table <- sum_fit@coef3
    
    SE_slopes <- coef_table[names(coeffs),     "Std. Error"]
    SE_inter  <- coef_table[names(intercepts), "Std. Error"]
    
    SE_vec <- c(SE_slopes, SE_inter)
    
    betas[i, ] <- all_params
    SEs[i, ]   <- SE_vec
    
    print(i)
  }
  
  list(betas = betas, SEs = SEs)
}


###############################################
### 1. RUN CCA FOR 5%, 15%, 30%
###############################################

res5  <- run_cca(NA5_list)
res15 <- run_cca(NA15_list)
res30 <- run_cca(NA30_list)


###############################################
### 2. FUNCTION TO CREATE LONG TABLE
###############################################

make_full_table <- function(betas, SEs) {
  
  params <- colnames(betas)
  
  out_list <- lapply(seq_along(params), function(j) {
    data.frame(
      iteration = 1:nrow(betas),
      parameter = params[j],
      beta      = betas[, j],
      SE        = SEs[, j]
    )
  })
  
  df <- do.call(rbind, out_list)
  
  df$poly_type <- dplyr::case_when(
    grepl("\\.L$", df$parameter) ~ "L",
    grepl("\\.Q$", df$parameter) ~ "Q",
    grepl("\\.C$", df$parameter) ~ "C",
    grepl("\\^4$", df$parameter) ~ "4",
    TRUE ~ "none"
  )
  
  df
}

###############################################
### 3. CREATE TABLES for raw data estimates
###############################################

table5_full_CCA  <- make_full_table(res5$betas,  res5$SEs)
table15_full_CCA <- make_full_table(res15$betas, res15$SEs)
table30_full_CCA <- make_full_table(res30$betas, res30$SEs)

write.csv(table5_full_CCA, "table5_full_CCA.csv", row.names = FALSE)
write.csv(table15_full_CCA, "table15_full_CCA.csv", row.names = FALSE)
write.csv(table30_full_CCA, "table30_full_CCA.csv", row.names = FALSE)


###############################################
### 4. CREATE TABLES for comparisons to True model
###############################################
## final model results
true_fit <- vglm(
  affairs_mod ~ children + religiousness + rating,
  data = data,
  family = cumulative(parallel = TRUE)
)

true_betas <- coef(true_fit)
true_intercepts <- true_fit@misc$zeta

true_params <- c(true_betas, true_intercepts)

true_SE_table <- summary(true_fit)@coef3
true_SEs <- true_SE_table[names(true_params), "Std. Error"]

# function for table creation
make_eval_table <- function(res, true_params, true_SEs) {
  
  betas <- res$betas
  SEs   <- res$SEs
  
  param_names <- colnames(betas)
  M <- nrow(betas)
  
  mean_beta <- colMeans(betas)
  true_beta <- true_params[param_names]
  
  out <- data.frame(
    parameter      = param_names,
    true_beta      = true_beta,
    mean_beta      = mean_beta,
    
    # absolute + relative bias
    abs_bias       = abs(mean_beta - true_beta),
    rel_bias       = (mean_beta - true_beta) / true_beta,
    
    # original bias
    bias           = mean_beta - true_beta,
    
    emp_sd         = apply(betas, 2, sd),
    true_se        = true_SEs[param_names],
    mean_se        = colMeans(SEs),
    coverage       = NA_real_
  )
  
  # Coverage calculation
  for (j in seq_along(param_names)) {
    beta_j <- betas[, j]
    se_j   <- SEs[, j]
    lower  <- beta_j - 1.96 * se_j
    upper  <- beta_j + 1.96 * se_j
    
    out$coverage[j] <- mean(
      out$true_beta[j] >= lower & out$true_beta[j] <= upper,
      na.rm = TRUE
    )
  }
  
  out
}

table5_eval  <- make_eval_table(res5,  true_params, true_SEs)
table15_eval <- make_eval_table(res15, true_params, true_SEs)
table30_eval <- make_eval_table(res30, true_params, true_SEs)

write.csv(table5_eval,  "table5_eval_CCA.csv",  row.names = FALSE)
write.csv(table15_eval, "table15_eval_CCA.csv", row.names = FALSE)
write.csv(table30_eval, "table30_eval_CCA.csv", row.names = FALSE)





