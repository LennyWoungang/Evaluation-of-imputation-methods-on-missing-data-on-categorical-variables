######################################### Model Selection
library(AER)
library(dplyr)
library(MASS)
library(glmnet)
library(nnet)
library(ordinalNet)
library(rms)
library(AER)
library(cat)
library(VGAM)
library(dplyr)

data("Affairs")
data <- Affairs

data <- data %>% 
  mutate(
    affairs_mod = case_when(
      affairs == 0 ~ 0,
      affairs %in% c(1, 2, 3) ~ 1,
      affairs %in% c(7, 12) ~ 2
    ),
    affairs_mod = factor(affairs_mod,
                         levels = 0:2,
                         labels = c("0 times", "1-3", "3+"),
                         ordered = T),
    
    religiousness = ordered(religiousness),
    rating = ordered(rating),
    education = factor(education),
    occupation = factor(occupation),
    yearsmarried = factor(yearsmarried),
    age = factor(age)
  ) %>% 
  dplyr::select(affairs_mod, gender, children, religiousness, rating, education, occupation, yearsmarried, age)

data_mod <- data


##### vglm function
model_F <- vglm(affairs_mod~ (gender + children + religiousness  + rating)^2, data = data, 
                family = cumulative(parallel = F ~ religiousness + ratin)) 
model_T <- vglm(affairs_mod~ (gender + children + religiousness  + rating)^2, data = data, 
                family = cumulative(parallel = T)) 

## test d'adequation du modele des cotes proportionnelles
L2= deviance(model_T) - deviance(model_F)
df= df.residual(model_T) - df.residual(model_F)
1-pchisq(L2, df)

##Model selection
model <- vglm(affairs_mod~ (gender + children + religiousness + rating)^2, data = data, 
              family = cumulative(parallel = T)) 
step_mod <- step4vglm(model, direction = "backward", trace = TRUE)
formula(step_mod)
Final_model <- vglm(affairs_mod~ children + religiousness + rating, data = data, 
                    family = cumulative(parallel = T)) 

## data creation for misssingness ##
data_for_missing <- data %>% 
  dplyr::select(affairs_mod, children, religiousness, rating)


## Check assumptions of vgam
check_hyp <- lrm(affairs_mod ~ children + religiousness + rating,
                 data = data,
                 x=T, y=T)

par(mfrow = c(3, 3))
residuals(check_hyp, type = "partial", pl = TRUE)









################################ Missing Data mechanism
# row-wise
create_NA <- function(data, p_row) {
  n <- nrow(data)
  p <- round(p_row * n)         
  cols <- colnames(data)
  rows_to_miss <- sample.int(n, p, replace = FALSE)
  out <- data
  for (r in rows_to_miss) {
    c <- sample(cols, 1)
    out[r, c] <- NA
  }
  out
}


### For different percentages of missingness
set.seed(6516)
M <- 500
NA5_list <- vector("list", M)
for(i in 1:M){
  NA5 <- create_NA(data = data_for_missing, p = 0.05)
  NA5_list[[i]] <- NA5
}


NA15_list <- vector("list", M)
for(i in 1:M){
  NA15 <- create_NA(data = data_for_missing, p = 0.15)
  NA15_list[[i]] <- NA15
}

NA30_list <- vector("list", M)
for(i in 1:M){
  NA30 <- create_NA(data = data_for_missing, p = 0.30)
  NA30_list[[i]] <- NA30
}

# check for missing data
count_missing_rows <- function(df) {
  sum(!complete.cases(df))
}
summary_df <- data.frame(
  five  = as.list(summary(sapply(NA5_list, count_missing_rows))),
  fifteen = as.list(summary(sapply(NA15_list, count_missing_rows))),
  thirty  = as.list(summary(sapply(NA30_list, count_missing_rows)))
)
count_total_NA <- function(df) sum(is.na(df))
total_NA_5  <- sapply(NA5_list, count_total_NA)
total_NA_15 <- sapply(NA15_list, count_total_NA)
total_NA_30 <- sapply(NA30_list, count_total_NA)









############################# FUNCTION TO RUN CCA FOR THE PARALLEL MODEL
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
  
  mean_beta <- colMeans(betas)
  true_beta <- true_params[param_names]
  
  out <- data.frame(
    parameter      = param_names,
    true_beta      = true_beta,
    mean_beta      = mean_beta,
    
    # absolute + relative bias of coefficients
    abs_bias       = abs(mean_beta - true_beta),
    rel_bias       = (mean_beta - true_beta) / true_beta,
    
    # signed bias
    bias           = mean_beta - true_beta,
    
    # empirical SD of beta estimates
    emp_sd         = apply(betas, 2, sd),
    
    # true model SE
    true_se        = true_SEs[param_names],
    
    # mean estimated SE across simulations
    mean_se        = colMeans(SEs),
    
    # empirical SD of SE estimates
    emp_sd_se      = apply(SEs, 2, sd),
    
    # SE bias
    se_bias        = colMeans(SEs) - true_SEs[param_names],
    
    # NEW: absolute SE bias
    abs_se_bias    = abs(colMeans(SEs) - true_SEs[param_names]),
    
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

table5_eval$rel_se_bias <- (table5_eval$mean_se - table5_eval$true_se) / table5_eval$true_se
table15_eval$rel_se_bias <- (table15_eval$mean_se - table15_eval$true_se) / table15_eval$true_se
table30_eval$rel_se_bias <- (table30_eval$mean_se - table30_eval$true_se) / table30_eval$true_se

write.csv(table5_eval,  "table5_eval_CCA.csv",  row.names = FALSE)
write.csv(table15_eval, "table15_eval_CCA.csv", row.names = FALSE)
write.csv(table30_eval, "table30_eval_CCA.csv", row.names = FALSE)












#################################### EM

############################################
### 1. Add extra variables to each dataset
############################################
cols_to_add <- Affairs[, c("age", "yearsmarried", "education", "occupation","gender")]

add_extra_vars <- function(df) {
  cbind(df, cols_to_add)
}

NA5_list_mod  <- lapply(NA5_list,  add_extra_vars)
NA15_list_mod <- lapply(NA15_list, add_extra_vars)
NA30_list_mod <- lapply(NA30_list, add_extra_vars)

############################################
### 2. Convert all variables to categorical integers
############################################
make_all_cat <- function(df) {
  df[] <- lapply(df, function(x) as.integer(factor(x)))
  df
}

############################################
### 3. EM imputation
############################################
######
run_EM_only <- function(df_list) {
  
  imputed_list <- vector("list", length(df_list))
  
  for (i in seq_along(df_list)) {
    cat("Running EM on dataset", i, "\n")
    
    df <- df_list[[i]]
    
    # convert everything to categorical integers for EM
    df_cat <- df
    df_cat[] <- lapply(df_cat, function(x) as.integer(factor(x)))
    
    M <- as.matrix(df_cat)
    pre <- prelim.cat(M)
    em_out <- em.cat(pre, maxits = 2000, eps = 1e-8)
    
    df_imp <- as.data.frame(imp.cat(pre, em_out))
    colnames(df_imp) <- colnames(df_cat)
    
    # store imputed dataset
    imputed_list[[i]] <- df_imp
  }
  
  imputed_list
}

EM5_imputed  <- run_EM_only(NA5_list_mod)
EM15_imputed <- run_EM_only(NA15_list_mod)
EM30_imputed <- run_EM_only(NA30_list_mod)

############################################
### 4. Model fit
############################################
fit_model_to_imputed <- function(imputed_list) {
  
  all_results <- vector("list", length(imputed_list))
  
  for (i in seq_along(imputed_list)) {
    cat("Fitting model on imputed dataset", i, "\n")
    
    df <- imputed_list[[i]]
  
    df$affairs_mod <- ordered(df$affairs)
    df$children   <- factor(df$children)  
    df$religiousness <- ordered(df$religiousness)
    df$rating <- ordered(df$rating)

    # Fit model
    fit <- vglm(
      affairs_mod ~ children + religiousness + rating,
      data = df,
      family = cumulative(parallel = TRUE)
    )
    
    # Extract coefficients
    est <- c(coef(fit), fit@misc$zeta)
    se_tab <- summary(fit)@coef3
    se <- se_tab[names(est), "Std. Error"]
    
    all_results[[i]] <- tibble(
      dataset  = i,
      term     = names(est),
      estimate = est,
      se       = se
    )
  }
  bind_rows(all_results)
}

EM5_results  <- fit_model_to_imputed(EM5_imputed)
EM15_results <- fit_model_to_imputed(EM15_imputed)
EM30_results <- fit_model_to_imputed(EM30_imputed)

############################################################
### 5.Clean results
############################################################
rename <- function(df) {
  df %>%
    mutate(term = ifelse(term == "children2", "childrenyes", term))
}

EM5_results  <- rename(EM5_results)
EM15_results <- rename(EM15_results)
EM30_results <- rename(EM30_results)

write.csv(EM5_results,  "table5_full_EM.csv",  row.names = FALSE)
write.csv(EM15_results, "table15_full_EM.csv", row.names = FALSE)
write.csv(EM30_results, "table30_full_EM.csv", row.names = FALSE)

reorder_EM_results <- function(EM_res) {
  df <- EM_res %>%
    dplyr::rename(
      iteration = dataset,
      parameter = term,
      beta      = estimate,
      SE        = se
    )
  df$iteration <- as.numeric(df$iteration)
  param_order <- c(
    "(Intercept):1",
    "(Intercept):2",
    "childrenyes",
    "religiousness.L",
    "religiousness.Q",
    "religiousness.C",
    "religiousness^4",
    "rating.L",
    "rating.Q",
    "rating.C",
    "rating^4"
  )
  param_order <- param_order[param_order %in% unique(df$parameter)]
  df2 <- df %>%
    dplyr::mutate(parameter = factor(parameter, levels = param_order)) %>%
    dplyr::arrange(parameter, iteration)
  df2$poly_type <- "none"
  df2
}

EM5_clean  <- reorder_EM_results(EM5_results)
EM15_clean <- reorder_EM_results(EM15_results)
EM30_clean <- reorder_EM_results(EM30_results)

prepare_EM_for_eval <- function(EM_df) {
  param_order <- EM_df$parameter[EM_df$iteration == 1]
  M <- max(EM_df$iteration)
  betas <- matrix(NA, nrow = M, ncol = length(param_order))
  SEs   <- matrix(NA, nrow = M, ncol = length(param_order))
  
  colnames(betas) <- param_order
  colnames(SEs)   <- param_order
  for (i in 1:M) {
    tmp <- EM_df[EM_df$iteration == i, ]
    
    tmp <- tmp[match(param_order, tmp$parameter), ]  # reorder
    
    betas[i, ] <- tmp$beta
    SEs[i, ]   <- tmp$SE
  }
  list(betas = betas, SEs = SEs)
}

resEM5  <- prepare_EM_for_eval(EM5_clean)
resEM15 <- prepare_EM_for_eval(EM15_clean)
resEM30 <- prepare_EM_for_eval(EM30_clean)

###############################################
### 6. CREATE TABLES for comparisons to True model
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

table5_eval_EM  <- make_eval_table(resEM5, true_params, true_SEs)
table15_eval_EM <- make_eval_table(resEM15, true_params, true_SEs)
table30_eval_EM <- make_eval_table(resEM30, true_params, true_SEs)

table5_eval_EM$rel_se_bias <- (table5_eval_EM$mean_se - table5_eval_EM$true_se) / table5_eval_EM$true_se
table15_eval_EM$rel_se_bias <- (table15_eval_EM$mean_se - table15_eval_EM$true_se) / table15_eval_EM$true_se
table30_eval_EM$rel_se_bias <- (table30_eval_EM$mean_se - table30_eval_EM$true_se) / table30_eval_EM$true_se

write.csv(table5_eval_EM,  "table5_eval_EM.csv",  row.names = FALSE)
write.csv(table15_eval_EM, "table15_eval_EM.csv", row.names = FALSE)
write.csv(table30_eval_EM, "table30_eval_EM.csv", row.names = FALSE)

################################## plots 
## load data

table5_eval_CCA <- read.csv("table5_eval_CCA.csv")
table15_eval_CCA <- read.csv("table15_eval_CCA.csv")
table30_eval_CCA <- read.csv("table30_eval_CCA.csv")
table5_eval_EM <- read.csv("table5_eval_EM.csv")
table15_eval_EM <- read.csv("table15_eval_EM.csv")
table30_eval_EM <- read.csv("table30_eval_EM.csv")
table5_eval_MI <- read.csv("table5_eval_MI.csv")
table15_eval_MI <- read.csv("table15_eval_MI.csv")
table30_eval_MI <- read.csv("table30_eval_MI.csv")

##### bias plots
plot_bias_with_cov <- function(df, title = "Bias Plot") {
  
  df$coverage_label <- ifelse(df$coverage != 1,
                              paste0("cov=", df$coverage),
                              "")
  
  df$label_y <- pmax(df$bias + df$emp_sd,
                     df$bias - df$emp_sd) + 0.01
  
  ggplot(df, aes(x = parameter, y = bias)) +
    
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    
    geom_point(aes(color = "Mean Bias"), size = 2) +
    
    geom_errorbar(
      aes(
        ymin = bias - emp_sd,
        ymax = bias + emp_sd,
        color = "Bias ± Empirical SD"
      ),
      width = 0.25
    ) +
    
  geom_point(aes(color = "Coverage label shown only when coverage ≠ 1"),
             alpha = 0) +
  
  geom_text(
    aes(y = label_y, label = coverage_label),
    size = 3
  ) +
    
    facet_grid(. ~ method) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.title = element_blank()
    ) +
    
    labs(
      title = title,
      y = "Bias",
      x = "Coefficients"
    ) +
    
    scale_color_manual(values = c(
      "Mean Bias" = "black",
      "Bias ± Empirical SD" = "blue",
      "Coverage label shown only when coverage ≠ 1" = "white"  # hidden color
    ))
}




##### SE plots
plot_se_bias_with_cov <- function(df, title = "SE Bias Plot") {
  
  ggplot(df, aes(x = parameter, y = se_bias)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    
    geom_point(aes(color = "Mean SE Bias"), size = 2) +
    
    geom_errorbar(
      aes(
        ymin = se_bias - emp_sd_se,
        ymax = se_bias + emp_sd_se,
        color = "SE Bias ± Empirical SD"
      ),
      width = 0.25
    )  +
    
    facet_grid(. ~ method) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.title = element_blank()
    ) +
    
    labs(
      title = title,
      y = "SE Bias",
      x = "Coefficients"
    ) +
    
    scale_color_manual(values = c(
      "Mean SE Bias" = "black",
      "SE Bias ± Empirical SD" = "blue"
    ))
}





###################### 5%
table5_eval_CCA  <- table5_eval_CCA  %>% mutate(missingness = "5%",  method = "CCA")
table5_eval_EM <- table5_eval_EM %>% mutate(missingness = "5%", method = "EM")
table5_eval_MI <- table5_eval_MI %>% mutate(missingness = "5%", method = "MI")

eval_df_5 <- bind_rows(table5_eval_CCA, table5_eval_EM, table5_eval_MI)

##### bias plots
plot_bias_with_cov(eval_df_5, "Bias at 5% Missingness for CCA, EM, MI")

##### SE plots
plot_se_bias_with_cov(eval_df_5, "SE Bias at 5% Missingness for CCA, EM, MI")



###################### 15%
table15_eval_CCA  <- table15_eval_CCA  %>% mutate(missingness = "15%",  method = "CCA")
table15_eval_EM <- table15_eval_EM %>% mutate(missingness = "15%", method = "EM")
table15_eval_MI <- table15_eval_MI %>% mutate(missingness = "15%", method = "MI")

eval_df_15 <- bind_rows(table15_eval_CCA, table15_eval_EM, table15_eval_MI)

##### bias plots
plot_bias_with_cov(eval_df_15, "Bias at 15% Missingness for CCA, EM, MI")

##### SE plots
plot_se_bias_with_cov(eval_df_15, "SE Bias at 15% Missingness for CCA, EM, MI")


###################### 30%
table30_eval_CCA  <- table30_eval_CCA  %>% mutate(missingness = "30%",  method = "CCA")
table30_eval_EM <- table30_eval_EM %>% mutate(missingness = "30%", method = "EM")
table30_eval_MI <- table30_eval_MI %>% mutate(missingness = "30%", method = "MI")

eval_df_30 <- bind_rows(table30_eval_CCA, table30_eval_EM, table30_eval_MI)

##### bias plots
plot_bias_with_cov(eval_df_30, "Bias at 30% Missingness for CCA, EM, MI")

##### SE plots
plot_se_bias_with_cov(eval_df_30, "SE Bias at 30% Missingness for CCA, EM, MI")




##################################### New EM function

# Suppose true_beta and true_se are vectors
# Make sure they are named according to the parameters
true_beta <- coef(Final_model)  # replace with your "true" coefficients
true_se   <- sqrt(diag(vcov(Final_model)))

# Function to summarize EM results
summarize_EM <- function(EM_list, true_beta, true_se) {
  M <- length(EM_list)
  n_param <- length(true_beta)
  
  # Extract coefficient and SE matrices
  beta_mat <- sapply(EM_list, function(x) x$coeff)
  se_mat   <- sapply(EM_list, function(x) x$SE)
  
  # Ensure matrices are M x n_param
  beta_mat <- t(beta_mat)
  se_mat   <- t(se_mat)
  
  # Summary statistics
  mean_beta <- colMeans(beta_mat)
  bias      <- mean_beta - true_beta
  abs_bias  <- abs(bias)
  rel_bias  <- abs_bias / abs(true_beta)
  emp_sd    <- apply(beta_mat, 2, sd)
  mean_se   <- colMeans(se_mat)
  emp_sd_se <- apply(se_mat, 2, sd)
  se_bias   <- mean_se - true_se
  abs_se_bias <- abs(se_bias)
  rel_se_bias <- abs_se_bias / true_se
  
  coverage <- sapply(seq_len(n_param), function(j) {
    mean((beta_mat[, j] - 1.96*se_mat[, j] <= true_beta[j]) &
           (true_beta[j] <= beta_mat[, j] + 1.96*se_mat[, j]))
  })
  
  # Compute RMSE
  rmse_beta <- sqrt(colMeans((beta_mat - matrix(true_beta, nrow=M, ncol=n_param, byrow=TRUE))^2))
  rmse_se   <- sqrt(colMeans((se_mat - matrix(true_se, nrow=M, ncol=n_param, byrow=TRUE))^2))
  
  # Combine into a dataframe
  df <- data.frame(
    parameter    = names(true_beta),
    true_beta    = true_beta,
    mean_beta    = mean_beta,
    abs_bias     = abs_bias,
    rel_bias     = rel_bias,
    bias         = bias,
    emp_sd       = emp_sd,
    true_se      = true_se,
    mean_se      = mean_se,
    emp_sd_se    = emp_sd_se,
    se_bias      = se_bias,
    abs_se_bias  = abs_se_bias,
    coverage     = coverage,
    rel_se_bias  = rel_se_bias,
    rmse_beta    = rmse_beta,
    rmse_se      = rmse_se
  )
  
  rownames(df) <- NULL  # Remove row names
  return(df)
}
# Apply for 5%, 15%, 30% missing
summary_5_EM  <- summarize_EM(EM_results_5, true_beta, true_se)
summary_15_EM <- summarize_EM(EM_results_15, true_beta, true_se)
summary_30_EM <- summarize_EM(EM_results_30, true_beta, true_se)

# Save summaries to CSV files
write.csv(summary_5_EM,  "summary_5_EM.csv",  row.names = FALSE)
write.csv(summary_15_EM, "summary_15_EM.csv", row.names = FALSE)
write.csv(summary_30_EM, "summary_30_EM.csv", row.names = FALSE)







