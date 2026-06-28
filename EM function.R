######################## Sample EM
library(AER)
library(cat)

data("Affairs")
df <- Affairs

# Data needs to be turned to categories with integers, Counts need to be converted to categories as well
df$affairs       <- as.integer(factor(df$affairs))
df$gender        <- as.integer(factor(df$gender))
df$age           <- as.integer(factor(df$age))
df$yearsmarried   <- as.integer(factor(df$yearsmarried ))
df$children      <- as.integer(factor(df$children))
df$religiousness <- as.integer(factor(df$religiousness))
df$education     <- as.integer(factor(df$education))
df$occupation     <- as.integer(factor(df$occupation))
df$rating        <- as.integer(factor(df$rating))

set.seed(123)
df$affairs[sample(1:nrow(df), 30)]         <- NA
df$gender[sample(1:nrow(df), 40)]          <- NA
df$age[sample(1:nrow(df), 35)]             <- NA
df$yearsmarried[sample(1:nrow(df), 35)]    <- NA
df$children[sample(1:nrow(df), 40)]        <- NA
df$religiousness[sample(1:nrow(df), 50)]   <- NA
df$education[sample(1:nrow(df), 30)]       <- NA
df$occupation[sample(1:nrow(df), 30)]      <- NA
df$rating[sample(1:nrow(df), 25)]          <- NA


df_mat <- as.matrix(df)
pre <- prelim.cat(df_mat) # describes how many categories each variable has and which values are missing

em_out <- em.cat(pre, maxits = 100, eps = 1e-08) 
# E step (Expectation Step) : it calculates the probability of each possible category for the missing entries, 
# based on the observed values in that same row and and based on the current estimate of the joint probability table for all variables with missing cells
# M step (Maximization Step): it updates joint probabilities of all categories, marginal distributions and conditional probabilities
# by using observed counts expected counts from the E-step

# eps = 1e-08: The algorithm stops when changes between iterations are less than 0.00000001.
# maxits = 2000: It will force-stop after 2000 iterations if it still hasn’t converged.
# em_out has :
# (1) estimated cell probabilities: full joint probability table for all categorical variables.
# (2)  the estimated mean vector (for categorical cells): Used for imputation.
# (3) posterior distributions: These are used by imp.cat() to randomly fill missing values.

df_imputed <- imp.cat(pre, em_out)
df_imputed <- as.data.frame(df_imputed)
# This uses the posterior probabilities from em_out to: sample a category for each missing cell and generate a fully completed dataset

colnames(df_imputed) <- c(
  "affairs",
  "gender",
  "age",
  "yearsmarried",
  "children",
  "religiousness",
  "education",
  "occupation",
  "rating"
)






############################################
##.    EM
############################################
library(AER)
library(cat)
library(VGAM)
library(dplyr)

############################################
### 1. Add extra variables to each dataset
############################################
cols_to_add <- Affairs[, c("age", "yearsmarried", "education", "occupation")]

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
run_EM_and_vglm <- function(df, i) {
  cat("  → Running dataset", i, "\n")
  
  df <- make_all_cat(df)
  
  M <- as.matrix(df)
  colnames(M) <- colnames(df)   # keep names
  
  pre <- prelim.cat(M)
  em_out <- em.cat(pre, maxits = 2000, eps = 1e-08)
  
  df_imp <- as.data.frame(imp.cat(pre, em_out))
  colnames(df_imp) <- colnames(df)   # restore names
  
  df_imp$affairs_mod <- ordered(df_imp$affairs)
  
  fit <- vglm(
    affairs_mod ~ children + religiousness + rating,
    data   = df_imp,
    family = cumulative(parallel = TRUE)
  )
  
  est <- c(coef(fit), fit@misc$zeta)
  
  se_tab <- summary(fit)@coef3
  se <- se_tab[names(est), "Std. Error"]
  
  tibble(
    dataset  = i,
    term     = names(est),
    estimate = est,
    se       = se
  )
}

############################################
### 4. Wrapper to run EM + VGLM for a list
############################################
run_for_list <- function(lst, label) {
  cat("\n===== Running", label, "datasets =====\n")
  bind_rows(
    lapply(seq_along(lst), function(i) run_EM_and_vglm(lst[[i]], i))
  )
}

############################################
### 5. Results
############################################
EM5_results  <- run_for_list(NA5_list_mod,  "5% Missingness")
EM15_results <- run_for_list(NA15_list_mod, "15% Missingness")
EM30_results <- run_for_list(NA30_list_mod, "30% Missingness")

# Converts EM results 
prepare_EM_results <- function(EM_res) {
  
  params <- unique(EM_res$term)
  
  M <- length(unique(EM_res$dataset))
  
  betas <- matrix(NA, nrow = M, ncol = length(params))
  SEs   <- matrix(NA, nrow = M, ncol = length(params))
  colnames(betas) <- params
  colnames(SEs)   <- params
  
  for (i in seq_len(M)) {
    subset_i <- EM_res[EM_res$dataset == i, ]
    betas[i, ] <- subset_i$estimate
    SEs[i, ]   <- subset_i$se
  }
  
  list(
    betas = betas,
    SEs   = SEs
  )
}

res5_EM  <- prepare_EM_results(EM5_results)
res15_EM <- prepare_EM_results(EM15_results)
res30_EM <- prepare_EM_results(EM30_results)



###############################################
### 6. CREATE TABLES for raw data estimates
###############################################
make_full_table <- function(betas, SEs) {
  data.frame(
    parameter   = colnames(betas),
    mean_beta   = colMeans(betas),
    emp_sd      = apply(betas, 2, sd),
    mean_se     = colMeans(SEs)
  )
}

table5_full_EM  <- make_full_table(res5_EM$betas,  res5_EM$SEs)
table15_full_EM <- make_full_table(res15_EM$betas, res15_EM$SEs)
table30_full_EM <- make_full_table(res30_EM$betas, res30_EM$SEs)

write.csv(table5_full_EM,  "table5_full_EM.csv",  row.names = FALSE)
write.csv(table15_full_EM, "table15_full_EM.csv", row.names = FALSE)
write.csv(table30_full_EM, "table30_full_EM.csv", row.names = FALSE)

###############################################
### 7. CREATE TABLES for comparisons to True model
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

table5_eval_EM  <- make_eval_table(res5_EM,  true_params, true_SEs)
table15_eval_EM <- make_eval_table(res15_EM, true_params, true_SEs)
table30_eval_EM <- make_eval_table(res30_EM, true_params, true_SEs)

write.csv(table5_eval_EM,  "table5_eval_EM.csv",  row.names = FALSE)
write.csv(table15_eval_EM, "table15_eval_EM.csv", row.names = FALSE)
write.csv(table30_eval_EM, "table30_eval_EM.csv", row.names = FALSE)













