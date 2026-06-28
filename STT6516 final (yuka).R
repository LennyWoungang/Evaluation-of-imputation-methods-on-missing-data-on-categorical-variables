library(AER)
library(dplyr)
library(MASS)
library(glmnet)
library(nnet)
library(ordinalNet)
library(rms)
library(VGAM)
library(future.apply)


create_NA <- function(data, p_row) {
  n <- nrow(data)
  p <- round(p_row * n)         
  cols <- colnames(data[, c("affairs_mod", "children", "religiousness", "rating")])
  rows_to_miss <- sample.int(n, p, replace = FALSE)
  out <- data
  for (r in rows_to_miss) {
    c <- sample(cols, 1)
    out[r, c] <- NA
  }
  out
}
############################################################################


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



model_po <- vglm(affairs_mod ~  (children + religiousness  + rating + gender)^2,
                 data = data, family = cumulative(parallel = T) )



model <- vglm(affairs_mod~ ( children + religiousness  + rating + gender)^2,
              data = data, family = cumulative(parallel = F ~ religiousness + rating)) 


dev <- deviance(model_po) - deviance(model)
df <- df.residual(model_po) - df.residual(model)
1-pchisq(dev, df)


step_mod <- step4vglm(model_po, direction = "backward", trace = TRUE)



# Check the hypothesis with the selected model

model_po <- vglm(affairs_mod ~ children + religiousness + rating,
                 data = data, family = cumulative(parallel = T) )



model <- vglm(affairs_mod ~ children + religiousness + rating,
              data = data, family = cumulative(parallel = F ~ religiousness + rating)) 


dev <- deviance(model_po) - deviance(model)
df <- df.residual(model_po) - df.residual(model)
1-pchisq(dev, df)

check_hyp <- lrm(affairs_mod ~ children + religiousness + rating,
                 data = data,
                 x=T, y=T)

par(mfrow = c(3, 3))
residuals(check_hyp, type = "partial", pl = TRUE)

par(mfrow = c(1, 1))



# lines are crossed, we should use partial proportional ordinal model 

final_mod <- vglm(affairs_mod ~ children + religiousness + rating,
                  data = data, family = cumulative(parallel = T) )
summary(final_mod)
SE <- sqrt(diag(vcov(final_mod)))


## --- Knock off some values --- ##


set.seed(6516)
M <- 500
NA5_list <- vector("list", M)
for(i in 1:M){
  NA5 <- create_NA(data = data, p = 0.05)
  NA5_list[[i]] <- NA5
}


NA15_list <- vector("list", M)

for(i in 1:M){
  NA15 <- create_NA(data = data, p = 0.15)
  NA15_list[[i]] <- NA15
}

NA30_list <- vector("list", M)
for(i in 1:M){
  NA30 <- create_NA(data = data, p = 0.3)
  NA30_list[[i]] <- NA30
}



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






true_beta <- coef(final_mod)
param_names <-  names(true_beta) 
#true_beta <- matrix(true_beta, nrow = M, ncol = length(true_beta), byrow = TRUE)







## --- CCA --- #



## function that does your whole simulation summary for ONE NA*_list
summarise_CCA <- function(NA_list, true_beta, param_names) {
  
  M <- length(NA_list)
  
  betas <- matrix(NA_real_, nrow = M, ncol = length(param_names))
  SEs   <- matrix(NA_real_, nrow = M, ncol = length(param_names))
  colnames(betas) <- param_names
  colnames(SEs)   <- param_names
 
  for (i in seq_len(M)) {
    data <- NA_list[[i]]
    data_cca <- na.omit(data)
    
    fit <- vglm(affairs_mod ~ children + religiousness + rating,
                data = data_cca, family = cumulative(parallel = TRUE))
    
    coeffs <- coef(fit)
    SE     <- sqrt(diag(vcov(fit)))  # <- use vcov(fit), not final_mod
    
    betas[i, ] <- coeffs
    SEs[i, ]   <- SE
  }
  
  ## true_beta should be a vector in the same order as param_names
  z <- qnorm(1 - 0.05/2)
  
  bias    <- colMeans(betas) - true_beta
  emp_sd  <- apply(betas, 2, sd)
  mean_se <- colMeans(SEs)
  
  low  <- betas - z * SEs
  high <- betas + z * SEs
  
  ## recycle true_beta across rows
  true_mat   <- matrix(true_beta, nrow = M, ncol = length(true_beta), byrow = TRUE)
  cove_check <- (true_mat >= low) & (true_mat <= high)
  coverage   <- colMeans(cove_check)
  
  data.frame(
    parameters = param_names,
    true_beta  = true_beta,
    bias       = bias,
    emp_sd     = emp_sd,
    mean_se    = mean_se,
    coverage   = coverage,
    row.names  = NULL
  )
}

summary5  <- summarise_CCA(NA5_list,  true_beta = true_beta, param_names = param_names)
summary15 <- summarise_CCA(NA15_list, true_beta = true_beta, param_names = param_names)
summary30 <- summarise_CCA(NA30_list, true_beta = true_beta, param_names = param_names)




CCA <- function(NA_list, true_beta, param_names) {
  
  M <- length(NA_list)
  
  betas <- matrix(NA_real_, nrow = M, ncol = length(param_names))
  SEs   <- matrix(NA_real_, nrow = M, ncol = length(param_names))
  colnames(betas) <- param_names
  colnames(SEs)   <- param_names
  
  for (i in seq_len(M)) {
    data <- NA_list[[i]]
    data_cca <- na.omit(data)
    
    fit <- vglm(affairs_mod ~ children + religiousness + rating,
                data = data_cca, family = cumulative(parallel = TRUE))
    
    coeffs <- coef(fit)
    SE     <- sqrt(diag(vcov(fit)))  # <- use vcov(fit), not final_mod
    
    betas[i, ] <- coeffs
    SEs[i, ]   <- SE
  }
  
 list(beta = betas, SEs = SEs)
}


CCA5  <- CCA(NA5_list,  true_beta = true_beta, param_names = param_names)
CCA15 <- CCA(NA15_list, true_beta = true_beta, param_names = param_names)
CCA30 <- CCA(NA30_list, true_beta = true_beta, param_names = param_names)





### --- MI --- ###

## Check the MI first 
## settings
m <- 20
set.seed(123)  # for choosing which 5 datasets to inspect

# randomly choose 5 out of your 500 datasets
i_test <- sample(seq_along(NA30_list), 5, replace = FALSE)

# PDF to store all diagnostics
pdf("MI_diagnostics_5datasets.pdf")

for (k in seq_along(i_test)) {
  idx <- i_test[k]
  cat("\n==============================\n")
  cat("Diagnostics for dataset", idx, "\n")
  cat("==============================\n")
  
  # 1) take one incomplete dataset
  data_test <- NA30_list[[idx]]
  
  # 2) run mice on it
  imp_test <- mice(
    data_test,
    m         = m,
    printFlag = FALSE,   # set TRUE if you want iteration logs
    seed      = 999
  )
  
  # quick textual summary in console
  print(summary(imp_test))
  
  ## --- 3) Convergence plots ---
  # mice::plot() returns a lattice object; use print() inside loops
  print(plot(imp_test))
  
  ## --- 4) Stripplot for key variables ---
  print(stripplot(
    imp_test,
    affairs_mod + children + religiousness + rating ~ .imp
  ))
  
  ## --- 5) Barplots: observed vs imputed for 4 categorical vars ---
  # long format with original + all imputations
  long <- complete(imp_test, "long", include = TRUE)
  long$source <- ifelse(long$.imp == 0, "observed", "imputed")
  
  # variables you want to check
  vars_cat <- c("affairs_mod", "children", "religiousness", "rating")
  
  # convert to character to avoid factor / ordered factor mixing issues
  long_cat <- long %>%
    mutate(across(all_of(vars_cat), as.character)) %>%
    dplyr::select(source, all_of(vars_cat)) %>%
    pivot_longer(
      cols      = all_of(vars_cat),
      names_to  = "variable",
      values_to = "level"
    )
  
  # one faceted barplot per dataset
  p_bar <- ggplot(long_cat, aes(x = level, fill = source)) +
    geom_bar(position = "fill") +
    facet_wrap(~ variable, scales = "free_x") +
    ylab("Proportion") +
    ggtitle(paste("Dataset", idx, "- observed vs imputed"))
  
  print(p_bar)
}

dev.off()






## Now analysis 
M <- 500
m <- 20
#true_beta <- matrix(true_beta, nrow = M, ncol = length(true_beta), byrow = TRUE)

plan(multisession, workers = 4)

run_MI_scenario <- function(NA_list,          # e.g. NA5_list
                            scenario_label,   # e.g. "5", "15", "30"
                            M, m, 
                            param_names, 
                            true_beta) {
  
  # true_beta can be a vector (length p) or M x p matrix with same rows
  if (is.matrix(true_beta) || is.data.frame(true_beta)) {
    beta_true_vec <- colMeans(as.matrix(true_beta))  # p-vector
  } else {
    beta_true_vec <- as.numeric(true_beta)
  }
  
  if (length(beta_true_vec) != length(param_names)) {
    stop("Length of true_beta does not match length of param_names.")
  }
  
  t0 <- Sys.time()
  
  MI_results <- future.apply::future_lapply(1:M, function(i){
    
    data <- NA_list[[i]]
    
    imp <- mice::mice(data = data,
                      m = m,
                      printFlag = FALSE, 
                      seed = 1212 + i)
    
    comp_list <- mice::complete(imp, action = "all")
    
    imp_betas <- matrix(NA_real_, nrow = m, ncol = length(param_names))
    imp_SEs   <- matrix(NA_real_, nrow = m, ncol = length(param_names))
    colnames(imp_betas) <- param_names
    colnames(imp_SEs)   <- param_names
    
    for (t in 1:m) {
      split_data <- comp_list[[t]] %>%
        dplyr::select(affairs_mod, children, religiousness, rating)
      
      fit_comp <- VGAM::vglm(
        affairs_mod ~ children + religiousness + rating,
        data   = split_data,
        family = VGAM::cumulative(parallel = TRUE)
      )
      
      coeffs_comp <- coef(fit_comp)
      SE_comp     <- sqrt(diag(vcov(fit_comp)))
      
      imp_betas[t, ] <- coeffs_comp
      imp_SEs[t,  ]  <- SE_comp
    }
    
    ## Rubin's Rules within each simulated dataset i
    beta_bar <- colMeans(imp_betas)           # pooled beta
    between  <- apply(imp_betas, 2, var)      # between-imputation var
    within   <- colMeans(imp_SEs^2)           # within-imputation var
    T_var    <- within + (1 + 1/m) * between  # total var
    se_pool  <- sqrt(T_var)
    
    list(beta_bar = beta_bar, se = se_pool, T_var = T_var)
    
  }, future.seed = TRUE)
  
  t_total <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  message("Scenario ", scenario_label, " finished in ", round(t_total, 1), " seconds.")
  
  ## Stack Rubin results across M datasets
  rubin_beta_bars <- do.call(rbind, lapply(MI_results, `[[`, "beta_bar"))  # M x p
  rubin_se        <- do.call(rbind, lapply(MI_results, `[[`, "se"))        # M x p
  rubin_T_var     <- do.call(rbind, lapply(MI_results, `[[`, "T_var"))     # M x p
  
  ## Monte Carlo summaries over M (sim repetitions)
  MI_bias   <- colMeans(rubin_beta_bars) - beta_true_vec
  emp_sd    <- apply(rubin_beta_bars, 2, sd)
  mean_Tvar <- colMeans(rubin_T_var)
  mean_se   <- sqrt(mean_Tvar)
  
  z    <- qnorm(1 - 0.05/2)
  low  <- rubin_beta_bars - z * sqrt(rubin_T_var)
  high <- rubin_beta_bars + z * sqrt(rubin_T_var)
  
  ## Coverage: for each parameter j, proportion of CIs that contain true beta_j
  inside   <- sweep(low,  2, beta_true_vec, "<=") &
    sweep(high, 2, beta_true_vec, ">=")
  coverage <- colMeans(inside)
  
  MI_summary <- data.frame(
    parameters = colnames(rubin_beta_bars),
    true_beta  = beta_true_vec,
    bias       = MI_bias,
    emp_sd     = emp_sd,
    mean_se    = mean_se,
    coverage   = coverage,
    row.names  = NULL
  )
  
  ## Save if you want
  saveRDS(MI_results, file = paste0("MI_results_", scenario_label, ".rds"))
  saveRDS(MI_summary, file = paste0("MI_summary_", scenario_label, ".rds"))
  
  ## Also return them to R
  list(
    MI_results      = MI_results,
    MI_summary      = MI_summary,
    rubin_beta_bars = rubin_beta_bars,
    rubin_T_var     = rubin_T_var,
    rubin_se        = rubin_se,
    runtime_sec     = t_total
  )
}



res5  <- run_MI_scenario(NA5_list,  scenario_label = "5",  M, m, param_names, true_beta)
res15 <- run_MI_scenario(NA15_list, scenario_label = "15", M, m, param_names, true_beta)
res30 <- run_MI_scenario(NA30_list, scenario_label = "30", M, m, param_names, true_beta)


MI_summary5 <- res5$MI_summary
MI_summary15 <- res15$MI_summary
MI_summary30 <- res30$MI_summary



relative_bias5 <- MI_summary5$bias/true_beta
relative_bias15 <- MI_summary15$bias/true_beta
relative_bias30 <- MI_summary30$bias/true_beta



##############################################################################################
library(tibble) 
library(flextable) 
library(gtsummary)




data2 <- data %>%
  mutate(
    children = factor(children,
                      levels = c("no", "yes"),
                      labels = c("No children", "Has children")),
    religiousness = factor(
      religiousness,
      levels = 1:5,
      labels = c("Anti", "Not at all", "Slightly", "Somewhat", "Very")
    ),
    rating = factor(
      rating,
      levels = 1:5,
      labels = c("Very unhappy", "Somewhat unhappy", "Average", "Happier than average", "Very happy")
    )
  )



library(gtsummary)
library(kableExtra)
library(dplyr)

## 1. Build the gtsummary table (no markdown in headers)
table1 <- data2 %>% 
  select(affairs_mod, children, religiousness, rating) %>% 
  tbl_summary(
    by = affairs_mod,
    statistic = list(all_categorical() ~ "{n} ({p}%)"),
    type = list(
      children      ~ "categorical",
      religiousness ~ "categorical",
      rating        ~ "categorical"
    )
  ) %>% 
  add_n() %>% 
  modify_header(
    label  ~ "Variable",
    stat_1 ~ "Affairs = 0 <br>(n = {n})",
    stat_2 ~ "Affairs = 1–3 <br>(n = {n})",
    stat_3 ~ "Affairs = 4+ <br>(n = {n})"
  ) %>% 
  modify_spanning_header(all_stat_cols() ~ "Affairs") %>%   
  modify_caption(
    "<span style='color:black'>
     Table 1. Descriptive statistics by categories of affairs
     </span>"
  ) %>%
  # remove the default gtsummary footnote "n (%)"
  remove_footnote_header(columns = everything()) %>%
  remove_footnote_body(columns = everything(), rows = TRUE)

## 2. Convert to kableExtra and style it
n_body <- nrow(table1$table_body)

table1 %>%
  as_kable_extra(
    booktabs    = TRUE,
    escape      = FALSE,   # allow HTML in caption / footnote
    addtl_fmt   = FALSE    # don't auto-escape our HTML
  ) %>%
  kable_styling(full_width = FALSE) %>%
  # header row: thick top and bottom border, bold text
  row_spec(
    0,
    extra_css = "border-top:2px solid black;
                 border-bottom:2px solid black;
                 font-weight:bold;"
  ) %>%
  # last body row: thick bottom border
  row_spec(
    n_body,
    extra_css = "border-bottom:2px solid black;"
  ) %>%
  # custom small footnote: "Note: n (%)"
  footnote(
    general = "<span style='font-size:0.85em;'>
                 <em>Note:</em> n (%)
               </span>",
    general_title  = "",      # put everything on one line
    threeparttable = TRUE,
    escape         = FALSE
  )




### Table for the analysis 
z <- true_beta/SE
alpha <- 0.05
z_crit <- qnorm(1-alpha/2)
CI_low <- true_beta - z_crit*SE
CI_high <- true_beta + z_crit*SE

tab_coef <- tibble(
  term      = names(true_beta),
  estimate  = true_beta,
  std_error = SE,                         
  z_value   = z,
  CI_lower  = CI_low,
  CI_upper  = CI_high
)

p_val <- 2 * pnorm(abs(z), lower.tail = FALSE)

library(knitr)
library(kableExtra)

tab_coef %>%
  mutate(
    estimate  = sprintf("%.3f", estimate),
    std_error = sprintf("%.3f", std_error),
    z_value   = sprintf("%.2f", z_value),
    CI        = sprintf("[%.3f, %.3f]", CI_lower, CI_upper)
  ) %>%
  dplyr::select(term, estimate, std_error, z_value, CI) %>%
  kable(booktabs = TRUE,
        col.names = c("Coefficient", "Estimate", "SE", 
                      "Z", "95% C.I."),
        caption   = "<span style='color:black'>Table 2. Proportional odds model for number of affairs</span>"
  ) %>%
kable_styling(full_width = FALSE) %>% 
  row_spec(0, extra_css = "border-top: 2px solid black; border-bottom: 2px solid black;") %>%
  # bottom horizontal line (use number of rows in tab_coef)
  row_spec(nrow(tab_coef), extra_css = "border-bottom: 2px solid black;") %>% 
  footnote(
    general = paste0(
      "<span style='font-size:0.85em;'><em>Note:</em> ",
      "Residual degrees of freedom: ", df.residual(final_mod),
      ". Estimates are log-odds.</span>"
    ),
    general_title  = "",     # <- no separate "Note:" cell
    threeparttable = TRUE,
    escape         = FALSE   # allow HTML in footnote
  )



library(dplyr)
library(kableExtra)

## 1) Extract coefficients & SEs
beta  <- coef(final_mod)
se    <- sqrt(diag(vcov(final_mod)))
zcrit <- qnorm(0.975)

## 2) Keep only the effects you want to present
keep <- c("childrenyes", "religiousness.L", "rating.L")

tab_or <- tibble(
  Predictor = c("Children: yes vs no",
                "Religiousness (L)",
                "Marital rating (L)"),
  beta  = beta[keep],
  se    = se[keep]
) %>%
  mutate(
    OR      = exp(beta),
    OR_low  = exp(beta - zcrit * se),
    OR_high = exp(beta + zcrit * se)
  )

## 3) Format for display
tab_or_disp <- tab_or %>%
  mutate(
    OR = sprintf("%.2f", OR),
    CI = sprintf("(%.2f, %.2f)", OR_low, OR_high)
  ) %>%
  select(Predictor, OR, CI)

n_rows <- nrow(tab_or_disp)

## 4) Publish-ready kable
tab_or_disp %>%
  kable(
    booktabs  = TRUE,
    col.names = c("Predictor", "Odds ratio", "95% C.I."),
    caption   = "<span style='color:black'>
                 Table 3. Selected odds ratios from proportional odds model
                 </span>",
    escape    = FALSE
  ) %>%
  kable_styling(full_width = FALSE) %>%
  # thick top + thick header-bottom line
  row_spec(
    0,
    extra_css = "border-top:2px solid black;
                 border-bottom:2px solid black;
                 font-weight:bold;"
  ) %>%
  # thick bottom line under last row
  row_spec(
    n_rows,
    extra_css = "border-bottom:2px solid black;"
  ) %>%
  footnote(
    general = "<span style='font-size:0.85em;'>
                 <em>Note:</em> Odds ratios are for cumulative odds of being
                 in lower-affair categories versus higher.
               </span>",
    general_title  = "",
    threeparttable = TRUE,
    escape         = FALSE
  )

##################################################################################
## --- RMSE calculatoin --- ##


setwd("C:\\Users\\yukas\\Documents\\School\\2025 Fall\\STT6516\\projet")
EM5 <- read.csv("table5_full_EM.csv")
EM15 <- read.csv("table15_full_EM.csv")
EM30 <- read.csv("table30_full_EM.csv")


library(tidyr)

EM5_coeff <-  EM5 %>%
  dplyr::select(dataset, term, estimate) %>%
  pivot_wider(
    id_cols   = dataset,
    names_from  = term,
    values_from = estimate
  )

EM15_coeff <-  EM15 %>%
  dplyr::select(dataset, term, estimate) %>%
  pivot_wider(
    id_cols   = dataset,
    names_from  = term,
    values_from = estimate
  )

EM30_coeff <-  EM30 %>%
  dplyr::select(dataset, term, estimate) %>%
  pivot_wider(
    id_cols   = dataset,
    names_from  = term,
    values_from = estimate
  ) 

EM5_SE <- EM5 %>%
  dplyr::select(dataset, term, se) %>%
  pivot_wider(
    id_cols   = dataset,
    names_from  = term,
    values_from = se
  )

EM15_SE <- EM15 %>%
  dplyr::select(dataset, term, se) %>%
  pivot_wider(
    id_cols   = dataset,
    names_from  = term,
    values_from = se
  )

EM30_SE <- EM30 %>%
  dplyr::select(dataset, term, se) %>%
  pivot_wider(
    id_cols   = dataset,
    names_from  = term,
    values_from = se
  )


EM5_coeff <- dplyr::select(EM5_coeff, -dataset)
EM15_coeff <- dplyr::select(EM15_coeff, -dataset)
EM30_coeff <- dplyr::select(EM30_coeff, -dataset)

EM5_SE <- dplyr::select(EM5_SE, -dataset)
EM15_SE <- dplyr::select(EM15_SE, -dataset)
EM30_SE <- dplyr::select(EM30_SE, -dataset)



MICE5_coeff <- do.call(rbind, lapply(MI_results_5, `[[`, "beta_bar")) 
MICE15_coeff <- do.call(rbind, lapply(MI_results_15, `[[`, "beta_bar")) 
MICE30_coeff <- do.call(rbind, lapply(MI_results_30, `[[`, "beta_bar")) 


MICE5_SE <- do.call(rbind, lapply(MI_results_5, `[[`, "se")) 
MICE15_SE <- do.call(rbind, lapply(MI_results_15, `[[`, "se")) 
MICE30_SE <- do.call(rbind, lapply(MI_results_30, `[[`, "se")) 


CCA5_coeff <- CCA5$beta
CCA15_coeff <- CCA15$beta
CCA30_coeff <- CCA30$beta

CCA5_SE <- CCA5$SEs
CCA15_SE <- CCA15$SEs
CCA30_SE <- CCA30$SEs

missing5_coeff <- cbind(CCA5_coeff, EM5_coeff, MICE5_coeff)
missing15_coeff <- cbind(CCA15_coeff, EM15_coeff, MICE15_coeff)
missing30_coeff <- cbind(CCA30_coeff, EM30_coeff, MICE30_coeff)

missing5_SE <- cbind(CCA5_SE, EM5_SE, MICE5_SE)
missing15_SE <- cbind(CCA15_SE, EM15_SE, MICE15_SE)
missing30_SE <- cbind(CCA30_SE, EM30_SE, MICE30_SE)


true_beta_mat <- matrix(rep(true_beta, times = 500),
                        nrow = 500, byrow = TRUE)
true_beta_mat <- cbind(true_beta_mat, true_beta_mat, true_beta_mat)



rmse5 <- sqrt(colMeans( (missing5_coeff - true_beta_mat)^2 ))
rmse15 <- sqrt(colMeans( (missing15_coeff - true_beta_mat)^2 ))
rmse30 <- sqrt(colMeans( (missing30_coeff - true_beta_mat)^2 ))



true_SE_mat <- matrix(rep(SE, times = 500),
                       nrow = 500, byrow = TRUE) 
true_SE_mat <- cbind(true_SE_mat, true_SE_mat, true_SE_mat)

rmse5_SE <- sqrt(colMeans( (missing5_SE - true_SE_mat)^2 ))
rmse15_SE <- sqrt(colMeans( (missing15_SE - true_SE_mat)^2 ))
rmse30_SE <- sqrt(colMeans( (missing30_SE - true_SE_mat)^2 ))





p <- length(param_names)

make_rmse_df <- function(rmse_vec, miss_lab) {
  data.frame(
    terms = param_names,
    CCA   = rmse_vec[1:p],
    EM    = rmse_vec[(p + 1):(2 * p)],
    MICE  = rmse_vec[(2 * p + 1):(3 * p)]
  ) |> 
    dplyr::rename_with(~ paste0(.x, "_", miss_lab), -terms)
}

rmse5_df  <- make_rmse_df(rmse5,  "5")
rmse15_df <- make_rmse_df(rmse15, "15")
rmse30_df <- make_rmse_df(rmse30, "30")


rmse5_df_SE  <- make_rmse_df(rmse5_SE,  "5")
rmse15_df_SE <- make_rmse_df(rmse15_SE, "15")
rmse30_df_SE <- make_rmse_df(rmse30_SE, "30")




library(dplyr)

rmse_coeff <- rmse5_df %>%
  left_join(rmse15_df, by = "terms") %>%
  left_join(rmse30_df, by = "terms")

rmse_SE <- rmse5_df_SE %>% 
  left_join(rmse15_df_SE, by = 'terms') %>% 
  left_join(rmse30_df_SE, by = 'terms')


rmse_coeff <- rmse_coeff %>%
  dplyr::select(
    terms,
    CCA_5,  EM_5,  MICE_5,
    CCA_15, EM_15, MICE_15,
    CCA_30, EM_30, MICE_30
  ) %>% 
  mutate(
    terms = case_when(
      terms == "(Intercept):1" ~ "Intercept 1",
      terms == "(Intercept):2" ~ "Intercept 2",
      TRUE ~ terms
    )
  )


rmse_SE <- rmse_SE %>%
  dplyr::select(
    terms,
    CCA_5,  EM_5,  MICE_5,
    CCA_15, EM_15, MICE_15,
    CCA_30, EM_30, MICE_30
  ) %>% 
  mutate(
    terms = case_when(
      terms == "(Intercept):1" ~ "Intercept 1",
      terms == "(Intercept):2" ~ "Intercept 2",
      TRUE ~ terms
    )
  )




write.csv(rmse_coeff, 
          file = "rmse_coeff.csv")


write.csv(rmse_SE, 
          file = "rmse_SE.csv")
