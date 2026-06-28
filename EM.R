
# Reference 
#https://github.com/arnabkrmaity/ProportionalOddsMissingResponse/blob/main/PropMissingResponse.sas
#Cumulative Logit Ordinal Regression with Proportional Odds
#under Nonignorable Missing Responses – Application to Phase III
#Trial

library(MASS)

#####################################################################################
## --- Helpers --- ##

## theta: vector (cutpoints first, then slopes) same order as coef(vglm)
## x_row: 1-row data.frame
## x_vars: predictor names
## y_levels: levels of Y (ordered)
prob_PO <- function(theta, x_row, x_vars, y_levels) {
  K <- length(y_levels)
  
  alpha <- theta[1:(K - 1)]
  beta  <- theta[K:length(theta)]
  
  # Build model matrix like in vglm fit
  X <- model.matrix(
    as.formula(paste("~", paste(x_vars, collapse = " + "))),
    data = x_row[, x_vars, drop = FALSE]
  )
  if ("(Intercept)" %in% colnames(X)) {
    X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
  }
  
  lp <- as.numeric(X %*% beta)
  
  eta   <- alpha + lp
  gamma <- plogis(eta)           # cum probabilities P(Y <= k | X)
  
  p <- numeric(K)
  p[1] <- gamma[1]
  if (K > 2) {
    for (k in 2:(K - 1)) {
      p[k] <- gamma[k] - gamma[k - 1]
    }
  }
  p[K] <- 1 - gamma[K - 1]
  
  p
}



##############################################################################################
## --- Single Imputation ---- ##


EM_handmade <- function(data, y_var, x_vars,
                        maxitr = 50, tol = 1e-6, verbose = TRUE) {
  
  N <- nrow(data)
  
  ## Preserve which were ordered BEFORE coercing, if you care
  is_ord <- sapply(c(y_var, x_vars), function(v) is.ordered(data[[v]]))
  
  
  data[[y_var]] <- factor(data[[y_var]])
  for (v in x_vars) {
    data[[v]] <- factor(data[[v]])
  }
  
  y_levels <- levels(data[[y_var]])
  K <- length(y_levels)
  
  x_levels_list <- lapply(x_vars, function(v) levels(data[[v]]))
  names(x_levels_list) <- x_vars
  
  ## Complete-cases index for initial fit
  cc_idx <- stats::complete.cases(data[, c(y_var, x_vars), drop = FALSE])
  
  ## Initial proportional odds fit
  mod <- as.formula(paste(y_var, "~", paste(x_vars, collapse = " + ")))
  fit <- vglm(mod, family = cumulative(parallel = TRUE),
              data = data[cc_idx, , drop = FALSE])
  coeffs <- coef(fit)
  
  ## Initial marginals for each X_j
  # j is an indicator of predictors, j = 1, 2, 3 (children, religiousness, rating)
  # l indicates levels in j 
  # pi_{jl}^{t=0} = Pr(X_j = x_{jl})
  pi_list <- lapply(x_vars, function(v) {
    tab <- table(data[[v]], useNA = "no")
    tab / sum(tab)
  })
  names(pi_list) <- x_vars
  
  get_param_vec <- function(coeffs, pi_list) {
    c(coeffs, unlist(lapply(pi_list, as.numeric)))
  }
  param_old <- get_param_vec(coeffs, pi_list)
  
  ## -------- EM loop -------- ##
  for (iter in seq_len(maxitr)) {
    
    if (verbose) cat("EM iteration:", iter, "\n")
    
    ## ---- E-step: build augmented dataset ---- ##
    aug_list <- vector("list", N)
    
    for (i in seq_len(N)) {
      row_i <- data[i, , drop = FALSE]
      
      mis_vec <- is.na(row_i[, c(y_var, x_vars), drop = TRUE])
      names(mis_vec) <- c(y_var, x_vars)
      n_mis <- sum(mis_vec)
      
      if (n_mis == 0) {
        ## no missingness is included in the observation
        row_i$wgt <- 1
        aug_list[[i]] <- row_i
        
      } else if (n_mis == 1) {
        
        miss_name <- names(mis_vec)[which(mis_vec)]
        
        if (miss_name == y_var) {
          ## Case 1: Y is missing
          # compute p_k = P(Y=k|X; beta)
          # we have k = 1(0times), 2(1-3 times), 3 (4+)
          # weight becomes its probability of k 
          
          newdata_x <- row_i
          p_y_given_x <- as.numeric(
            predict(fit, newdata = newdata_x, type = "response")
          )
          
          rows_k <- lapply(seq_len(K), function(k) {
            r <- row_i
            r[[y_var]] <- y_levels[k]
            r$wgt <- p_y_given_x[k]
            r
          })
          
          aug_list[[i]] <- do.call(rbind, rows_k)
          
        } else {
          ## Case 2: one of the predictors is missing, Y observed
          # fitting a joint model: 
          # f(Y, X; theta, pi) = f(Y|X; theta) \prod_{l=1}^p f(X_k:pi_k)
          
          
          j <- miss_name
          x_levels_j <- x_levels_list[[j]]
          L_j <- length(x_levels_j)
          pi_j <- pi_list[[j]]
          
          augment <- lapply(seq_len(L_j), function(l) {
            r <- row_i
            r[[j]] <- factor(x_levels_j[l], levels = x_levels_j)
            r
          })
          aug_df <- do.call(rbind, augment)
          
          for (v in x_vars) {
            aug_df[[v]] <- factor(aug_df[[v]], levels = x_levels_list[[v]])
          }
          
          pred_mat <- predict(fit, newdata = aug_df, type = "response")
          
          y_i <- as.character(row_i[[y_var]])
          y_idx <- match(y_i, y_levels)
          
          # conditional prob. of the observed outcome for i if X_j were at level l
          # p_l = P(Y=y|X = x_{jl}, X_{-j}; \theta^t)
          p_y_given_x <- pred_mat[, y_idx]
          
          
          # calculate posterior prob. = w
          for_post_prob <- as.numeric(pi_j) * p_y_given_x
          w <- for_post_prob / sum(for_post_prob)
          
          rows_j <- lapply(seq_len(L_j), function(l) {
            r <- aug_df[l, , drop = FALSE]
            r$wgt <- w[l]
            r
          })
          
          aug_list[[i]] <- do.call(rbind, rows_j)
        }
        
      } else {
        stop("Row ", i, " has ", n_mis,
             " missing values among {", y_var, ", ",
             paste(x_vars, collapse = ", "), "}.")
      }
    } 
    
    ## build full augmented dataset for this EM iteration
    aug_data <- do.call(rbind, aug_list)
    
    ## Imposing the ordering to aug_data
    for (v in c(y_var, x_vars)) {
      lev <- if (v == y_var) y_levels else x_levels_list[[v]]
      if (is_ord[v]) {
        aug_data[[v]] <- ordered(aug_data[[v]], levels = lev)
      } else {
        aug_data[[v]] <- factor(aug_data[[v]], levels = lev)
      }
    }
    
    ## ---- M-step ---- ##
    fit_new <- vglm(mod,
                    family = cumulative(parallel = TRUE),
                    data   = aug_data,
                    weights = aug_data$wgt)
    coef_new <- coef(fit_new)
    
    pi_list_new <- list()
    for (v in x_vars) {
      tab_w <- tapply(aug_data$wgt, aug_data[[v]], sum)
      tab_w[is.na(tab_w)] <- 0
      pi_list_new[[v]] <- tab_w / sum(tab_w)
    }
    
    param_new <- get_param_vec(coef_new, pi_list_new)
    eps <- max(abs(param_new - param_old))
    
    if (verbose) cat("  max parameter change =", eps, "\n")
    
    fit      <- fit_new
    coeffs   <- coef_new
    pi_list  <- pi_list_new
    param_old <- param_new
    
    if (eps < tol) {
      if (verbose) cat("Converged at iteration", iter, "\n")
      break
    }
  }
  
  coeff <- coef(fit)
  SE <- sqrt(diag(vcov(fit)))
  
  list(
    fit        = fit,
    coeff      = coeff,
    SE         = SE,
    pi_list    = pi_list,
    y_levels   = y_levels,
    x_levels   = x_levels_list,
    is_ord     = is_ord,
    iterations = iter,
    converged  = (iter < maxitr && eps < tol)
  )
}

### ------ Example ------ ##

t0 <- Sys.time()
fit_try <- EM_handmade(
  data   = NA30_list[[1]],
  y_var  = "affairs_mod",
  x_vars = c("children", "religiousness", "rating"),
  maxitr = 100, tol= 1e-8, verbose = TRUE
)

t1 <- Sys.time()
cat("EM run took", round(as.numeric(t1 - t0), 2), "seconds\n")



x_vars <- c("children", "religiousness", "rating")
y_var <- "affairs_mod"

# Function to run EM on a list of datasets with progress printing
run_EM_progress <- function(NA_list, perc_missing) {
  EM_results <- vector("list", length(NA_list))
  
  for (i in seq_along(NA_list)) {
    print(paste0("Running EM for ", perc_missing, "% missing: dataset ", i, "/", length(NA_list)))
    
    EM_results[[i]] <- EM_handmade(
      data = NA_list[[i]],
      y_var = y_var,
      x_vars = x_vars,
      maxitr = 100,
      tol = 1e-8,
      verbose = FALSE  # set to TRUE if you want EM iteration messages per dataset
    )
    
    if (i %% 50 == 0) print(paste0("Completed ", i, " datasets"))
  }
  
  EM_results
}

# Run sequentially for 5%, 15%, 30% missing
EM_results_5  <- run_EM_progress(NA5_list, 5)
EM_results_15 <- run_EM_progress(NA15_list, 15)
EM_results_30 <- run_EM_progress(NA30_list, 30)












##################################################################################
## ---- Multiple Imputation ---- ##



EM_MI <- function(data, y_var, x_vars,
                  Imp_M = 20, maxitr = 50, tol = 1e-6, verbose = TRUE) {
  
  ## 1. Run EM once
  em_out <- EM_handmade(data, y_var, x_vars,
                        maxitr = maxitr, tol = tol, verbose = verbose)
  
  fit       <- em_out$fit
  theta_hat <- em_out$coeff
  V_hat     <- vcov(fit)
  
  pi_list_hat <- em_out$pi_list
  y_levels    <- em_out$y_levels
  x_levels    <- em_out$x_levels
  is_ord      <- em_out$is_ord
  
  N <- nrow(data)
  p <- length(theta_hat)  # number of PO parameters
  
  ## analysis model
  mod <- as.formula(paste(y_var, "~", paste(x_vars, collapse = " + ")))
  
  ## storage
  imp_coef <- matrix(NA, nrow = Imp_M, ncol = p)
  imp_SE   <- matrix(NA, nrow = Imp_M, ncol = p)
  colnames(imp_coef) <- names(theta_hat)
  colnames(imp_SE)   <- names(theta_hat)
  
  for (m in seq_len(Imp_M)) {
    if (verbose) cat("Imputation", m, "...\n")
    
    ## 2(a) draw theta^(m)
    theta_m <- as.numeric(MASS::mvrnorm(1, mu = theta_hat, Sigma = V_hat))
    pi_list_m <- pi_list_hat
    
    ## 2(b) copy data and impute
    data_imp <- data
    
    for (i in seq_len(N)) {
      row_i <- data_imp[i, , drop = TRUE]
      
      mis_vec <- is.na(row_i[c(y_var, x_vars)])
      names(mis_vec) <- c(y_var, x_vars)
      n_mis <- sum(mis_vec)
      
      if (n_mis == 0) next
      
      if (n_mis > 1) {
        stop("Row ", i, " has >1 missing among {", y_var, ", ",
             paste(x_vars, collapse = ", "), "}.")
      }
      
      miss_name <- names(mis_vec)[which(mis_vec)]
      
      if (miss_name == y_var) {
        ## --- Impute missing Y | X, theta^(m) ---
        row_df <- data_imp[i, , drop = FALSE]
        p_vec  <- prob_PO(theta_m, row_df, x_vars, y_levels)
        data_imp[i, y_var] <- sample(y_levels, size = 1, prob = p_vec)
        
      } else {
        ## --- Impute missing X_j | Y, X_-j; theta^(m), pi_j^(m) ---
        j <- miss_name
        
        x_levels_j <- x_levels[[j]]
        L_j        <- length(x_levels_j)
        pi_j_m     <- as.numeric(pi_list_m[[j]])
        
        y_i   <- as.character(data_imp[i, y_var])
        y_idx <- match(y_i, y_levels)
        
        post_w <- numeric(L_j)
        
        for (l in seq_len(L_j)) {
          row_i_tmp <- data_imp[i, , drop = FALSE]
          row_i_tmp[[j]] <- factor(x_levels_j[l], levels = x_levels_j)
          
          ## Ensure all X have the correct levels
          for (v in x_vars) {
            row_i_tmp[[v]] <- factor(row_i_tmp[[v]], levels = x_levels[[v]])
          }
          
          p_vec_l    <- prob_PO(theta_m, row_i_tmp, x_vars, y_levels)
          post_w[l]  <- pi_j_m[l] * p_vec_l[y_idx]
        }
        
        post_w <- post_w / sum(post_w)
        sampled_level <- sample(x_levels_j, size = 1, prob = post_w)
        data_imp[i, j] <- sampled_level
      }
    }
    
    ## 2(c) reimpose factor / ordered structure
    data_imp[[y_var]] <- if (is_ord[y_var]) {
      ordered(data_imp[[y_var]], levels = y_levels)
    } else {
      factor(data_imp[[y_var]], levels = y_levels)
    }
    
    for (v in x_vars) {
      lev <- x_levels[[v]]
      data_imp[[v]] <- if (is_ord[v]) {
        ordered(data_imp[[v]], levels = lev)
      } else {
        factor(data_imp[[v]], levels = lev)
      }
    }
    
    ## 2(d) fit PO model on imputed data
    fit_m <- vglm(mod, family = cumulative(parallel = TRUE),
                  data = data_imp)
    
    beta_m <- coef(fit_m)
    se_m <- sqrt(diag(vcov(fit_m)))
    
    imp_coef[m, ] <- beta_m
    imp_SE[m, ]   <- se_m
  }
  
  list(
    beta_mat = imp_coef,  # M x p
    SE_mat   = imp_SE,    # M x p
    theta_hat = theta_hat,
    V_hat     = V_hat
  )
}




## ---- Example ---- ##
t0 <- Sys.time()
check <- EM_MI(NA30_list[[1]], 
               y_var = "affairs_mod",
               x_vars = c("children", "religiousness", "rating"),
               Imp_M = 2,
               maxitr = 50,
               tol = 1e-6,
               verbose = T)
t1 <- Sys.time()
cat("MI: EM run took", round(as.numeric(t1 - t0), 2), "seconds\n")
#EM run took 31.61 seconds



####################################################################

