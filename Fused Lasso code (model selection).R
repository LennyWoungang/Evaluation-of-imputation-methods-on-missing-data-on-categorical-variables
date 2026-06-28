library(AER)
library(dplyr)
library(MASS)
library(glmnet)
library(nnet)
library(ordinalNet)
library(rms)
set.seed(123)

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

data_mod$religiousness1 <- as.integer(data_mod$religiousness >= 2)
data_mod$religiousness2 <- as.integer(data_mod$religiousness >= 3)
data_mod$religiousness3 <- as.integer(data_mod$religiousness >= 4)
data_mod$religiousness4 <- as.integer(data_mod$religiousness >= 5)

data_mod$rating1 <- as.integer(data_mod$rating >= 2)
data_mod$rating2 <- as.integer(data_mod$rating >= 3)
data_mod$rating3 <- as.integer(data_mod$rating >= 4)
data_mod$rating4 <- as.integer(data_mod$rating >= 5)


step_vars <- c(
  paste0("religiousness", 1:4),
  paste0("rating", 1:4)
)


y <- data_mod$affairs_mod
child_num <- as.numeric(ifelse(data_mod$children == "yes", 1, 0))
X <- as.matrix(data_mod[, c(step_vars)])
X <- cbind(X, child_num)


## --- Fused Lasso --- ##
fit <- ordinalNet(X, y, alpha = 1, family = "acat", link = "logit")

beta_hat <- coef(fit)
round(beta_hat, 3)


## --- data cleaning --- ##

data_mod <- data_mod %>% mutate(
  reli_mod = case_when(
    religiousness == 1 ~ 1,
    religiousness == 2 ~ 2, 
    religiousness == 3 ~ 3, 
    religiousness %in% c(4, 5) ~ 4
  ),
  reli_mod = ordered(reli_mod),
  
  rate_mod = case_when(
    rating %in% c(1, 2) ~ 1,
    rating == 3 ~ 2,
    rating == 4 ~ 3,
    rating == 5 ~ 4
  ),
  rate_mod = ordered(rate_mod)
)


######################################### Model Selection

##### polr function
str(data)
data$religiousness <- ordered(data$religiousness)
data$rating <- ordered(data$rating)
model <- polr(affairs_mod ~ (gender + children + religiousness  + rating)^2,
              data = data_polr, method = "logistic", Hess = T)
step_mod <- stepAIC(model, direction = "backward", trace = TRUE)
formula(step_mod)


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






