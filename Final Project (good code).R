################################ Original dataset
set.seed(123)
library(AER)
library(glmnet)
library(dplyr)
library(ordinalNet)
library(VGAM)
data("Affairs")
data <- Affairs