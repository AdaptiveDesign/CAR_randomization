library(MASS)
library(pwr)
library(dplyr)
source('CAR_functions.R')

## Generate Data
nsample <- 210
sigma <- 15.4
power <- 0.8
d <- pwr.t.test(n = nsample/2, power = power)$d # effect size
beta_trt <- d*sigma # Difference between the means

# generate covariates
cor_ZH <- 0.12
cor_ZL <- 0.16
cor_ZD <- 0.27
cor_matrix <- matrix(c(1, cor_ZH, cor_ZL, cor_ZD,
                       cor_ZH, 1, 0, 0, 
                       cor_ZL, 0, 1, 0, 
                       cor_ZD, 0, 0, 1), ncol = 4)

data <- mvrnorm(nsample, mu = rep(0,4), Sigma = cor_matrix)
## Z
Z <- data[,1]
## Hospital-acquired infection
breakpoints <- quantile(data[,2], probs = c(0, 0.12, 1))
X_H <- cut(data[,2], breaks = breakpoints, labels = c(1,0), include.lowest = TRUE)
## Large tube size
breakpoints <- quantile(data[,3], probs = c(0, 0.13, 1))
X_L <- cut(data[,3], breaks = breakpoints, labels = c(1,0), include.lowest = TRUE)
## Drain present
breakpoints <- quantile(data[,4], probs = c(0, 0.58, 1))
X_D <- cut(data[,4], breaks = breakpoints, labels = c(1,0), include.lowest = TRUE)
covmat <- cbind(X_H, X_L, X_D)
ntrt <- 3
trtseq <- c(1,2,3)
block_size <- 6

CRD(nsample, ntrt, trtseq)
SPB(covmat, ntrt, trtseq, block_size)
SBSD(covmat, ntrt, trtseq, block_size)
Minimization(covmat, ntrt, trtseq, p_formula = 1, p = 0.7)
HDBR(covmat, ntrt, trtseq, priority = c('X_H','X_L','X_D'), limit = c(2,3,4,5))
