# Modifications to existing functions from:
# ------------------------------------------------
# The 'Minirand' R package version 0.1.3, originally developed by Man Jin, Adam Polis, 
# and Jonathan Hartzel, was used as a starting point for further development. 
# Original source code can be found at: https://github.com/cran/Minirand
# This software is distributed under the MIT License.
#
# Changes implemented in this version:
# - Implemented probability allocation methods as described by Pocock, S. J., & Simon, R. 
#   in their paper "Sequential treatment assignment with balancing for prognostic factors 
#   in the controlled clinical trial".
# - Modifications to the original 'Minirand' code to accommodate scenarios 
#   with a single factor randomization, enhancing its applicability to a broader 
#   range of clinical trial designs.
#
# Hierarchical Dynamic Balancing Randomization:
# ------------------------------------------------
# The 'Hierarchical Dynamic Balancing' randomization algorithm implemented in our code
# is inspired by the design developed in the paper by Heritier, S., Gebski, V., & Pillai, A.,
# titled 'Dynamic balancing randomization in controlled clinical trials'.
# 
# We have also developed 'Complete Randomization', 'Stratified Permuted Block 
# Randomization', and 'Stratified Big Stick Design' for comparison.
#
# The foundational work by the aforementioned authors is acknowledged as a significant contribution
# to the field. Our code represents an academic endeavor to operationalize their proposed designs
# into practical, executable algorithms, while also introducing novel adaptations.
#
# References:
# Pocock, S. J., & Simon, R. (1975). SEQUENTIAL TREATMENT ASSIGNMENT WITH BALANCING FOR 
# PROGNOSTIC FACTORS IN THE CONTROLLED CLINICAL TRIAL. Biometrics, 31, 103-115.
# Heritier, S., Gebski, V., & Pillai, A. (2005). Dynamic balancing randomization in controlled clinical trials.
# Statistics in Medicine, 24, 3729–3741. DOI: 10.1002/sim.2421.
#
# Note on Licensing:
# Users are responsible for adhering to licensing agreements and use policies of the 'Minirand' 
# package and the publications by Pocock & Simon, as well as Heritier et al. when applying 
# these modified functions in research.

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
