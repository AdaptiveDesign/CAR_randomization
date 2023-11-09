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

Minirand <- function(covmat = covmat, j, covwt = covwt, ratio = ratio, ntrt = ntrt, trtseq = trtseq, method = "Range", result = res, p)
{
  if (j > 1) {
    matchx = apply(covmat[1:(j - 1), , drop = FALSE], 1, 
                   function(x, xrow) {
                     as.numeric(x == xrow)
                   }, covmat[j, ])    
    n_matchx <- matrix(0, ncol(covmat), ntrt)
    for (k in 1:ntrt)
    {
      n_matchx[,k] <- apply(as.matrix(matchx[, result[1:(j-1)]==trtseq[k]]), 1, sum)
    }   
    
    if (method == "Range")
    {
      imbalance <- rep(0, ntrt)
      for (i in 1:ntrt)
      {
        temp <- n_matchx
        temp[,i] <- temp[, i]+1 
        num_level <- temp %*% diag(1/ratio)   
        range_level <- apply(num_level, 1, range)
        imb_margin <- range_level[2, ] - range_level[1, ]
        imbalance[i] <- sum(covwt %*% imb_margin)
      }
    }  
    
    if (method == "SD")
    {
      imbalance <- rep(0, ntrt)
      for (i in 1:ntrt)
      {
        temp <- n_matchx
        temp[, i] <- temp[, i]+1 
        num_level <- temp %*% diag(1/ratio)   
        sd_level <- apply(num_level, 1, sd)
        imb_margin <- sd_level
        imbalance[i] <- sum(covwt %*% imb_margin)
      }
    }
    
    if (method == "Var") 
    {
      imbalance <- rep(0, ntrt)
      for (i in 1:ntrt)
      {
        temp <- n_matchx
        temp[,i] <- temp[, i]+1 
        num_level <- temp %*% diag(1/ratio)   
        var_level <- apply(num_level, 1, var)
        imb_margin <- var_level
        imbalance[i] <- sum(covwt %*% imb_margin)
      }
    }
    
    trt.mini <- trtseq[imbalance == min(imbalance)]  
    trt.highprob <- trt.mini 
    trt.lowprob <- trtseq[-trt.mini]
    res <- ifelse(length(trt.highprob) < ntrt, sample(c(trt.highprob, trt.lowprob), 1, replace = TRUE, prob = c(rep(p/length(trt.highprob), length(trt.highprob)),
                                                                                                                rep((1-p)/length(trt.lowprob), length(trt.lowprob)))), sample(trtseq, 1, replace = TRUE, prob = rep(1/ntrt, ntrt)))
  }
  return(res)
}

randbalance <- function(trt, covmat, ntrt, trtseq)
{
  balance <- vector(length = ncol(covmat), "list")
  names(balance) = colnames(covmat)
  for (i in 1:ncol(covmat))
  {
    trt <- factor(trt, levels = c(1:ntrt))
    balance[[i]] <- table(trt, covmat[, i])
    
  }
  return(balance)
}

totimbal <- function (trt = trt, covmat = covmat, covwt = covwt, ratio = ratio, ntrt = ntrt, trtseq = trtseq, method = "Range")
{
  balance <- randbalance(trt, covmat, ntrt, trtseq)
  if (method == "Range")
  {
    imbalance<-rep(0, ncol(covmat)) 
    for (i in 1:ncol(covmat))
    {
      temp <- balance[[i]]
      num_level <- t(temp) %*% diag(1/ratio)   
      range_level <- apply(num_level, 1, range)
      imb_margin <- range_level[2,] - range_level[1,]
      imbalance[i] <- sum(imb_margin)
    }
    total_imb <- sum(imbalance*covwt)
  }  
  
  if (method == "SD")
  {
    imbalance <- rep(0, ncol(covmat))
    for (i in 1:ncol(covmat))
    {
      temp <- balance[[i]]
      num_level <- t(temp)%*%diag(1/ratio)   
      sd_level <- apply(num_level, 1, sd)
      imb_margin <- sd_level
      imbalance[i] <- sum(imb_margin)
    }
    total_imb <- sum(imbalance*covwt)
  }
  
  if (method == "Var")
  {
    imbalance <- rep(0, ncol(covmat))
    for (i in 1:ncol(covmat))
    {
      temp <- balance[[i]]
      num_level <- t(temp) %*% diag(1/ratio)   
      var_level <- apply(num_level, 1, var)
      imb_margin <- var_level
      imbalance[i] <- sum(imb_margin)
    }
    total_imb <- sum(imbalance*covwt)
  }
  return(total_imb)
}

Minirand_newAdd <- function(covmat = covmat, 
                            j, 
                            covwt = covwt, 
                            ratio = ratio, 
                            ntrt = ntrt, 
                            trtseq = trtseq, 
                            method = "Range", 
                            result = res, 
                            p_formula = 1,
                            p = NA,
                            q = NA,
                            t = NA)
{ 
  if (j > 1) {
    n_matchx <- matrix(0, ncol(covmat), ntrt)
    matchx = apply(covmat[1:(j - 1), , drop = FALSE], 1, 
                   function(x, xrow) {
                     as.numeric(x == xrow)
                   }, covmat[j, ]) 
    for (k in 1:ntrt)
    {
      n_matchx[,k] <- apply(as.matrix(matchx[, result[1:(j-1)]==trtseq[k]]), 1, sum)
    }   
    
    
    if (method == "Range")
    {
      imbalance <- rep(0, ntrt)
      for (i in 1:ntrt)
      {
        temp <- n_matchx
        temp[,i] <- temp[, i]+1 
        num_level <- temp %*% diag(1/ratio)   
        range_level <- apply(num_level, 1, range)
        imb_margin <- range_level[2, ] - range_level[1, ]
        imbalance[i] <- sum(covwt %*% imb_margin)
      }
    }  
    
    if (method == "SD")
    {
      imbalance <- rep(0, ntrt)
      for (i in 1:ntrt)
      {
        temp <- n_matchx
        temp[, i] <- temp[, i]+1 
        num_level <- temp %*% diag(1/ratio)   
        sd_level <- apply(num_level, 1, sd)
        imb_margin <- sd_level
        imbalance[i] <- sum(covwt %*% imb_margin)
      }
    }
    
    if (method == "Var") 
    {
      imbalance <- rep(0, ntrt)
      for (i in 1:ntrt)
      {
        temp <- n_matchx
        temp[,i] <- temp[, i]+1 
        num_level <- temp %*% diag(1/ratio)   
        var_level <- apply(num_level, 1, var)
        imb_margin <- var_level
        imbalance[i] <- sum(covwt %*% imb_margin)
      }
    }
    
    if(p_formula == 1){
      trt.mini <- trtseq[imbalance == min(imbalance)]  
      trt.highprob <- trt.mini 
      trt.lowprob <- trtseq[-trt.mini]
      res <- ifelse(length(trt.highprob) < ntrt, sample(c(trt.highprob, trt.lowprob), 1, replace = TRUE, prob = c(rep(p/length(trt.highprob), length(trt.highprob)),
                                                                                                                  rep((1-p)/length(trt.lowprob), length(trt.lowprob)))), sample(trtseq, 1, replace = TRUE, prob = rep(1/ntrt, ntrt)))
      #  print(c(trt.lowprob,trt.highprob,p))
    }
    if(p_formula == 2){
      prob <- q - 2*(ntrt*q - 1)/ (ntrt*(ntrt + 1)) * rank(imbalance)
      res <- sample(trtseq, 1, prob = prob)
    }
    if(p_formula == 3){
      prob <- 1 / (ntrt - t) * (1 - t*imbalance/sum(imbalance))
      res <- sample(trtseq, 1, prob = prob)
    }
  }
  
  ## if the assigned group has the lowest imbalance score, then set Gj = 1; otherwise, set Gj = 0
  if(imbalance[res] == min(imbalance)){
    Gj <- 1
  } else {
    Gj <- 0
  }
  
  ## Dj = |N_max - N_min|
  result[j] <- res
  res_f<- factor(result[1:j], levels = c(1:ntrt))
  Dja <- diff(range(table(res_f)))
  Dj=Dja^2
  resA<-res
  
  res<-c(resA, imbalance, Gj,  Dj)
  
  return(res)
}

Minirand_onefactor <- function(covmat = covmat, 
                               j, 
                               covwt = covwt, 
                               ratio = ratio, 
                               ntrt = ntrt, 
                               trtseq = trtseq, 
                               method = "Range", 
                               result = res, 
                               p_formula = 1,
                               p = NA,
                               q = NA,
                               t = NA)
{ 
  if (j > 1) {
    n_matchx <- matrix(0, 1, ntrt)
    matchx = as.numeric(covmat[1:j-1] == covmat[j]) 
    
    for (k in 1:ntrt)
    { 
      n_matchx[,k] <- sum(matchx[result[1:(j-1)]==trtseq[k]])
    } 
    
    
    if (method == "Range")
    {
      imbalance <- rep(0, ntrt)
      for (i in 1:ntrt)
      {
        temp <- n_matchx
        temp[,i] <- temp[, i]+1 
        num_level <- temp %*% diag(1/ratio)   
        range_level <- apply(num_level, 1, range)
        imb_margin <- range_level[2, ] - range_level[1, ]
        imbalance[i] <- sum(covwt %*% imb_margin)
      }
    }  
    
    if (method == "SD")
    {
      imbalance <- rep(0, ntrt)
      for (i in 1:ntrt)
      {
        temp <- n_matchx
        temp[, i] <- temp[, i]+1 
        num_level <- temp %*% diag(1/ratio)   
        sd_level <- apply(num_level, 1, sd)
        imb_margin <- sd_level
        imbalance[i] <- sum(covwt %*% imb_margin)
      }
    }
    
    if (method == "Var") 
    {
      imbalance <- rep(0, ntrt)
      for (i in 1:ntrt)
      {
        temp <- n_matchx
        temp[,i] <- temp[, i]+1 
        num_level <- temp %*% diag(1/ratio)   
        var_level <- apply(num_level, 1, var)
        imb_margin <- var_level
        imbalance[i] <- sum(covwt %*% imb_margin)
      }
    }
    
    if(p_formula == 1){
      trt.mini <- trtseq[imbalance == min(imbalance)]  
      trt.highprob <- trt.mini 
      trt.lowprob <- trtseq[-trt.mini]
      res <- ifelse(length(trt.highprob) < ntrt, sample(c(trt.highprob, trt.lowprob), 1, replace = TRUE, prob = c(rep(p/length(trt.highprob), length(trt.highprob)),
                                                                                                                  rep((1-p)/length(trt.lowprob), length(trt.lowprob)))), sample(trtseq, 1, replace = TRUE, prob = rep(1/ntrt, ntrt)))
      #  print(c(trt.lowprob,trt.highprob,p))
    }
    if(p_formula == 2){
      prob <- q - 2*(ntrt*q - 1)/ (ntrt*(ntrt + 1)) * rank(imbalance)
      res <- sample(trtseq, 1, prob = prob)
    }
    if(p_formula == 3){
      prob <- 1 / (ntrt - t) * (1 - t*imbalance/sum(imbalance))
      res <- sample(trtseq, 1, prob = prob)
    }
  }
  
  ## if the assigned group has the lowest imbalance score, then set Gj = 1; otherwise, set Gj = 0
  if(imbalance[res] == min(imbalance)){
    Gj <- 1
  } else {
    Gj <- 0
  }
  
  ## Dj = |N_max - N_min|
  result[j] <- res
  res_f<- factor(result[1:j], levels = c(1:ntrt))
  Dja <- diff(range(table(res_f)))
  Dj=Dja^2
  resA<-res
  
  res<-c(resA, imbalance, Gj,  Dj)
  
  return(res)
}

CRD <- function(nsample, ntrt, trtseq){
  res <- sample(x = trtseq, size = nsample, replace = TRUE)
  resguess <- sample(x = trtseq, size = nsample, replace = TRUE)
  Gidiff <- c()
  for(k in 1:nsample){
    Gidiff[k] <- diff(range(table(res[1:k])))
  }
  abImba<-mean(Gidiff^2/(1:nsample))
  abRand<-(mean(res == resguess) - 1/ntrt) * (ntrt/(ntrt - 1))
  
  return(list(res = factor(res),
              abImba = abImba,
              abRand = abRand))
}

Minimization <- function(covmat,
                         ntrt,
                         trtseq,
                         method = 'Range',
                         p_formula = 1,
                         p,
                         q,
                         t)
{
  if(is.null(ncol(covmat))){
    nsample <- length(covmat)
    PS_df <- data.frame(matrix(nrow = nsample, ncol = ntrt + 3))
    colnames(PS_df) <- c('res', paste0('IB', 1:ntrt), 'Gi', 'Di')
    PS_df$res[1] <- sample(trtseq, 1)
    PS_df$Di[1] <- 1
    Gi1 <- sample(trtseq, 1)
    PS_df$Gi[1]<-0
    if(Gi1 == PS_df$res[1]) {PS_df$Gi[1]<-1}
    for(j in 2:nsample){
      PS_df[j,] <- Minirand_onefactor(covmat=covmat, 
                                      j, 
                                      covwt=1, 
                                      ratio=rep(1/ntrt, ntrt), 
                                      ntrt=ntrt, 
                                      trtseq=trtseq, 
                                      method=method, 
                                      result=PS_df$res,
                                      p_formula = p_formula,
                                      p = p,
                                      q = q,
                                      t = t)
    }
  } else {
    nsample <- nrow(covmat)
    PS_df <- data.frame(matrix(nrow = nsample, ncol = ntrt + 3))
    colnames(PS_df) <- c('res', paste0('IB', 1:ntrt), 'Gi', 'Di')
    PS_df$res[1] <- sample(trtseq, 1)
    PS_df$Di[1] <- 1
    Gi1 <- sample(trtseq, 1)
    PS_df$Gi[1]<-0
    if(Gi1 == PS_df$res[1]) {PS_df$Gi[1]<-1}
    for(j in 2:nsample){
      PS_df[j,] <- Minirand_newAdd(covmat=covmat, 
                                   j, 
                                   covwt=rep(1/ncol(covmat), ncol(covmat)), 
                                   ratio=rep(1/ntrt, ntrt), 
                                   ntrt=ntrt, 
                                   trtseq=trtseq, 
                                   method=method, 
                                   result=PS_df$res,
                                   p_formula = p_formula,
                                   p = p,
                                   q = q,
                                   t = t)
    }
  }
  abImba <- mean(PS_df$Di/(1:nsample))
  abRand <- (ntrt * mean(PS_df$Gi) - 1) / (ntrt - 1)
  return(list(res = factor(PS_df$res),
              abImba = abImba,
              abRand = abRand))
}

DBR <- function(covmat,
                j,
                ntrt,
                trtseq,
                result,
                priority,
                limit)
{
  if(j == 1){
    res <- sample(trtseq, 1) # The result
    guess <- sample(trtseq, 1) # The guess
    return(list(res = res, guess = guess))
  } else {
    for (k in 1:length(priority)) {
      level_k <- covmat[j, priority[k]] 
      result_k <- result[1:j-1][covmat[1:(j-1), priority[k]] == level_k]
      result_k <- factor(result_k, levels = trtseq)
      count_k <- table(result_k)
      D_k <- max(count_k) - min(count_k)
      if (D_k >= limit[k]){
        ## Find the treatment groups which have the minimum number of subjects allocated to
        mini_group <- names(count_k[count_k == min(count_k)])
        res <- sample(mini_group, 1) # The result
        guess <- sample(mini_group, 1) # The guess
        return(list(res = res, guess = guess))
      }
    }
    result_overall <- result[1:j-1]
    result_overall <- factor(result_overall, levels = trtseq)
    count_overall <- table(result_overall)
    D_overall <- max(count_overall) - min(count_overall)
    if(D_overall >= limit[length(limit)]) {
      mini_group <- names(count_overall[count_overall == min(count_overall)])
      res <- sample(mini_group, 1)
      guess <- sample(mini_group, 1)
      return(list(res = res, guess = guess))
    } else {
      res <- sample(trtseq, 1)
      guess <- sample(trtseq, 1)
      return(list(res = res, guess = guess))
    }
  }
}

HDBR <- function(covmat, 
                 ntrt, 
                 trtseq,
                 priority,
                 limit){
  if(is.null(nrow(covmat))) {
    covmat <- as.data.frame(covmat)
    colnames(covmat) <- priority 
  }
  nsample <- nrow(covmat)
  res <- c()
  resguess <- c()
  for(j in 1:nsample){
    res_list <- DBR(covmat, j, ntrt, trtseq, res, priority, limit)
    res <- c(res, res_list$res)
    resguess <- c(resguess, res_list$guess)
  }
  Gidiff <- c()
  for(k in 1:nsample){
    Gidiff[k] <- diff(range(table(res[1:k])))
  }
  abImba<-mean(Gidiff^2/(1:nsample))
  abRand<-(mean(res == resguess) - 1/ntrt) * (ntrt/(ntrt - 1))
  return(list(res = factor(res),
              abImba = abImba,
              abRand = abRand))
}

# FUN for GUESS

GUESS <- function(pre.allo, i, b, npre, ntrt, n = 1) {
  P.allo <- numeric(ntrt)
  for (t in 1:ntrt) {
    P.allo[t] <- ((b/sum(pre.allo == t)) * (1 + floor((i-1)/b)) - npre[t]) / (b * (1 + floor((i-1)/b)) - (i-1))
  }
  
  max_indices <- which(P.allo == max(P.allo))
  if (length(max_indices) > 1) {
    max_index <- sample(max_indices, 1)
  } else {
    max_index <- max_indices
  }
  
  guess <- rep(0, ntrt)
  guess[max_index] <- max_index
  
  return(list(P.allo = P.allo, guess = guess, pre.allo = pre.allo, npre = npre))
}



RESguess <- function(nsample, pre.allo.true, b, npre, ntrt) {
  resguess <- vector("list", nsample)
  for (i in 1:nsample) {
    if (i == 1) {
      pre.allo <- integer(0)
      npre <- rep(0, ntrt)
    }
    if (i > 1) {
      pre.allo <- pre.allo.true[1:(i - 1)]
      npre <- sapply(1:ntrt, function(x) sum(pre.allo == x))
    }
    GUESSinfo <- GUESS(pre.allo, i, b, npre, ntrt)
    resguess[[i]] <- c(which(GUESSinfo$guess != 0))
  }
  resguess <- unlist(resguess)
  return(resguess)
}

# Function for Stratified Permuted Block Randomization
stratified_permuted_block_randomization <- function(dataset, strata_var, block_size, num_treatments) {
  # Create a function to perform permuted block randomization within each stratum
  permuted_block_randomization <- function(stratum_data) {
    # Randomly permute the treatment assignment within blocks
    num_obs <- nrow(stratum_data)
    num_blocks <- ceiling(num_obs / block_size)
    block_assignment <- rep(1:num_blocks, each = block_size)[1:num_obs]
    # block_assignment <- sample(block_assignment)
    treatment_order <- rep(NA, num_obs)
    resguess <- rep(NA, num_obs)
    for (i in 1:num_blocks) {
      block_indices <- which(block_assignment == i)
      if(length(block_indices) == block_size){
        block_treatments <- sample(rep(1:num_treatments, each = block_size / num_treatments))
        #       guess <- RESguess(block_size, pre.allo=c(), b=6, npre=c(0,0,0), ntrt=3)
        pre.allo.true<-block_treatments
        guess<-RESguess(block_size,pre.allo.true,b=block_size,npre,num_treatments)
      } else {
        block_treatments <- sample(num_treatments, length(block_indices), replace = TRUE)
        #       guess<-RESguess(length(block_indices), pre.allo=c(), b=length(block_indices), npre=c(0,0,0), ntrt=3)
        pre.allo.true<-block_treatments
        guess<-RESguess(length(block_indices),pre.allo.true,b=length(block_indices),npre,num_treatments)
      }
      treatment_order[block_indices] <- block_treatments
      resguess[block_indices]<-guess
    }
    # Create a treatment vector based on the permutation
    stratum_data$treatment <- treatment_order
    stratum_data$resguess <- resguess
    return(stratum_data)
  }
  
  # Perform stratified permuted block randomization
  randomized_dataset <- dataset %>%
    group_by(across(all_of(strata_var))) %>%
    do(permuted_block_randomization(.))
  
  tmt <- randomized_dataset[order(randomized_dataset$participant_id),'treatment']
  resguess <- randomized_dataset[order(randomized_dataset$participant_id), 'resguess']
  
  return(list("tmt" = unlist(tmt), "resguess" = unlist(resguess)))
  #return(unlist(tmt))
}

SPB <- function(covmat, ntrt, trtseq, block_size){
  covmat_df <- as.data.frame(covmat)
  covmat_df$stratum <- interaction(covmat_df)
  nsample <- nrow(covmat_df)
  participant_id <- 1:nsample
  dataset <- data.frame(participant_id, covmat_df)
  num_treatments <- ntrt
  res_list <- stratified_permuted_block_randomization(dataset, "stratum", block_size, num_treatments)
  res <- res_list$tmt
  names(res) <- NULL
  resguess <- res_list$resguess
  Gidiff <- c()
  for(k in 1:nsample){
    Gidiff[k] <- diff(range(table(res[1:k])))
  }
  abImba<-mean(Gidiff^2/(1:nsample))
  abRand<-(mean(res == resguess) - 1/ntrt) * (ntrt/(ntrt - 1))
  res <- factor(res, levels = 1:ntrt, labels = trtseq)
  return(list(res = res,
              abImba = abImba,
              abRand = abRand))
}

# Function for stratified BSD

check_pair_available <- function(urn) {
  count1 <- sum(urn == 1)
  count2 <- sum(urn == 2)
  count3 <- sum(urn == 3)
  if (count1 >= 1 && count2 >= 1 && count3 >= 1) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

stratified_BSD <- function(dataset, strata_var, block_size, num_treatments) {
  # Perform big stick randomization within each stratum
  BSD_randomization <- function(stratum_data) {
    n_obs <- nrow(stratum_data)
    trt_assign <- numeric(n_obs)
    active_urn <- vector("list", n_obs+1)
    inactive_urn <- vector("list", n_obs+1)
    resguess <- rep(NA, n_obs)
    
    for (i in 1:(n_obs+1)) {
      if (i == 1) {
        #       num_treatments <- 3
        #       block_size <- 6
        active_urn[[i]] <- sample(rep(1:num_treatments, each = block_size/num_treatments))
        inactive_urn[[i]] <- c()
      } else {
        if (check_pair_available(active_urn[[i - 1]])) {
          available_indices <- c(which(active_urn[[i - 1]] == 1)[1], which(active_urn[[i - 1]] == 2)[1], which(active_urn[[i - 1]] == 3)[1])
          selected_index <- sample(available_indices, 1)
          trt_assign[i-1] <- active_urn[[i - 1]][selected_index]
          active_urn[[i]] <- active_urn[[i - 1]][-selected_index]
          inactive_urn[[i]] <- c(inactive_urn[[i - 1]], trt_assign[i-1])
          
          if (check_pair_available(inactive_urn[[i]])) {
            pair_indices <- c(which(inactive_urn[[i]] == 1)[1], which(inactive_urn[[i]] == 2)[1], which(inactive_urn[[i]] == 3)[1])
            pair <- inactive_urn[[i]][pair_indices]
            inactive_urn[[i]] <- inactive_urn[[i]][-pair_indices]
            active_urn[[i]] <- c(active_urn[[i]], pair)
          } else {
            inactive_urn[[i]] <- inactive_urn[[i]]
            active_urn[[i]] <- active_urn[[i]]
          }
          
        } else {
          available_indices <- c(which(active_urn[[i - 1]] == 1)[1], which(active_urn[[i - 1]] == 2)[1], which(active_urn[[i - 1]] == 3)[1])
          selected_index <- sample(available_indices[!is.na(available_indices)], 1)
          trt_assign[i-1] <- active_urn[[i - 1]][selected_index]
          active_urn[[i]] <- active_urn[[i - 1]][-selected_index]
          inactive_urn[[i]] <- c(inactive_urn[[i - 1]], trt_assign[i-1])
          
          if (check_pair_available(inactive_urn[[i]])) {
            pair_indices <- c(which(inactive_urn[[i]] == 1)[1], which(inactive_urn[[i]] == 2)[1], which(inactive_urn[[i]] == 3)[1])
            pair <- inactive_urn[[i]][pair_indices]
            inactive_urn[[i]] <- inactive_urn[[i]][-pair_indices]
            active_urn[[i]] <- c(active_urn[[i]], pair)
          } else {
            inactive_urn[[i]] <- inactive_urn[[i]]
            active_urn[[i]] <- active_urn[[i]]
          }
        }
      }
      
      
      pre.allo.true<-trt_assign
      resguess<-RESguess(n_obs,pre.allo.true,b=block_size,npre,3)
    }
    stratum_data$treatment <- trt_assign
    stratum_data$resguess <- resguess
    return(stratum_data)
  }
  # Perform stratified big stick randomization
  randomized_dataset <- dataset %>%
    group_by(across(all_of(strata_var))) %>%
    do(BSD_randomization(.))
  
  tmt <- randomized_dataset[order(randomized_dataset$participant_id),'treatment']
  resguess <- randomized_dataset[order(randomized_dataset$participant_id), 'resguess']
  
  return(list("tmt" = unlist(tmt), "resguess" = unlist(resguess)))
}

SBSD <- function(covmat, ntrt, trtseq, block_size){
  covmat_df <- as.data.frame(covmat)
  covmat_df$stratum <- interaction(covmat_df)
  nsample <- nrow(covmat_df)
  participant_id <- 1:nsample
  dataset <- data.frame(participant_id, covmat_df)
  num_treatments <- ntrt
  res_list <- stratified_BSD(dataset, "stratum", block_size, num_treatments)
  res <- res_list$tmt
  names(res) <- NULL
  resguess <- res_list$resguess
  Gidiff <- c()
  for(k in 1:nsample){
    Gidiff[k] <- diff(range(table(res[1:k])))
  }
  abImba<-mean(Gidiff^2/(1:nsample))
  abRand<-(mean(res == resguess) - 1/ntrt) * (ntrt/(ntrt - 1))
  res <- factor(res, levels = 1:ntrt, labels = trtseq)
  return(list(res = res,
              abImba = abImba,
              abRand = abRand))
}
