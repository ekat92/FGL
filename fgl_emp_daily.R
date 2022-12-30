# This script produces Tables 1-2 of the empirical application (daily data)

######################################################################################
# load packages
if (!require('huge')) install.packages('huge'); library(huge)
if (!require('Matrix')) install.packages('Matrix'); library(Matrix)
if (!require('zoo')) install.packages('zoo'); library(zoo)
if (!require('glmnet')) install.packages('glmnet'); library(glmnet)
if (!require('MASS')) install.packages('MASS'); library(MASS)
if (!require('dplyr')) install.packages('dplyr'); library(dplyr)
if (!require('CVXR')) install.packages('CVXR'); library(CVXR)
if (!require('nlshrink')) install.packages('nlshrink'); library(nlshrink)
if (!require('tidyverse')) install.packages('tidyverse'); library(tidyverse)
if (!require('tidyr')) install.packages('tidyr'); library(tidyr)
if (!require('POET')) install.packages('POET'); library(POET)
if (!require('fastclime')) install.packages('fastclime'); library(fastclime)
if (!require('glassoFast')) install.packages('glassoFast'); library(glassoFast)
if (!require('igraph')) install.packages('igraph'); library(igraph)

##############################FUNCTIONS###############################################
portfolios <- function(X,estsigm){
  mu = colMeans(X)
  # Compute preliminary values
  p <- length(mu)
  one <- matrix(1,p,1)
  phi <- as.single(t(one) %*% estsigm %*% one)
  psi <- as.single(t(one) %*% estsigm %*% mu)
  Phi <- as.single(t(mu) %*% estsigm %*% mu)
  
  # Global min variance
  gmv <- (estsigm %*% one) / phi 
  gmvvar <- 1 / phi
  
  # Markowitz 
  denom  <- phi*Phi - psi^2
  mkw    <- (Phi-target*psi)/denom * estsigm %*% one + (target*phi-psi)/denom * estsigm%*%mu 
  mkwvar <- (phi*target^2 - 2*psi*target + Phi)/denom
  
  #Fan
  fanw <- targetsigma*(Phi)^(-0.5) * (estsigm %*% mu)
  
  ptfret <- list(cbind(gmv),cbind(mkw),cbind(fanw),c(gmvvar),c(mkwvar))
  
  return(ptfret)
}

rolling <- function(excess,weight) {
  rets <- excess[(pred+1):nrow(excess),]
  portf_ret <- NULL
  for (k in 1:R) {
    set11 <- zoo(rets[k,StockInd[k,]])
    set12 <- na.locf(set11, fromLast = FALSE)
    set12 <- as.matrix(set12)
    rets[k,StockInd[k,]] = set12
    portf_ret[k] = t(as.matrix(na.omit(weight[k,])))%*%as.matrix(rets[k,StockInd[k,]])
  }
  mur <- mean(portf_ret)
  # out-of-sample variance 
  sigma2  <- var(portf_ret) 
  # out-of-sample sharpe ratio 
  sharpe  <- mur/sqrt(sigma2)
  #compute portfolio total return
  rf_new <- rf[(pred+1):nrow(excess)]
  rf_new <- as.matrix(rf_new)
  wexret1 <- matrix(c(NA), R, ncol(excess))  #R*p matrix
  for (k in 1:R) {
    # set11 <- zoo(rets[k,StockInd[k,]])
    # set12 <- na.locf(set11, fromLast = FALSE)
    # rets[k,StockInd[k,]] <- set12
    wexret1[k,1:length(rets[k,StockInd[k,]])] = as.matrix(t(rets[k,StockInd[k,]]))
    
    for (i in 1:length(wexret1[k,])) {
      wexret1[k,i]=wexret1[k,i]+rf_new[k]
    }
  }
  portret <- NULL
  for (k in 1:R) {
    portret[k] = t(as.matrix(na.omit(weight[k,])))%*%as.matrix(na.omit(wexret1[k,]))
  }
  
  retra2 <- rets #R*p matrix
  for (k in 1:R) {
    set11 <- zoo(rets[k,StockInd[k,]])
    set12 <- na.locf(set11, fromLast = FALSE)
    rets[k,StockInd[k,]] <- set12
    retra2[k,1:length(rets[k,StockInd[k,]])] = 1+as.matrix(t(rets[k,StockInd[k,]]))
    
    for (i in 1:length(retra2[k,])) {
      retra2[k,i]=retra2[k,i]/(1+portret[k])
    }
  }
  # portfolio weight before rebalancing
  wmins <- weight * retra2 #elementwise product, R*p
  # portfolio weight after rebalancing
  wtpl  <- weight[2:R,]
  wmin  <- wmins[1:(R-1),]
  
  turn = matrix((rep(NA, (R-1)*1)), nrow=(R-1))
  # compute turnover rate vector
  for (i in 1:(R-1)) {
    turn[i] = sum(abs(wtpl[i,] - wmin[i,]), na.rm = TRUE)  # colSums(abs(wtpl - wmin))
  }
  # turnover rate for all assets
  turnover <- mean(turn)
  
  # return with transaction cost
  reventc        <- matrix((rep(NA, (R-1))))
  for (i in 1:length(reventc)){
    reventc[i]     <- portf_ret[i] - (c * (1+portf_ret[i])) * (sum(abs(wtpl[i,] - wmin[i,]), na.rm = TRUE))
  }
  # oos return with transaction cost
  mutc = mean(reventc)
  # oos variance with transaction cost
  sigtc = var(reventc)
  # oos sharpe with transaction cost
  shartc = mutc/sqrt(sigtc)
  
  
  perest <- list(c(mur,sigma2,sharpe,turnover),c(mutc,sigtc,shartc,NA))
  
  return(perest)
}

SYMPOSDEF <- function(A, r){ #r is number of rounding digits
  for (i in 1:nrow(A)) {
    for (j in 1:ncol(A)) {
      if (i != j){
        if (abs(A[i,j])<=abs(A[j,i])){
          A[i,j] = A[i,j]
        } else {
          A[i,j] = A[j,i]
        }
      }
    }
  }
  svd <- eigen(A)
  D = diag(svd[["values"]], ncol(A))
  V = svd[["vectors"]]
  d = svd[["values"]]
  for (k in 1:length(d)) {
    if(d[k]<=0){
      d[k]=d[k-1]
    }    
  }
  D = diag(d, ncol(A))
  A <- V%*%D%*%t(V)
  A = round(A,r)
  return(A)
}

######################################################################################
####Daily Data#####
FF5_DAILY <- FF5_DAILY[11:5050,]
FF5_DAILY <- as.matrix(FF5_DAILY)
FF1 <- FF5_DAILY[,1]
FF1 <- as.matrix(FF1)
FF3 <- FF5_DAILY[,1:3]
FF5 <- FF5_DAILY[,1:5]
rf <- FF5_DAILY[,6]
adj_returns2 <- as.matrix(daily_excess)
#####Make time period of daily data divisible by 21###
snp_excess <- snp_excess[11:5050,]
snp_excess <- as.matrix(snp_excess)

c=0.001 #transaction cost per stock
##DAILY TARGETS##
target <- 0.000379 #((1+0.000379)^(252)-1=0.1)
targetsigma <- 0.013

R  <-4536 #18 years (daily data)
# R <- 63 #3 months
nvec = nrow(adj_returns2)
pred <- nvec-R 
Rseq <- seq(1,R,21)

StockInd = matrix(c(0), R, ncol(adj_returns2)) # each row saves indexes of selected columns in adj_reutrns2
weightFF1NS_global  <- matrix(c(NA), R, ncol(adj_returns2))
weightFF1NS_mkw  <- matrix(c(NA), R, ncol(adj_returns2))
weightFF1NS_fan  <- matrix(c(NA), R, ncol(adj_returns2))

weightFF3NS_global  <- matrix(c(NA), R, ncol(adj_returns2))
weightFF3NS_mkw  <- matrix(c(NA), R, ncol(adj_returns2))
weightFF3NS_fan  <- matrix(c(NA), R, ncol(adj_returns2))

weightFF5NS_global  <- matrix(c(NA), R, ncol(adj_returns2))
weightFF5NS_mkw  <- matrix(c(NA), R, ncol(adj_returns2))
weightFF5NS_fan  <- matrix(c(NA), R, ncol(adj_returns2))

weightFF1_global  <- matrix(c(NA), R, ncol(adj_returns2))
weightFF1_mkw  <- matrix(c(NA), R, ncol(adj_returns2))
weightFF1_fan  <- matrix(c(NA), R, ncol(adj_returns2))
# 

weightFF3_global  <- matrix(c(NA), R, ncol(adj_returns2))
weightFF3_mkw  <- matrix(c(NA), R, ncol(adj_returns2))
weightFF3_fan  <- matrix(c(NA), R, ncol(adj_returns2))

weightFF5_global  <- matrix(c(NA), R, ncol(adj_returns2))
weightFF5_mkw  <- matrix(c(NA), R, ncol(adj_returns2))
weightFF5_fan  <- matrix(c(NA), R, ncol(adj_returns2))

weightFGL_global  <- matrix(c(NA), R, ncol(adj_returns2))
weightFGL_mkw  <- matrix(c(NA), R, ncol(adj_returns2))
weightFGL_fan  <- matrix(c(NA), R, ncol(adj_returns2))

weightLW_global  <- matrix(c(NA), R, ncol(adj_returns2))
weightLW_mkw  <- matrix(c(NA), R, ncol(adj_returns2))
weightLW_fan  <- matrix(c(NA), R, ncol(adj_returns2))

weightNLW_global<- matrix(c(NA), R, ncol(adj_returns2))
weightNLW_mkw<- matrix(c(NA), R, ncol(adj_returns2))
weightNLW_fan<- matrix(c(NA), R, ncol(adj_returns2))

weightPOETproj_global<- matrix(c(NA), R, ncol(adj_returns2))
weightPOETproj_mkw<- matrix(c(NA), R, ncol(adj_returns2))
weightPOETproj_fan<- matrix(c(NA), R, ncol(adj_returns2))

weightPOET_global  <- matrix(c(NA), R, ncol(adj_returns2))
weightPOET_mkw  <- matrix(c(NA), R, ncol(adj_returns2))
weightPOET_fan  <- matrix(c(NA), R, ncol(adj_returns2))

weightCLIME_global  <- matrix(c(NA), R, ncol(adj_returns2))
weightCLIME_mkw  <- matrix(c(NA), R, ncol(adj_returns2))
weightCLIME_fan  <- matrix(c(NA), R, ncol(adj_returns2))
###########################
counter1 =0 
for(j in Rseq)  {
  print(j)
  counter1 = counter1+1
  set1 = adj_returns2
  # set1 = adj_returns2
  list1 = c() #saves indexes of not-NA columns
  for (i in 1:ncol(adj_returns2)){
    if (!is.na(set1[j,i]) & !is.na(set1[(j+503),i]) ){
      list1 = c(list1,i)
    }
  }
  StockInd[counter1,1:length(list1)] = list1
  set1 = set1[,list1]
  set1 = set1[j:(j+503),]
  set11 <- zoo(set1)
  set12 <- na.locf(set11, fromLast = FALSE)
  set1 <- as.matrix(set12)
  
  n <- nrow(set1)
  p <- ncol(set1)
  #########Estimated Factors##########
  k=POETKhat(t(set1))
  k= k[["K1BN"]]
  betas <- eigen(t(set1)%*%set1)[["vectors"]][,1:k]
  covariate <- as.matrix(prcomp(set1)[["x"]][,(1:k)]) ##these are factors like Fhat
  betas = as.matrix(betas) # p*K
  betas = t(betas) #K*p
  residuals <- set1 - covariate%*%betas
  residuals <- as.matrix(residuals)
  
  Y <- residuals
  n <- nrow(Y)
  p <- ncol(Y)
  covSET1 = linshrink_cov(set1, k = 0)
  covU = linshrink_cov(Y, k = 0)

  ##################################
  #######################################
  #####When factors are known###########
  ###FF1
  covariate1 <- FF1[j:(pred+j-1),]
  covariate1 <- as.matrix(covariate1) #needed if use 1 factor
  FactorModel1 <- lm(formula = set1~covariate1) #with intercept
  betas1 <- FactorModel1$coefficients #k*pred
  betas1 = betas1[2,]
  betas1 = as.matrix(betas1) # p*1 if 1 factor, or K*p if more than 1 factor
  betas1 = t(betas1) #1*p If use more than 1 factor dont need this step
  residuals1 <- FactorModel1$residuals
  # 
  ###FF3
  covariate2 <- FF3[j:(pred+j-1),]
  FactorModel2 <- lm(formula = set1~covariate2) #with intercept
  betas2 <- FactorModel2$coefficients #k*pred
  betas2 = betas2[2:4,]
  betas2 = as.matrix(betas2) # p*1 if 1 factor, or K*p if more than 1 factor
  residuals2 <- FactorModel2$residuals
  # 
  ###FF5
  covariate3 <- FF5[j:(pred+j-1),]
  FactorModel3 <- lm(formula = set1~covariate3) #with intercept
  betas3 <- FactorModel3$coefficients #k*pred
  betas3 = betas3[2:6,]
  betas3 = as.matrix(betas3) # p*1 if 1 factor, or K*p if more than 1 factor
  residuals3 <- FactorModel3$residuals
  # 
  ##########FFs WITHOUT SPARSITY RESTRICTION
  ###FF1
  prec1 = ginv(cov(residuals1))
  
  Theta1 = prec1 - prec1%*%t(betas1)%*%solve( (var(covariate1))^(-1)+
                                                betas1%*%prec1%*%t(betas1))%*%betas1%*%prec1
  w1 <- portfolios(set1, estsigm = Theta1)
  
  ###FF3
  prec2 = ginv(cov(residuals2))
  Theta2 = prec2 - prec2%*%t(betas2)%*%solve( solve(cov(covariate2))+
                                                betas2%*%prec2%*%t(betas2))%*%betas2%*%prec2
  w2 <- portfolios(set1, estsigm = Theta2)
  
  ###FF5
  prec3 = ginv(cov(residuals3))
  Theta3 = prec3 - prec3%*%t(betas3)%*%solve( solve(cov(covariate3))+
                                                betas3%*%prec3%*%t(betas3))%*%betas3%*%prec3
  w3 <- portfolios(set1, estsigm = Theta3)
  
 ##########################################
  
  ###FF1
  Y1 <- residuals1
  covU1 = linshrink_cov(Y1, k = 0)
  fwgl1 <- rglasso(covU1,weight = TRUE, N=nrow(set1)-1)
  ThetaFWGL1_1 = fwgl1[[2]]
  
  ThetaFWGL2_1 = ThetaFWGL1_1 - ThetaFWGL1_1%*%t(betas1)%*%solve( (var(covariate1))^(-1)+
                                                                    betas1%*%ThetaFWGL1_1%*%t(betas1))%*%betas1%*%ThetaFWGL1_1
  wFGL1 <- portfolios(set1, estsigm = ThetaFWGL2_1)
  
  ###FF3
  Y2 <- residuals2
  covU2 = linshrink_cov(Y2, k = 0)
  fwgl2 <- rglasso(covU2,weight = TRUE, N=nrow(set1)-1)
  ThetaFWGL1_2 = fwgl2[[2]]
  
  ThetaFWGL2_2 = ThetaFWGL1_2 - ThetaFWGL1_2%*%t(betas2)%*%solve( solve(cov(covariate2))+
                                                                    betas2%*%ThetaFWGL1_2%*%t(betas2))%*%betas2%*%ThetaFWGL1_2
  wFGL2 <- portfolios(set1, estsigm = ThetaFWGL2_2)
  ###FF5
  Y3 <- residuals3
  covU3 = linshrink_cov(Y3, k = 0)
  fwgl3 <- rglasso(covU3,weight = TRUE, N=nrow(set1)-1)
  ThetaFWGL1_3 = fwgl3[[2]]
  
  ThetaFWGL2_3 = ThetaFWGL1_3 - ThetaFWGL1_3%*%t(betas3)%*%solve( solve(cov(covariate3))+
                                                                    betas3%*%ThetaFWGL1_3%*%t(betas3))%*%betas3%*%ThetaFWGL1_3
  wFGL3 <- portfolios(set1, estsigm = ThetaFWGL2_3)
  
  ########WEIGHTED GLASSO###########
  
  ###If only need factor-based models remove three lines below###
  wgl <- rglasso(covSET1, lambda = "bic", N=nrow(set1)-1)
  ThetaWGL = wgl[[2]]
  
  fwgl <- rglasso(covU,weight = TRUE, N=nrow(set1)-1)
  ThetaFWGL1 = fwgl[[2]]
 
  ##If number of factors =1
  if (k==1){
    ThetaFWGL2 = ThetaFWGL1 - ThetaFWGL1%*%t(betas)%*%solve( (var(covariate))^(-1)+
                                                               betas%*%ThetaFWGL1%*%t(betas))%*%betas%*%ThetaFWGL1
  } else {
    ThetaFWGL2 = ThetaFWGL1 - ThetaFWGL1%*%t(betas)%*%solve( solve(cov(covariate))+
                                                               betas%*%ThetaFWGL1%*%t(betas))%*%betas%*%ThetaFWGL1
  }
  ########LW#############
  if (k==1){
    LWTheta = solve(covU) -  solve(covU)%*%t(betas)%*%solve( (var(covariate))^(-1)+
                                                               betas%*% solve(covU)%*%t(betas))%*%betas%*% solve(covU)
  } else {
    LWTheta = solve(covU) -  solve(covU)%*%t(betas)%*%solve( solve(cov(covariate))+
                                                               betas%*% solve(covU)%*%t(betas))%*%betas%*% solve(covU)
  }
  
  covU_NL=nlshrink_cov(Y, k = 0)
  ########NLW#############
  NLWTheta = ginv(covU_NL) -  ginv(covU_NL)%*%t(betas)%*%solve( solve(cov(covariate))+
                                                                  betas%*% ginv(covU_NL)%*%t(betas))%*%betas%*% ginv(covU_NL)
  
  ############PROJECTED POET##################
  
  covPOET2=POET(t(set1),K=k,0.5,"soft","vad")$SigmaY
  ThetaPOET_nd=solve(covPOET2)#INVERSE MATRIX PRECISION ESTIMATE
  ThetaPOET2 = SYMPOSDEF(ThetaPOET_nd, 5)
  ############POET##################
  
  covPOET=POET(t(set1),K=k,0.5,"soft","vad")$SigmaY
  #covY are sample covariancematrix of returns by POET,0.5,SOFT,VAD ARE DEFAULTS
  ThetaPOET=solve(covPOET)#INVERSE MATRIX PRECISION ESTIMATE
  
  ##########CLIME########
  clime_est = fastclime(Y)
  out2 = fastclime.selector(clime_est$lambdamtx, clime_est$icovlist,0.2)
  ClimeTheta = out2$icov
  if (k==1){
    ThetaFClime = ClimeTheta - ClimeTheta%*%t(betas)%*%solve( (var(covariate))^(-1)+
                                                                betas%*%ClimeTheta%*%t(betas))%*%betas%*%ClimeTheta
  } else {
    ThetaFClime = ClimeTheta - ClimeTheta%*%t(betas)%*%solve( solve(cov(covariate))+
                                                                betas%*%ClimeTheta%*%t(betas))%*%betas%*%ClimeTheta
  }
  
  # Weight vector (R*p)
  wFGL <- portfolios(set1, estsigm = ThetaFWGL2)
  wLW <- portfolios(set1, estsigm = LWTheta)
  wPOET <- portfolios(set1, estsigm = ThetaPOET)
  wCLIME <- portfolios(set1, estsigm = ThetaFClime)
  wNLW <- portfolios(set1, estsigm = NLWTheta)
  wPOET_proj <- portfolios(set1, estsigm = ThetaPOET2)
  
  # ####ADD IT BACK LATER###
  
  weightFGL_global[(counter1),1:length(list1)]  <- t(wFGL[[1]])
  weightFGL_mkw[(counter1),1:length(list1)]  <- t(wFGL[[2]])
  weightFGL_fan[(counter1),1:length(list1)]  <- t(wFGL[[3]])
  
  
  weightCLIME_global[(counter1),1:length(list1)]  <- t(wCLIME[[1]])
  weightCLIME_mkw[(counter1),1:length(list1)]  <- t(wCLIME[[2]])
  weightCLIME_fan[(counter1),1:length(list1)]  <- t(wCLIME[[3]])
  
  weightLW_global[(counter1),1:length(list1)]  <- t(wLW[[1]])
  weightLW_mkw[(counter1),1:length(list1)]  <- t(wLW[[2]])
  weightLW_fan[(counter1),1:length(list1)]  <- t(wLW[[3]])
  
  weightPOET_global[(counter1),1:length(list1)]  <- t(wPOET[[1]])
  weightPOET_mkw[(counter1),1:length(list1)]  <- t(wPOET[[2]])
  weightPOET_fan[(counter1),1:length(list1)]  <- t(wPOET[[3]])
  
  weightNLW_global[(counter1),1:length(list1)]  <- t(wNLW[[1]])
  weightNLW_mkw[(counter1),1:length(list1)]  <- t(wNLW[[2]])
  weightNLW_fan[(counter1),1:length(list1)]  <- t(wNLW[[3]])
  
  weightPOETproj_global[(counter1),1:length(list1)]  <- t(wPOET_proj[[1]])
  weightPOETproj_mkw[(counter1),1:length(list1)]  <- t(wPOET_proj[[2]])
  weightPOETproj_fan[(counter1),1:length(list1)]  <- t(wPOET_proj[[3]])
  
  weightFF1_global[(counter1),1:length(list1)]  <- t(wFGL1[[1]])
  weightFF1_mkw[(counter1),1:length(list1)]  <- t(wFGL1[[2]])
  weightFF1_fan[(counter1),1:length(list1)]  <- t(wFGL1[[3]])
  
  weightFF3_global[(counter1),1:length(list1)]  <- t(wFGL2[[1]])
  weightFF3_mkw[(counter1),1:length(list1)]  <- t(wFGL2[[2]])
  weightFF3_fan[(counter1),1:length(list1)]  <- t(wFGL2[[3]])
  
  weightFF5_global[(counter1),1:length(list1)]  <- t(wFGL3[[1]])
  weightFF5_mkw[(counter1),1:length(list1)]  <- t(wFGL3[[2]])
  weightFF5_fan[(counter1),1:length(list1)]  <- t(wFGL3[[3]])
  
  ###NO SPARSITY IN FACTORS
  weightFF1NS_global[(counter1),1:length(list1)]  <- t(w1[[1]])
  weightFF1NS_mkw[(counter1),1:length(list1)]  <- t(w1[[2]])
  weightFF1NS_fan[(counter1),1:length(list1)]  <- t(w1[[3]])
  
  weightFF3NS_global[(counter1),1:length(list1)]  <- t(w2[[1]])
  weightFF3NS_mkw[(counter1),1:length(list1)]  <- t(w2[[2]])
  weightFF3NS_fan[(counter1),1:length(list1)]  <- t(w2[[3]])
  
  weightFF5NS_global[(counter1),1:length(list1)]  <- t(w3[[1]])
  weightFF5NS_mkw[(counter1),1:length(list1)]  <- t(w3[[2]])
  weightFF5NS_fan[(counter1),1:length(list1)]  <- t(w3[[3]])
  
}

StockInd <- matrix(rep(StockInd[1:216,],each=21),nrow=R)

weightFGL_global <- matrix(rep(weightFGL_global[1:216,],each=21),nrow=R)
weightFGL_mkw <- matrix(rep(weightFGL_mkw[1:216,],each=21),nrow=R)
weightFGL_fan <- matrix(rep(weightFGL_fan[1:216,],each=21),nrow=R)

weightLW_global <- matrix(rep(weightLW_global[1:216,],each=21),nrow=R)
weightLW_mkw <- matrix(rep(weightLW_mkw[1:216,],each=21),nrow=R)
weightLW_fan <- matrix(rep(weightLW_fan[1:216,],each=21),nrow=R)

weightPOET_global <- matrix(rep(weightPOET_global[1:216,],each=21),nrow=R)
weightPOET_mkw <- matrix(rep(weightPOET_mkw[1:216,],each=21),nrow=R)
weightPOET_fan <- matrix(rep(weightPOET_fan[1:216,],each=21),nrow=R)

weightCLIME_global <- matrix(rep(weightCLIME_global[1:216,],each=21),nrow=R)
weightCLIME_mkw <- matrix(rep(weightCLIME_mkw[1:216,],each=21),nrow=R)
weightCLIME_fan <- matrix(rep(weightCLIME_fan[1:216,],each=21),nrow=R)

weightNLW_global<- matrix(rep(weightNLW_global[1:216,],each=21),nrow=R)
weightNLW_mkw<- matrix(rep(weightNLW_mkw[1:216,],each=21),nrow=R)
weightNLW_fan<- matrix(rep(weightNLW_fan[1:216,],each=21),nrow=R)

weightPOETproj_global<- matrix(rep(weightPOETproj_global[1:216,],each=21),nrow=R)
weightPOETproj_mkw<- matrix(rep(weightPOETproj_mkw[1:216,],each=21),nrow=R)
weightPOETproj_fan<- matrix(rep(weightPOETproj_fan[1:216,],each=21),nrow=R)

weightFF1_global<- matrix(rep(weightFF1_global[1:216,],each=21),nrow=R)
weightFF1_mkw<- matrix(rep(weightFF1_mkw[1:216,],each=21),nrow=R)
weightFF1_fan<- matrix(rep(weightFF1_fan[1:216,],each=21),nrow=R)

weightFF3_global<- matrix(rep(weightFF3_global[1:216,],each=21),nrow=R)
weightFF3_mkw<- matrix(rep(weightFF3_mkw[1:216,],each=21),nrow=R)
weightFF3_fan<- matrix(rep(weightFF3_fan[1:216,],each=21),nrow=R)

weightFF5_global<- matrix(rep(weightFF5_global[1:216,],each=21),nrow=R)
weightFF5_mkw<- matrix(rep(weightFF5_mkw[1:216,],each=21),nrow=R)
weightFF5_fan<- matrix(rep(weightFF5_fan[1:216,],each=21),nrow=R)

###NO SPARSITY IN FACTORS
weightFF1NS_global<- matrix(rep(weightFF1NS_global[1:216,],each=21),nrow=R)
weightFF1NS_mkw<- matrix(rep(weightFF1NS_mkw[1:216,],each=21),nrow=R)
weightFF1NS_fan<- matrix(rep(weightFF1NS_fan[1:216,],each=21),nrow=R)

weightFF3NS_global<- matrix(rep(weightFF3NS_global[1:216,],each=21),nrow=R)
weightFF3NS_mkw<- matrix(rep(weightFF3NS_mkw[1:216,],each=21),nrow=R)
weightFF3NS_fan<- matrix(rep(weightFF3NS_fan[1:216,],each=21),nrow=R)

weightFF5NS_global<- matrix(rep(weightFF5NS_global[1:216,],each=21),nrow=R)
weightFF5NS_mkw<- matrix(rep(weightFF5NS_mkw[1:216,],each=21),nrow=R)
weightFF5NS_fan<- matrix(rep(weightFF5NS_fan[1:216,],each=21),nrow=R)

########################################
portfFGL_global  <- rolling(adj_returns2, weight = weightFGL_global)
portfFGL_mkw  <- rolling(adj_returns2, weight = weightFGL_mkw)
portfFGL_fan  <- rolling(adj_returns2, weight = weightFGL_fan)

portfLW_global  <- rolling(adj_returns2, weight = weightLW_global)
portfLW_mkw  <- rolling(adj_returns2, weight = weightLW_mkw)
portfLW_fan  <- rolling(adj_returns2, weight = weightLW_fan)

portfPOET_global  <- rolling(adj_returns2, weight = weightPOET_global)
portfPOET_mkw  <- rolling(adj_returns2, weight = weightPOET_mkw)
portfPOET_fan  <- rolling(adj_returns2, weight = weightPOET_fan)

portfCLIME_global  <- rolling(adj_returns2, weight = weightCLIME_global)
portfCLIME_mkw  <- rolling(adj_returns2, weight = weightCLIME_mkw)
portfCLIME_fan  <- rolling(adj_returns2, weight = weightCLIME_fan)

portfNLW_global<- rolling(adj_returns2, weight = weightNLW_global)
portfNLW_mkw<- rolling(adj_returns2, weight = weightNLW_mkw)
portfNLW_fan<- rolling(adj_returns2, weight = weightNLW_fan)

portfPOETproj_global<- rolling(adj_returns2, weight = weightPOETproj_global)
portfPOETproj_mkw<- rolling(adj_returns2, weight = weightPOETproj_mkw)
portfPOETproj_fan<- rolling(adj_returns2, weight = weightPOETproj_fan)

portfFF1_global<- rolling(adj_returns2, weight = weightFF1_global)
portfFF1_mkw<- rolling(adj_returns2, weight = weightFF1_mkw)
portfFF1_fan<- rolling(adj_returns2, weight = weightFF1_fan)

portfFF3_global<- rolling(adj_returns2, weight = weightFF3_global)
portfFF3_mkw<- rolling(adj_returns2, weight = weightFF3_mkw)
portfFF3_fan<- rolling(adj_returns2, weight = weightFF3_fan)

portfFF5_global<- rolling(adj_returns2, weight = weightFF5_global)
portfFF5_mkw<- rolling(adj_returns2, weight = weightFF5_mkw)
portfFF5_fan<- rolling(adj_returns2, weight = weightFF5_fan)

portfFF1NS_global  <- rolling(adj_returns2, weight = weightFF1NS_global)
portfFF1NS_mkw  <- rolling(adj_returns2, weight = weightFF1NS_mkw)
portfFF1NS_fan  <- rolling(adj_returns2, weight = weightFF1NS_fan)

portfFF3NS_global  <- rolling(adj_returns2, weight = weightFF3NS_global)
portfFF3NS_mkw  <- rolling(adj_returns2, weight = weightFF3NS_mkw)
portfFF3NS_fan  <- rolling(adj_returns2, weight = weightFF3NS_fan)

portfFF5NS_global  <- rolling(adj_returns2, weight = weightFF5NS_global)
portfFF5NS_mkw  <- rolling(adj_returns2, weight = weightFF5NS_mkw)
portfFF5NS_fan  <- rolling(adj_returns2, weight = weightFF5NS_fan)

#####EW#####
weight_EW  <- matrix(c(NA), R, ncol(adj_returns2))
for (k in 1:R) {
  set11 <- zoo(rets[k,StockInd[k,]])
  set12 <- na.locf(set11, fromLast = FALSE)
  rets[k,StockInd[k,]] <- set12
  weight_EW[k,1:length(rets[k,StockInd[k,]])] <- rep(1/length(rets[k,StockInd[k,]]),length(rets[k,StockInd[k,]]))
}
portf_EW <- rolling(adj_returns2, weight = weight_EW)
#####INDEX#####
adj_composite_EXCESS <- as.matrix(spx_excess_daily)
snp_test <- snp_excess[(pred+1):nrow(excess),]
snp_test <- as.matrix(snp_test)
mu_snp_test <- mean(snp_test)
var_snp_test <- var(snp_test)
sharpe_snp_test <- mu_snp_test/sqrt(var_snp_test)
portf_Index <- c(mu_snp_test,var_snp_test,sharpe_snp_test)
############################################
#########################################################
A1 = matrix(NA,12,4)
list1 = list(portfFGL_global,portfLW_global,portfPOET_global,portfCLIME_global,portfNLW_global,portfPOETproj_global, weightFF1_global, weightFF3_global, weightFF5_global, portfFF1NS_global, portfFF3NS_global, portfFF5NS_global)
for(i in seq_along(list1)){
  A1[i,1:3] = list1[i][[1]][[1]][1:3]
  A1[i,2] = sqrt(A1[i,2])}

B1 = matrix(NA,12,4)
list1 = list(portfFGL_global,portfLW_global,portfPOET_global,portfCLIME_global,portfNLW_global,portfPOETproj_global, weightFF1_global, weightFF3_global, weightFF5_global, portfFF1NS_global, portfFF3NS_global, portfFF5NS_global)
for(i in seq_along(list1)){
  B1[i,1:3] = list1[i][[1]][[2]][1:3]
  B1[i,4] = list1[i][[1]][[1]][4]
  B1[i,2] = sqrt(B1[i,2])}

A2 = matrix(NA,12,4)
list1 = list(portfFGL_mkw,portfLW_mkw,portfPOET_mkw,portfCLIME_mkw,portfNLW_mkw,portfPOETproj_mkw, weightFF1_mkw, weightFF3_mkw, weightFF5_mkw, portfFF1NS_mkw, portfFF3NS_mkw, portfFF5NS_mkw)
for(i in seq_along(list1)){
  A2[i,1:3] = list1[i][[1]][[1]][1:3]
  A2[i,2] = sqrt(A2[i,2])}

B2 = matrix(NA,12,4)
list1 = list(portfFGL_mkw,portfLW_mkw,portfPOET_mkw,portfCLIME_mkw,portfNLW_mkw,portfPOETproj_mkw, weightFF1_mkw, weightFF3_mkw, weightFF5_mkw, portfFF1NS_mkw, portfFF3NS_mkw, portfFF5NS_mkw)
for(i in seq_along(list1)){
  B2[i,1:3] = list1[i][[1]][[2]][1:3]
  B2[i,4] = list1[i][[1]][[1]][4]
  B2[i,2] = sqrt(B2[i,2])}

A3 = matrix(NA,11,4)
list1 = list(portfFGL_fan,portfLW_fan,portfCLIME_fan,portfNLW_fan,portfPOETproj_fan, weightFF1_fan, weightFF3_fan, weightFF5_fan, portfFF1NS_fan, portfFF3NS_fan, portfFF5NS_fan)
for(i in seq_along(list1)){
  A3[i,1:3] = list1[i][[1]][[1]][1:3]
  A3[i,2] = sqrt(A3[i,2])}

B3 = matrix(NA,11,4)
list1 = list(portfFGL_fan,portfLW_fan,portfCLIME_fan,portfNLW_fan,portfPOETproj_fan, weightFF1_fan, weightFF3_fan, weightFF5_fan, portfFF1NS_fan, portfFF3NS_fan, portfFF5NS_fan)
for(i in seq_along(list1)){
  B3[i,1:3] = list1[i][[1]][[2]][1:3]
  B3[i,4] = list1[i][[1]][[1]][4]
  B3[i,2] = sqrt(B1[i,2])}


A = rbind(A1,B1,A2,B2,A3,B3)
write.csv(A,file="daily.csv")


