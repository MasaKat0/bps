library(openxlsx)
library(tidyverse)

#Target data (quarterly inflation)
yI <- read.xlsx("data.xlsx",sheet="Inflation")
yI <- matrix(as.numeric(as.matrix(yI)[,2:ncol(yI)]), nrow(yI), ncol(yI)-1)
#Agent mean
a <- read.xlsx("data.xlsx",sheet="Agent Mean")
a <- matrix(as.numeric(as.matrix(a)[,2:ncol(a)]), nrow(a), ncol(a)-1)
#Agent variance
A <- read.xlsx('data.xlsx',sheet='Agent Var')
A <- matrix(as.numeric(as.matrix(A)[,2:ncol(A)]), nrow(A), ncol(A)-1)
#Agent degrees of freedom
n  <- read.xlsx('data.xlsx',sheet='Agent dof')
n  <- matrix(as.numeric(as.matrix(n)[,2:ncol(n)]), nrow(n), ncol(n)-1)


#Total time in analysis
T <- length(yI)
#Number of agents
J = dim(a)[2]
#Number of agents + intercept
p = J+1
#Number of forecast periods
K = 100


## Set priors
#　Discount factor [state observation] variance
disc_rate = c(0.95, 0.99)
#　Prior on mean of BPS coefficients
prior_mean = c(0, rep(1, J)*1/J)
#　Prior on covariance of BPS coefficients
prior_var = diag(p)*1
#　Prior on BPS degrees of freedom
prior_df = 1/(1-disc_rate[2])
#　Prior on BPS observation variance
prior_obs_var = 0.01



## Burn_in and MCMC
burn_in <- 3000;
mcmc_iter <- 5000;


## Run BPA
# Posterior BPS mean
E_BPS <-　matrix(0, mcmc_iter, K)
# Posterior BPS mean
E_BPS = matrix(0, mcmc_iter,K); 
# Posterior BPS variance
V_BPS = matrix(0, mcmc_iter,K); 
error = matrix(0, mcmc_iter,K);
mlike = matrix(0, mcmc_iter,K);
#　Posterior BPS forecast coefficient mean
ak_results = rep(0, K);
# Posterior BPS forecast coefficient variance
Rk_results = rep(0, K);
# Posterior BPS forecast observation variance
vt_results = rep(0, K);
# Posterior BPS forecast degrees of freedom
nt_results = rep(0, K);
nu = matrix(0, K, 1)
std_var <- function(x) {
  return ((x + t(x)) / 2)
}


source("BPS.R")

for (t in 50:(T-1)) {
  print(t)
  y = yI[1:t]
  mean_agent = a[1:t,];
  var_agent = A[1:t,];
  df_agent = n[1:t,];
  
  results <- BPS(y, mean_agent, var_agent, df_agent, disc_rate, prior_mean, prior_var, prior_df, prior_obs_var, burn_in, mcmc_iter)
  
  a_k <- results$a_k
  R_k <- results$R_k
  v_k <- results$v_k
  n_k <- results$n_k
  theta_post <- results$theta_post
  X_post <- results$X_post
  
  #ak_results[t] = a_k;
  #Rk_results[t] = R_k;
  #vt_results[t] = v_k;
  #nt_results[t] = n_k;
  nu[t - 49] = n_k;
  
  for (i in 1:mcmc_iter) {
    rand_gamma <- rgamma(n = 1, shape =disc_rate[2]*n[t+1]/2)
    lambda <- sqrt(0.5*disc_rate[2]*n[t+1]/rand_gamma)
    
    chol_variance <- length(a[t+1,])*chol(std_var(diag(A[t+1,])))
    x_t <- mvrnorm(1, matrix(0, nrow(chol_variance)), chol_variance)
    x_t <- t(matrix(c(1, a[t+1,] + x_t)))
    
    # compute aggregated mean and variance
    E_BPS[i,t-49] = x_t %*% t(t(a_k[i,]))
    V_BPS[i,t-49] = x_t %*% R_k[(p*i-(p-1)):(p*i),] %*% t(x_t) +v_k[i]
    
    error[i,t -49] = yI[t+1]-E_BPS[i,t-49]
    mlike[i,t-49] = exp(log(gamma(0.5*(nu[t-49]+1)))-log(gamma(0.5*nu[t-49]))-0.5*log(pi*nu[t-49]*V_BPS[i,t-49])-(0.5*(nu[t-49]+1))*log(1+1/(nu[t-49]*V_BPS[i,t-49])*(yI[t+1]-E_BPS[i,t-49])^2))
    
  }
}


# Average error of BPS
BPS_error = apply(error[,(ncol(error)-(K-1)):ncol(error)],2, mean); 
w = (1:length(BPS_error));
#　RMSFE of BPS
BPS_rmse = ((cumsum(BPS_error[1:length(BPS_error)]^2))/w)^(1/2);
# Marginal likelihood of BPS
BPS_mlike = cumsum(log(mean(mlike[(ncol(mlike)-(K-1)):ncol(mlike)])))

