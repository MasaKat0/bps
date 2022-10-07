library(MASS)
library(dplyr)

BPS <- function(y, mean_agent, var_agent, df_agent, disc_rate, prior_mean, prior_var, prior_df, prior_obs_var, burn_in, mcmc_iter){
  # Bayesian Predictive Synthesis (McAlinn & West, 2019, Journal of Econometrics 210: 155-169):
  #
  #  Synthesis Function:      
  #           y_t = x_t'\theta_t + \nu_t
  #      \theta_t = \theta_{t-1} + \omega_t
  #         \nu_t \sim Normal(0,v_t)
  #      \omega_t \sim Normal(0,v_tW_t)
  #
  #  Forecasts: 
  #       x_{j,t} \sim t(a_{j,t},A_{j,t},n_{j,t})
  #
  #  inputs:
  #    y: target time series data from t=1:T
  #    mean_agent: TxJ matrix of mean of agent forecast t=1:T
  #    var_agent: TxJ matrix of variance of agent forecast t=1:T
  #    df_agent: TxJ matrix of d.o.f. of agent forecast t=1:T
  #    disc_rate: discount rate for [state obs_var]
  #    m0 --  1xJ vector prior mean for state
  #    C0 --  JxJ prior var matrix
  #    n0 --  prior d.o.f.
  #    s0 --  prior of obs var
  #    burn_in, mcmc_iter: number of burn-in/MCMC iterations
  #  outputs:
  #    a_k: forecast mean of \theta_{T+1}
  #    R_k: forecast variance of \theta_{T+1}
  #    v_k: forecast obs var
  #    n_k: forecast d.o.f.
  #    theta_post: posterior \theta_{1:T}
  #    X_post: posterior sample of agent forecast
  #
  #  ? 2017, Kenichiro McAlinn, All rights reserved.
  
  if ((length(mean_agent) == 1) && (mean_agent > 0)) {
    data <-  y
    y <- y[(mean_agent + 1):length(y)]
    T <- length(y)
    
    x_t <- matrix(0, length(y), mean_agent)
    for (i in 1:x) {
      i <- 1
      x_t[,i] <- data[(mean_agent - i + 1):(length(data) - i)]
    }
    
    X_t <- matrix(0, T, p_x*(mcmc_iter+1))
    p <- x + 1
  } else if ((length(mean_agent) == 1) && (mean_agent == 0)) {
    T <- length(y)
    X_t <- matrix(0, T, p_x*(mcmc_iter+1))
    p <- 1
  } else {
    T <- length(y)
    p_x <- ncol(mean_agent)
    p <- p_x + 1
    X_t <- matrix(0, T, p_x*(mcmc_iter+1))
  }
  
  
  mcmc_iter = burn_in + mcmc_iter
  
  m_t <- matrix(0, T + 1, p)
  C_t <- matrix(0, (T + 1)*p, p)
  n_t <- matrix(0, T + 1, 1)
  s_t <- matrix(0, T + 1, mcmc_iter)
  a_t <- matrix(0, T, p)
  R_t <- matrix(0, T*p, p)
  f_t <- matrix(0, T, 1)
  q_t <- matrix(0, T, 1)
  
  er_t <- matrix(0, T, 1)
  m_like <- matrix(0, T, 1)
  v_t <- matrix(0, T, mcmc_iter)
  
  phi_t <- matrix(0, T, p_x)
  X_t <- matrix(0, T, p_x*(mcmc_iter + 1))
  theta_t <- matrix(0, T, p*mcmc_iter)
  
  a_k <- matrix(0, mcmc_iter, p)
  R_k <- matrix(0, p*mcmc_iter, p)
  v_k <- matrix(0, mcmc_iter, 1)
  n_k <- matrix(0, 1, 1)
  
  m_t[1,] <- prior_mean
  C_t[1:p,] <- prior_var
  n_t[1,] <- prior_df
  s_t[1,] <- prior_obs_var
  v_t[1,] <- s_t[1,]
  
  d <- disc_rate[1]
  beta <- disc_rate[2]
  
  std_var <- function(x) {
    return ((x + t(x)) / 2)
  }
  
  
  for (t in 1:T) {
    for (n in 1:p_x) {
      var_gamma <- rgamma(1, beta*df_agent[t,n]/2)
      phi_t[t,n] = 0.5*beta*df_agent[t,n]/var_gamma
    }
    
    chol_variance <- chol(std_var(diag(phi_t[t,]*var_agent[t,])))
    var_normal <- mvrnorm(1, matrix(0, p_x), chol_variance)
    X_t[t,1:p_x] = mean_agent[t,]+var_normal
  }
  
  
  for (i in 1:mcmc_iter) {
    print("mcmc iter:")
    print(i)
    # Forward Filter
    for (t in 1:T) {
      #F_t <- X_t[t, (p*i - (p_x - 1)):(p*i)] %>% t()
      F_t = c(1, X_t[t,(p_x*i-(p_x-1)):(p_x*i)]) %>% t()
      
      # prior for time t (1期先予測分布)mu_t|D_{t-1}
      a_t[t,] <- m_t[t,] #状態変数の平均
      R_t[(p*t - (p - 1)):(p*t),] <- C_t[(p*t - (p - 1)):(p*t),]/d #状態変数の分散共分散行列
      
      # predict time t (1期先予測尤度)Y_t|D_{t-1}
      f_t[t,1] <- F_t %*% a_t[t,]
      q_t[t,1] <- F_t %*% (C_t[(p*t - (p - 1)):(p*t),] %*% t(F_t)/d) + s_t[t,i]
      # compute forecast error and adaptive vector
      e_t <- y[t] - f_t[t,1] #1期先予測誤差(イノベーション)
      er_t[t,1] <- y[t] - f_t[t,1]
      nu <- beta*n_t[t,]
      m_like[t,1] <- log(gamma(0.5*(nu + 1))) - log(gamma(0.5*nu)) - 0.5*log(pi*nu*q_t[t,1]) - (0.5*(nu + 1))*log(1 + 1/(nu*q_t[t, 1])*(y[t] - f_t[t,1])^2)
      # m_like[t,1] <- dnorm(y[t],f_t[t,1],sqrt(q_t[t,1]))
      A_t <- R_t[(p*t - (p - 1)):(p*t),] %*% t(F_t)/q_t[t,1]
      
      # posterior for time t
      n_t[t+1,] <- beta*n_t[t,] + 1
      r_t <- (beta*n_t[t,] + e_t^2/q_t[t,1])/n_t[t+1,]
      s_t[t+1,i] <- r_t*s_t[t,i]
      v_t[t,1] <- n_t[t+1,]*s_t[t+1,1]/rchisq(1, n_t[t+1,])
      # v_t(t+1,1) = (beta*n_t(t+1,1))*s_t(t+1,1)/chi2rnd(beta*n_t(t+1,1));
      m_t[t+1,] <- a_t[t,] + t(A_t*e_t)
      C_t[(p*(t + 1) - (p - 1)):(p*(t + 1)),] <- std_var(r_t*(R_t[(p*t - (p - 1)):(p*t),] - q_t[t,1]*(A_t %*% t(A_t))))
    }
    # sample theta at T
    rand_gamma <- rgamma(n = 1, shape = n_t[nrow(n_t)]/2, rate = 2/(n_t[nrow(n_t)]*s_t[nrow(s_t),i]))
    v_t[nrow(v_t),i] = 1/rand_gamma
    
    chol_variance <- length(m_t[nrow(m_t),])*chol(std_var(C_t[(nrow(C_t)-(p-1)):nrow(C_t),]*v_t[nrow(v_t),i]/s_t[nrow(s_t),i]))
    theta_t[nrow(theta_t), (p*i-(p-1)):(p*i)] = mvrnorm(1, m_t[nrow(m_t),], chol_variance)
    
    # theta at T+1
    n_k = beta*n_t[nrow(n_t)]+1;
    rand_gamma <- rgamma(n = 1, shape = beta*n_t[nrow(n_t)]/2, rate = 2/(beta*n_t[nrow(n_t)]*s_t[nrow(s_t),i]))
    v_k[i,1] <- 1/rand_gamma
    a_k[i,] = m_t[nrow(m_t),]
    R_k[(p*i-(p-1)):(p*i),] = (C_t[(nrow(C_t)-(p-1)):nrow(C_t),]/d)*(v_k[i,1]/s_t[nrow(s_t),i])
    # backward-sampler
    for (t in rev(1:(T-1))) {
      rand_gamma <- rgamma(n = 1, shape = (1-beta)*n_t[t+1]/2, rate = 2/(n_t[t+1]*s_t[t+1,i]))
      v_t[t,i] = 1/(1/v_t[t+1,i]*beta+rand_gamma)
      m_star_t = t(m_t[t+1,]) + d*(t(theta_t[t+1,(p*i-(p-1)):(p*i)])-t(a_t[t+1,]))
      C_star_t = C_t[(p*(t+1)-(p-1)):(p*(t+1)),]*(1-d)*(v_t[t,i]/s_t[t+1,i])
      chol_variance <- length(m_star_t)*chol(std_var(C_star_t))
      theta_t[t,(p*i-(p-1)):(p*i)] <- mvrnorm(1, m_star_t[1,], length(m_star_t)*chol(std_var(C_star_t)))
    }
        
    for (t in 1:T) {
      A_st = diag(phi_t[t,] * var_agent[t,]);
      a_st = t(matrix(mean_agent[t,]))
      theta_p = t(matrix(theta_t[t,(p*i-(p-2)):(p*i)]))
      theta_1 = theta_t[t,p*i-(p-1)]
      sigma = (theta_p %*% A_st)/((v_t[t,i] + theta_p %*% A_st %*% t(theta_p) )[1.1])
      a_star = a_st+sigma*(y[t]-(theta_1+(theta_p %*% t(a_st))[1,1]));
      A_star = std_var(A_st-A_st %*% t(theta_p) %*% sigma);
      #A_star = std_var(A_st);
      X_t[t,(p_x*(i+1)-(p_x-1)):(p_x*(i+1))] = mvrnorm(1, a_star, length(a_star)*chol(std_var(A_star))) 
      scale_values = (df_agent[t,]+(X_t[t,(p_x*(i+1)-(p_x-1)):(p_x*(i+1))]-a_st)^2/var_agent[t,]) / 2
      for (n in 1:length(scale_values)) {
        phi_t[t,n] = 0.5*(df_agent[t,n]+1)/rgamma(n = 1, shape = scale_values[n])
      }
    }
    
  }
  
  theta_post = theta_t[,(burn_in*p_x+1):ncol(theta_t)]
  X_post = X_t[,(burn_in*p_x+1):ncol(X_t)]
  a_k = a_k[(burn_in+1):nrow(a_k),];
  R_k = R_k[(p*burn_in+1):nrow(R_k),];
  v_k = v_k[(burn_in+1):nrow(v_k),]
  
  list("theta_post" = theta_post, "X_post" = X_post, "a_k" = a_k, "R_k" = R_k, "v_k" = v_k, "n_k" = n_k) %>% 
    return()
}
