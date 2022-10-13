#### User-Defined Function ####
makeK<-function(M,difford){
  if(nargs()<2){
    difford=2;
  }
  P = diag(M)
  if(difford>0){
    for(d in 1:difford){
      P = diff(P);
    }
  }
  P=t(P) %*% P
  return(P)  
}

llh_theta_int<-function(mu_vec){
  # mu_j: scalar
  intTerm<-matrix(rep(mu_vec,each=n),nrow=n)
  return(intTerm)
}

llh_theta_lin<-function(eta_mat,linearCov_mat, mix_normal_lin){
  # linearCov_mat: n by q covariate matrix
  # eta_mat: q by p vector
  linTerm<-linearCov_mat %*% (eta_mat*mix_normal_lin)
  return(linTerm)
}

llh_theta_nonlin<-function(nonLinearCov_mat, modified_gamma_mat, mix_normal_nonlin){
  # nonLinearCov_mat: n by sum_k b^*_k matrix
  # gamma_vec: sum_k b^*_k by p vetor

  nonlinTerm<-nonLinearCov_mat %*% (modified_gamma_mat*mix_normal_nonlin)
  return(nonlinTerm)
}


ALD_log_likelihood <- function(temp_errors, quantile_level)
{
  
  temp_sum_loglike = as.numeric(temp_errors>0)*(-1*quantile_level*temp_errors)+
    as.numeric(temp_errors<0)*((1-quantile_level)*temp_errors)
  
  
  return(sum(temp_sum_loglike))
}

Update_t_Ga_NEW<-function(t_old,node_mat,response_vec,
                          tau_quantile,
                          t_std, a_t, b_t,
                          intTerm, linTerm, nonlinTerm, llh_old)
{
  
  t_prime = t_old + rnorm(1,0,t_std)
  
  # llh
  theta_mat<- intTerm + linTerm + nonlinTerm
  beta_mat_prime<- theta_mat * as.numeric(abs(theta_mat)>t_prime)


  temp_errors = response_vec-apply(node_mat * beta_mat_prime,1,sum)
  llh_prime<- ALD_log_likelihood(  temp_errors, tau_quantile)

  
  log_AR= llh_prime - llh_old + log(dgamma(t_prime, shape = a_t, scale = b_t)) -
    log(dgamma(t_old, shape = a_t, scale = b_t))
  u = runif(1,0,1)
  if(log_AR>=log(u)){
    t_new=t_prime
    accept=1
    llh=llh_prime
  }else{
    t_new=t_old
    accept=0
    llh=llh_old
  }
  res<-list(t_new=t_new,accept=accept,llh_out=llh)
  return(res)
}


Update_nonLin_gamma<-function(t,
                              gamma_mat,nonLinearCov_mat,node_mat,response_vec,sigma_gamma,tau_vec,lambda_mat,
                              tau_quantile,
                              mix_normal_nonlin,
                              intTerm, linTerm, modified_gamma_mat,
                              llh_old){
  # linearCov_mat: n by q covariate matrix
  # eta_mat: q by p matrix
  # nonLinearCov_mat: n by sum_k b^*_k matrix
  # gamma_mat: sum_k b^*_k by p vetor
  # lambda_mat: same size as gamma mat (local shrinkage for gamma)
  # tau_vec: p by 1 vector (global shrinkage for gamma)
  
  # symmetric normal proposal function
  gamma_mat_prop<-matrix(rnorm(nrow(gamma_mat)*ncol(gamma_mat),gamma_mat,sd=sigma_gamma),nrow=nrow(gamma_mat))
  gamma_mat_old<-gamma_mat

  #prior
  tau_aug_tmp<-matrix(rep(tau_vec,each=q),ncol=p)
  prior_den_prime<-apply((gamma_mat_prop/lambda_mat)^2/(-2*tau_aug_tmp^2),2,sum)
  prior_den_old<-apply((gamma_mat/lambda_mat)^2/(-2*tau_aug_tmp^2),2,sum)
  
  accept_vec<-rep(0,p)
  for(j in 1:p){
    #llh: nonlinear term
    gamma_mat_prime<-gamma_mat_old
    modified_gamma_mat_prime = modified_gamma_mat
    gamma_j_old<-gamma_mat[,j]
    gamma_j_prime<-gamma_mat_prop[,j]
    
    gamma_mat_prime[,j]<-gamma_j_prime
    modified_gamma_mat_prime[,j] = rep(gamma_j_prime, b_star)
    nonlinTerm_prime<-llh_theta_nonlin(nonLinearCov_mat, modified_gamma_mat_prime, mix_normal_nonlin)
    theta_mat_prime<- intTerm + linTerm + nonlinTerm_prime
    beta_mat_prime<- theta_mat_prime * as.numeric(abs(theta_mat_prime)>t)
    

    temp_errors = response_vec-apply(node_mat * beta_mat_prime,1,sum)

    llh_prime<- ALD_log_likelihood(  temp_errors, tau_quantile)
    ######################################################################
    
    log_AR_j<- (llh_prime + prior_den_prime[j]) - (llh_old+prior_den_old[j])
    u<-runif(1,0,1)
    
    if(log_AR_j>=log(u)){
      gamma_mat_old[,j]=gamma_j_prime
      modified_gamma_mat[,j] = rep(gamma_j_prime, b_star)
      accept_vec[j]=1
      llh_old=llh_prime
    }else{
      accept_vec[j]=0
    }
  }  
  
  res<-list(gamma_mat_out=gamma_mat_old,modified_gamma_mat_out = modified_gamma_mat,
            accept=accept_vec, llh_out = llh_old)
  return(res)
}

Update_nonLin_tau2<-function(xi_vec,lambda_mat,gamma_mat,old_tau2){
  sum_b<-nrow(gamma_mat)
  tau2_vec<-1/rgamma(length(xi_vec),shape=(1+sum_b)/2,rate=((1/xi_vec)+apply((gamma_mat/lambda_mat)^2,2,sum)/2))
  return(tau2_vec)
}

Update_nonLin_xi<-function(tau){
  return(1/rgamma(length(tau),shape=1,rate=(1+1/tau^2)))
}

Update_nonLin_lambda2<-function(zeta,tau_augment,gamma_mat){
  tmp_df<-((gamma_mat/tau_augment)^2)/2
  new_rate=1/zeta + tmp_df
  return(1/matrix(rgamma(nrow(zeta)*ncol(zeta),shape=(1+rep(rep(1,q),p))/2,rate=new_rate),nrow=nrow(zeta)))
}

Update_nonLin_zeta<-function(lambda_mat){
  return(1/matrix(rgamma(nrow(lambda_mat)*ncol(lambda_mat),shape=1,rate=(1+1/lambda_mat^2)),ncol=ncol(lambda_mat)))
}

Update_lin_eta<-function(t,
                         eta_mat,linearCov_mat,node_mat,response_vec,
                         sigma_eta,alpha_vec,kappa_mat,
                         tau_quantile, mix_normal_lin, 
                         intTerm,nonlinTerm,
                         llh_old){
  # eta_mat: q by p matrix
  # alpha_vec: p by 1 vector (global shrinkage)
  # kappa_mat: q by p matrix (local shrinkage)
  
  eta_mat_prop<-matrix(rnorm(nrow(eta_mat)*ncol(eta_mat),mean=eta_mat,sd=sigma_eta),nrow=nrow(eta_mat))
  eta_mat_old<-eta_mat

  #prior
  alpha_aug<-matrix(rep(alpha_vec,each=q),ncol=p)
  prior_den_old<-apply((eta_mat/kappa_mat)^2/(-2*alpha_aug^2),2,sum)
  prior_den_prime<-apply((eta_mat_prop/kappa_mat)^2/(-2*alpha_aug^2),2,sum)
  
  
  accept_vec<-rep(NA,p)
  for(j in 1:p){
    #llh: linear term
    eta_mat_prime<-eta_mat_old
    eta_j_old<-eta_mat[,j]
    eta_j_prime<-eta_mat_prop[,j]
    
    eta_mat_prime[,j]<-eta_j_prime
    linTerm_prime<-llh_theta_lin(eta_mat_prime,linearCov_mat,mix_normal_lin)
    theta_mat_prime<- intTerm + linTerm_prime + nonlinTerm
    beta_mat_prime<- theta_mat_prime * as.numeric(abs(theta_mat_prime)>t)
    
    temp_errors = response_vec-apply(node_mat * beta_mat_prime,1,sum)

    llh_prime<- ALD_log_likelihood(  temp_errors, tau_quantile)
    ######################################################################
    
    log_AR_j<-(llh_prime + prior_den_prime[j]) - (llh_old+prior_den_old[j])
    u<-runif(1,0,1)
    
    if(log_AR_j>=log(u)){
      eta_mat_old[,j]=eta_j_prime
      accept_vec[j]=1
      llh_old=llh_prime
    }else{
      accept_vec[j]=0
    }
  }
  res<-list(eta_mat_out=eta_mat_old,accept=accept_vec,
            llh_out = llh_old)
  return(res)
}

Update_lin_alpha2<-function(rho_vec,kappa_mat,eta_mat,old_alpha2){
  alpha2_vec<-1/rgamma(length(rho_vec),shape=(1+nrow(kappa_mat))/2 , rate=(1/rho_vec + apply((eta_mat/kappa_mat)^2,2,sum)/2))

  return(alpha2_vec)
}

Update_lin_rho<-function(alpha){
  return(1/rgamma(length(alpha),shape=1,rate=(1+1/alpha^2)))
}

Update_lin_kappa2<-function(phi_mat,alpha_augment,eta_mat){
  return(1/matrix(rgamma(nrow(eta_mat)*ncol(eta_mat),shape=1,rate=(1/phi_mat+(eta_mat/alpha_augment)^2/2)),ncol=ncol(eta_mat)))
}

Update_lin_phi<-function(kappa_mat){
  return(1/matrix(rgamma(nrow(kappa_mat)*ncol(kappa_mat),shape=1,rate=(1+1/kappa_mat^2)),ncol=ncol(kappa_mat)))
}

Update_Int<-function(t,mu_vec,
                     node_mat,response_vec,sigma_int,sigma_mu,
                     tau_quantile, 
                     linTerm, nonlinTerm, llh_old){
  #symmetric normal proposal
  mu_prop<-rnorm(length(mu_vec),mean=mu_vec,sd=sigma_int)
  mu_vec_old<-mu_vec
  
  #prior
  prior_den_old<- mu_vec^2/(-2*sigma_mu^2)
  prior_den_prime<- mu_prop^2/(-2*sigma_mu^2)
  
  accept_vec<-rep(NA,p)
  for(j in 1:p){
    #llh: intercept term
    mu_vec_prime<-mu_vec_old
    mu_j_old<-mu_vec[j]
    mu_j_prime<-mu_prop[j]
    
    mu_vec_prime[j]<-mu_j_prime
    intTerm_prime<-llh_theta_int(mu_vec_prime)
    theta_mat_prime<- intTerm_prime + linTerm + nonlinTerm
    beta_mat_prime<- theta_mat_prime * as.numeric(abs(theta_mat_prime)>t)

    temp_errors = response_vec-apply(node_mat * beta_mat_prime,1,sum)
    llh_prime<- ALD_log_likelihood(  temp_errors, tau_quantile)
    ######################################################################
    
    log_AR_j<-(llh_prime + prior_den_prime[j]) - (llh_old+prior_den_old[j])
    u<-runif(1,0,1)
    
    if(log_AR_j>=log(u)){
      mu_vec_old[j]=mu_j_prime
      accept_vec[j]=1
      llh_old=llh_prime
    }else{
      accept_vec[j]=0
    }
  }
  res<-list(mu_vec_out=mu_vec_old,accept=accept_vec, llh_out = llh_old)
  return(res)
}

mix_normal_MH <- function(t,
                          eta_mat,linearCov_mat,nonLinearCov_mat,
                          node_mat,response_vec,
                          tau_quantile, mix_normal_lin, mix_normal_nonlin,
                          mix_nor_std, 
                          intTerm, modified_gamma_mat,llh_old)
{
  accept_vec = rep(0,p)
  
  mix_normal_lin_old = mix_normal_lin
  mix_normal_nonlin_old = mix_normal_nonlin
  
  for(j in 1:p)
  {
    binom_lin = 2*rbinom(q,1,1/(1+exp(-2*mix_normal_lin_old[,j]))) -1
    binom_nonlin = 2*rbinom(sum(b_star),1,1/(1+exp(-2*mix_normal_nonlin_old[,j]))) -1
    
    prop_mix_normal_lin = mix_normal_lin_old[,j] + mix_nor_std*rnorm(q,0,1)
    prop_mix_normal_nonlin = mix_normal_nonlin_old[,j] + mix_nor_std*rnorm(sum(b_star),0,1)
    
    mix_normal_lin_new = mix_normal_lin_old
    mix_normal_lin_new[,j] = prop_mix_normal_lin
    
    mix_normal_nonlin_new = mix_normal_nonlin_old
    mix_normal_nonlin_new[,j] = prop_mix_normal_nonlin
    
    theta_mat_new <-intTerm + 
      llh_theta_lin(eta_mat,linearCov_mat, mix_normal_lin_new) + 
      llh_theta_nonlin(nonLinearCov_mat,modified_gamma_mat , mix_normal_nonlin_new)
    
    beta_mat_new<- theta_mat_new * as.numeric(abs(theta_mat_new)>t)
    
    temp_errors = response_vec-apply(node_mat * beta_mat_new,1,sum)
    llh_new<- ALD_log_likelihood(  temp_errors, tau_quantile)
    
    ######################################################################
    
    llh_diff = llh_new - llh_old - 0.5*(sum((prop_mix_normal_lin - binom_lin)^2))-
      0.5*(sum((prop_mix_normal_nonlin - binom_nonlin)^2))+
      0.5*(sum((mix_normal_lin_old[,j] - binom_lin)^2))+
      0.5*(sum((mix_normal_nonlin_old[,j] - binom_nonlin)^2))
    
    if(llh_diff>log(runif(1)))
    {
      mix_normal_lin_old = mix_normal_lin_new
      mix_normal_nonlin_old = mix_normal_nonlin_new
      llh_old = llh_new
      accept_vec[j] = 1
    }
  }
  
  res = list(lin_out = mix_normal_lin_old,
             nonlin_out = mix_normal_nonlin_old,
             accept_out = accept_vec,
             llh_out = llh_old)
}
