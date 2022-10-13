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
#########################################################################
llh_theta_int<-function(mu_vec){
  # mu_j: scalar
  intTerm<-matrix(rep(mu_vec,each=n),nrow=n)
  return(intTerm)
}
#########################################################################
llh_theta_lin<-function(eta_mat,linearCov_mat, mix_normal_lin){
  # linearCov_mat: n by q covariate matrix
  # eta_mat: q by p vector
  # mix_normal_lin: q by p vector
  
  linTerm<-linearCov_mat %*% (eta_mat*mix_normal_lin)
  return(linTerm)
}
#########################################################################
llh_theta_nonlin<-function(nonLinearCov_mat, modified_gamma_mat, mix_normal_nonlin){
  # nonLinearCov_mat: n by sum_k b^*_k matrix
  # modified_gamma_mat: sum_k b^*_k by p matrix
  # min_normal_nonlin b^*_k by p matrix
  
  nonlinTerm<-nonLinearCov_mat %*% (modified_gamma_mat*mix_normal_nonlin)
  return(nonlinTerm)
}
#########################################################################
ALD_log_likelihood <- function(temp_errors, quantile_level, indicator_th)
{
  
  temp_sum_loglike = as.numeric(temp_errors>indicator_th)*(-1*quantile_level*temp_errors)+
    as.numeric(temp_errors>indicator_th)*(quantile_level*indicator_th)+
    as.numeric(temp_errors<indicator_th)*((1-quantile_level)*temp_errors)+
    as.numeric(temp_errors<indicator_th)*(-1*(1-quantile_level)*indicator_th)
  
  return(sum(temp_sum_loglike))
}
#########################################################################
Update_Int<-function(t,mu_vec,
                     node_mat,response_vec,sigma_int,sigma_mu,
                     tau_quantile, 
                     linTerm, nonlinTerm, llh_old, indicator_th){
  #symmetric normal proposal
  mu_prop<-rnorm(1,mean=mu_vec,sd=sigma_int)
  mu_vec_old<-mu_vec
  
  #prior
  prior_den_old<- mu_vec^2/(-2*sigma_mu^2)
  prior_den_prime<- mu_prop^2/(-2*sigma_mu^2)
  
  accept_vec<-NA
  
  #llh: intercept term
  mu_vec_prime<-mu_vec_old
  mu_j_old<-mu_vec
  mu_j_prime<-mu_prop
  
  mu_vec_prime<-mu_j_prime
  intTerm_prime<-llh_theta_int(mu_vec_prime)
  theta_mat_prime<- intTerm_prime + linTerm + nonlinTerm
  beta_mat_prime<- theta_mat_prime * as.numeric(abs(theta_mat_prime)>t)
  
  temp_errors = response_vec- node_mat * beta_mat_prime
  
  llh_prime<- ALD_log_likelihood(  temp_errors, tau_quantile, indicator_th)
  ######################################################################
  
  log_AR <-(llh_prime + prior_den_prime) - (llh_old+prior_den_old)
  u<-runif(1,0,1)
  
  if(log_AR > log(u)){
    mu_vec_old=mu_j_prime
    accept_vec=1
    llh_old=llh_prime
  }else{
    accept_vec=0
  }
  
  res<-list(mu_vec_out=mu_vec_old,accept=accept_vec, llh_out = llh_old)
  return(res)
}
#########################################################################
Update_nonLin_gamma<-function(t,
                              gamma_mat,nonLinearCov_mat,node_mat,response_vec,sigma_gamma,tau_vec,lambda_mat,
                              tau_quantile,
                              mix_normal_nonlin,
                              intTerm, linTerm, modified_gamma_mat,
                              llh_old, indicator_th){
  # linearCov_mat: n by q covariate matrix
  # eta_mat: q by p matrix
  # nonLinearCov_mat: n by sum_k b^*_k matrix
  # gamma_mat: sum_k b^*_k by p vetor
  # lambda_mat: same size as gamma mat (local shrinkage for gamma)
  # tau_vec: p by 1 vector (global shrinkage for gamma)
  
  # symmetric normal proposal function
  gamma_mat_prop<-rnorm(length(gamma_mat),gamma_mat,sd=sigma_gamma)
  gamma_mat_old<-gamma_mat
  
  ######################################################################
  
  #prior
  tau_aug_tmp<-rep(tau_vec,each=q)
  prior_den_prime<-sum((gamma_mat_prop/lambda_mat)^2/(-2*tau_aug_tmp^2))
  prior_den_old<-sum((gamma_mat/lambda_mat)^2/(-2*tau_aug_tmp^2))
  
  accept_vec<-NA
  
  #llh: nonlinear term
  gamma_mat_prime<-gamma_mat_old
  modified_gamma_mat_prime = modified_gamma_mat
  gamma_j_old<-gamma_mat
  gamma_j_prime<-gamma_mat_prop
  
  gamma_mat_prime<-gamma_j_prime
  modified_gamma_mat_prime = rep(gamma_j_prime, b_star)
  nonlinTerm_prime<-llh_theta_nonlin(nonLinearCov_mat, modified_gamma_mat_prime, mix_normal_nonlin)
  theta_mat_prime<- intTerm + linTerm + nonlinTerm_prime
  beta_mat_prime<- theta_mat_prime * as.numeric(abs(theta_mat_prime)>t)
  
  
  temp_errors = response_vec- node_mat * beta_mat_prime
  
  llh_prime<- ALD_log_likelihood(  temp_errors, tau_quantile, indicator_th = NZV_FII_LLH)
  ######################################################################
  
  log_AR <- (llh_prime + prior_den_prime) - (llh_old+prior_den_old)
  
  u<-runif(1,0,1)
  
  if(log_AR>=log(u)){
    gamma_mat_old=gamma_j_prime
    modified_gamma_mat = rep(gamma_j_prime, b_star)
    accept_vec=1
    llh_old=llh_prime
  }else{
    accept_vec=0
  }
  
  
  res<-list(gamma_mat_out=gamma_mat_old,modified_gamma_mat_out = modified_gamma_mat,
            accept=accept_vec, llh_out = llh_old)
  return(res)
}
##################################################################################
Update_nonLin_tau2<-function(xi_vec,lambda_mat,gamma_mat,old_tau2){
  # xi_vec: p by 1 vector
  # lambda_mat: same size as gamma mat (local shrinkage for gamma)
  # tau_mat: p by 1 matrix (global shrinkage for gamma)
  sum_b<-length(gamma_mat)
  tau2_vec<-1/rgamma(length(xi_vec),shape=(1+sum_b)/2,rate=((1/xi_vec)+sum((gamma_mat/lambda_mat)^2)/2))
  
  return(tau2_vec)
}
##################################################################################
Update_nonLin_xi<-function(tau){
  return(1/rgamma(length(tau),shape=1,rate=(1+1/tau^2)))
}
##################################################################################
Update_nonLin_lambda2<-function(zeta,tau_augment,gamma_mat){
  
  tmp_df<-((gamma_mat/tau_augment)^2)/2
  new_rate=1/zeta + tmp_df
  return(1/rgamma(length(zeta),shape=1,rate=new_rate))
}
##################################################################################
Update_nonLin_zeta<-function(lambda_mat){
  return(1/rgamma(length(lambda_mat),shape=1,rate=(1+1/lambda_mat^2)))
}
##################################################################################
Update_lin_eta<-function(t,
                         eta_mat,linearCov_mat,node_mat,response_vec,
                         sigma_eta,alpha_vec,kappa_mat,
                         tau_quantile, mix_normal_lin, 
                         intTerm,nonlinTerm,
                         llh_old, indicator_th){
  # eta_mat: q by p matrix
  # alpha_vec: p by 1 vector (global shrinkage)
  # kappa_mat: q by p matrix (local shrinkage)
  
  eta_mat_prop<-rnorm(length(eta_mat),mean=eta_mat,sd=sigma_eta)
  eta_mat_old<-eta_mat
  
  
  #prior
  alpha_aug<-rep(alpha_vec,each=q)
  prior_den_old<-sum((eta_mat/kappa_mat)^2/(-2*alpha_aug^2))
  prior_den_prime<-sum((eta_mat_prop/kappa_mat)^2/(-2*alpha_aug^2))
  
  
  accept_vec<-NA
  
  #llh: linear term
  eta_mat_prime<-eta_mat_old
  eta_j_old<-eta_mat
  eta_j_prime<-eta_mat_prop
  
  eta_mat_prime<-eta_j_prime
  linTerm_prime<-llh_theta_lin(eta_mat_prime,linearCov_mat,mix_normal_lin)
  theta_mat_prime<- intTerm + linTerm_prime + nonlinTerm
  beta_mat_prime<- theta_mat_prime * as.numeric(abs(theta_mat_prime)>t)
  
  
  temp_errors = response_vec- node_mat * beta_mat_prime
  
  llh_prime<- ALD_log_likelihood(  temp_errors, tau_quantile, indicator_th)
  ######################################################################
  
  log_AR_j<-(llh_prime + prior_den_prime) - (llh_old+prior_den_old)
  u<-runif(1,0,1)
  
  if(log_AR_j>=log(u)){
    eta_mat_old=eta_j_prime
    accept_vec=1
    llh_old=llh_prime
  }else{
    accept_vec=0
  }
  
  res<-list(eta_mat_out=eta_mat_old,accept=accept_vec,
            llh_out = llh_old)
  return(res)
}
##################################################################################
Update_lin_alpha2<-function(rho_vec,kappa_mat,eta_mat,old_alpha2){
  
  alpha2_vec<-1/rgamma(length(rho_vec),shape=(1+length(kappa_mat))/2 , rate=(1/rho_vec + sum((eta_mat/kappa_mat)^2)/2))
  return(alpha2_vec)
}
##################################################################################
Update_lin_rho<-function(alpha){
  return(1/rgamma(length(alpha),shape=1,rate=(1+1/alpha^2)))
}
##################################################################################
Update_lin_kappa2<-function(phi_mat,alpha_augment,eta_mat){
  return(1/rgamma(length(eta_mat),shape=1,rate=(1/phi_mat+(eta_mat/alpha_augment)^2/2)))
}
##################################################################################
Update_lin_phi<-function(kappa_mat){
  return(1/rgamma(length(kappa_mat),shape=1,rate=(1+1/kappa_mat^2)))
}
##################################################################################
mix_normal_MH <- function(t,
                          eta_mat,linearCov_mat,nonLinearCov_mat,
                          node_mat,response_vec,
                          tau_quantile, mix_normal_lin, mix_normal_nonlin,
                          mix_nor_std, 
                          intTerm, modified_gamma_mat,llh_old, indicator_th)
{
  accept_vec = NA
  
  mix_normal_lin_old = mix_normal_lin
  mix_normal_nonlin_old = mix_normal_nonlin
  
  
  binom_lin = 2*rbinom(q,1,1/(1+exp(-2*mix_normal_lin_old))) -1
  binom_nonlin = 2*rbinom(sum(b_star),1,1/(1+exp(-2*mix_normal_nonlin_old))) -1
  
  prop_mix_normal_lin = mix_normal_lin_old + mix_nor_std*rnorm(q,0,1)
  prop_mix_normal_nonlin = mix_normal_nonlin_old + mix_nor_std*rnorm(sum(b_star),0,1)
  
  mix_normal_lin_new = mix_normal_lin_old
  mix_normal_lin_new = prop_mix_normal_lin
  
  mix_normal_nonlin_new = mix_normal_nonlin_old
  mix_normal_nonlin_new = prop_mix_normal_nonlin
  
  theta_mat_new <-intTerm + 
    llh_theta_lin(eta_mat,linearCov_mat, mix_normal_lin_new) + 
    llh_theta_nonlin(nonLinearCov_mat,modified_gamma_mat , mix_normal_nonlin_new)
  
  beta_mat_new<- theta_mat_new * as.numeric(abs(theta_mat_new)>t)
  
  temp_errors = response_vec - node_mat*beta_mat_new
  
  llh_new<- ALD_log_likelihood(  temp_errors, tau_quantile,indicator_th )
  
  ######################################################################
  
  llh_diff = llh_new - llh_old - 0.5*(sum((prop_mix_normal_lin - binom_lin)^2))-
    0.5*(sum((prop_mix_normal_nonlin - binom_nonlin)^2))+
    0.5*(sum((mix_normal_lin_old - binom_lin)^2))+
    0.5*(sum((mix_normal_nonlin_old - binom_nonlin)^2))
  
  if(llh_diff>log(runif(1)))
  {
    mix_normal_lin_old = mix_normal_lin_new
    mix_normal_nonlin_old = mix_normal_nonlin_new
    llh_old = llh_new
    accept_vec = 1
  }else{
    accept_vec = 0
  }
  
  
  res = list(lin_out = mix_normal_lin_old,
             nonlin_out = mix_normal_nonlin_old,
             accept_out = accept_vec,
             llh_out = llh_old)
}
##################################################################################
Update_t_Ga_NEW<-function(t_old,
                          node_mat,response_vec,
                          tau_quantile,
                          t_std, a_t, b_t,
                          intTerm, linTerm, nonlinTerm, llh_old, indicator_th)
{
  # independent proposal function
  #t_prime<-rgamma(1,shape=2.5,scale = 1)
  
  t_prime = t_old + rnorm(1,0,t_std)
  
  # llh
  
  theta_mat<- intTerm + linTerm + nonlinTerm
  beta_mat_prime<- theta_mat * as.numeric(abs(theta_mat)>t_prime)
  
  temp_errors = response_vec-node_mat*beta_mat_prime
  
  llh_prime<- ALD_log_likelihood(temp_errors, tau_quantile, indicator_th)
  
  ######################################################################
  
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
##################################################################################
check_for_DAGs <- function(DAG_all_patients_this_iter)
{
  # DAG_index = rep(FALSE, sample_size)
  # for(per_id in 1:sample_size)
  # {
  #    # DAG_this_person = DAG_all_patients_this_iter[per_id,,]
  #    # adj_mat_this_person = graph_from_adjacency_matrix(DAG_this_person)
  #    # DAG_index[per_id] = is_dag(adj_mat_this_person)
  #    DAG_index[per_id] = is.DAG(DAG_all_patients_this_iter[per_id,,])
  # }
  # 
  # sum_DAG_index = sum(DAG_index)
  # if(sum_DAG_index == sample_size)
  # {
  #   return(TRUE)
  # }else{
  #   return(FALSE)
  # }
  
  #union_DAG = apply(DAG_all_patients_this_iter, c(2,3),sum)
  #union_DAG[union_DAG!=0] =1
  ##
  #union_DAG = t(union_DAG)
  ##
  return(is.DAG(DAG_all_patients_this_iter))
  
}
##################################################################################
################# Update functions with DAG check ################################
##################################################################################
Update_t_Ga_NEW_WDC<-function(t_old,
                              node_mat,response_vec,
                              tau_quantile,
                              t_std, a_t, b_t,
                              intTerm, linTerm, nonlinTerm, llh_old, indicator_th, 
                              DAG_all_patients_this_iter,
                              num_beta_nnz)
{
  # independent proposal function
  #t_prime<-rgamma(1,shape=2.5,scale = 1)
  
  t_prime = t_old + rnorm(1,0,t_std)
  
  if(t_prime>t_old)
  {
    theta_mat<- intTerm + linTerm + nonlinTerm
    beta_mat_prime<- theta_mat * as.numeric(abs(theta_mat)>t_prime)
    
    temp_errors = response_vec-node_mat*beta_mat_prime
    
    llh_prime<- ALD_log_likelihood(temp_errors, tau_quantile, indicator_th)
    
    ######################################################################
    
    log_AR= llh_prime - llh_old + log(dgamma(t_prime, shape = a_t, scale = b_t)) -
      log(dgamma(t_old, shape = a_t, scale = b_t))
    u = runif(1,0,1)
    if(log_AR>=log(u)){
      t_new=t_prime
      accept=1
      llh=llh_prime
      #DAG_all_patients_this_iter[,l,p_reg] = as.numeric(beta_mat_prime!=0)
      sum_num_changes = sum(as.numeric(beta_mat_prime!=0))
    }else{
      t_new=t_old
      accept=0
      llh=llh_old
      #sum_num_changes = 0
      sum_num_changes = num_beta_nnz
    }
    res<-list(t_new=t_new,accept=accept,llh_out=llh, dag_check = TRUE, 
              #dag_all = DAG_all_patients_this_iter
              sum_num_changes = sum_num_changes)
    return(res)
    
  }else{
    theta_mat<- intTerm + linTerm + nonlinTerm
    beta_mat_prime<- theta_mat * as.numeric(abs(theta_mat)>t_prime)
    
    temp_DAG_all_patients_this_iter = DAG_all_patients_this_iter
    sum_num_changes = sum(as.numeric(beta_mat_prime!=0))
    temp_DAG_all_patients_this_iter[p_reg,l] = as.numeric(sum_num_changes!=0)
    
    if(check_for_DAGs(temp_DAG_all_patients_this_iter))
    {
      temp_errors = response_vec-node_mat*beta_mat_prime
      
      llh_prime<- ALD_log_likelihood(temp_errors, tau_quantile, indicator_th)
      
      ######################################################################
      
      log_AR= llh_prime - llh_old + log(dgamma(t_prime, shape = a_t, scale = b_t)) -
        log(dgamma(t_old, shape = a_t, scale = b_t))
      u = runif(1,0,1)
      if(log_AR>=log(u)){
        t_new=t_prime
        accept=1
        llh=llh_prime
        #DAG_all_patients_this_iter[,l,p_reg] = as.numeric(beta_mat_prime!=0)
      }else{
        t_new=t_old
        accept=0
        llh=llh_old
        #sum_num_changes = 0
        sum_num_changes = num_beta_nnz
      }
      res<-list(t_new=t_new,accept=accept,llh_out=llh, dag_check = TRUE,
                #dag_all = DAG_all_patients_this_iter
                sum_num_changes = sum_num_changes)
      return(res)
    }else{
      t_new=t_old
      accept=0
      llh=llh_old
      res<-list(t_new=t_new,accept=accept,llh_out=llh, dag_check = FALSE,
                #dag_all = DAG_all_patients_this_iter
                #sum_num_changes = 0
                sum_num_changes = num_beta_nnz)
      return(res)
    }
  }
}
##################################################################################
Update_Int_WDC<-function(t,mu_vec,
                         node_mat,response_vec,sigma_int,sigma_mu,
                         tau_quantile, 
                         linTerm, nonlinTerm, llh_old, indicator_th, 
                         DAG_all_patients_this_iter,
                         num_beta_nnz){
  #symmetric normal proposal
  mu_prop<-rnorm(1,mean=mu_vec,sd=sigma_int)
  mu_vec_old<-mu_vec
  
  #prior
  prior_den_old<- mu_vec^2/(-2*sigma_mu^2)
  prior_den_prime<- mu_prop^2/(-2*sigma_mu^2)
  
  accept_vec<-NA
  
  #llh: intercept term
  mu_vec_prime<-mu_vec_old
  mu_j_old<-mu_vec
  mu_j_prime<-mu_prop
  
  mu_vec_prime<-mu_j_prime
  intTerm_prime<-llh_theta_int(mu_vec_prime)
  theta_mat_prime<- intTerm_prime + linTerm + nonlinTerm
  beta_mat_prime<- theta_mat_prime * as.numeric(abs(theta_mat_prime)>t)
  
  temp_DAG_all_patients_this_iter = DAG_all_patients_this_iter
  sum_num_changes = sum(as.numeric(beta_mat_prime!=0))
  temp_DAG_all_patients_this_iter[p_reg,l] = as.numeric(sum_num_changes!=0)
  
  if(check_for_DAGs(temp_DAG_all_patients_this_iter))
  {
    temp_errors = response_vec- node_mat * beta_mat_prime
    
    llh_prime<- ALD_log_likelihood(  temp_errors, tau_quantile, indicator_th)
    ######################################################################
    
    log_AR <-(llh_prime + prior_den_prime) - (llh_old+prior_den_old)
    u<-runif(1,0,1)
    
    if(log_AR > log(u)){
      mu_vec_old=mu_j_prime
      accept_vec=1
      llh_old=llh_prime
      #DAG_all_patients_this_iter[,l,p_reg] = as.numeric(beta_mat_prime!=0)
    }else{
      accept_vec=0
      #sum_num_changes = 0
      sum_num_changes = num_beta_nnz
    }
  }else{
    accept_vec=0
    #sum_num_changes = 0
    sum_num_changes = num_beta_nnz
  }
  
  res<-list(mu_vec_out=mu_vec_old,accept=accept_vec, llh_out = llh_old,
            #dag_all = DAG_all_patients_this_iter
            sum_num_changes = sum_num_changes)
  return(res)
  
  
}
#########################################################################
Update_nonLin_gamma_WDC<-function(t,
                                  gamma_mat,nonLinearCov_mat,node_mat,response_vec,sigma_gamma,tau_vec,lambda_mat,
                                  tau_quantile,
                                  mix_normal_nonlin,
                                  intTerm, linTerm, modified_gamma_mat,
                                  llh_old, indicator_th,DAG_all_patients_this_iter,
                                  num_beta_nnz){
  # linearCov_mat: n by q covariate matrix
  # eta_mat: q by p matrix
  # nonLinearCov_mat: n by sum_k b^*_k matrix
  # gamma_mat: sum_k b^*_k by p vetor
  # lambda_mat: same size as gamma mat (local shrinkage for gamma)
  # tau_vec: p by 1 vector (global shrinkage for gamma)
  
  # symmetric normal proposal function
  gamma_mat_prop<-rnorm(length(gamma_mat),gamma_mat,sd=sigma_gamma)
  gamma_mat_old<-gamma_mat
  
  ######################################################################
  
  #prior
  tau_aug_tmp<-rep(tau_vec,each=q)
  prior_den_prime<-sum((gamma_mat_prop/lambda_mat)^2/(-2*tau_aug_tmp^2))
  prior_den_old<-sum((gamma_mat/lambda_mat)^2/(-2*tau_aug_tmp^2))
  
  accept_vec<-NA
  
  #llh: nonlinear term
  gamma_mat_prime<-gamma_mat_old
  modified_gamma_mat_prime = modified_gamma_mat
  gamma_j_old<-gamma_mat
  gamma_j_prime<-gamma_mat_prop
  
  gamma_mat_prime<-gamma_j_prime
  modified_gamma_mat_prime = rep(gamma_j_prime, b_star)
  nonlinTerm_prime<-llh_theta_nonlin(nonLinearCov_mat, modified_gamma_mat_prime, mix_normal_nonlin)
  theta_mat_prime<- intTerm + linTerm + nonlinTerm_prime
  beta_mat_prime<- theta_mat_prime * as.numeric(abs(theta_mat_prime)>t)
  
  temp_DAG_all_patients_this_iter = DAG_all_patients_this_iter
  sum_num_changes = sum(as.numeric(beta_mat_prime!=0))
  temp_DAG_all_patients_this_iter[p_reg,l] = as.numeric(sum_num_changes!=0)
  
  if(check_for_DAGs(temp_DAG_all_patients_this_iter))
  {
    temp_errors = response_vec- node_mat * beta_mat_prime
    
    llh_prime<- ALD_log_likelihood(  temp_errors, tau_quantile, indicator_th = NZV_FII_LLH)
    ######################################################################
    
    log_AR <- (llh_prime + prior_den_prime) - (llh_old+prior_den_old)
    
    u<-runif(1,0,1)
    
    if(log_AR>=log(u)){
      gamma_mat_old=gamma_j_prime
      modified_gamma_mat = rep(gamma_j_prime, b_star)
      accept_vec=1
      llh_old=llh_prime
      #DAG_all_patients_this_iter[,l,p_reg] = as.numeric(beta_mat_prime!=0)
    }else{
      accept_vec=0
      #sum_num_changes = 0
      sum_num_changes = num_beta_nnz
    }
  }else{
    accept_vec=0
    #sum_num_changes = 0
    sum_num_changes = num_beta_nnz
  }
  
  res<-list(gamma_mat_out=gamma_mat_old,modified_gamma_mat_out = modified_gamma_mat,
            accept=accept_vec, llh_out = llh_old, 
            #dag_all = DAG_all_patients_this_iter
            sum_num_changes = sum_num_changes)
  return(res)
}
##################################################################################
Update_lin_eta_WDC<-function(t,
                             eta_mat,linearCov_mat,node_mat,response_vec,
                             sigma_eta,alpha_vec,kappa_mat,
                             tau_quantile, mix_normal_lin, 
                             intTerm,nonlinTerm,
                             llh_old, indicator_th, 
                             DAG_all_patients_this_iter,
                             num_beta_nnz){
  # eta_mat: q by p matrix
  # alpha_vec: p by 1 vector (global shrinkage)
  # kappa_mat: q by p matrix (local shrinkage)
  
  eta_mat_prop<-rnorm(length(eta_mat),mean=eta_mat,sd=sigma_eta)
  eta_mat_old<-eta_mat
  
  
  #prior
  alpha_aug<-rep(alpha_vec,each=q)
  prior_den_old<-sum((eta_mat/kappa_mat)^2/(-2*alpha_aug^2))
  prior_den_prime<-sum((eta_mat_prop/kappa_mat)^2/(-2*alpha_aug^2))
  
  
  accept_vec<-NA
  
  #llh: linear term
  eta_mat_prime<-eta_mat_old
  eta_j_old<-eta_mat
  eta_j_prime<-eta_mat_prop
  
  eta_mat_prime<-eta_j_prime
  linTerm_prime<-llh_theta_lin(eta_mat_prime,linearCov_mat,mix_normal_lin)
  theta_mat_prime<- intTerm + linTerm_prime + nonlinTerm
  beta_mat_prime<- theta_mat_prime * as.numeric(abs(theta_mat_prime)>t)
  
  temp_DAG_all_patients_this_iter = DAG_all_patients_this_iter
  sum_num_changes = sum(as.numeric(beta_mat_prime!=0))
  temp_DAG_all_patients_this_iter[p_reg,l] = as.numeric(sum_num_changes!=0)

  if(check_for_DAGs(temp_DAG_all_patients_this_iter))
  {
    temp_errors = response_vec- node_mat * beta_mat_prime
    
    llh_prime<- ALD_log_likelihood(  temp_errors, tau_quantile, indicator_th)
    ######################################################################
    
    log_AR_j<-(llh_prime + prior_den_prime) - (llh_old+prior_den_old)
    u<-runif(1,0,1)
    
    if(log_AR_j>=log(u)){
      eta_mat_old=eta_j_prime
      accept_vec=1
      llh_old=llh_prime
      #DAG_all_patients_this_iter[,l,p_reg] = as.numeric(beta_mat_prime!=0)
    }else{
      accept_vec=0
      #sum_num_changes = 0
      sum_num_changes = num_beta_nnz
    }
    
  }else{
    accept_vec=0
    #sum_num_changes = 0
    sum_num_changes = num_beta_nnz
  }
  res<-list(eta_mat_out=eta_mat_old,accept=accept_vec,
            llh_out = llh_old, 
            #dag_all = DAG_all_patients_this_iter
            sum_num_changes = sum_num_changes)
  return(res)
}
##################################################################################
mix_normal_MH_WDC <- function(t,
                              eta_mat,linearCov_mat,nonLinearCov_mat,
                              node_mat,response_vec,
                              tau_quantile, mix_normal_lin, mix_normal_nonlin,
                              mix_nor_std, 
                              intTerm, modified_gamma_mat,llh_old, indicator_th,
                              DAG_all_patients_this_iter,
                              num_beta_nnz)
{
  accept_vec = NA
  
  mix_normal_lin_old = mix_normal_lin
  mix_normal_nonlin_old = mix_normal_nonlin
  
  
  binom_lin = 2*rbinom(q,1,1/(1+exp(-2*mix_normal_lin_old))) -1
  binom_nonlin = 2*rbinom(sum(b_star),1,1/(1+exp(-2*mix_normal_nonlin_old))) -1
  
  prop_mix_normal_lin = mix_normal_lin_old + mix_nor_std*rnorm(q,0,1)
  prop_mix_normal_nonlin = mix_normal_nonlin_old + mix_nor_std*rnorm(sum(b_star),0,1)
  
  mix_normal_lin_new = mix_normal_lin_old
  mix_normal_lin_new = prop_mix_normal_lin
  
  mix_normal_nonlin_new = mix_normal_nonlin_old
  mix_normal_nonlin_new = prop_mix_normal_nonlin
  
  theta_mat_new <-intTerm + 
    llh_theta_lin(eta_mat,linearCov_mat, mix_normal_lin_new) + 
    llh_theta_nonlin(nonLinearCov_mat,modified_gamma_mat , mix_normal_nonlin_new)
  
  beta_mat_new<- theta_mat_new * as.numeric(abs(theta_mat_new)>t)
  
  temp_DAG_all_patients_this_iter = DAG_all_patients_this_iter
  sum_num_changes = sum(as.numeric(beta_mat_new!=0))
  temp_DAG_all_patients_this_iter[p_reg,l] = as.numeric(sum_num_changes!=0)

  if(check_for_DAGs(temp_DAG_all_patients_this_iter))
  {
    temp_errors = response_vec - node_mat*beta_mat_new
    
    llh_new<- ALD_log_likelihood(  temp_errors, tau_quantile,indicator_th )
    
    ######################################################################
    
    llh_diff = llh_new - llh_old - 0.5*(sum((prop_mix_normal_lin - binom_lin)^2))-
      0.5*(sum((prop_mix_normal_nonlin - binom_nonlin)^2))+
      0.5*(sum((mix_normal_lin_old - binom_lin)^2))+
      0.5*(sum((mix_normal_nonlin_old - binom_nonlin)^2))
    
    if(llh_diff>log(runif(1)))
    {
      mix_normal_lin_old = mix_normal_lin_new
      mix_normal_nonlin_old = mix_normal_nonlin_new
      llh_old = llh_new
      accept_vec = 1
      #DAG_all_patients_this_iter[,l,p_reg] = as.numeric(beta_mat_new!=0)
    }else{
      accept_vec = 0
      #sum_num_changes = 0
      sum_num_changes = num_beta_nnz
    }
    
  }else{
    accept_vec=0
    #sum_num_changes = 0
    sum_num_changes = num_beta_nnz
  }
  
  res = list(lin_out = mix_normal_lin_old,
             nonlin_out = mix_normal_nonlin_old,
             accept_out = accept_vec,
             llh_out = llh_old, 
             #dag_all = DAG_all_patients_this_iter
             sum_num_changes = sum_num_changes)
}
##################################################################################
SWAP_function <- function()
{
  for(p_reg in 1:num_of_p)
  {
    for(l in setdiff(1:num_of_p,p_reg))
    {
      if(union_DAG[p_reg, l] == 1)
      {
        temp_union_DAG = union_DAG
        temp_union_DAG[p_reg, l] = 0
        temp_union_DAG[l, p_reg] = 1
        
        if(is.DAG(temp_union_DAG))
        {
          # We will proceed with the swap iff the graph with flipped edge is a DAG
          initial_LLH_1 = LLH_MC[p_reg,l,iter_num  + 1]
          initial_LLH_2 = LLH_MC[l,p_reg,iter_num  + 1]
          
          beta_vec_2b_swap_1 =  BETA_MAT_MC[,l,p_reg,iter_num + 1]
          beta_vec_2b_swap_2 =  BETA_MAT_MC[,p_reg,l,iter_num + 1]
          
          temp_beta_vec_swap = beta_vec_2b_swap_1
          beta_vec_2b_swap_1 = beta_vec_2b_swap_2
          beta_vec_2b_swap_2 = temp_beta_vec_swap
          
          x_dot_beta_swap_1 = beta_vec_2b_swap_1 * rawX[,l,p_reg]
          x_dot_beta_swap_2 = beta_vec_2b_swap_2 * rawX[,p_reg,l]
          
          X_DOT_BETA_MAT_SWAP_1 = X_DOT_BETA_MC[,,p_reg]
          X_DOT_BETA_MAT_SWAP_2 = X_DOT_BETA_MC[,,l]
          
          X_DOT_BETA_MAT_SWAP_1[,l] = x_dot_beta_swap_1
          X_DOT_BETA_MAT_SWAP_2[,p_reg] = x_dot_beta_swap_2
          
          nzv_fii_llh_1 = apply(X_DOT_BETA_MAT_SWAP_1[,setdiff(1:num_of_p,l)],1,sum)
          nzv_fii_llh_2 = apply(X_DOT_BETA_MAT_SWAP_2[,setdiff(1:num_of_p,p_reg)],1,sum)
          
          temp_error_1 = Y[,p_reg] - x_dot_beta_swap_1
          temp_error_2 = Y[,l] - x_dot_beta_swap_2
          
          swap_LLH_1 = ALD_log_likelihood(temp_error_1, tau_index, nzv_fii_llh_1)
          swap_LLH_2 = ALD_log_likelihood(temp_error_2, tau_index, nzv_fii_llh_2)
          
          LLH_diff_af_and_bf_swap = (swap_LLH_1 + swap_LLH_2) - (initial_LLH_1 + initial_LLH_2)
          
          if(LLH_diff_af_and_bf_swap > log(runif(1,0,1)))
          {
            # This means nodes should be swapped and all associated parameters have to swapped
            union_DAG = temp_union_DAG
            swap_iter = iter_num + 1
            ######################################
            var_type_1 = c("THRESHOLD_MC", "MU_MC", "TAU_MC", "XI_MC", "ALPHA_MC", "RHO_MC")
            for(var_2b_swapped in var_type_1)
            {
              command = paste0('temp_var = ',var_2b_swapped,'[p_reg,l,swap_iter];',
                               var_2b_swapped,'[p_reg,l,swap_iter] =',var_2b_swapped,'[l,p_reg,swap_iter];',
                               var_2b_swapped,'[l,p_reg,swap_iter] = temp_var;')
              #print(command)
              eval(parse(text= command))
            }
            ######################################
            var_type_2 = c("GAMMA_MC", "LAMBDA_MC", "ZETA_MC", "ETA_MC", "KAPPA_MC", "PHI_MC",
                           "MIX_NOR_ETA_MC", "MIX_NOR_GAMMA_MC", "THETA_MAT_MC", "INT_TERM_MC",
                           "LIN_TERM_MC", "NONLIN_TERM_MC", "MODIFIED_GAMMA_MAT_MC")
            for(var_2b_swapped in var_type_2)
            {
              command = paste0('temp_var = ',var_2b_swapped,'[,l,p_reg,swap_iter];',
                               var_2b_swapped,'[,l,p_reg,swap_iter] =',var_2b_swapped,'[,p_reg,l,swap_iter];',
                               var_2b_swapped,'[,p_reg,l,swap_iter] = temp_var;')
              #print(command)
              eval(parse(text= command))
            }
            ######################################
            LLH_MC[p_reg,l,swap_iter] = swap_LLH_1
            LLH_MC[l,p_reg,swap_iter] = swap_LLH_2
            ######################################
            X_DOT_BETA_MC[,l,p_reg] = x_dot_beta_swap_1
            X_DOT_BETA_MC[,p_reg,l] = x_dot_beta_swap_2
            ######################################
            BETA_MAT_MC[,l,p_reg,swap_iter] = beta_vec_2b_swap_1 
            BETA_MAT_MC[,p_reg,l,swap_iter] = beta_vec_2b_swap_2 
            ######################################
            Total_edge_swap[p_reg,l] = Total_edge_swap[p_reg,l] + 1
            print("Edge Swapped")
            ######################################
          }
        }
      }
    }
  }
  
}



















