if(!require("splines",quietly = TRUE)){install.packages("splines")}
if(!require("R.matlab", quietly = TRUE)){install.packages("R.matlab")}

library(splines)
library(R.matlab)
source("pxHS_Functions.R")
sample_size = 10
num_of_p = 5
num_of_q = 2
KT_val = 0.5
number_of_simulations = 5

if(num_of_q == 2)
{
  true_th = "0p5"
}else if(num_of_q == 5)
{
  true_th = "1"
}

if(KT_val == 0.5)
{
  kendall_value = "0p5"
}else if(KT_val == 0.25)
{
  kendall_value = "0p25"
}

Results_directory_name = paste0("Misspecified_HS_Results_KT_",kendall_value,"_with_th_",true_th,"_n_",sample_size,"_p_",num_of_p,"_q_",num_of_q)
command_mkdir = paste0("mkdir ",Results_directory_name)
system(command_mkdir)

misspecified_order = read.csv(paste0("Misspecified_order_with_KT_",kendall_value,"_th_",true_th,"_n_",sample_size,
                                     "_p_",num_of_p,"_q_",num_of_q,".csv"), header = FALSE)
misspecified_order = as.matrix(misspecified_order)
misspecified_order = as.vector(misspecified_order)


for(simulation_index in c(1:number_of_simulations)){
  
  for(p_regression in c(1:(num_of_p-1))){
    
    for(TAU_INDEX in c(1:9)){
      
      misspecified_p_regression = misspecified_order[p_regression]
      tau_index =TAU_INDEX/10
      print_satement = paste0("Dataset-",simulation_index,". Starting QR at node ",misspecified_p_regression,", quantile level ",tau_index,"\n")
      cat(print_satement)
      ###########################################################
      file_name = paste0("data_with_th_",true_th,"_n_",sample_size,"_p_",num_of_p,"_q_",num_of_q,
                         "/QR_X_sim_n_", sample_size, "_seed_",simulation_index,".csv", sep="")
      
      
      rawU = read.csv(file_name, header = FALSE)
      rawU = as.matrix(rawU)
      ###########################################################
      filename_output = paste0("data_with_th_",true_th,"_n_",sample_size,"_p_",num_of_p,"_q_",num_of_q,
                               "/QR_Y_sim_n_", sample_size, "_seed_",simulation_index,".csv", sep="")
      temp = read.csv(filename_output, header = FALSE)
      temp = as.matrix(temp)
      QR_P_sim = temp[,misspecified_p_regression]
      Y = QR_P_sim
      ###########################################################
      QR_Simple_P_sim =cbind(rep(1,sample_size), temp[,misspecified_order[(p_regression+1):num_of_p]]);
      rawX = QR_Simple_P_sim
      rawX = as.matrix(rawX)
      q=ncol(rawU)
      q_input = q;
      n=length(Y)
      p=ncol(rawX)
      ############################################################
      
      #### Initailization:B-spline ####
      order = 4 #order of splines
      nknots = 16 #number of interior knots
      M = order + nknots
      basis_obj<-list()
      basisMat = matrix(NA,nrow=n,ncol=0)
      
      for(i in 1:q){
        basis_obj[[i]]<-bs(rawU[,i],knots=quantile(unique(rawU[,i]),seq(0,1,1/(nknots+1)))[2:(nknots+1)],intercept=T)
        basisMat=cbind(basisMat,basis_obj[[i]])
      }
      
      #### Normalization  ####
      K = makeK(M);
      Kinv = pracma::pinv(K);
      U_nonLin = matrix(NA,nrow=n,ncol=0)
      U_lin=matrix(NA,nrow=nrow(rawU),ncol=ncol(rawU))
      b_star =vector()
      
      
      for(i in 1:q){
        #print(i)
        Utmp = basisMat[,((i-1)*M+1):(i*M)]
        tmp_svd=svd(Utmp %*% Kinv %*% t(Utmp))
        svd_U=tmp_svd$u
        svd_S=tmp_svd$d
        nullvals = svd_S < 10^-10;
        d = max(3, min(which( cumsum(svd_S[!nullvals])/sum(svd_S[!nullvals]) > .995 )))
        d = min(d, sum(!nullvals))
        Utmp = svd_U[,1:d] %*% diag(sqrt(svd_S[1:d]))
        Utmp2 = cbind(1,rawU[,i])
        Utmp = Utmp - Utmp2 %*% qr.solve(Utmp2,Utmp)
        Utmp = Utmp/norm(Utmp,'F')
        U_lin[,i] = rawU[,i]-mean(rawU[,i])
        U_lin[,i] = U_lin[,i]/norm(as.matrix(U_lin[,i]),'F')
        U_nonLin<-cbind(U_nonLin,Utmp)
        b_star = c(b_star,dim(Utmp)[2])
      }
      
      
      #### Initialization: Coefficient Value ####
      set.seed(123456+simulation_index*1000+q_input*10)
      
      zeta_<-abs(matrix(rnorm(q*p,mean=1,sd=0.5),nrow=q,ncol=p))
      lambda_<-sqrt(abs(matrix(rnorm(q*p,mean=1,sd=0.5),nrow=q,ncol=p)))
      xi_<-abs(rnorm(p,mean=1,sd=0.5))
      tau_<-sqrt(abs(rnorm(p,mean=1,sd=0.5)))
      
      rho_<- abs(rnorm(p,mean=1,sd=0.5))
      alpha_<-sqrt(abs(rnorm(p,mean=1,sd=0.5)))
      phi_<-matrix(abs(rnorm(q*p,mean=1,sd=0.5)),nrow=q,ncol=p)
      kappa_<-sqrt(abs(matrix(rnorm(q*p,mean=1,sd=0.5),nrow=q,ncol=p)))
      
      gamma_<-matrix(rnorm(q*p,1,0.5),nrow=q,ncol=p)
      eta_<-matrix(rnorm(q*p,1,0.5),nrow=q,ncol=p)
      
      temp_binom = rbinom(sum(b_star)*p,1,0.5)
      temp_binom = 2*temp_binom - 1
      mix_nor_gamma_<-matrix(rnorm(sum(b_star)*p,temp_binom,0.5),nrow=sum(b_star),ncol=p)
      
      temp_binom = rbinom(q*p,1,0.5)
      temp_binom = 2*temp_binom - 1
      mix_nor_eta_<-matrix(rnorm(q*p,temp_binom,0.5),nrow=q,ncol=p)
      
      #########################
      sigma_mu<-0.5
      mu_<-rnorm(p,1,0.5)
      
      a_epi<-1
      b_epi<-1
      
      t_<-rgamma(1,shape=10,scale = 0.1) #choosing such that mean of gamma is 1
      
      ##### some additional params for t update wwith normal proposal ######
      
      t_sd_vec = 0.1*c(1/16,1/8,1/4,1/2,1,2,4,8,16) 
      t_iter = 5
      t_std = t_sd_vec[t_iter]
      
      mix_nor_std = 0.1
      eta_gamma_std = 0.1
      sigma_int<-0.5
      sigma_t<-0.5
      
      
      #### MCMC ####
      set.seed(123456789+simulation_index*100000+q_input*100)
      iteNum<-20000
      burnIn<-10000
      thinning = 10
      
      post_size = (iteNum - burnIn)/thinning
      
      coef_mat<-matrix(NA,nrow=(iteNum - burnIn)/thinning,
                       ncol=p+q*p+q*p+q*p+sum(b_star)*p+1)
      
      colnames(coef_mat)<-c(paste0("mu_",1:p),paste0("eta_",1:(q*p)),
                            paste0("gamma_",1:(q*p)),"t",
                            paste0("mix_norm_eta_",1:(q*p)),
                            paste0("mix_norm_gamma_",1:(sum(b_star)*p)))
      llh_vec<-rep(NA,iteNum)
      hyperParam_mat<-matrix(NA,nrow=(iteNum - burnIn)/thinning,ncol=2*(p+q*p)+2*(p+q*p))
      colnames(hyperParam_mat)<-c(paste0("alpha_",1:p),paste0("rho_",1:p),paste0("kappa_",1:(q*p)),paste0("phi_",1:(q*p)),
                                  paste0("tau_",1:p),paste0("xi_",1:p),paste0("lambda_",1:(q*p)),paste0("zeta_",1:(q*p)))
      
      MH_acceptance<-rep(0,(1+p+p+p+p))
      names(MH_acceptance)<-c(paste0("mu_",1:p),paste0("eta_",1:p),paste0("gamma_",1:p),"t",
                              paste0("mix_nor_",1:p))
      
      record_flag<-F
      to_be_saved_idx = 0
      t_saves = rep(0, iteNum)
      
      ##############################################################
      modified_gamma_mat = matrix(0, nrow = sum(b_star), ncol = p)
      
      for(mod_gamma_col in 1:p)
      {
        modified_gamma_mat[,mod_gamma_col] = rep(gamma_[,mod_gamma_col], b_star) 
      }
      
      intTerm = llh_theta_int(mu_)
      linTerm = llh_theta_lin(eta_,U_lin , mix_nor_eta_)
      nonlinTerm = llh_theta_nonlin(U_nonLin, modified_gamma_mat, mix_nor_gamma_)
      
      theta_mat<- intTerm + linTerm + nonlinTerm
      
      beta_mat<- theta_mat * as.numeric(abs(theta_mat)>t_)
      temp_errors = Y - apply(rawX * beta_mat,1,sum)
      
      llh_old<- ALD_log_likelihood(temp_errors, tau_index)
      v = array(0, dim = c(sample_size, ncol(rawX), post_size ))
      v_Q = array(0, dim = c(num_of_q, ncol(rawX)-1,  post_size ))
      ##############################################################
      time_rec<-system.time({
        for(ite in 1:iteNum){
          
          if(ite%%5000==0){
            print_satement = paste0(ite," number of iterations complete","\n")
            cat(print_satement)
          }
          
          if(ite %in% seq(burnIn+1, iteNum, thinning))
          {
            record_flag = T
            to_be_saved_idx = to_be_saved_idx + 1
          }
          #print("update nonlinear")
          #update nonlinear coefficient: MH
          
          gamma_MH<-Update_nonLin_gamma(t = t_, 
                                        gamma_mat =  gamma_, 
                                        nonLinearCov_mat = U_nonLin, node_mat = rawX, response_vec = Y,
                                        sigma_gamma = eta_gamma_std, tau_vec = tau_, lambda_mat =lambda_,
                                        tau_quantile = tau_index,
                                        mix_normal_nonlin = mix_nor_gamma_,
                                        intTerm = intTerm,
                                        linTerm = linTerm,
                                        modified_gamma_mat = modified_gamma_mat,
                                        llh_old = llh_old
          )
          if(record_flag){coef_mat[to_be_saved_idx,paste0("gamma_",1:(q*p))]<-as.vector(gamma_MH$gamma_mat_out)}
          gamma_ <- gamma_MH$gamma_mat_out
          modified_gamma_mat = gamma_MH$modified_gamma_mat_out
          llh_old = gamma_MH$llh_out
          nonlinTerm = llh_theta_nonlin(U_nonLin, modified_gamma_mat, mix_nor_gamma_)
          MH_acceptance[paste0("gamma_",1:p)]<- MH_acceptance[paste0("gamma_",1:p)] + gamma_MH$accept
          
          #update nonlinear coefficient: hyperParam
          tau_Gibbs<-sqrt(Update_nonLin_tau2(xi_vec = xi_, lambda_mat = lambda_, gamma_mat = gamma_, old_tau2 = tau_^2))
          if(record_flag){hyperParam_mat[to_be_saved_idx,paste0("tau_",1:p)]<-tau_Gibbs}
          tau_<-tau_Gibbs
          
          xi_Gibbs<-Update_nonLin_xi(tau = tau_)
          if(record_flag){hyperParam_mat[to_be_saved_idx,paste0("xi_",1:p)]<-xi_Gibbs}
          xi_<-xi_Gibbs
          
          tau_aug<-matrix(rep(tau_,each=q),ncol=p)
          lambda_Gibbs<-sqrt(Update_nonLin_lambda2(zeta = zeta_, tau_augment = tau_aug, gamma_mat=gamma_))
          if(record_flag){hyperParam_mat[to_be_saved_idx,paste0("lambda_",1:(q*p))]<-as.vector(lambda_Gibbs)}
          lambda_<-lambda_Gibbs
          
          zeta_Gibbs<-Update_nonLin_zeta(lambda_mat = lambda_)
          if(record_flag){hyperParam_mat[to_be_saved_idx,paste0("zeta_",1:(q*p))]<-zeta_Gibbs}
          zeta_<-zeta_Gibbs
          
          #print("update linear")
          #update linear coefficient: MH
          eta_MH<-Update_lin_eta(t = t_ ,
                                 eta_mat = eta_ ,linearCov_mat = U_lin,node_mat = rawX,response_vec = Y,
                                 sigma_eta = eta_gamma_std,alpha_vec = alpha_,kappa_mat = kappa_, 
                                 tau_quantile = tau_index,
                                 mix_normal_lin = mix_nor_eta_,
                                 intTerm = intTerm,
                                 nonlinTerm = nonlinTerm,
                                 llh_old = llh_old
          )
          if(record_flag){coef_mat[to_be_saved_idx,paste0("eta_",1:(q*p))]<-eta_MH$eta_mat_out}
          eta_<-eta_MH$eta_mat_out
          linTerm = llh_theta_lin(eta_,U_lin , mix_nor_eta_)
          llh_old = eta_MH$llh_out
          MH_acceptance[paste0("eta_",1:p)]<- MH_acceptance[paste0("eta_",1:p)] + eta_MH$accept
          
          #update linear coefficient: hyperParam
          alpha_Gibbs<-sqrt(Update_lin_alpha2(rho_vec = rho_, kappa_mat = kappa_, eta_mat = eta_, old_alpha2 = alpha_^2))
          if(record_flag){hyperParam_mat[to_be_saved_idx,paste0("alpha_",1:p)]<-alpha_Gibbs}
          alpha_<-alpha_Gibbs
          
          rho_Gibbs<-Update_lin_rho(alpha_)
          if(record_flag){hyperParam_mat[to_be_saved_idx,paste0("rho_",1:p)]<-rho_Gibbs}
          rho_<-rho_Gibbs
          
          alpha_aug<-matrix(rep(alpha_,each=q),ncol=p)
          kappa_Gibbs<-sqrt(Update_lin_kappa2(phi_mat = phi_, alpha_augment = alpha_aug, eta_mat = eta_))
          if(record_flag){hyperParam_mat[to_be_saved_idx,paste0("kappa_",1:(q*p))]<-kappa_Gibbs}
          kappa_<-kappa_Gibbs
          
          phi_Gibbs<-Update_lin_phi(kappa_)
          if(record_flag){hyperParam_mat[to_be_saved_idx,paste0("phi_",1:(q*p))]<-phi_Gibbs}
          phi_<-phi_Gibbs
          
          #print("update int")
          #update intercept: MH
          mu_MH<-Update_Int(t = t_, mu_vec = mu_, 
                            node_mat = rawX, response_vec = Y, 
                            sigma_int = sigma_int, sigma_mu = sigma_mu,
                            tau_quantile = tau_index,
                            linTerm = linTerm,
                            nonlinTerm = nonlinTerm,
                            llh_old = llh_old
          )
          if(record_flag){coef_mat[to_be_saved_idx,paste0("mu_",1:p)]<-mu_MH$mu_vec_out}
          mu_ <- mu_MH$mu_vec_out
          intTerm = llh_theta_int(mu_)
          llh_old = mu_MH$llh_out
          MH_acceptance[paste0("mu_",1:p)]<- MH_acceptance[paste0("mu_",1:p)] + mu_MH$accept
          
          t_MH<-Update_t_Ga_NEW(t_old = t_, 
                                node_mat = rawX, response_vec = Y,
                                tau_quantile = tau_index,
                                t_std = t_std,
                                a_t = 10,
                                b_t = 0.1,
                                intTerm = intTerm,
                                linTerm = linTerm,
                                nonlinTerm = nonlinTerm,
                                llh_old = llh_old
          )
          ######
          t_saves[ite] = t_MH$t_new
          ######
          
          if(record_flag){coef_mat[to_be_saved_idx,"t"]<- t_MH$t_new}
          t_ <- t_MH$t_new
          llh_old = t_MH$llh_out
          MH_acceptance["t"]<-MH_acceptance["t"] + t_MH$accept
          
          ############### Steps to decide t_std #############
          if(ite%%100 ==0 && ite<burnIn)
          {
            ac_t_monitor = as.numeric(MH_acceptance["t"])/(ite -1)
            
            if(ac_t_monitor > 0.5 && t_iter < 9)
            {
              t_iter = t_iter +1
              t_std = t_sd_vec[t_iter]
            }else if(ac_t_monitor < 0.2 && t_iter > 1)
            {
              t_iter = t_iter -1
              t_std = t_sd_vec[t_iter]
            }
          }
          ###################################################
          
          mix_normal_updates = mix_normal_MH(t = t_ ,
                                             eta_mat = eta_ ,
                                             linearCov_mat = U_lin,
                                             nonLinearCov_mat = U_nonLin,
                                             node_mat = rawX,response_vec = Y,
                                             tau_quantile = tau_index, 
                                             mix_normal_lin = mix_nor_eta_, 
                                             mix_normal_nonlin = mix_nor_gamma_,
                                             mix_nor_std = mix_nor_std,
                                             intTerm = intTerm,
                                             modified_gamma_mat = modified_gamma_mat,
                                             llh_old = llh_old)
          
          mix_nor_eta_ = mix_normal_updates$lin_out
          mix_nor_gamma_ = mix_normal_updates$nonlin_out
          llh_old = mix_normal_updates$llh_out
          linTerm = llh_theta_lin(eta_,U_lin , mix_nor_eta_)
          nonlinTerm = llh_theta_nonlin(U_nonLin, modified_gamma_mat, mix_nor_gamma_)
          
          #########################
          
          MH_acceptance[paste0("mix_nor_",1:p)]<- MH_acceptance[paste0("mix_nor_",1:p)] + mix_normal_updates$accept_out
          llh_vec[ite]<- (mix_normal_updates$llh_out)
          
          if(record_flag)
          {
            coef_mat[to_be_saved_idx,paste0("mix_norm_gamma_",1:(sum(b_star)*p))]<-
              as.vector(mix_nor_gamma_)
          }
          
          if(record_flag)
          {
            coef_mat[to_be_saved_idx,paste0("mix_norm_eta_",1:(q*p))]<-
              as.vector(mix_nor_eta_)
          }
          
          if(record_flag)
          {
            theta_mat = intTerm + linTerm + nonlinTerm
            beta_mat = theta_mat * as.numeric(abs(theta_mat)>t_)
            temp_intercept_column = theta_mat[,1]
            beta_mat[,1]  = temp_intercept_column
            v[,,to_be_saved_idx] = beta_mat
            #################################################
            tmp_mat = tau_aug*lambda_
            nlinrate = 1-1/(1+tmp_mat^2)
            
            tmp_mat = alpha_aug*kappa_
            linrate = 1-1/(1+tmp_mat^2)
            
            lnrate<-pmax(linrate,nlinrate)
            lnrate = lnrate[,-1]
            v_Q[,,to_be_saved_idx] = lnrate
          }
          ###################################################
        }
      })
      
      Time_Taken = as.numeric(time_rec["elapsed"])
      tau_index_save = tau_index *10
      print(time_rec)
      
      
      ###################################################
      ############Post Processing #######################
      
      hs_mean_coef = as.matrix(apply(coef_mat,2,mean))
      hs_mean_coef = t(hs_mean_coef)
      
      hs_mu_vec<-hs_mean_coef[,paste0("mu_",1:p)]
      hs_eta_mat<-matrix(hs_mean_coef[,paste0("eta_",1:(q*p))],ncol=p)
      hs_gamma_mat<-matrix(hs_mean_coef[,paste0("gamma_",1:(q*p))],ncol=p)
      hs_mix_norm_eta_mat<-matrix(hs_mean_coef[,paste0("mix_norm_eta_",1:(q*p))],ncol=p)
      hs_mix_norm_gamma_mat<-matrix(hs_mean_coef[,paste0("mix_norm_gamma_",1:(sum(b_star)*p))],ncol=p)
      
      modified_gamma_mat = matrix(0, nrow = sum(b_star), ncol = p)
      
      for(mod_gamma_col in 1:p)
      {
        modified_gamma_mat[,mod_gamma_col] = rep(hs_gamma_mat[,mod_gamma_col], b_star) 
      }
      
      hs_t_est<-hs_mean_coef[,"t"]
      EST_threshold = hs_t_est
      
      hs_theta_est<-llh_theta_int(hs_mu_vec) + llh_theta_lin(hs_eta_mat,U_lin,hs_mix_norm_eta_mat) + llh_theta_nonlin(U_nonLin,modified_gamma_mat ,hs_mix_norm_gamma_mat)
      temp_intercept_column = hs_theta_est[,1]
      hs_beta_est<- hs_theta_est * as.numeric(abs(hs_theta_est)>hs_t_est)
      hs_beta_est[,1] = temp_intercept_column
      #########################################################
      hs_Y_est<- apply(rawX * hs_beta_est,1,sum)
      
      general_est_beta = matrix(0, sample_size, num_of_p)
      general_est_theta = matrix(0, sample_size, num_of_p)
      
      general_est_beta[,misspecified_order[(p_regression+1):num_of_p]] = hs_beta_est[,2:ncol(rawX)]
      general_est_theta[,misspecified_order[(p_regression+1):num_of_p]] = hs_theta_est[,2:ncol(rawX)]
      
      general_est_beta = general_est_beta[,-c(misspecified_p_regression)]
      general_est_theta = general_est_theta[,-c(misspecified_p_regression)]
      
      
      if (misspecified_p_regression != num_of_p)
      {
        file_name = paste0("data_with_th_",true_th,"_n_",sample_size,"_p_",num_of_p,"_q_",num_of_q,
                           "/True_betas/True_beta_matrix_n",sample_size,"_p_",num_of_p,"_q_",num_of_q,"_p_reg_",misspecified_p_regression,"_tau_",tau_index_save,"_simu_",simulation_index,".csv")
        True_beta_matrix = read.csv(file_name, header = FALSE)
        True_beta_matrix = as.matrix(True_beta_matrix)
        True_beta_matrix = cbind(matrix(0, sample_size, misspecified_p_regression), True_beta_matrix)
        True_beta_matrix = True_beta_matrix[,-misspecified_p_regression]
        
        file_name = paste0("data_with_th_",true_th,"_n_",sample_size,"_p_",num_of_p,"_q_",num_of_q,
                           "/True_thetas/True_theta_matrix_n",sample_size,"_p_",num_of_p,"_q_",num_of_q,"_p_reg_",misspecified_p_regression,"_tau_",tau_index_save,"_simu_",simulation_index,".csv")
        True_theta_matrix = read.csv(file_name, header = FALSE)
        True_theta_matrix = as.matrix(True_theta_matrix)
        True_theta_matrix = cbind(matrix(0, sample_size, misspecified_p_regression), True_theta_matrix)
        True_theta_matrix = True_theta_matrix[,-misspecified_p_regression]
        
      }else{
        True_beta_matrix = matrix(0, sample_size, num_of_p -1)
        True_theta_matrix = matrix(0, sample_size, num_of_p -1)
      }
      
      FROBENIOUS_norm_beta = norm(True_beta_matrix - general_est_beta, "F")
      FROBENIOUS_norm_theta = norm(True_theta_matrix - general_est_theta, "F")
      
      #########################################
      if(misspecified_p_regression != num_of_p)
      {
        filename_output = paste0("data_with_th_",true_th,"_n_",sample_size,"_p_",num_of_p,"_q_",num_of_q,"/QR_Y_sim_TRUEs_n_",sample_size,
                                 "_seed_",simulation_index,".csv");
        temp = read.csv(filename_output, header = FALSE);
        QR_P_sim_TRUEs = temp[,((misspecified_p_regression-1)*9+1):(misspecified_p_regression*9)];
        QR_P_sim_TRUEs = as.matrix(QR_P_sim_TRUEs)
        true_Y = QR_P_sim_TRUEs[,tau_index_save]
      }else{
        true_quantiles = quantile(QR_P_sim, seq(0.1,0.9,0.1))
        QR_P_sim_TRUEs = matrix(rep(as.numeric(true_quantiles),sample_size),ncol = 9, byrow = TRUE)
        true_Y = QR_P_sim_TRUEs[,tau_index_save]
      }
      
      MSE = sum((true_Y - hs_Y_est)^2)/sample_size
      
      True_T_List = readMat(paste0("data_with_th_",true_th,"_n_",sample_size,"_p_",num_of_p,"_q_",num_of_q,
                                   "/Y_order_n_",sample_size,"_p_",num_of_p,"_q_",num_of_q,"_v6.mat"))
      True_G_List = readMat(paste0("data_with_th_",true_th,"_n_",sample_size,"_p_",num_of_p,"_q_",num_of_q,
                                   "/X_order_n_",sample_size,"_p_",num_of_p,"_q_",num_of_q,"_v6.mat"))
      
      adjusted_true_cords = list(NULL)
      adjusted_null_cords = list(NULL)
      for(P_REG in 1:(num_of_p-1))
      {
        adjusted_true_cords[[P_REG]] = as.vector(True_T_List$which.Ys.TP[[P_REG]][[1]]) - P_REG +1
        adjusted_null_cords[[P_REG]] = setdiff(2:(num_of_p - P_REG+1), adjusted_true_cords[[P_REG]])
      }
      
      if(misspecified_p_regression != num_of_p)
      {
        temp_boolean_matrix = matrix(FALSE, nrow = num_of_q, ncol = num_of_p - misspecified_p_regression)
        TP_ts_this_reg_mod  = adjusted_true_cords[[misspecified_p_regression]] - 1
        for(bool_col in 1:length(TP_ts_this_reg_mod))
        {
          temp_G_array = True_G_List$which.Xs.TP[[misspecified_p_regression]][[1]][[bool_col]][[1]]
          temp_G_array = unlist(temp_G_array)
          temp_G_array = as.numeric(temp_G_array)
          
          if(min(temp_G_array)!= -999)
          {
            for(temp_row_num in temp_G_array)
            {
              temp_boolean_matrix[temp_row_num, 
                                  TP_ts_this_reg_mod[bool_col]] = TRUE
            }
          }
        }
        temp_boolean_matrix = cbind(matrix(FALSE, num_of_q, misspecified_p_regression), temp_boolean_matrix)
        temp_boolean_matrix = temp_boolean_matrix[,-misspecified_p_regression]
      }else{
        temp_boolean_matrix = matrix(FALSE, num_of_q, num_of_p -1)
      }
      ###################################################
      th_for_FDR = 0.1
      source("FDR_Control_qDAGx_m.R")
      ###################################################
      summary = c(FROBENIOUS_norm_theta,
                  FROBENIOUS_norm_beta,
                  as.numeric(EST_threshold),
                  TPR_P, FPR_P, AUC_P, TPR_Q, FPR_Q, AUC_Q,MSE,
                  Time_Taken)
      
      file_name = paste0(Results_directory_name,"/Post_process_MC_p_reg_",misspecified_p_regression,
                         "_simulation_index_",simulation_index,
                         "_tau_index_",tau_index_save,".csv")
      
      write.table(summary,file_name,sep = ",",row.names = FALSE,col.names = FALSE )
      ###################################################
    }
  }
}


  