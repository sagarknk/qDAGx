if(!require("splines",quietly = TRUE)){install.packages("splines")}
if(!require("R.matlab", quietly = TRUE)){install.packages("R.matlab")}
if(!require("gRbase", quietly = TRUE)){
  install.packages("BiocManager")
  require("BiocManager")
  BiocManager::install("gRbase")
}

library(splines)
library(R.matlab)
library(gRbase)
source("qDAGx_functions_pxHS.R")
sample_size = 10
num_of_p = 5
num_of_q = 2
th_for_FDR = 0.1
if(num_of_q == 2)
{
  true_th = "0p5"
}else if(num_of_q == 5)
{
  true_th = "1"
}

Results_directory_name = paste0("DAG_Results_with_th_",true_th,"_n_",sample_size,"_p_",num_of_p,"_q_",num_of_q)
command_mkdir = paste0("mkdir ",Results_directory_name)
system(command_mkdir)

number_of_simulations = 5

for(simulation_index in c(1:number_of_simulations)){
  
    for(TAU_INDEX in c(1:9)){
      
      tau_index =TAU_INDEX/10
      print_satement = paste0("Dataset-",simulation_index,". Starting qDAGx at quantile level ",tau_index,"\n")
      cat(print_satement)
      
      ##########################################################
      MSE_Scales = read.csv(paste0("data_with_th_",true_th,"_n_",sample_size,"_p_",num_of_p,"_q_",num_of_q,"/MSE_Scales.csv"), header = FALSE)
      MSE_Scales = as.numeric(as.matrix(MSE_Scales))
      ###########################################################
      file_name = paste0("data_with_th_",true_th,"_n_",sample_size,"_p_",num_of_p,"_q_",num_of_q,
                         "/QR_X_sim_n_", sample_size, "_seed_",simulation_index,".csv", sep="")
      
      
      rawU = read.csv(file_name, header = FALSE)
      rawU = as.matrix(rawU)
      ###########################################################
      filename_output = paste0("data_with_th_",true_th,"_n_",sample_size,"_p_",num_of_p,"_q_",num_of_q,
                               "/QR_Y_sim_n_", sample_size, "_seed_",simulation_index,".csv", sep="")
      temp = read.csv(filename_output, header = FALSE);
      temp = as.matrix(temp)
      rawX = array(0, c(sample_size, num_of_p, num_of_p))
      for(p_reg in 1:num_of_p)
      {
        augumented_temp = temp
        augumented_temp[,p_reg] = matrix(1,sample_size,1)
        rawX[,,p_reg] = augumented_temp
      }
      ###########################################################
      filename_output = paste0("data_with_th_",true_th,"_n_",sample_size,"_p_",num_of_p,"_q_",num_of_q,
                               "/QR_Y_sim_n_", sample_size, "_seed_",simulation_index,".csv", sep="")
      Y = read.csv(filename_output, header = FALSE);
      Y = as.matrix(Y)
      q=ncol(rawU)
      q_input = q;
      n=sample_size
      p=num_of_p
      
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
      
      total_MC_iters = 5000
      burnIn = 2500
      #################################################################
      ZETA_MC = array(0, c(q,p,p,total_MC_iters))
      for(p_reg in 1:num_of_p)
      {
        zeta_<-abs(matrix(rnorm(q*p,mean=1,sd=0.5),nrow=q,ncol=p))
        ZETA_MC[,,p_reg,1] = zeta_
      }
      #################################################################
      LAMBDA_MC = array(0,c(q,p,p,total_MC_iters))
      for(p_reg in 1:num_of_p)
      {
        lambda_<-sqrt(abs(matrix(rnorm(q*p,mean=1,sd=0.5),nrow=q,ncol=p)))
        LAMBDA_MC[,,p_reg,1] = lambda_
      }
      #################################################################
      XI_MC = array(0,c(p,p,total_MC_iters))
      for(p_reg in 1:num_of_p)
      {
        xi_<-abs(rnorm(p,mean=1,sd=0.5))
        XI_MC[p_reg,,1] = xi_  
      }
      #################################################################
      TAU_MC = array(0,c(p,p,total_MC_iters))
      for(p_reg in 1:num_of_p)
      {
        tau_<-sqrt(abs(rnorm(p,mean=1,sd=0.5)))
        TAU_MC[p_reg,,1] = tau_
      }
      #################################################################
      RHO_MC = array(0,c(p,p,total_MC_iters))
      for(p_reg in 1:num_of_p)
      {
        rho_<- abs(rnorm(p,mean=1,sd=0.5))
        RHO_MC[p_reg,,1] = rho_  
      }
      ################################################################
      ALPHA_MC = array(0,c(p,p,total_MC_iters))
      for(p_reg in 1:num_of_p)
      {
        alpha_<-sqrt(abs(rnorm(p,mean=1,sd=0.5)))
        ALPHA_MC[p_reg,,1] = alpha_
      }
      ################################################################
      PHI_MC = array(0,c(q,p,p,total_MC_iters))
      for(p_reg in 1:num_of_p)
      {
        phi_<-matrix(abs(rnorm(q*p,mean=1,sd=0.5)),nrow=q,ncol=p)
        PHI_MC[,,p_reg,1] = phi_  
      }
      ################################################################
      KAPPA_MC = array(0,c(q,p,p,total_MC_iters))
      for(p_reg in 1:num_of_p)
      {
        kappa_<-sqrt(abs(matrix(rnorm(q*p,mean=1,sd=0.5),nrow=q,ncol=p)))
        KAPPA_MC[,,p_reg,1] = kappa_ 
      }
      ################################################################
      GAMMA_MC = array(0,c(q,p,p,total_MC_iters))
      for(p_reg in 1:num_of_p)
      {
        gamma_<-matrix(rnorm(q*p,1,0.5),nrow=q,ncol=p)
        GAMMA_MC[,,p_reg,1] = gamma_  
      }
      ################################################################
      ETA_MC = array(0,c(q,p,p,total_MC_iters))
      for(p_reg in 1:num_of_p)
      {
        eta_<-matrix(rnorm(q*p,1,0.5),nrow=q,ncol=p)
        ETA_MC[,,p_reg,1] = eta_ 
      }
      ################################################################
      MIX_NOR_GAMMA_MC = array(0,c(sum(b_star),p,p,total_MC_iters))
      for(p_reg in 1:num_of_p)
      {
        temp_binom = rbinom(sum(b_star)*p,1,0.5)
        temp_binom = 2*temp_binom - 1
        mix_nor_gamma_<-matrix(rnorm(sum(b_star)*p,temp_binom,0.5),nrow=sum(b_star),ncol=p)
        MIX_NOR_GAMMA_MC[,,p_reg,1] = mix_nor_gamma_ 
      }
      ###############################################################
      MIX_NOR_ETA_MC = array(0,c(q,p,p,total_MC_iters))
      for(p_reg in 1:num_of_p)
      {
        temp_binom = rbinom(q*p,1,0.5)
        temp_binom = 2*temp_binom - 1
        mix_nor_eta_<-matrix(rnorm(q*p,temp_binom,0.5),nrow=q,ncol=p)
        MIX_NOR_ETA_MC[,,p_reg,1] = mix_nor_eta_
      }
      ###############################################################
      MU_MC = array(0,c(p,p,total_MC_iters))
      for(p_reg in 1:num_of_p)
      {
        sigma_mu<-0.5
        mu_<-rnorm(p,1,sigma_mu)
        MU_MC[p_reg,,1] = mu_ 
      }
      ###############################################################
      MODIFIED_GAMMA_MAT_MC = array(0,c(sum(b_star),p,p,total_MC_iters))
      for(p_reg in 1:num_of_p)
      {
        modified_gamma_mat = matrix(0, nrow = sum(b_star), ncol = p)
        
        for(mod_gamma_col in 1:p)
        {
          modified_gamma_mat[,mod_gamma_col] = rep(GAMMA_MC[,mod_gamma_col,p_reg,1], b_star) 
        }
        MODIFIED_GAMMA_MAT_MC[,,p_reg,1] = modified_gamma_mat
      }
      ##############################################################
      INT_TERM_MC = array(0,c(n,p,p,total_MC_iters))
      for(p_reg in 1:num_of_p)
      {
        mu_ = MU_MC[p_reg,,1]
        intTerm = llh_theta_int(mu_)
        INT_TERM_MC[,,p_reg,1] = intTerm
      }
      #############################################################
      LIN_TERM_MC = array(0,c(n,p,p,total_MC_iters))
      for(p_reg in 1:num_of_p)
      {
        linTerm = llh_theta_lin(ETA_MC[,,p_reg,1],U_lin , MIX_NOR_ETA_MC[,,p_reg,1])
        LIN_TERM_MC[,,p_reg,1] = linTerm
      }
      #############################################################
      NONLIN_TERM_MC = array(0,c(n,p,p,total_MC_iters))
      for(p_reg in 1:num_of_p)
      {
        nonlinTerm = llh_theta_nonlin(U_nonLin, MODIFIED_GAMMA_MAT_MC[,,p_reg,1], MIX_NOR_GAMMA_MC[,,p_reg,1])
        NONLIN_TERM_MC[,,p_reg,1] = nonlinTerm
      }
      #############################################################
      THETA_MAT_MC = array(0,c(n,p,p,total_MC_iters))
      THETA_MAT_MC = INT_TERM_MC + LIN_TERM_MC + NONLIN_TERM_MC
      #############################################################
      THRESHOLD_MC = array(0,c(p,p,total_MC_iters))
      for(p_reg in 1:num_of_p)
      {
        for(l in 1:num_of_p)
        {
          temp_max = max(abs(THETA_MAT_MC[,l,p_reg,1]))
          jitter_th = runif(1,0,0.01)
          THRESHOLD_MC[p_reg,l,1] = temp_max + jitter_th
        }
      }
      ##############################################################
      BETA_MAT_MC = array(0,c(n,p,p,total_MC_iters))
      for(per_id in 1:sample_size)
      {
        BETA_MAT_MC[per_id,,,1] = THETA_MAT_MC[per_id,,,1]*as.numeric(abs(t(THETA_MAT_MC[per_id,,,1]))>THRESHOLD_MC[,,1])
      }
      ##############################################################
      X_DOT_BETA_MC = array(0,c(n,p,p))
      X_DOT_BETA_MC = rawX*BETA_MAT_MC[,,,1]
      ##############################################################
      LLH_MC = array(0,c(p,p,total_MC_iters))
      for(p_reg in 1:num_of_p)
      {
        for(l in 1:num_of_p)
        {
          errors = Y[,p_reg] - X_DOT_BETA_MC[,l,p_reg]
          NZV_FII_LLH = apply(X_DOT_BETA_MC[,setdiff(1:num_of_p,l),p_reg], 1, sum)
          LLH_MC[p_reg,l,1] = ALD_log_likelihood(errors,tau_index,indicator_th = NZV_FII_LLH)
        }
      }
      ##############################################################
      union_DAG = matrix(0,num_of_p, num_of_p)
      ###############################################################
      ACCEPTANCE_MC = array(0,c(p,p,5))
      # 1 for mu
      # 2 for eta
      # 3 for eta*
      # 4 for xi and xi*
      # 5 for threshold
      UPDATE_INDEX_MC = array(0,c(p,p,total_MC_iters-1))
      v_Q_MC = array(0, dim = c(q,p,p,total_MC_iters ))
      for(p_reg in 1:num_of_p)
      {
        for(l in 1:num_of_p)
        {
          tau_aug<-rep(TAU_MC[p_reg,l, 1],each=q)
          lambda_ = LAMBDA_MC[,l,p_reg, 1]
          
          tmp_mat = tau_aug*lambda_
          nlinrate = 1-1/(1+tmp_mat^2)
          
          alpha_aug<-rep(ALPHA_MC[p_reg,l, 1],each=q)
          kappa_ = KAPPA_MC[,l,p_reg, 1]
          
          tmp_mat = alpha_aug*kappa_
          linrate = 1-1/(1+tmp_mat^2)
          
          lnrate<-pmax(linrate,nlinrate)
          v_Q_MC[,l,p_reg, 1] = lnrate
        }
      }
      ###############################################################
      mix_nor_std = 0.1
      eta_gamma_std = 0.1
      sigma_int<-0.5
      sigma_t<-0.5
      a_t = 10
      b_t = 0.1
      ##some additional params for t update with normal proposal#####
      t_sd_vec = 0.1*c(1/16,1/8,1/4,1/2,1,2,4,8,16) 
      t_iter = matrix(5,p,p)
      t_std = matrix(0.1,p,p)
      ##############################################################
      ############### MCMC Begins ##################################
      ##############################################################
      set.seed(123456789+simulation_index*100000+q_input*100)
      
      #profvis
      print(paste("MCMC begins"))
      time_rec = system.time({
        for(iter in 2:total_MC_iters)
        {
          iter_num = iter - 1
          if(iter %% 1000 ==0){print(iter)}
          
          ###################
          ## Update mu_jj ###
          ###################
          for(p_reg in 1:num_of_p)
          {
            t_ = THRESHOLD_MC[p_reg,p_reg,iter_num]
            mu_ = MU_MC[p_reg,p_reg,iter_num]
            llh_old = LLH_MC[p_reg,p_reg,iter_num]
            NZV_FII_LLH = apply(X_DOT_BETA_MC[,setdiff(1:num_of_p,p_reg),p_reg], 1, sum)
            
            ### While evaluating log likelihood, we cannot consider zero for the indicator threshold.
            ### Because for a local update we need to compute \Sum_{j !=l}\beta_{jl}^{(i)}*x_{il}
            ### Which forms the threshold in the indicator. Hence we evaluate NZV_FII_LLH
            ### which stands for "Non-Zero value For Indicator In Log-LikeliHood" 
            
            mu_MH<-Update_Int(t = t_, mu_vec = mu_, node_mat = rawX[,p_reg,p_reg],response_vec = Y[,p_reg],
                              sigma_int = sigma_int, sigma_mu = sigma_mu, tau_quantile = tau_index,
                              linTerm = LIN_TERM_MC[,p_reg,p_reg,iter_num],nonlinTerm = NONLIN_TERM_MC[,p_reg,p_reg,iter_num],
                              llh_old = llh_old, indicator_th = NZV_FII_LLH)
            
            MU_MC[p_reg,p_reg,iter_num + 1] = mu_MH$mu_vec_out
            LLH_MC[p_reg,p_reg, iter_num + 1] = mu_MH$llh_out
            ACCEPTANCE_MC[p_reg,p_reg,1] = ACCEPTANCE_MC[p_reg,p_reg,1] + mu_MH$accept
            INT_TERM_MC[,p_reg,p_reg,iter_num + 1] = llh_theta_int(mu_MH$mu_vec_out)
          }
          #####################
          ## Update eta*_jj ###
          #####################
          for(p_reg in 1:num_of_p)
          {
            t_ = THRESHOLD_MC[p_reg,p_reg,iter_num]
            gamma_ = GAMMA_MC[,p_reg, p_reg,iter_num]
            tau_ = TAU_MC[p_reg,p_reg,iter_num]
            lambda_ = LAMBDA_MC[,p_reg,p_reg,iter_num]
            mix_nor_gamma_  = MIX_NOR_GAMMA_MC[,p_reg,p_reg,iter_num]
            modified_gamma_mat = MODIFIED_GAMMA_MAT_MC[,p_reg,p_reg,iter_num]
            llh_old = LLH_MC[p_reg,p_reg, iter_num + 1]
            NZV_FII_LLH = apply(X_DOT_BETA_MC[,setdiff(1:num_of_p,p_reg),p_reg], 1, sum)
            
            gamma_MH<-Update_nonLin_gamma(t = t_, 
                                          gamma_mat =  gamma_, 
                                          nonLinearCov_mat = U_nonLin, node_mat = rawX[,p_reg,p_reg], response_vec =  Y[,p_reg],
                                          sigma_gamma = eta_gamma_std, tau_vec = tau_, lambda_mat =lambda_,
                                          tau_quantile = tau_index,
                                          mix_normal_nonlin = mix_nor_gamma_,
                                          intTerm = INT_TERM_MC[,p_reg,p_reg,iter_num + 1],
                                          linTerm = LIN_TERM_MC[,p_reg,p_reg,iter_num],
                                          modified_gamma_mat = modified_gamma_mat,
                                          llh_old = llh_old, indicator_th = NZV_FII_LLH)
            
            
            GAMMA_MC[,p_reg,p_reg,iter_num + 1] <- gamma_MH$gamma_mat_out
            MODIFIED_GAMMA_MAT_MC[,p_reg,p_reg,iter_num + 1] = gamma_MH$modified_gamma_mat_out
            LLH_MC[p_reg,p_reg,iter_num + 1] = gamma_MH$llh_out
            NONLIN_TERM_MC[,p_reg,p_reg,iter_num + 1] = llh_theta_nonlin(U_nonLin, 
                                                                         MODIFIED_GAMMA_MAT_MC[,p_reg,p_reg,iter_num + 1], 
                                                                         mix_nor_gamma_)
            ACCEPTANCE_MC[p_reg,p_reg,3] = ACCEPTANCE_MC[p_reg,p_reg,3] + gamma_MH$accept
          }
          #####################
          ## Update T_jj ######
          #####################
          for(p_reg in 1:num_of_p)
          {
            xi_ = XI_MC[p_reg,p_reg,iter_num]
            lambda_ = LAMBDA_MC[,p_reg,p_reg,iter_num]
            gamma_ = GAMMA_MC[,p_reg,p_reg,iter_num + 1]
            tau_ = TAU_MC[p_reg,p_reg,iter_num]
            
            tau_Gibbs<-sqrt(Update_nonLin_tau2(xi_vec = xi_, lambda_mat = lambda_, gamma_mat = gamma_, old_tau2 = tau_^2))
            
            TAU_MC[p_reg,p_reg,iter_num + 1]<-tau_Gibbs
          }
          #####################
          ## Update c_jj ######
          #####################
          for(p_reg in 1:num_of_p)
          {
            xi_Gibbs<-Update_nonLin_xi(tau = TAU_MC[p_reg,p_reg,iter_num + 1])
            XI_MC[p_reg,p_reg,iter_num + 1]<-xi_Gibbs
          }
          #####################
          ## Update L_jj ######
          #####################
          for(p_reg in 1:num_of_p)
          {
            
            tau_aug<-rep(TAU_MC[p_reg,p_reg,iter_num + 1],each=q)
            lambda_Gibbs<-sqrt(Update_nonLin_lambda2(zeta = ZETA_MC[,p_reg,p_reg,iter_num], 
                                                     tau_augment = tau_aug, gamma_mat= GAMMA_MC[,p_reg,p_reg,iter_num + 1]))
            
            LAMBDA_MC[,p_reg,p_reg,iter_num + 1]<-lambda_Gibbs
          }
          ########################
          ## Update Zeta_jj ######
          ########################
          for(p_reg in 1:num_of_p)
          {
            zeta_Gibbs<-Update_nonLin_zeta(lambda_mat = LAMBDA_MC[,p_reg,p_reg,iter_num + 1])
            
            ZETA_MC[,p_reg,p_reg,iter_num + 1]<-zeta_Gibbs
          }
          #####################
          ## Update eta_jj ####
          #####################
          for(p_reg in 1:num_of_p)
          {
            t_ = THRESHOLD_MC[p_reg,p_reg,iter_num]
            eta_ = ETA_MC[,p_reg,p_reg,iter_num]
            alpha_ = ALPHA_MC[p_reg,p_reg,iter_num]
            kappa_ = KAPPA_MC[,p_reg,p_reg,iter_num]
            mix_nor_eta_ = MIX_NOR_ETA_MC[,p_reg,p_reg,iter_num]
            NZV_FII_LLH = apply(X_DOT_BETA_MC[,setdiff(1:num_of_p,p_reg),p_reg], 1, sum)
            
            eta_MH<-Update_lin_eta(t = t_ ,
                                   eta_mat = eta_ ,linearCov_mat = U_lin,node_mat = rawX[,p_reg,p_reg],response_vec = Y[,p_reg],
                                   sigma_eta = eta_gamma_std,alpha_vec = alpha_,kappa_mat = kappa_, 
                                   tau_quantile = tau_index,
                                   mix_normal_lin = mix_nor_eta_,
                                   intTerm = INT_TERM_MC[,p_reg,p_reg,iter_num + 1],
                                   nonlinTerm = NONLIN_TERM_MC[,p_reg,p_reg,iter_num + 1],
                                   llh_old = LLH_MC[p_reg,p_reg,iter_num + 1], indicator_th = NZV_FII_LLH)
            
            ETA_MC[,p_reg,p_reg,iter_num + 1] <-eta_MH$eta_mat_out
            LIN_TERM_MC[,p_reg,p_reg, iter_num + 1] = llh_theta_lin(ETA_MC[,p_reg,p_reg,iter_num + 1] ,U_lin , mix_nor_eta_)
            LLH_MC[p_reg,p_reg,iter_num + 1] = eta_MH$llh_out
            ACCEPTANCE_MC[p_reg,p_reg,2]<- ACCEPTANCE_MC[p_reg,p_reg,2] + eta_MH$accept
          }
          #####################
          ## Update A_jj ######
          #####################
          for(p_reg in 1:num_of_p)
          {
            rho_ = RHO_MC[p_reg,p_reg,iter_num]
            kappa_ = KAPPA_MC[,p_reg,p_reg,iter_num]
            eta_ = ETA_MC[,p_reg,p_reg,iter_num + 1]
            alpha_ = ALPHA_MC[p_reg,p_reg,iter_num]
            
            alpha_Gibbs<-sqrt(Update_lin_alpha2(rho_vec = rho_, kappa_mat = kappa_, eta_mat = eta_, old_alpha2 = alpha_^2))
            
            ALPHA_MC[p_reg,p_reg,iter_num + 1]<-alpha_Gibbs
          }
          #######################
          ## Update rho_jj ######
          #######################
          for(p_reg in 1:num_of_p)
          {
            alpha_ = ALPHA_MC[p_reg,p_reg,iter_num + 1]
            rho_Gibbs<-Update_lin_rho(alpha_)
            
            RHO_MC[p_reg,p_reg, iter_num + 1]<-rho_Gibbs
          }
          #####################
          ## Update K_jj ######
          #####################
          for(p_reg in 1:num_of_p)
          {
            alpha_aug<-rep(ALPHA_MC[p_reg,p_reg,iter_num + 1],each=q)
            eta_ = ETA_MC[,p_reg,p_reg,iter_num + 1]
            phi_ = PHI_MC[,p_reg,p_reg,iter_num]
            
            kappa_Gibbs<-sqrt(Update_lin_kappa2(phi_mat = phi_, alpha_augment = alpha_aug, eta_mat = eta_))
            
            KAPPA_MC[,p_reg,p_reg,iter_num + 1]<-kappa_Gibbs
          }
          #######################
          ## Update phi_jj ######
          #######################
          for(p_reg in 1:num_of_p)
          {
            phi_Gibbs<-Update_lin_phi(KAPPA_MC[,p_reg,p_reg,iter_num + 1])
            PHI_MC[,p_reg,p_reg,iter_num + 1]<-phi_Gibbs
          }
          ################################
          ## Update xi_jj and xi*_jj######
          ################################
          for(p_reg in 1:num_of_p)
          {
            t_ = THRESHOLD_MC[p_reg,p_reg,iter_num]
            eta_ = ETA_MC[,p_reg,p_reg,iter_num + 1]
            mix_nor_eta_ = MIX_NOR_ETA_MC[,p_reg,p_reg,iter_num]
            mix_nor_gamma_ = MIX_NOR_GAMMA_MC[,p_reg,p_reg,iter_num]
            intTerm = INT_TERM_MC[,p_reg,p_reg,iter_num + 1]
            modified_gamma_mat = MODIFIED_GAMMA_MAT_MC[,p_reg,p_reg,iter_num + 1]
            llh_old = LLH_MC[p_reg,p_reg,iter_num  + 1]
            NZV_FII_LLH = apply(X_DOT_BETA_MC[,setdiff(1:num_of_p,p_reg),p_reg], 1, sum)
            
            mix_normal_updates = mix_normal_MH(t = t_ ,
                                               eta_mat = eta_ ,
                                               linearCov_mat = U_lin,
                                               nonLinearCov_mat = U_nonLin,
                                               node_mat = rawX[,p_reg,p_reg],response_vec = Y[,p_reg],
                                               tau_quantile = tau_index, 
                                               mix_normal_lin = mix_nor_eta_, 
                                               mix_normal_nonlin = mix_nor_gamma_,
                                               mix_nor_std = mix_nor_std,
                                               intTerm = intTerm,
                                               modified_gamma_mat = modified_gamma_mat,
                                               llh_old = llh_old, indicator_th = NZV_FII_LLH)
            
            MIX_NOR_ETA_MC[,p_reg,p_reg,iter_num + 1] = mix_normal_updates$lin_out
            MIX_NOR_GAMMA_MC[,p_reg,p_reg,iter_num + 1] = mix_normal_updates$nonlin_out
            LLH_MC[p_reg,p_reg,iter_num  + 1] = mix_normal_updates$llh_out
            LIN_TERM_MC[,p_reg,p_reg,iter_num + 1] = llh_theta_lin(eta_,U_lin , MIX_NOR_ETA_MC[,p_reg,p_reg,iter_num + 1])
            NONLIN_TERM_MC[,p_reg,p_reg,iter_num + 1]  = llh_theta_nonlin(U_nonLin, modified_gamma_mat,
                                                                          MIX_NOR_GAMMA_MC[,p_reg,p_reg,iter_num + 1])
            
            ACCEPTANCE_MC[p_reg,p_reg,4]<- ACCEPTANCE_MC[p_reg,p_reg,4]+ mix_normal_updates$accept_out
          }
          #####################
          ## Update t_jj ######
          #####################
          for(p_reg in 1:num_of_p)
          {
            NZV_FII_LLH = apply(X_DOT_BETA_MC[,setdiff(1:num_of_p,p_reg),p_reg], 1, sum)
            
            t_MH<-Update_t_Ga_NEW(t_old = THRESHOLD_MC[p_reg,p_reg,iter_num], 
                                  node_mat = rawX[,p_reg,p_reg], response_vec = Y[,p_reg],
                                  tau_quantile = tau_index,
                                  t_std = t_std,
                                  a_t = 10,
                                  b_t = 0.1,
                                  intTerm = INT_TERM_MC[,p_reg,p_reg,iter_num + 1],
                                  linTerm = LIN_TERM_MC[,p_reg,p_reg,iter_num + 1],
                                  nonlinTerm = NONLIN_TERM_MC[,p_reg,p_reg,iter_num + 1],
                                  llh_old = LLH_MC[p_reg,p_reg,iter_num  + 1], indicator_th = NZV_FII_LLH)
            
            THRESHOLD_MC[p_reg,p_reg,iter_num + 1] <- t_MH$t_new
            LLH_MC[p_reg,p_reg,iter_num  + 1] = t_MH$llh_out
            ACCEPTANCE_MC[p_reg,p_reg,5]<-ACCEPTANCE_MC[p_reg,p_reg,5] + t_MH$accept
            
            ############### Steps to decide t_std #############
            if(iter%%100 ==0 && iter <burnIn)
            {
              ac_t_monitor = as.numeric(ACCEPTANCE_MC[p_reg,p_reg,5])/(iter -1)
              
              if(ac_t_monitor > 0.5 && t_iter[p_reg,p_reg] < 9)
              {
                t_iter[p_reg,p_reg] = t_iter[p_reg,p_reg] +1
                t_std[p_reg,p_reg] = t_sd_vec[t_iter[p_reg,p_reg]]
                
              }else if(ac_t_monitor < 0.2 && t_iter[p_reg,p_reg] > 1)
              {
                t_iter[p_reg,p_reg] = t_iter[p_reg,p_reg] -1
                t_std[p_reg,p_reg] = t_sd_vec[t_iter[p_reg,p_reg]]
              }
            }
            THETA_MAT_MC[,p_reg,p_reg,iter_num + 1] = INT_TERM_MC[,p_reg,p_reg,iter_num + 1] + 
              LIN_TERM_MC[,p_reg,p_reg,iter_num + 1] + 
              NONLIN_TERM_MC[,p_reg,p_reg,iter_num + 1]
            
            BETA_MAT_MC[,p_reg,p_reg,iter_num + 1] = THETA_MAT_MC[,p_reg,p_reg,iter_num + 1]*
              as.numeric(abs(THETA_MAT_MC[,p_reg,p_reg,iter_num + 1]) >  THRESHOLD_MC[p_reg,p_reg,iter_num + 1])
            
            X_DOT_BETA_MC[,p_reg,p_reg] = rawX[,p_reg,p_reg]*BETA_MAT_MC[,p_reg,p_reg,iter_num + 1]
            
            # Threshold update Ends
          }
          #####################
          ## Update t_jl ######
          #####################
          
          update_index = UPDATE_INDEX_MC[,,iter_num]
          for(p_reg in 1:num_of_p)
          {
            for(l in setdiff(1:num_of_p,p_reg))
            {
              NZV_FII_LLH = apply(X_DOT_BETA_MC[,setdiff(1:num_of_p,l),p_reg], 1, sum)
              
              ###############################################################
              num_beta_nnz = sum(BETA_MAT_MC[,l,p_reg,iter_num] !=0)
              ###############################################################
              
              t_MH<-Update_t_Ga_NEW_WDC(t_old = THRESHOLD_MC[p_reg,l,iter_num], 
                                        node_mat = rawX[,l,p_reg], response_vec = Y[,p_reg],
                                        tau_quantile = tau_index,
                                        t_std = t_std,
                                        a_t = 10,
                                        b_t = 0.1,
                                        intTerm = INT_TERM_MC[,l,p_reg,iter_num ],
                                        linTerm = LIN_TERM_MC[,l,p_reg,iter_num ],
                                        nonlinTerm = NONLIN_TERM_MC[,l,p_reg,iter_num ],
                                        llh_old = LLH_MC[p_reg,l,iter_num  ], indicator_th = NZV_FII_LLH,
                                        DAG_all_patients_this_iter = union_DAG,
                                        num_beta_nnz)
              
              THRESHOLD_MC[p_reg,l,iter_num + 1] <- t_MH$t_new
              LLH_MC[p_reg,l,iter_num  + 1] = t_MH$llh_out
              ACCEPTANCE_MC[p_reg,l,5]<-ACCEPTANCE_MC[p_reg,l,5] + t_MH$accept
              
              SUM_NUM_CHANGES = t_MH$sum_num_changes
              union_DAG[p_reg, l] = as.numeric(SUM_NUM_CHANGES != 0)
              
              if(t_MH$dag_check)
              {
                update_index[p_reg,l] = 1
              }
              ############### Steps to decide t_std #############
              if(iter%%100 ==0 && iter <burnIn)
              {
                ac_t_monitor = as.numeric(ACCEPTANCE_MC[p_reg,l,5])/(iter -1)
                
                if(ac_t_monitor > 0.5 && t_iter[p_reg,l] < 9)
                {
                  t_iter[p_reg,l] = t_iter[p_reg,l] +1
                  t_std[p_reg,l] = t_sd_vec[t_iter[p_reg,l]]
                  
                }else if(ac_t_monitor < 0.2 && t_iter[p_reg,l] > 1)
                {
                  t_iter[p_reg,l] = t_iter[p_reg,l] -1
                  t_std[p_reg,l] = t_sd_vec[t_iter[p_reg,l]]
                }
              }
              ##################################################
              THETA_MAT_MC[,l,p_reg,iter_num + 1] = INT_TERM_MC[,l,p_reg,iter_num ] + 
                LIN_TERM_MC[,l,p_reg,iter_num ] + 
                NONLIN_TERM_MC[,l,p_reg,iter_num ]
              
              BETA_MAT_MC[,l,p_reg,iter_num + 1] = THETA_MAT_MC[,l,p_reg,iter_num + 1]*
                as.numeric(abs(THETA_MAT_MC[,l,p_reg,iter_num + 1]) >  THRESHOLD_MC[p_reg,l,iter_num + 1])
              
              X_DOT_BETA_MC[,l,p_reg] = rawX[,l,p_reg]*BETA_MAT_MC[,l,p_reg,iter_num + 1]
              
            }
          }
          ###################
          UPDATE_INDEX_MC[,,iter_num] = update_index
          ###################
          ## Update mu_jl ###
          ###################
          for(p_reg in 1:num_of_p)
          {
            for(l in setdiff(1:num_of_p, p_reg))
            {
              t_ = THRESHOLD_MC[p_reg,l,iter_num + 1]
              mu_ = MU_MC[p_reg,l,iter_num]
              llh_old = LLH_MC[p_reg,l,iter_num + 1]
              NZV_FII_LLH = apply(X_DOT_BETA_MC[,setdiff(1:num_of_p,l),p_reg], 1, sum)
              
              ###############################################################
              num_beta_nnz = sum(BETA_MAT_MC[,l,p_reg,iter_num + 1] !=0)
              ###############################################################
              
              mu_MH<-Update_Int_WDC(t = t_, mu_vec = mu_, node_mat = rawX[,l,p_reg],response_vec = Y[,p_reg],
                                    sigma_int = sigma_int, sigma_mu = sigma_mu, tau_quantile = tau_index,
                                    linTerm = LIN_TERM_MC[,l,p_reg,iter_num],nonlinTerm = NONLIN_TERM_MC[,l,p_reg,iter_num],
                                    llh_old = llh_old, indicator_th = NZV_FII_LLH, 
                                    DAG_all_patients_this_iter = union_DAG,
                                    num_beta_nnz)
              
              MU_MC[p_reg,l,iter_num + 1] = mu_MH$mu_vec_out
              LLH_MC[p_reg,l, iter_num + 1] = mu_MH$llh_out
              ACCEPTANCE_MC[p_reg,l,1] = ACCEPTANCE_MC[p_reg,l,1] + mu_MH$accept
              INT_TERM_MC[,l,p_reg,iter_num + 1] = llh_theta_int(mu_MH$mu_vec_out)
              
              SUM_NUM_CHANGES = mu_MH$sum_num_changes
              union_DAG[p_reg, l] = as.numeric(SUM_NUM_CHANGES != 0)
              
              THETA_MAT_MC[,l,p_reg,iter_num + 1] = INT_TERM_MC[,l,p_reg,iter_num + 1] + 
                LIN_TERM_MC[,l,p_reg,iter_num ] + 
                NONLIN_TERM_MC[,l,p_reg,iter_num ]
              
              BETA_MAT_MC[,l,p_reg,iter_num + 1] = THETA_MAT_MC[,l,p_reg,iter_num + 1]*
                as.numeric(abs(THETA_MAT_MC[,l,p_reg,iter_num + 1]) >  THRESHOLD_MC[p_reg,l,iter_num + 1])
              
              X_DOT_BETA_MC[,l,p_reg] = rawX[,l,p_reg]*BETA_MAT_MC[,l,p_reg,iter_num + 1]
            }
          }
          
          #####################
          ## Update eta*_jl ###
          #####################
          for(p_reg in 1:num_of_p)
          {
            for(l in setdiff(1:num_of_p, p_reg))
            {
              t_ = THRESHOLD_MC[p_reg,l,iter_num + 1]
              gamma_ = GAMMA_MC[,l, p_reg,iter_num]
              tau_ = TAU_MC[p_reg,l,iter_num]
              lambda_ = LAMBDA_MC[,l,p_reg,iter_num]
              mix_nor_gamma_  = MIX_NOR_GAMMA_MC[,l,p_reg,iter_num]
              modified_gamma_mat = MODIFIED_GAMMA_MAT_MC[,l,p_reg,iter_num]
              llh_old = LLH_MC[p_reg,l, iter_num + 1]
              NZV_FII_LLH = apply(X_DOT_BETA_MC[,setdiff(1:num_of_p,l),p_reg], 1, sum)
              
              ###############################################################
              num_beta_nnz = sum(BETA_MAT_MC[,l,p_reg,iter_num + 1] !=0)
              ###############################################################
              
              gamma_MH<-Update_nonLin_gamma_WDC(t = t_, 
                                                gamma_mat =  gamma_, 
                                                nonLinearCov_mat = U_nonLin, node_mat = rawX[,l,p_reg], response_vec =  Y[,p_reg],
                                                sigma_gamma = eta_gamma_std, tau_vec = tau_, lambda_mat =lambda_,
                                                tau_quantile = tau_index,
                                                mix_normal_nonlin = mix_nor_gamma_,
                                                intTerm = INT_TERM_MC[,l,p_reg,iter_num + 1],
                                                linTerm = LIN_TERM_MC[,l,p_reg,iter_num],
                                                modified_gamma_mat = modified_gamma_mat,
                                                llh_old = llh_old, indicator_th = NZV_FII_LLH,
                                                DAG_all_patients_this_iter = union_DAG,
                                                num_beta_nnz)
              
              GAMMA_MC[,l,p_reg,iter_num + 1] <- gamma_MH$gamma_mat_out
              MODIFIED_GAMMA_MAT_MC[,l,p_reg,iter_num + 1] = gamma_MH$modified_gamma_mat_out
              LLH_MC[p_reg,l,iter_num + 1] = gamma_MH$llh_out
              NONLIN_TERM_MC[,l,p_reg,iter_num + 1] = llh_theta_nonlin(U_nonLin, 
                                                                       MODIFIED_GAMMA_MAT_MC[,l,p_reg,iter_num + 1], 
                                                                       mix_nor_gamma_)
              ACCEPTANCE_MC[p_reg,l,3] = ACCEPTANCE_MC[p_reg,l,3] + gamma_MH$accept
              
              SUM_NUM_CHANGES = gamma_MH$sum_num_changes
              union_DAG[p_reg, l] = as.numeric(SUM_NUM_CHANGES != 0)
              
              THETA_MAT_MC[,l,p_reg,iter_num + 1] = INT_TERM_MC[,l,p_reg,iter_num + 1] + 
                LIN_TERM_MC[,l,p_reg,iter_num ] + 
                NONLIN_TERM_MC[,l,p_reg,iter_num + 1 ]
              
              BETA_MAT_MC[,l,p_reg,iter_num + 1] = THETA_MAT_MC[,l,p_reg,iter_num + 1]*
                as.numeric(abs(THETA_MAT_MC[,l,p_reg,iter_num + 1]) >  THRESHOLD_MC[p_reg,l,iter_num + 1])
              
              X_DOT_BETA_MC[,l,p_reg] = rawX[,l,p_reg]*BETA_MAT_MC[,l,p_reg,iter_num + 1]
            }
          }
          #####################
          ## Update T_jl ######
          #####################
          for(p_reg in 1:num_of_p)
          {
            for(l in setdiff(1:num_of_p, p_reg))
            {
              xi_ = XI_MC[p_reg,l,iter_num]
              lambda_ = LAMBDA_MC[,l,p_reg,iter_num]
              gamma_ = GAMMA_MC[,l,p_reg,iter_num + 1]
              tau_ = TAU_MC[p_reg,l,iter_num]
              
              tau_Gibbs<-sqrt(Update_nonLin_tau2(xi_vec = xi_, lambda_mat = lambda_, gamma_mat = gamma_, old_tau2 = tau_^2))
              
              TAU_MC[p_reg,l,iter_num + 1]<-tau_Gibbs
            }
          }
          #####################
          ## Update c_jl ######
          #####################
          for(p_reg in 1:num_of_p)
          {
            for(l in setdiff(1:num_of_p,p_reg))
            {
              xi_Gibbs<-Update_nonLin_xi(tau = TAU_MC[p_reg,l,iter_num + 1])
              XI_MC[p_reg,l,iter_num + 1]<-xi_Gibbs 
            }
          }
          #####################
          ## Update L_jl ######
          #####################
          for(p_reg in 1:num_of_p)
          {
            for(l in setdiff(1:num_of_p,p_reg))
            {
              tau_aug<-rep(TAU_MC[p_reg,l,iter_num + 1],each=q)
              lambda_Gibbs<-sqrt(Update_nonLin_lambda2(zeta = ZETA_MC[,l,p_reg,iter_num], 
                                                       tau_augment = tau_aug, gamma_mat= GAMMA_MC[,l,p_reg,iter_num + 1]))
              
              LAMBDA_MC[,l,p_reg,iter_num + 1]<-lambda_Gibbs 
            }
          }
          ########################
          ## Update Zeta_jl ######
          ########################
          for(p_reg in 1:num_of_p)
          {
            for(l in setdiff(1:num_of_p,p_reg))
            {
              zeta_Gibbs<-Update_nonLin_zeta(lambda_mat = LAMBDA_MC[,l,p_reg,iter_num + 1])
              
              ZETA_MC[,l,p_reg,iter_num + 1]<-zeta_Gibbs
            }
          }
          #####################
          ## Update eta_jl ####
          #####################
          for(p_reg in 1:num_of_p)
          {
            for(l in setdiff(1:num_of_p,p_reg))
            {
              t_ = THRESHOLD_MC[p_reg,l,iter_num + 1]
              eta_ = ETA_MC[,l,p_reg,iter_num]
              alpha_ = ALPHA_MC[p_reg,l,iter_num]
              kappa_ = KAPPA_MC[,l,p_reg,iter_num]
              mix_nor_eta_ = MIX_NOR_ETA_MC[,l,p_reg,iter_num]
              NZV_FII_LLH = apply(X_DOT_BETA_MC[,setdiff(1:num_of_p,l),p_reg], 1, sum)
              
              ###############################################################
              num_beta_nnz = sum(BETA_MAT_MC[,l,p_reg,iter_num + 1] !=0)
              ###############################################################
              
              eta_MH<-Update_lin_eta_WDC(t = t_ ,
                                         eta_mat = eta_ ,linearCov_mat = U_lin,node_mat = rawX[,l,p_reg],response_vec = Y[,p_reg],
                                         sigma_eta = eta_gamma_std,alpha_vec = alpha_,kappa_mat = kappa_, 
                                         tau_quantile = tau_index,
                                         mix_normal_lin = mix_nor_eta_,
                                         intTerm = INT_TERM_MC[,l,p_reg,iter_num + 1],
                                         nonlinTerm = NONLIN_TERM_MC[,l,p_reg,iter_num + 1],
                                         llh_old = LLH_MC[p_reg,l,iter_num + 1], indicator_th = NZV_FII_LLH,
                                         DAG_all_patients_this_iter = union_DAG,
                                         num_beta_nnz)
              
              ETA_MC[,l,p_reg,iter_num + 1] <-eta_MH$eta_mat_out
              LIN_TERM_MC[,l,p_reg, iter_num + 1] = llh_theta_lin(ETA_MC[,l,p_reg,iter_num + 1] ,U_lin , mix_nor_eta_)
              LLH_MC[p_reg,l,iter_num + 1] = eta_MH$llh_out
              ACCEPTANCE_MC[p_reg,l,2]<- ACCEPTANCE_MC[p_reg,l,2] + eta_MH$accept
              
              SUM_NUM_CHANGES = eta_MH$sum_num_changes
              union_DAG[p_reg, l] = as.numeric(SUM_NUM_CHANGES != 0)
              
              THETA_MAT_MC[,l,p_reg,iter_num + 1] = INT_TERM_MC[,l,p_reg,iter_num + 1] + 
                LIN_TERM_MC[,l,p_reg,iter_num + 1 ] + 
                NONLIN_TERM_MC[,l,p_reg,iter_num + 1 ]
              
              BETA_MAT_MC[,l,p_reg,iter_num + 1] = THETA_MAT_MC[,l,p_reg,iter_num + 1]*
                as.numeric(abs(THETA_MAT_MC[,l,p_reg,iter_num + 1]) >  THRESHOLD_MC[p_reg,l,iter_num + 1])
              
              X_DOT_BETA_MC[,l,p_reg] = rawX[,l,p_reg]*BETA_MAT_MC[,l,p_reg,iter_num + 1]
            }
          }
          
          #####################
          ## Update A_jl ######
          #####################
          for(p_reg in 1:num_of_p)
          {
            for(l in setdiff(1:num_of_p,p_reg))
            {
              rho_ = RHO_MC[p_reg,l,iter_num]
              kappa_ = KAPPA_MC[,l,p_reg,iter_num]
              eta_ = ETA_MC[,l,p_reg,iter_num + 1]
              alpha_ = ALPHA_MC[p_reg,l,iter_num]
              
              alpha_Gibbs<-sqrt(Update_lin_alpha2(rho_vec = rho_, kappa_mat = kappa_, eta_mat = eta_, old_alpha2 = alpha_^2))
              
              ALPHA_MC[p_reg,l,iter_num + 1]<-alpha_Gibbs 
            }
          }
          #######################
          ## Update rho_jl ######
          #######################
          for(p_reg in 1:num_of_p)
          {
            for(l in setdiff(1:num_of_p,p_reg))
            {
              alpha_ = ALPHA_MC[p_reg,l,iter_num + 1]
              rho_Gibbs<-Update_lin_rho(alpha_)
              
              RHO_MC[p_reg,l, iter_num + 1]<-rho_Gibbs 
            }
          }
          #####################
          ## Update K_jl ######
          #####################
          for(p_reg in 1:num_of_p)
          {
            for(l in setdiff(1:num_of_p,p_reg))
            {
              alpha_aug<-rep(ALPHA_MC[p_reg,l,iter_num + 1],each=q)
              eta_ = ETA_MC[,l,p_reg,iter_num + 1]
              phi_ = PHI_MC[,l,p_reg,iter_num]
              
              kappa_Gibbs<-sqrt(Update_lin_kappa2(phi_mat = phi_, alpha_augment = alpha_aug, eta_mat = eta_))
              
              KAPPA_MC[,l,p_reg,iter_num + 1]<-kappa_Gibbs 
            }
          }
          #######################
          ## Update phi_jl ######
          #######################
          for(p_reg in 1:num_of_p)
          {
            for(l in setdiff(1:num_of_p,p_reg))
            {
              phi_Gibbs<-Update_lin_phi(KAPPA_MC[,l,p_reg,iter_num + 1])
              PHI_MC[,l,p_reg,iter_num + 1]<-phi_Gibbs 
            }
          }
          ################################
          ## Update xi_jl and xi*_jl######
          ################################
          for(p_reg in 1:num_of_p)
          {
            for(l in setdiff(1:num_of_p,p_reg))
            {
              t_ = THRESHOLD_MC[p_reg,l,iter_num + 1]
              eta_ = ETA_MC[,l,p_reg,iter_num + 1]
              mix_nor_eta_ = MIX_NOR_ETA_MC[,l,p_reg,iter_num]
              mix_nor_gamma_ = MIX_NOR_GAMMA_MC[,l,p_reg,iter_num]
              intTerm = INT_TERM_MC[,l,p_reg,iter_num + 1]
              modified_gamma_mat = MODIFIED_GAMMA_MAT_MC[,l,p_reg,iter_num + 1]
              llh_old = LLH_MC[p_reg,l,iter_num  + 1]
              NZV_FII_LLH = apply(X_DOT_BETA_MC[,setdiff(1:num_of_p,l),p_reg], 1, sum)
              
              ###############################################################
              num_beta_nnz = sum(BETA_MAT_MC[,l,p_reg,iter_num + 1] !=0)
              ###############################################################
              
              mix_normal_updates = mix_normal_MH_WDC(t = t_ ,
                                                     eta_mat = eta_ ,
                                                     linearCov_mat = U_lin,
                                                     nonLinearCov_mat = U_nonLin,
                                                     node_mat = rawX[,l,p_reg],response_vec = Y[,p_reg],
                                                     tau_quantile = tau_index, 
                                                     mix_normal_lin = mix_nor_eta_, 
                                                     mix_normal_nonlin = mix_nor_gamma_,
                                                     mix_nor_std = mix_nor_std,
                                                     intTerm = intTerm,
                                                     modified_gamma_mat = modified_gamma_mat,
                                                     llh_old = llh_old, indicator_th = NZV_FII_LLH,
                                                     DAG_all_patients_this_iter = union_DAG,
                                                     num_beta_nnz)
              
              MIX_NOR_ETA_MC[,l,p_reg,iter_num + 1] = mix_normal_updates$lin_out
              MIX_NOR_GAMMA_MC[,l,p_reg,iter_num + 1] = mix_normal_updates$nonlin_out
              LLH_MC[p_reg,l,iter_num  + 1] = mix_normal_updates$llh_out
              LIN_TERM_MC[,l,p_reg,iter_num + 1] = llh_theta_lin(eta_,U_lin , MIX_NOR_ETA_MC[,l,p_reg,iter_num + 1])
              NONLIN_TERM_MC[,l,p_reg,iter_num + 1]  = llh_theta_nonlin(U_nonLin, modified_gamma_mat,
                                                                        MIX_NOR_GAMMA_MC[,l,p_reg,iter_num + 1])
              
              ACCEPTANCE_MC[p_reg,l,4]<- ACCEPTANCE_MC[p_reg,l,4]+ mix_normal_updates$accept_out
              
              SUM_NUM_CHANGES = mix_normal_updates$sum_num_changes
              union_DAG[p_reg, l] = as.numeric(SUM_NUM_CHANGES != 0)
              
              THETA_MAT_MC[,l,p_reg,iter_num + 1] = INT_TERM_MC[,l,p_reg,iter_num + 1] + 
                LIN_TERM_MC[,l,p_reg,iter_num + 1 ] + 
                NONLIN_TERM_MC[,l,p_reg,iter_num + 1 ]
              
              BETA_MAT_MC[,l,p_reg,iter_num + 1] = THETA_MAT_MC[,l,p_reg,iter_num + 1]*
                as.numeric(abs(THETA_MAT_MC[,l,p_reg,iter_num + 1]) >  THRESHOLD_MC[p_reg,l,iter_num + 1])
              
              X_DOT_BETA_MC[,l,p_reg] = rawX[,l,p_reg]*BETA_MAT_MC[,l,p_reg,iter_num + 1]
            }
          }
          ######################################
          for(p_reg in 1:num_of_p)
          {
            for(l in 1:num_of_p)
            {
              tau_aug<-rep(TAU_MC[p_reg,l,iter_num + 1],each=q)
              lambda_ = LAMBDA_MC[,l,p_reg,iter_num + 1]
              
              tmp_mat = tau_aug*lambda_
              nlinrate = 1-1/(1+tmp_mat^2)
              
              alpha_aug<-rep(ALPHA_MC[p_reg,l,iter_num + 1],each=q)
              kappa_ = KAPPA_MC[,l,p_reg,iter_num + 1]
              
              tmp_mat = alpha_aug*kappa_
              linrate = 1-1/(1+tmp_mat^2)
              
              lnrate<-pmax(linrate,nlinrate)
              v_Q_MC[,l,p_reg,iter_num + 1] = lnrate
            }
          }
          ######################################
        }
        
      })
      
      ############################# Posterior Inference #####################################
      tau_index = tau_index*10
      thinning = 10
      thin_iter_seq = seq(burnIn+thinning,total_MC_iters,thinning)
      post_size = (total_MC_iters - burnIn)/thinning
      #######################################################################################
      
      POST_MEAN_MU = apply(MU_MC[,,thin_iter_seq],c(1,2),mean)
      POST_MEAN_ETA = apply(ETA_MC[,,,thin_iter_seq], c(1,2,3), mean)
      POST_MEAN_GAMMA = apply(GAMMA_MC[,,,thin_iter_seq], c(1,2,3), mean)
      POST_MEAN_MIX_NORM_ETA = apply(MIX_NOR_ETA_MC[,,,thin_iter_seq], c(1,2,3), mean)
      POST_MEAN_MIX_NORM_GAMMA = apply(MIX_NOR_GAMMA_MC[,,,thin_iter_seq], c(1,2,3), mean)
      POST_MEAN_THRESHOLD = apply(THRESHOLD_MC[,,thin_iter_seq], c(1,2),mean)
      
      POST_THETA_EST = array(0,c(n,p,p))
      POST_BETA_EST = array(0,c(n,p,p))
      
      for(p_reg in 1:num_of_p)
      {
        hs_mu_vec = POST_MEAN_MU[p_reg,]
        hs_eta_mat = POST_MEAN_ETA[,,p_reg]
        hs_gamma_mat = POST_MEAN_GAMMA[,,p_reg]
        hs_mix_norm_eta_mat = POST_MEAN_MIX_NORM_ETA[,,p_reg]
        hs_mix_norm_gamma_mat = POST_MEAN_MIX_NORM_GAMMA[,,p_reg]
        hs_threshold_est = POST_MEAN_THRESHOLD[p_reg,]
        
        modified_gamma_mat = matrix(0, nrow = sum(b_star), ncol = p)
        for(mod_gamma_col in 1:p)
        {
          modified_gamma_mat[,mod_gamma_col] = rep(hs_gamma_mat[,mod_gamma_col], b_star) 
        }
        
        hs_theta_est = llh_theta_int(hs_mu_vec) + 
          llh_theta_lin(hs_eta_mat,U_lin,hs_mix_norm_eta_mat) + 
          llh_theta_nonlin(U_nonLin,modified_gamma_mat ,hs_mix_norm_gamma_mat)
        
        POST_THETA_EST[,,p_reg] = hs_theta_est
        
        for(p_reg_col in 1:p)
        {
          POST_BETA_EST[,p_reg_col,p_reg] = POST_THETA_EST[,p_reg_col,p_reg]*as.numeric(abs(POST_THETA_EST[,p_reg_col,p_reg])>hs_threshold_est[p_reg_col])
        }
      }
      ######################################################################
      TRUE_BETA_MAT = array(0,c(n,p,p-1))
      TRUE_THETA_MAT = array(0,c(n,p,p-1))
      
      for(p_reg in 1:(num_of_p - 1))
      {
        file_name = paste0("data_with_th_",true_th,"_n_",sample_size,"_p_",num_of_p,"_q_",num_of_q,"/True_betas/True_beta_matrix_n",sample_size,"_p_",num_of_p,"_q_",num_of_q,"_p_reg_",p_reg,"_tau_",tau_index,"_simu_",simulation_index,".csv")
        True_beta_matrix = read.csv(file_name, header = FALSE)
        True_beta_matrix = cbind(matrix(0,sample_size,p_reg), True_beta_matrix)
        
        TRUE_BETA_MAT[,,p_reg] = as.matrix(True_beta_matrix)
        #########################################
        file_name = paste0("data_with_th_",true_th,"_n_",sample_size,"_p_",num_of_p,"_q_",num_of_q,"/True_thetas/True_theta_matrix_n",sample_size,"_p_",num_of_p,"_q_",num_of_q,"_p_reg_",p_reg,"_tau_",tau_index,"_simu_",simulation_index,".csv")
        True_theta_matrix = read.csv(file_name, header = FALSE)
        True_theta_matrix = cbind(matrix(0,sample_size,p_reg), True_theta_matrix)
        
        TRUE_THETA_MAT[,,p_reg] = as.matrix(True_theta_matrix)
      }
      ######################################################################
      TRUE_Y = matrix(0,sample_size,num_of_p - 1)
      
      for(p_reg in 1:(num_of_p - 1))
      {
        filename_output = paste0("data_with_th_",true_th,"_n_",sample_size,"_p_",num_of_p,"_q_",num_of_q,"/QR_Y_sim_TRUEs_n_",sample_size,
                                 "_seed_",simulation_index,'.csv');
        temp = read.csv(filename_output, header = FALSE);
        temp = temp[,((p_reg-1)*9+1):(p_reg*9)];
        temp = as.matrix(temp)
        
        TRUE_Y[,p_reg] = temp[,tau_index]
      }
      
      EST_Y = matrix(0,nrow = sample_size, ncol = num_of_p -1)
      MSE= rep(0,num_of_p - 1)
      
      for(p_reg in 1:(num_of_p - 1))
      {
        EST_Y[,p_reg] = apply(rawX[,,p_reg]*POST_BETA_EST[,,p_reg],1,sum)
        MSE[p_reg] =  sum((TRUE_Y[,p_reg] - EST_Y[,p_reg])^2)/sample_size
        
        MSE[p_reg] = MSE[p_reg]/(MSE_Scales[p_reg]+1)
      }
      ######################################################################
      MEAN_TBS_MSE = mean(MSE)
      ######################################################################
      MEAN_TBS_THETA_NORM = 0.5*sqrt(sum((TRUE_BETA_MAT - POST_BETA_EST[,,1:(num_of_p-1)])*(TRUE_BETA_MAT - POST_BETA_EST[,,1:(num_of_p-1)])))
      MEAN_TBS_BETA_NORM = 0.5*sqrt(sum((TRUE_THETA_MAT - POST_THETA_EST[,,1:(num_of_p-1)])*(TRUE_THETA_MAT - POST_THETA_EST[,,1:(num_of_p-1)])))
      
      ######################################################################
      #### Processing again for TPR and FPR calculations ###################
      
      TRUE_BETA_MAT = array(0,c(n,p-1,p-1))
      TRUE_THETA_MAT = array(0,c(n,p-1,p-1))
      
      ########################################################################
      for(p_reg in 1:(num_of_p - 1))
      {
        file_name = paste0("data_with_th_",true_th,"_n_",sample_size,"_p_",num_of_p,"_q_",num_of_q,"/True_betas/True_beta_matrix_n",sample_size,"_p_",num_of_p,"_q_",num_of_q,"_p_reg_",p_reg,"_tau_",tau_index,"_simu_",simulation_index,".csv")
        True_beta_matrix = read.csv(file_name, header = FALSE)
        if(p_reg != 1)
        {
          True_beta_matrix = cbind(matrix(0,sample_size,p_reg -1), True_beta_matrix)
        }
        TRUE_BETA_MAT[,,p_reg] = as.matrix(True_beta_matrix)
        #########################################
        file_name = paste0("data_with_th_",true_th,"_n_",sample_size,"_p_",num_of_p,"_q_",num_of_q,"/True_thetas/True_theta_matrix_n",sample_size,"_p_",num_of_p,"_q_",num_of_q,"_p_reg_",p_reg,"_tau_",tau_index,"_simu_",simulation_index,".csv")
        True_theta_matrix = read.csv(file_name, header = FALSE)
        if(p_reg != 1)
        {
          True_theta_matrix = cbind(matrix(0,sample_size,p_reg -1), True_theta_matrix)
        }
        TRUE_THETA_MAT[,,p_reg] = as.matrix(True_theta_matrix)
      }
      ##################################################################
      TPR_FPR_AUC_P_ALL_PREGS = matrix(0,3, num_of_p - 1)
      BETA_MAT_MC_RED = BETA_MAT_MC[,,,seq(burnIn+thinning,total_MC_iters,thinning)]
      THETA_MAT_MC_RED = THETA_MAT_MC[,,,seq(burnIn+thinning,total_MC_iters,thinning)]
      BOOL_TP_BETA_MAT_MC_RED = array(FALSE, c(sample_size, num_of_p - 1, num_of_p - 1, post_size))
      BOOL_FP_BETA_MAT_MC_RED = array(FALSE, c(sample_size, num_of_p - 1, num_of_p - 1, post_size))
      
      for(p_reg in 1:(num_of_p - 1))
      {
        T_selec_3k_iters_1 = array(0, c(sample_size, num_of_p - 1, post_size))
        T_selec_3k_iters_2 = array(0, c(sample_size, num_of_p - 1, post_size))
        T_selec_no_TP_FP =   array(0, c(sample_size, num_of_p - 1, post_size))
        
        for (iter in 1:post_size)
        {
          temp_beta = BETA_MAT_MC_RED[,,p_reg,iter]
          temp_beta = temp_beta[,setdiff(1:num_of_p,p_reg)]
          T_selec_3k_iters_1[,,iter] = (temp_beta !=0) &  (TRUE_BETA_MAT[,,p_reg]!=0)
          BOOL_TP_BETA_MAT_MC_RED[,,p_reg,iter] = (temp_beta !=0) &  (TRUE_BETA_MAT[,,p_reg]!=0)
          
          T_selec_3k_iters_2[,,iter] = (temp_beta !=0) &  (TRUE_BETA_MAT[,,p_reg]==0)
          BOOL_FP_BETA_MAT_MC_RED[,,p_reg,iter] = (temp_beta !=0) &  (TRUE_BETA_MAT[,,p_reg]==0)
          
          T_selec_no_TP_FP[,,iter]   = (temp_beta!=0)
        }
        
        T_post_prob_selec_1 = apply(T_selec_3k_iters_1, c(1,2), mean)
        T_post_prob_selec_2 = apply(T_selec_3k_iters_2, c(1,2), mean)
        T_post_prob_selec_no_TP_FP = apply(T_selec_no_TP_FP, c(1,2), mean)
        
        ################### FDR addition ########################################
        th_seq = seq(0.01,1,0.01)
        TP_for_var_th = rep(0, length(th_seq))
        FP_for_var_th = rep(0, length(th_seq))
        
        for(iter in 1:length(th_seq))
        {
          TP_for_var_th[iter] = sum((T_post_prob_selec_1 > th_seq[iter]) & (TRUE_BETA_MAT[,,p_reg]!=0))
          FP_for_var_th[iter] = sum((T_post_prob_selec_2 > th_seq[iter]) & (TRUE_BETA_MAT[,,p_reg]==0))
        }
        
        FDR = FP_for_var_th/(FP_for_var_th + TP_for_var_th)
        FDR = FDR[!is.na(FDR)]
        diff_with_req_th = abs(FDR - th_for_FDR)
        
        desired_th_for_P = c()
        tryCatch(
          {
            min_indices = which(diff_with_req_th == min(diff_with_req_th))
            desired_th_for_P = th_seq[max(min_indices)]
          },
          error = function(e){
          },
          warning = function(w){
          }
        )
        if(is.null(desired_th_for_P))
        {
          desired_th_for_P = 0.5
        }
        TPR_P = sum((T_post_prob_selec_1 > desired_th_for_P) & (TRUE_BETA_MAT[,,p_reg]!=0))/sum(TRUE_BETA_MAT[,,p_reg]!=0)
        FPR_P = sum((T_post_prob_selec_2 > desired_th_for_P) & (TRUE_BETA_MAT[,,p_reg]==0))/sum(TRUE_BETA_MAT[,,p_reg]==0)
        #########################################################################
        threshold_seq = seq(0,1,1/1000)
        tpr = rep(0,length(threshold_seq))
        fpr = rep(0,length(threshold_seq))
        
        for(iter in 1:length(threshold_seq))
        {
          tpr[iter] = sum(sum((T_post_prob_selec_no_TP_FP>=threshold_seq[iter]) & (TRUE_BETA_MAT[,,p_reg]!=0)))/sum(sum(TRUE_BETA_MAT[,,p_reg] !=0));
          fpr[iter] = sum(sum((T_post_prob_selec_no_TP_FP>=threshold_seq[iter]) & (TRUE_BETA_MAT[,,p_reg]==0)))/sum(sum(TRUE_BETA_MAT[,,p_reg] ==0));
        }
        auc = 0;
        for (iter in 1:(length(threshold_seq)-1)){
          auc = auc + tpr[iter]*(fpr[iter]-fpr[iter+1]);
        }
        if(auc>1) {AUC_P = 1}
        AUC_P = auc
        
        TPR_FPR_AUC_P_ALL_PREGS[1,p_reg] = TPR_P
        TPR_FPR_AUC_P_ALL_PREGS[2,p_reg] = FPR_P
        TPR_FPR_AUC_P_ALL_PREGS[3,p_reg] = AUC_P
      }
      ######################################################################
      MEAN_TBS_TPR_FPR_AUC_P_ALL_PREGS = apply(TPR_FPR_AUC_P_ALL_PREGS,1,mean)
      ######################################################################
      True_T_List = readMat(paste0("data_with_th_",true_th,"_n_",sample_size,"_p_",num_of_p,"_q_",num_of_q,"/Y_order_n_",sample_size,"_p_",num_of_p,"_q_",num_of_q,"_v6.mat"))
      True_G_List = readMat(paste0("data_with_th_",true_th,"_n_",sample_size,"_p_",num_of_p,"_q_",num_of_q,"/X_order_n_",sample_size,"_p_",num_of_p,"_q_",num_of_q,"_v6.mat"))
      
      adjusted_true_cords = list(NULL)
      adjusted_null_cords = list(NULL)
      for(P_REG in 1:(num_of_p-1))
      {
        adjusted_true_cords[[P_REG]] = as.vector(True_T_List$which.Ys.TP[[P_REG]][[1]]) - P_REG +1
        adjusted_null_cords[[P_REG]] = setdiff(2:(num_of_p - P_REG+1), adjusted_true_cords[[P_REG]])
      }
      
      BOOLEAN_Q_ALL_PREGS = array(FALSE, c(q,p,num_of_p - 1))
      for(p_reg in 1:(num_of_p- 1))
      {
        temp_boolean_matrix = matrix(FALSE, nrow = num_of_q, ncol = num_of_p -p_reg)
        
        TP_ts_this_reg_mod  = adjusted_true_cords[[p_reg]] - 1
        
        for(bool_col in 1:length(TP_ts_this_reg_mod))
        {
          temp_G_array = True_G_List$which.Xs.TP[[p_reg]][[1]][[bool_col]][[1]]
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
        
        BOOLEAN_Q_ALL_PREGS[,(p_reg+1):num_of_p, p_reg] = temp_boolean_matrix
      }
      
      TPR_FPR_AUC_Q_ALL_PREGS = matrix(0,3, num_of_p - 1)
      v_Q_RED = v_Q_MC[,,,seq(burnIn+thinning,total_MC_iters,thinning)]
      
      for(p_reg in 1:(num_of_p - 1))
      {
        ######################## FDR Addition ####################################
        T_selec_3k_iters_1 = array(0, c(num_of_q, num_of_p , post_size))
        T_selec_3k_iters_2 = array(0, c(num_of_q, num_of_p , post_size))
        T_selec_no_TP_FP =   array(0, c(num_of_q, num_of_p , post_size))
        th_seq = seq(0.01,1,0.01)
        TP_for_var_th = rep(0, length(th_seq))
        FP_for_var_th = rep(0, length(th_seq))
        FDR = rep(0, length(th_seq))
        
        for(th_iter in 1:length(th_seq))
        {
          for (iter in 1:post_size)
          {
            temp_v_q = v_Q_RED[,,p_reg,iter]
            T_selec_3k_iters_1[,,iter] = (temp_v_q >th_seq[th_iter]) &  (BOOLEAN_Q_ALL_PREGS[,,p_reg]!=0)
            T_selec_3k_iters_2[,,iter] = (temp_v_q >th_seq[th_iter]) &  (BOOLEAN_Q_ALL_PREGS[,,p_reg]==0)
            T_selec_no_TP_FP[,,iter]   = (temp_v_q >th_seq[th_iter])
          }
          
          T_post_prob_selec_1 = apply(T_selec_3k_iters_1, c(1,2), mean)
          T_post_prob_selec_2 = apply(T_selec_3k_iters_2, c(1,2), mean)
          T_post_prob_selec_no_TP_FP = apply(T_selec_no_TP_FP, c(1,2), mean)
          
          TP_for_var_th = sum((T_post_prob_selec_1 > th_seq[th_iter]) & (BOOLEAN_Q_ALL_PREGS[,,p_reg]!=0))
          FP_for_var_th = sum((T_post_prob_selec_2 > th_seq[th_iter]) & (BOOLEAN_Q_ALL_PREGS[,,p_reg]==0))
          
          FDR[th_iter] = FP_for_var_th/(FP_for_var_th + TP_for_var_th)
          
        }
        FDR = FDR[!is.na(FDR)]
        diff_with_req_th = abs(FDR - th_for_FDR)
        
        desired_th_for_G = c()
        tryCatch(
          {
            min_indices = which(diff_with_req_th == min(diff_with_req_th))
            desired_th_for_G = th_seq[min(min_indices)]
          },
          error = function(e){
          },
          warning = function(w){
          }
        )
        if(is.null(desired_th_for_G))
        {
          desired_th_for_G = 0.5
        }
        
        T_selec_3k_iters_1 = array(0, c(num_of_q, num_of_p , post_size))
        T_selec_3k_iters_2 = array(0, c(num_of_q, num_of_p , post_size))
        T_selec_no_TP_FP =   array(0, c(num_of_q, num_of_p , post_size))
        
        for (iter in 1:post_size)
        {
          temp_v_q = v_Q_RED[,,p_reg,iter]
          T_selec_3k_iters_1[,,iter] = (temp_v_q >desired_th_for_G) &  (BOOLEAN_Q_ALL_PREGS[,,p_reg]!=0)
          T_selec_3k_iters_2[,,iter] = (temp_v_q >desired_th_for_G) &  (BOOLEAN_Q_ALL_PREGS[,,p_reg]==0)
          T_selec_no_TP_FP[,,iter]   = (temp_v_q >desired_th_for_G)
        }
        
        T_post_prob_selec_1 = apply(T_selec_3k_iters_1, c(1,2), mean)
        T_post_prob_selec_2 = apply(T_selec_3k_iters_2, c(1,2), mean)
        T_post_prob_selec_no_TP_FP = apply(T_selec_no_TP_FP, c(1,2), mean)
        
        TPR_Q = sum((T_post_prob_selec_1 > desired_th_for_G) & (BOOLEAN_Q_ALL_PREGS[,,p_reg]!=0))/sum(BOOLEAN_Q_ALL_PREGS[,,p_reg]!=0)
        FPR_Q = sum((T_post_prob_selec_2 > desired_th_for_G) & (BOOLEAN_Q_ALL_PREGS[,,p_reg]==0))/sum(BOOLEAN_Q_ALL_PREGS[,,p_reg]==0)
        
        ##########################################################################
        threshold_seq = seq(0,1,1/1000)
        tpr = rep(0,length(threshold_seq))
        fpr = rep(0,length(threshold_seq))
        
        for(iter in 1:length(threshold_seq))
        {
          tpr[iter] = sum(sum((T_post_prob_selec_no_TP_FP>=threshold_seq[iter]) & (BOOLEAN_Q_ALL_PREGS[,,p_reg]!=0)))/sum(sum(BOOLEAN_Q_ALL_PREGS[,,p_reg] !=0));
          fpr[iter] = sum(sum((T_post_prob_selec_no_TP_FP>=threshold_seq[iter]) & (BOOLEAN_Q_ALL_PREGS[,,p_reg]==0)))/sum(sum(BOOLEAN_Q_ALL_PREGS[,,p_reg] ==0));
        }
        auc = 0;
        for (iter in 1:(length(threshold_seq)-1)){
          auc = auc + tpr[iter]*(fpr[iter]-fpr[iter+1]);
        }
        if(auc>1) {AUC_Q = 1}
        AUC_Q = auc
        
        TPR_FPR_AUC_Q_ALL_PREGS[1,p_reg] = TPR_Q
        TPR_FPR_AUC_Q_ALL_PREGS[2,p_reg] = FPR_Q
        TPR_FPR_AUC_Q_ALL_PREGS[3,p_reg] = AUC_Q
      }
      ######################################################################
      MEAN_TBS_TPR_FPR_AUC_Q_ALL_PREGS = apply(TPR_FPR_AUC_Q_ALL_PREGS,1,mean)
      ######################################################################
      THRESHOLD_MC_RED = THRESHOLD_MC[,,seq(burnIn+thinning,total_MC_iters,thinning)]
      for(p_reg in 1:(num_of_p))
      {
        THRESHOLD_MC_RED[p_reg,p_reg,] = 0
      }
      THRESHOLD_MC_RED = apply(THRESHOLD_MC_RED, c(1,2), mean)
      for(p_reg in 1:num_of_p)
      {
        THRESHOLD_MC_RED[p_reg,] = c(THRESHOLD_MC_RED[p_reg,c(THRESHOLD_MC_RED[p_reg,]!=0)],0)
      }
      THRESHOLD_MC_RED = THRESHOLD_MC_RED[1:(num_of_p-1),1:(num_of_p-1)]
      ######################################################################
      MEAN_TBS_THRESHOLD = mean(apply(THRESHOLD_MC_RED,1,mean))
      ######################################################################
      
      summary = c(MEAN_TBS_TPR_FPR_AUC_P_ALL_PREGS,MEAN_TBS_TPR_FPR_AUC_Q_ALL_PREGS,
                  MEAN_TBS_MSE, MEAN_TBS_THRESHOLD, MEAN_TBS_THETA_NORM, MEAN_TBS_BETA_NORM, as.numeric(time_rec[3])/60)
      
      file_name = paste0(Results_directory_name,"/Post_process_MC_",
                         "_simulation_index_",simulation_index,
                         "_tau_index_",tau_index,".csv")
      
      write.table(summary,file_name,sep = ",",row.names = FALSE,col.names = FALSE )
      
      ######################################################################
    }
}

