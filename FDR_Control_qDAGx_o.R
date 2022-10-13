if(sum(True_beta_matrix == 0) !=0)
{
  T_selec_3k_iters_1 = array(0, c(sample_size, num_of_p - p_regression, post_size))
  T_selec_3k_iters_2 = array(0, c(sample_size, num_of_p - p_regression, post_size))
  
  for (iter in 1:post_size)
  {
    T_selec_3k_iters_1[,,iter] = (v[,2:(num_of_p - p_regression + 1),iter] !=0) &  (True_beta_matrix !=0)
    T_selec_3k_iters_2[,,iter] = (v[,2:(num_of_p - p_regression + 1),iter] !=0) &  (True_beta_matrix ==0)
  }
  
  T_post_prob_selec_1 = apply(T_selec_3k_iters_1, c(1,2), mean)
  T_post_prob_selec_2 = apply(T_selec_3k_iters_2, c(1,2), mean)
  
  th_seq = seq(0.01,1,0.01)
  
  TP_for_var_th = rep(0, length(th_seq))
  FP_for_var_th = rep(0, length(th_seq))
  
  for(iter in 1:length(th_seq))
  {
    TP_for_var_th[iter] = sum((T_post_prob_selec_1 > th_seq[iter]) & (True_beta_matrix!=0))
    FP_for_var_th[iter] = sum((T_post_prob_selec_2 > th_seq[iter]) & (True_beta_matrix==0))
  }
  
  FDR = FP_for_var_th/(FP_for_var_th + TP_for_var_th)
  FDR = FDR[!is.na(FDR)]
  diff_with_req_th = abs(FDR - th_for_FDR)
  min_indices = which(diff_with_req_th == min(diff_with_req_th))
  
  desired_th_for_P = th_seq[max(min_indices)]
  #######################
  TPR_P = sum((T_post_prob_selec_1 > desired_th_for_P) & (True_beta_matrix!=0))/sum(True_beta_matrix!=0)
  FPR_P = sum((T_post_prob_selec_2 > desired_th_for_P) & (True_beta_matrix==0))/sum(True_beta_matrix==0)
  ########################
  
  T_selec_no_TP_FP = array(0, c(sample_size, num_of_p - p_regression, post_size))
  for (iter in 1:post_size)
  {
    T_selec_no_TP_FP[,,iter] = (v[,2:(num_of_p - p_regression + 1),iter] !=0)
  }
  
  T_post_prob_selec_no_TP_FP = apply(T_selec_no_TP_FP, c(1,2), mean)
  
  ######################################
  threshold_seq = seq(0,1,1/1000);
  
  tpr = rep(0,length(threshold_seq));
  fpr = rep(0,length(threshold_seq));
  
  for(iter in 1:length(threshold_seq))
  {
    tpr[iter] = sum(sum((T_post_prob_selec_no_TP_FP>=threshold_seq[iter]) & (True_beta_matrix!=0)))/sum(sum(True_beta_matrix !=0));
    fpr[iter] = sum(sum((T_post_prob_selec_no_TP_FP>=threshold_seq[iter]) & (True_beta_matrix==0)))/sum(sum(True_beta_matrix ==0));
  }
  
  auc = 0;
  
  for (iter in 1:(length(threshold_seq)-1)){
    auc = auc + tpr[iter]*(fpr[iter]-fpr[iter+1]);
  }
  
  AUC_P = auc
  
  if(AUC_P>1)
  {
    AUC_P =1
  }
}else{
  #######################
  T_selec_3k_iters = array(0, c(sample_size, num_of_p - p_regression, post_size))
   
  for (iter in 1:post_size)
  {
     T_selec_3k_iters[,,iter] = (v[,2:(num_of_p - p_regression + 1),iter] !=0) &  (True_beta_matrix !=0)
  }
  
  T_post_prob_selec = apply(T_selec_3k_iters, c(1,2), mean)
  
  TPR_P = sum((T_post_prob_selec > 0.5) & (True_beta_matrix!=0))/sum(True_beta_matrix!=0)
  ########################
  
  FPR_P = 0
  AUC_P = TPR_P
}



#######################

if(sum(temp_boolean_matrix==0)!=0)
{
  T_selec_3k_iters_1 = array(0, c(num_of_q, num_of_p - p_regression, post_size))
  T_selec_3k_iters_2 = array(0, c(num_of_q, num_of_p - p_regression, post_size))
  
  th_seq = seq(0.01,1,0.01)
  
  TP_for_var_th = rep(0, length(th_seq))
  FP_for_var_th = rep(0, length(th_seq))
  FDR = rep(0, length(th_seq))
  
  for(th_iter in 1:length(th_seq))
  {
    for (iter in 1:post_size)
    {
      T_selec_3k_iters_1[,,iter] = (v_Q[,,iter] > th_seq[th_iter]) &  (temp_boolean_matrix !=0)
      T_selec_3k_iters_2[,,iter] = (v_Q[,,iter] > th_seq[th_iter]) &  (temp_boolean_matrix ==0)
    }
    
    T_post_prob_selec_1 = apply(T_selec_3k_iters_1, c(1,2), mean)
    T_post_prob_selec_2 = apply(T_selec_3k_iters_2, c(1,2), mean)
    
    TP_for_var_th = sum((T_post_prob_selec_1 > th_seq[th_iter]) & (temp_boolean_matrix!=0))
    FP_for_var_th = sum((T_post_prob_selec_2 > th_seq[th_iter]) & (temp_boolean_matrix==0))
    
    FDR[th_iter] = FP_for_var_th/(FP_for_var_th + TP_for_var_th)
    
  }
  
  FDR = FDR[!is.na(FDR)]
  diff_with_req_th = abs(FDR - th_for_FDR)
  min_indices = which(diff_with_req_th == min(diff_with_req_th))
  desired_th_for_G = th_seq[min(min_indices)]
  
  #######################
  T_selec_3k_iters = array(0, c(num_of_q, num_of_p - p_regression, post_size))
  
  for (iter in 1:post_size)
  {
    T_selec_3k_iters[,,iter] = (v_Q[,,iter] >desired_th_for_G) &  (temp_boolean_matrix !=0)
  }
  
  T_post_prob_selec = apply(T_selec_3k_iters, c(1,2), mean)
  
  TPR_Q = sum((T_post_prob_selec > desired_th_for_G) & (temp_boolean_matrix!=0))/sum(temp_boolean_matrix!=0)
  ##############################
  T_selec_3k_iters = array(0, c(num_of_q, num_of_p - p_regression, post_size))
  
  for (iter in 1:post_size)
  {
    T_selec_3k_iters[,,iter] = (v_Q[,,iter] >desired_th_for_G) &  (temp_boolean_matrix ==0)
  }
  
  T_post_prob_selec = apply(T_selec_3k_iters, c(1,2), mean)
  
  FPR_Q = sum((T_post_prob_selec > desired_th_for_G) & (temp_boolean_matrix==0))/sum(temp_boolean_matrix==0)
  ########################
  
  T_selec_no_TP_FP = array(0, c(num_of_q, num_of_p - p_regression,  post_size))
  for (iter in 1:post_size)
  {
    T_selec_no_TP_FP[,,iter] = (v_Q[,,iter] > desired_th_for_G)
  }
  
  T_post_prob_selec_no_TP_FP = apply(T_selec_no_TP_FP, c(1,2), mean)
  #############
  
  threshold_seq = seq(0,1,1/1000);
  
  tpr = rep(0,length(threshold_seq));
  fpr = rep(0,length(threshold_seq));
  
  
  for(iter in 1:length(threshold_seq))
  {
    tpr[iter] = sum(sum((T_post_prob_selec_no_TP_FP>=threshold_seq[iter]) & (temp_boolean_matrix!=0)))/sum(sum(temp_boolean_matrix !=0));
    fpr[iter] = sum(sum((T_post_prob_selec_no_TP_FP>=threshold_seq[iter]) & (temp_boolean_matrix==0)))/sum(sum(temp_boolean_matrix ==0));
  }
  
  
  auc = 0;
  
  for (iter in 1:(length(threshold_seq)-1)){
    auc = auc + tpr[iter]*(fpr[iter]-fpr[iter+1]);
  }
  
  AUC_Q = auc
  if(AUC_Q > 1)
  {
    AUC_Q = 1
  }
  ###################################################
  
}else{
  ###################################################
  T_selec_3k_iters = array(0, c(num_of_q, num_of_p - p_regression, post_size))
  
  for (iter in 1:post_size)
  {
    T_selec_3k_iters[,,iter] = (v_Q[,,iter] >0.5) &  (temp_boolean_matrix !=0)
  }
  
  T_post_prob_selec = apply(T_selec_3k_iters, c(1,2), mean)
  
  TPR_Q = sum((T_post_prob_selec > 0.5) & (temp_boolean_matrix!=0))/sum(temp_boolean_matrix!=0)
  FPR_Q = 0
  AUC_Q = TPR_Q
  
}
