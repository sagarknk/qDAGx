This code respository has codes to reproduce plots similar to Fig. 2 in the paper. To produce the plot comparing the three estimation procedures, qDAGx with known ordering (qDAGx_0), qDAGx with misspecified ordering (qDAGx_m) and qDAGx with unknown ordering (qDAGx), follow the steps below:

#####################
## DATA Generation ##
#####################

--> Generate synthetic data by running the MATLAB script "DATA_GEN_with_threshold.m". Number of response variables, number of covariates, threshold, sample size and number of simulations can be changed by setting the variables "num_of_p", "num_of_q", "threshold_val", "sample_size" and "num_of_simulations" respectively. Varying sparsity at the response and covariates level can be changed by setting the variables "sparsity_factor_y" and "sparsity_factor_x" respectively. Running this script generates a directory which consists of required data files for desired number of replications ("num_of_simulations").

#####################
#### qDAGx_0 ########
#####################

--> To run the esimation procedure, qDAGx with known ordering (qDAGx_0), source the R script "qDAGx_o.R". This script uses "pxHS_Functions.R" and "FDR_Control_qDAGx_o.R" for required MCMC samplers and variable selection. After running "qDAGX_o.R" is complete, run the MATLAB script "qDAGX_o_Result_compilation.m". This generates 2 .csv files containing the mean and sd of the nine performance measures, across 9 quantiles, which will be used to generate the line corresponding to "qDAGX_0" in the plot (similar to Fig. 2). 

#####################
#### qDAGx_m ########
#####################

--> First generate a misspecified order with the desired rank correlation coefficient (Kendall's Tau) by running the MATLAB script "misspecified_order.m". The value of Kendall's Tau can be set to 0.5 or 0.25 (as in the paper), by setting the variable "KT_val" in "misspecified_order.m".

--> To run the esimation procedure, qDAGx with misspecified ordering (qDAGx_m), source the R script "qDAGx_m.R". This script uses "pxHS_Functions.R" and "FDR_Control_qDAGx_m.R" for required MCMC samplers and variable selection. After running "qDAGX_m.R" is complete, run the MATLAB script "qDAGX_m_Result_compilation.m". This generates 2 .csv files containing the mean and sd of the nine performance measures, across 9 quantiles, which will be used to generate the line corresponding to "qDAGX_m" in the plot (similar to Fig. 2).

#####################
#### qDAGx ##########
#####################

--> To run the esimation procedure, qDAGx with unknown ordering (qDAGx), source the R script "qDAGx.R". This script uses "qDAGx_functions_pxHS.R" for required MCMC samplers. After running "qDAGX.R" is complete, run the MATLAB script "qDAGX_Result_compilation.m". This generates 2 .csv files containing the mean and sd of the nine performance measures, across 9 quantiles, which will be used to generate the line corresponding to "qDAGX" in the plot (similar to Fig. 2). 

#####################
## Generating plot ##
#####################

--> After the required result files are generated, source the file "qDAGxo_vs_qDAGxm_vs_qDAGx_plot.R" to generate the plot, comparing the three estimation procedures qDAGx_0, qDAGx_m and qDAGx.  
