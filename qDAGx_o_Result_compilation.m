sample_size = 10;
num_of_p = 5;
num_of_q = 2;
number_of_simulations = 5;

if num_of_q == 2
    true_th = '0p5';
elseif num_of_q == 5
    true_th = '1';
end

temp_skip_total = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MSE_scales = csvread(['data_with_th_',true_th,'_n_',num2str(sample_size),'_p_',...
    num2str(num_of_p),'_q_',num2str(num_of_q),'/MSE_Scales.csv']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
grand_result_mean_summary = zeros(11,9,25);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for data_set = 1:number_of_simulations
    for tau = 1:9
        temp_mat = [];
        for p_reg = 1:(num_of_p-1)
            try
                filename = ['HS_Results_with_th_',true_th,'_n_',num2str(sample_size),'_p_',...
                    num2str(num_of_p),'_q_',num2str(num_of_q),'/Post_process_MC_p_reg_',num2str(p_reg),...
                    '_simulation_index_',num2str(data_set),'_tau_index_',num2str(tau),'.csv'];
                
                out_from_file = csvread(filename);
                out_from_file = out_from_file';
                
                % scaling MSE % +1 is added because of intercept
                out_from_file(end -1) = out_from_file(end - 1)/(MSE_scales(p_reg)+1);
                
                temp_mat=[temp_mat ; out_from_file];
                
            catch
                fprintf("%d is the total missed as of now\n", temp_skip_total);
                temp_skip_total = temp_skip_total + 1;
            end
        end
        
        theta_norm_adj = sqrt(sum(temp_mat(:,1).*temp_mat(:,1)))/sample_size;
        beta_norm_adj = sqrt(sum(temp_mat(:,2).*temp_mat(:,2)))/sample_size;
        max_time_taken = max(temp_mat(:,end))/60; %% time in minutes
        
        mean_result = mean(temp_mat);
        
        mean_result(1) = theta_norm_adj;
        mean_result(2) = beta_norm_adj;
        mean_result(end) = max_time_taken;
        
        grand_result_mean_summary(:,tau,data_set) = mean_result;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
temp_mean = mean(grand_result_mean_summary,3);
temp_std = std(grand_result_mean_summary,0,3);

file_name = ['HS_All_in_one_result_mean_with_th_',true_th,'_',num2str(num_of_p),'_q_',num2str(num_of_q),'_n_',num2str(sample_size),'.csv'];
csvwrite(file_name, temp_mean)

file_name = ['HS_All_in_one_result_std_with_th_',true_th,'_',num2str(num_of_p),'_q_',num2str(num_of_q),'_n_',num2str(sample_size),'.csv'];
csvwrite(file_name, temp_std)
