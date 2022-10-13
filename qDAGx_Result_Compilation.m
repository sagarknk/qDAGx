num_of_q = 2;
num_of_p = 5;
sample_size = 10;
temp_skip_total = 1;
number_of_simulations = 5;

if num_of_q == 2
    true_th = '0p5';
elseif num_of_q == 5
    true_th = '1';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
grand_result_mean_summary = zeros(11,9,number_of_simulations);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for tau = 1:9
    temp_mat = zeros(11,number_of_simulations);
    for data_set = 1:number_of_simulations  
            try
                filename = ['DAG_Results_with_th_',true_th,'_n_',num2str(sample_size),'_p_',...
                    num2str(num_of_p),'_q_',num2str(num_of_q),'/Post_process_MC_',...
                    '_simulation_index_',num2str(data_set),'_tau_index_',num2str(tau),'.csv'];
                               
                out_from_file = csvread(filename);
                %%%% Scaling theta and beta diff norms
                out_from_file(9,1) = out_from_file(9,1)/sample_size;
                out_from_file(10,1) = out_from_file(10,1)/sample_size;
                temp_mat(:,data_set) = out_from_file;
                
            catch
                fprintf("%d is the total missed as of now\n", temp_skip_total);
                temp_skip_total = temp_skip_total + 1;
            end 
   
    end 
    grand_result_mean_summary(:,tau,:) = temp_mat;
   
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
temp_mean = mean(grand_result_mean_summary,3);
temp_std = std(grand_result_mean_summary,0,3);

file_name = ['Union_DAG_All_in_one_result_mean_with_th_',true_th,'_',num2str(num_of_p),'_q_',num2str(num_of_q),'_n_',num2str(sample_size),'.csv'];
csvwrite(file_name, temp_mean)

file_name = ['Union_DAG_All_in_one_result_std_with_th_',true_th,'_',num2str(num_of_p),'_q_',num2str(num_of_q),'_n_',num2str(sample_size),'.csv'];
csvwrite(file_name, temp_std)
