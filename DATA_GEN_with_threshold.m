% clear all;
threshold_val = 0.5;
num_of_p = 5;
num_of_q = 2;
sample_size = 10;
sparsity_factor_y = 5;
sparsity_factor_X = 2;
num_of_simulations = 5;

if threshold_val == 0.5
    dir_name_create= ['mkdir data_with_th_0p5_n_',num2str(sample_size),'_p_',num2str(num_of_p),'_q_',num2str(num_of_q)];
    dir_name = ['data_with_th_0p5_n_',num2str(sample_size),'_p_',num2str(num_of_p),'_q_',num2str(num_of_q)];
    system(dir_name_create);
elseif threshold_val == 1
    dir_name_create= ['mkdir data_with_th_1_n_',num2str(sample_size),'_p_',num2str(num_of_p),'_q_',num2str(num_of_q)];
    dir_name = ['data_with_th_1_n_',num2str(sample_size),'_p_',num2str(num_of_p),'_q_',num2str(num_of_q)];
    system(dir_name_create);
end
ts = linspace(0.1,0.9,9); %%% 9 equally spaced quantile levels
which_Ys_TP = cell(1, num_of_p);
which_Xs_TP = cell(1, num_of_p);

for p_regression = 1:(num_of_p-1)
    temp_array = (p_regression+1):num_of_p;
    temp_sparsity_num_Y = int16(length(temp_array)/sparsity_factor_y);
    
    if temp_sparsity_num_Y<1
        temp_sparsity_num_Y = 1;
    end
   
    temp_Ys_TP = datasample(temp_array, temp_sparsity_num_Y, 'Replace', false);
    which_Ys_TP{1,p_regression} = temp_Ys_TP;
    which_Xs_TP{1,p_regression} = cell(1, temp_sparsity_num_Y);
    
    zero_counter = 1;
    for x_index = 1:temp_sparsity_num_Y
        
        temp_sparsity_num_X = datasample(0:sparsity_factor_X, 1);
        
        if temp_sparsity_num_Y ==1 &&  temp_sparsity_num_X == 0
            temp_sparsity_num_X =1;
        end
        
        if zero_counter == temp_sparsity_num_Y
            temp_sparsity_num_X = 1;
        end 
        
        if temp_sparsity_num_X>0
            temp_Gs_TP = datasample(1:num_of_q, temp_sparsity_num_X, 'Replace', false);
            which_Xs_TP{1, p_regression}{1,x_index} = temp_Gs_TP;
        else 
            % dummy value
            which_Xs_TP{1, p_regression}{1,x_index} = -999;
        end 
        zero_counter = zero_counter+1;
    end 
end 

file_name = [dir_name,'/Y_order_n_',num2str(sample_size),'_p_',num2str(num_of_p),'_q_'...
    ,num2str(num_of_q),'.mat'];
save(file_name, 'which_Ys_TP');

file_name = [dir_name,'/X_order_n_',num2str(sample_size),'_p_',num2str(num_of_p),'_q_'...
    ,num2str(num_of_q),'.mat'];
save(file_name, 'which_Xs_TP');

dir_name_create = ['mkdir ',dir_name,'/True_betas/'];
system(dir_name_create);
dir_name_create = ['mkdir ',dir_name,'/True_thetas/'];
system(dir_name_create);

for num_simu = 1:num_of_simulations
    
    rng(num_simu*100);
    
    QR_X_sim = normrnd(0,1,sample_size, num_of_q);    
    QR_Y_sim = zeros(sample_size, num_of_p);
    QR_Y_sim(:,num_of_p) = normrnd(0,1,sample_size,1);
    
    for p_regression = flip(1:(num_of_p-1))
        taus = rand(sample_size,1);
        corresponding_Ys = which_Ys_TP{1, p_regression};
        corresponding_Xs = which_Xs_TP{1, p_regression};
        
        for jjj = 1:length(corresponding_Ys)
            
            temp_X_factor_for_Y = zeros(sample_size, length(corresponding_Xs{1,jjj}));
            temp_X_array = corresponding_Xs{1,jjj};
            
            if min(temp_X_array) ~= -999
                temp_X_factor_for_Y = QR_X_sim(:,temp_X_array);
                if length(temp_X_array) == 1
                    temp_X_factor_for_Y(:,1) = temp_X_factor_for_Y(:,1).*temp_X_factor_for_Y(:,1)+...
                        log(taus.*taus + 1);
                    
                elseif (length(temp_X_array) ==2)
                    temp_X_factor_for_Y(:,1) = temp_X_factor_for_Y(:,1).*temp_X_factor_for_Y(:,1)+...
                        log(taus.*taus + 1);
                    temp_X_factor_for_Y(:,2) = exp(temp_X_factor_for_Y(:,2));
                    
                elseif (length(temp_X_array) ==3)
                    temp_X_factor_for_Y(:,1) = temp_X_factor_for_Y(:,1).*temp_X_factor_for_Y(:,1)+...
                        log(taus.*taus + 1);
                    temp_X_factor_for_Y(:,2) = exp(temp_X_factor_for_Y(:,2));
                    temp_X_factor_for_Y(:,3) = log(abs(temp_X_factor_for_Y(:,3)));
                end
            else
                temp_X_factor_for_Y(:,1)  = 1+taus.*taus;
            end 
            
            % sum along rows
            QR_Y_sim(:, p_regression) = QR_Y_sim(:,p_regression) + ...
                (sum(temp_X_factor_for_Y,2).*(abs(sum(temp_X_factor_for_Y,2))>threshold_val))...
                .*QR_Y_sim(:,corresponding_Ys(jjj));
        end 
        
        QR_Y_sim(:,p_regression) = QR_Y_sim(:,p_regression) + 2*taus;
        
    end 
    
    QR_Y_sim_TRUEs = zeros(sample_size, length(ts),num_of_p-1);
    
    for tau_index = 1:length(ts)
        true_tau_vec = ts(tau_index)*ones(sample_size,1);
        
        for p_regression = 1:(num_of_p-1)
            corresponding_Ys = which_Ys_TP{1, p_regression};
            corresponding_Xs = which_Xs_TP{1, p_regression};
            true_beta_vec = zeros(sample_size, num_of_p - p_regression);
            true_theta_vec = zeros(sample_size, num_of_p - p_regression);
            
            for jjj = 1:length(corresponding_Ys)

                temp_X_factor_for_Y = zeros(sample_size, length(corresponding_Xs{1,jjj}));
                temp_X_array = corresponding_Xs{1,jjj};

                if min(temp_X_array) ~= -999
                    temp_X_factor_for_Y = QR_X_sim(:,temp_X_array);
                    if length(temp_X_array) == 1
                        temp_X_factor_for_Y(:,1) = temp_X_factor_for_Y(:,1).*temp_X_factor_for_Y(:,1)+...
                            log(true_tau_vec.*true_tau_vec + 1);
                        
                    elseif (length(temp_X_array) ==2)
                        temp_X_factor_for_Y(:,1) = temp_X_factor_for_Y(:,1).*temp_X_factor_for_Y(:,1)+...
                            log(true_tau_vec.*true_tau_vec + 1);
                        temp_X_factor_for_Y(:,2) = exp(temp_X_factor_for_Y(:,2));
                        
                    elseif (length(temp_X_array) ==3)
                        temp_X_factor_for_Y(:,1) = temp_X_factor_for_Y(:,1).*temp_X_factor_for_Y(:,1)+...
                            log(true_tau_vec.*true_tau_vec + 1);
                        temp_X_factor_for_Y(:,2) = exp(temp_X_factor_for_Y(:,2));
                        temp_X_factor_for_Y(:,3) = log(abs(temp_X_factor_for_Y(:,3)));
                    end
                else
                    temp_X_factor_for_Y(:,1)  = 1+true_tau_vec.*true_tau_vec;
                end

                % sum along rows
                QR_Y_sim_TRUEs(:, tau_index ,p_regression) = QR_Y_sim_TRUEs(:,tau_index,p_regression) + ...
                    (sum(temp_X_factor_for_Y,2).*(abs(sum(temp_X_factor_for_Y,2))>threshold_val))...
                    .*QR_Y_sim(:,corresponding_Ys(jjj));
                
                true_beta_vec(:, corresponding_Ys(jjj) - p_regression) = ...
                    (sum(temp_X_factor_for_Y,2).*(abs(sum(temp_X_factor_for_Y,2))>threshold_val));
                
                true_theta_vec(:, corresponding_Ys(jjj) - p_regression) = ...
                    (sum(temp_X_factor_for_Y,2));
                
            end 

            QR_Y_sim_TRUEs(:,tau_index, p_regression) = QR_Y_sim_TRUEs(:,tau_index , p_regression) + 2*true_tau_vec;
            
            file_name = [dir_name,'/True_betas/True_beta_matrix_n',num2str(sample_size), '_p_',num2str(num_of_p),'_q_', ...
                num2str(num_of_q),'_p_reg_',num2str(p_regression), '_tau_',num2str(tau_index),...
                '_simu_',num2str(num_simu),'.csv'];
            csvwrite(file_name, true_beta_vec);
      
            file_name = [dir_name,'/True_thetas/True_theta_matrix_n',num2str(sample_size), '_p_',num2str(num_of_p),'_q_', ...
                num2str(num_of_q),'_p_reg_',num2str(p_regression), '_tau_',num2str(tau_index),...
                '_simu_',num2str(num_simu),'.csv'];
            csvwrite(file_name, true_theta_vec);
            
        end
    end 
    
    filename_output = [dir_name,'/QR_Y_sim_TRUEs_n_',num2str(sample_size),'_seed_',num2str(num_simu),'.csv'];
    csvwrite(filename_output, QR_Y_sim_TRUEs);
    
    filename_output = [dir_name,'/QR_X_sim_n_',num2str(sample_size),'_seed_',num2str(num_simu),'.csv'];
    csvwrite(filename_output, QR_X_sim);
    
    filename_output = [dir_name,'/QR_Y_sim_n_',num2str(sample_size),'_seed_',num2str(num_simu),'.csv'];
    csvwrite(filename_output, QR_Y_sim);
    
end

file_name = [dir_name,'/Y_order_n_',num2str(sample_size),'_p_',num2str(num_of_p),'_q_'...
    ,num2str(num_of_q),'_v6.mat'];
save(file_name, 'which_Ys_TP', '-v6');

file_name = [dir_name,'/X_order_n_',num2str(sample_size),'_p_',num2str(num_of_p),'_q_'...
    ,num2str(num_of_q),'_v6.mat'];
save(file_name, 'which_Xs_TP','-v6');

for i = 1:(num_of_p -1)
a(i) = length(which_Ys_TP{i});
end

csvwrite([dir_name,'/MSE_Scales.csv'],a);
