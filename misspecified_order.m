n = 10;
p = 5;
q = 2;
KT_val = 0.5;

if q == 2
    th = '0p5';
elseif q == 5
    th = '1';
end 

if KT_val == 0.5
    kendall_threshold = '0p5';
elseif KT_val == 0.25
    kendall_threshold = '0p25';
end 

rng(12345*(n+p+q))
load(['data_with_th_',th,...
    '_n_',num2str(n),'_p_',num2str(p),'_q_',num2str(q),...
    '/Y_order_n_',num2str(n),'_p_',num2str(p),'_q_',num2str(q),'.mat'])

a = [1:p]';
order_store = zeros(1e5, p);
kendall_tau_value = zeros(1,1e5);

for i = 1:1e5
b = a(randperm(length(a)));
order_store(i,:) = b;
kendall_tau_value(i) = corr(a,b,'Type','Kendall', 'rows','all');
end

diff_with_KT_val = abs(KT_val - kendall_tau_value);
required_ordernigs = order_store(diff_with_KT_val == min(diff_with_KT_val),:);

[num_req_ord,~] = size(required_ordernigs);

if(num_req_ord == 1)
    final_ordering = required_ordernigs;
else
    num_correct_edges = zeros(1, num_req_ord);
    for i = 1:num_req_ord
        for j = 1:p
            cov = a(required_ordernigs(i,j));
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if cov ~=p
                true_cov = which_Ys_TP{cov};
            else
                true_cov = [];
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if j<=p-1
                    match_cov = intersect(true_cov, required_ordernigs(i,(j+1):p));
            else
                match_cov = [];
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            num_correct_edges(i) = num_correct_edges(i) + length(match_cov);
        end 
    end
    index_min_correct = find(num_correct_edges == min(num_correct_edges));
    if length(index_min_correct)>1
        final_ordering = required_ordernigs(randsample(index_min_correct,1),:);
    else
        final_ordering = required_ordernigs(index_min_correct,:);
    end
end 

csvwrite(['Misspecified_order_with_KT_',kendall_threshold,'_th_',th,...
    '_n_',num2str(n),'_p_',num2str(p),'_q_',num2str(q),'.csv'], final_ordering)

