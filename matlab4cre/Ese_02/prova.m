T_s = 1:10;
lambda_s = zeros(3, length(T_s));
n_s = zeros(45, 6, length(T_s));

for ii = 1:length(T_s)
    
    T = T_s(ii);
    Nome_Script
    lambda_s(:, ii) = lambda;
    n_s(:, :, ii) = n_final;
    
end
    

clearvars -except T T_s n_s lambda_s ii
