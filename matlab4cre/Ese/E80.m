% kmol, m^3, h, K

% A, B, C
C0 = [1, 0, 0];

r = @(C, T) exp(4.5 - 2500/T)*(C(1)^2 - C(2)*C(3)/exp(28.8 - 0.037*T - 5178/T));
DT = 1;

T_s = 550:0.01:555;
tau_s = 0.4:0.01:0.65;
DP = zeros(length(T_s), length(tau_s));

for ii = 1:length(T_s)
    for jj = 1:length(tau_s)
        
        T = T_s(ii);
        tau = tau_s(jj);
        
        [t, C] = ode23s(@(t, C)function_1(t, C, T, r), [0, tau], C0);
        
        DP(ii, jj) = C(end, 2)/(t(end)+1);
        
    end
end

[i, j] = find(DP == max(max(DP)));

T = T_s(i)
tau = tau_s(j)


function YY = function_1(~, C, T, r)
    
    YY = [-r(C, T), r(C, T), r(C, T)]';

end