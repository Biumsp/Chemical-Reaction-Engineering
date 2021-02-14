% Template fot the tanks in series model

close all
clear variables
clc

% -------------------------------------------------------------------------
% Data
% -------------------------------------------------------------------------

% Time [s]
t = [];

% RTD [1/s]
E = [];

% CDF [-]
% F = [];
% E(1) = (F(2)-F(1))/(t(2) - t(1));
% for ii = 2:length(F)-1
%     E(ii) = (F(ii+1)-F(ii-1))/(t(ii+1) - t(ii-1));
% end
% E(ii+1) = (F(end)-F(end-1))/(t(end) - t(end-1));

% Other data
tau = ;

% -------------------------------------------------------------------------
% Solution
% -------------------------------------------------------------------------

sigmaTheta2 = ZeroParametersModels(t, E);

n = 1/sigmaTheta2;

% Closest Upper Value
n_up = ceil(n);
tau_up = tau/n_up;
for ii = 1:n_up
    
    %TODO
    
end


% Closest Lower Value
n_low = floor(n);
tau_low = tau/n_low;
for ii = 1:n_low
    
    %TODO
    
end

% Interpolation
% quantity = quantity_low + (quantity_up - quantity_low)*factor
factor = n - n_low;

% TODO



