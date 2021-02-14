% Ese 01
% Parte 3/3
% Parallel reactions in a CSTR

close all;  clear variables;
clc

%% Input data

% Stoichiometric coefficients [A, B, C]
sc = [-1 -1; 1 0; 0 1];

k1 = 0.5;       % [1/min]
k2 = 0.1;       % [1/min]
k = [k1, k2];

CAin = 2;       % [mol/l]
CBin = 0;       % [mol/l]
CCin = 0;       % [mol/l]

% Conversion target
XTarget = 0.95;         % [-]

% Inlet conditions
ICs = [CAin, CBin, CCin];

%% Non linear system solution

% Change tau untill the conversion is the right one
tau_s = 1:60;

for tau = tau_s
    
    % Matrix of coefficients (Linear System)
    Z = zeros(3);
    A = -tau*(sc*k');
    Z(:,1) = A;
    A = eye(3) + Z;
    
    % Solution 
    C = A\ICs';
    
    % Conversion 
    X = (CAin - C(1))/CAin;
    
    if X >= XTarget
        break
    end
    
end

fprintf('X = %f @ tau = %f \n', X, tau)
