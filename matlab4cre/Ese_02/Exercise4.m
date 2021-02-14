% Ese 02
% Parte 4/3
% Isothermal CSTR with parallel reactions

close all; clear variables;
clc

% -------------------------------------------------------------------------
% Data
% -------------------------------------------------------------------------

R = 8.314462;       % [J/mol/K]
T = 353.15;         % [K]

E(1) = 2e4*4.186;   % [J/mol]
E(2) = 1.8e4*4.186; % [J/mol]
A(1) = 3e14;        % [1/s]
A(2) = 2e13;        % [1/s]
k = A.*exp(-E/R/T); % [1/s]

% A, B, C
Cin(1) = 55e3;      % [mol/m^3]
Cin(2) = 0;         % [mol/m^3]
Cin(3) = 0;         % [mol/m^3]

% Residence time
tau = 60;           % [s]

% Stoichiometry
sc = [-1, -1; 1, 0; 0, 1];

% -------------------------------------------------------------------------
% Solution
% -------------------------------------------------------------------------

% Matrix of coefficients (Linear System)
Z = zeros(3);
A = -tau*(sc*k');
Z(:,1) = A;
A = eye(3) + Z;

% Solution 
C = A\Cin';

% Conversion 
X = (Cin(1) - C(1))/Cin(1);

fprintf("Ca = %f [mol/m^3], Cb = %f [mol/m^3], Cc = %f [mol/m^3]\n", C);