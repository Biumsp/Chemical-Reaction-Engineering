% Ese 01
% Parte 4/3
% Optimization of batch operations

close all; clear variables;
options = [];
clc

% -------------------------------------------------------------------------
% Data
% -------------------------------------------------------------------------

R = 8.314462;   % [J/mol/K]
T = 500;        % [K]
V = 0.5;        % [m^3]

% Stoichiometry (Nc x Nr)
sc = [-1, -1; 1, 0; 0, 1];

Ea(1) = 9000*4.186;  % [J/mol]
Ea(2) = 19000*4.186; % [J/mol]
A(1) = 1.5e4/3600;   % [1/s]
A(2) = 6e6/3600;     % [1/s]
k = A.*exp(-Ea/R/T); % [1/s]

Cin = zeros(3, 1);
Cin(1) = 20e3/V;  % [mol/m^3]
Cin(2) = 0;       % [mol/m^3]
Cin(3) = 0;       % [mol/m^3]

% Downtimes per load
DT = 3600;        % [s]

% -------------------------------------------------------------------------
% Optimization
% -------------------------------------------------------------------------

tau_s = 0.001*3600:30:3*3600;

% Yeld of B
Y = zeros(1, length(tau_s));     % [-]

% Conversion of A
X = zeros(1, length(tau_s));     % [-]

% Cycle Length
CL = zeros(1, length(tau_s));    % [s]

% Cycles Per Day
CPD = zeros(1, length(tau_s));   % [-]

% Production of B
P = zeros(1, length(tau_s));     % [mol]

for ii = 1:length(tau_s)
    
    tau = tau_s(ii);
    [t, C] = ode45(@(t, C)sc*(k.*C(1:2)')', [0, tau], Cin, options);
    
    X(ii) = (Cin(1)-C(length(t), 1))/(Cin(1));
    Y(ii) = C(length(t), 2)/Cin(1);
    CL(ii) = tau + DT;
    CPD(ii) = 1440*60./CL(ii);
    P(ii) = CPD(ii)*C(length(t), 2)*V;
    
end

yyaxis left
plot(tau_s, P/1000)
xlabel("Residence Time [s]")
ylabel("Productivity [kmol/d]")
yyaxis right
plot(tau_s, Y)
ylabel("Yeld of B [-]")


[~, index] = max(P);
t_max = tau_s(index);
h = floor(t_max/3600);
min = floor((t_max - h*3600)/60);
fprintf("Residence time for maximum productivity: %f h %f min\n", h, min)

