% Ese 01
% Parte 2/3
% Parallel reactions in a Batch Reactor

close all;  clear variables;

%% Input data

% Stoichiometric coefficients [A, B, C]
sc = [-1 -1; 1 0; 0 1];

T = 100+273.15; % [K]
V = 1;          % [m3]
TauD = 1;       % [h]
rho = 800;      % [kg/m3]
MW = 35;        % [kg/kmol]
XTarget = 0.98; % [-]
PriceB = 15;    % [$/kg]
C1 = 100;       % [$/h]
C2 = 25;        % [$/h]
C3 = 8000;      % [$/h]

% Kinetic constants
k1 = 8e4*exp(-8000/1.987/T);    % [1/h]
k2 = 1e5*exp(-10000/1.987/T);   % [1/h]
k = [k1, k2];

% Initial conditions
Ctot = rho/MW;                  % [kmol/m3]
CA0 = Ctot;                     % [kmol/m3]
CB0 = 0;                        % [kmol/m3]
CC0 = 0;                        % [kmol/m3]

% Maximum simulation time 
tau = 3;                        % [h]

%% ODE solution
ICs = [CA0, CB0, CC0]';
[t, C] = ode45(@(t, C)isothermalBatch(t,C,k, sc), [0, tau], ICs, []);

%% Post processing
CA = C(:, 1);
CB = C(:, 2);
CC = C(:, 3);
X = (CA0 - CA)/CA0;

% Conversion
figure;
subplot(121)
plot(t, X)

% Find residence time (target conditions)
iTarget = find(X> XTarget, 1);
tauTarget = t(iTarget);
fprintf('tau | X = 0.98 : %f\n', tauTarget);

n = 24/(tauTarget + TauD);
PB = (CB(iTarget)*V)*n;     % [kmolB/day]
fprintf('n: %f PB: %f\n', n, PB);

% Optimal conditions
% I calculate the production of B at all residence times: 
% I get an income for every possible residence time
Income = (CB*V*MW)*24./(t + TauD)*PriceB;    
Costs = (C1*TauD + C2*t + C3)*24./(t + TauD);
Margin = Income - Costs;

% Economic analysis
subplot(122)
plot(t, Income, t, Costs, t, Margin);
legend("I", "C", "M");

[maxM, imaxM] = max(Margin);
taumaxM = t(imaxM);
fprintf('Max Margin = %f @ tau = %f \n', maxM, taumaxM);

%% ODE system

function dC = isothermalBatch(~,C,k, sc)

    % Reaction rates
    r = [k(1)*C(1), k(2)*C(1)]';

    % Calculate the derivatives of the concentrations 
    dC = sc*r;

end
