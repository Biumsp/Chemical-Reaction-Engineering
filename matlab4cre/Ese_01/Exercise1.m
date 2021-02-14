% Ese 01
% Parte 1/3
% Reactions in series in a Plug Flow Reactor

close all;  clear variables;


%% Input data

% Stoichiometric coefficients [A, B, C]
sc = [-1 0; 1 -1; 0 1];

T = 750+273.15; % [K]
P = 3*1e5;      % [Pa]
MW  = 25;       % [kg/kmol]
FA0 = 20;       % [kmol/h]
FB0 = 0;        % [kmol/h]
FC0 = 0;        % [kmol/h]
L = 100;        % [m]
d = 0.08;       % [m]

A1 = 1.2e8;     % [1/s]
A2 = 4.e8;      % [1/s]
E1 = 37000;     % [cal/mol]
E2 = 39000;     % [cal/mol]

% Kinetic constants
k1 = A1*exp(-E1/1.987/T);   % [1/s]
k2 = A2*exp(-E2/1.987/T);   % [1/s]

% Initial conditions
Ctot = P/8314/T;            % [kmol/m3]
Ftot = FA0 + FB0 + FC0;     % [kmol/h]
CA0 = FA0/Ftot*Ctot;        % [kmol/m3]
CB0 = FB0/Ftot*Ctot;        % [kmol/m3]
CC0 = FC0/Ftot*Ctot;        % [kmol/m3]

Atot = pi/4*d^2;            % [m^2]
Qtot = Ftot/Ctot;           % [m3/h]
v = Qtot/3600/Atot;         % [m/s]
tau = L/v;                  % [s]

%% ODE solution

% Initial Conditions
ICs = [CA0, CB0, CC0]';

% Solution of the ODEs
[t, C] = ode45(@(t, C)isothermalPFR(t, C, [k1, k2], sc), [0, tau], ICs, []);


%% Post processing

% Concentrations
figure;
subplot(221)
plot(t, C, "-o")
subplot(222)
% Compare with the analitical solution
plot(t, C(:, 1), '-', t, CA0*exp(-k1*t), 'o')
title("Compare with analitical solution")

% Yields
YB = C(:, 2)/CA0;
YC = C(:, 3)/CA0;
subplot(223)
plot(t, YB, t, YC)
legend('YB', 'YC')

% Selectivities
sB = C(:, 2)./(CA0 - C(:, 1));
sC = C(:, 3)./(CA0 - C(:, 1));
subplot(224)
plot(t, sB, t, sC)
legend('sB', 'sC')

% Max B concentration and yelds
[maxB , imaxB] = max(C(:, 2));
tmaxB = t(imaxB);
fprintf('Maximum YB = %f, Maximum CB = %f @ t = %f s\n', YB(imaxB), maxB, tmaxB);

%% ODE system

function dC = isothermalPFR(~, C, k, sc)

    % Reaction rates
    r = [k(1)*C(1), k(2)*C(2)]';

    % Calculate the derivatives of the concentrations 
    dC = sc*r;
    
end
