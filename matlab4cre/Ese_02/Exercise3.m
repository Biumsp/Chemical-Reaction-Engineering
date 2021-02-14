close all; clear variables;

global A;       % [1/s]
global Ea;      % [cal/mol]
global P;       % [Pa]
global S;       % [m2]
global Te;      % [K]
global U;       % [cal/m2/h/K]
global D;       % [m]
global deltaHR; % [cal/mol]
global Cp;      % [cal/mol/K]


%% User data
A = 2e8;            % frequency factor [1/s]
Ea = 24000;         % activation energy [cal/mol]
deltaHR = -5000;    % reaction heat [cal/mol]
Cp = 30;            % cp [cal/mol/K]

Te = 300 + 273.15;  % external temperature [K]
P = 3*101325;       % pressure [Pa]

FAin = 50;          % inlet flow rate [kmol/h]
FBin = 0;           % inlet flow rate [kmol/h]
Tin = 300 + 273.15; % inlet temperature [K]

U = 50;             % global heat exchange coefficient [kcal/m2/h/K]
D = 8e-2;           % internal diameter [m]
L = 150;            % reactor length [m]


%% Calculations

% Cross section area
S = pi/4*D^2;   % cross section area [m2]

% ODE solution
% [TODO]


%% Post processing
% [TODO]

%Plotting results
% [TODO]


%% Parametric analysis
% [TODO]


%% ODE system (defined locally, requires at least MATLAB-2016b)
function dY = PFR(z,Y)

    global A;       % [1/s]
    global Ea;      % [cal/mol]
    global P;       % [Pa]
    global S;       % [m2]
    global Te;      % [K]
    global U;       % [cal/m2/h/K]
    global D;       % [m]
    global deltaHR; % [cal/mol]
    global Cp;      % [cal/mol/K]

    % [TODO]

end
