% close all; 
clear variables;

global A;       % [1/s]
global Ea;      % [cal/mol]
global P;       % [Pa]
global S;       % [m2]
global U;       % [cal/m2/h/K]
global D;       % [m]
global deltaHR; % [cal/mol]
global Cp;      % [cal/mol/K]
global Fe;      % [kmol/h]
global Cpe;     % [cal/mol/K]
global Se;      % [m2]


%% User data
A = 2e8;            % frequency factor [1/s]
Ea = 24000;         % activation energy [cal/mol]

% Independent from T
deltaHR = -5000;    % reaction heat [cal/mol]
Cp = 30;            % cp [cal/mol/K]

% External Fluid 
Fe = 30;        % [kmol/h]
Cpe = 30;       % [cal/mol/K]
Se = 0.005;     % [m2]

% Independent from T
Tein = 200+273.15;  % external temperature [K]
P = 3*101325;       % pressure [Pa]

FAin = 50;          % inlet flow rate [kmol/h]
FBin = 0;           % inlet flow rate [kmol/h]
Tin = 300 + 273.15; % inlet temperature [K]

U = 20;             % global heat exchange coefficient [kcal/m2/h/K]
D = 8e-2;           % internal diameter [m]
L = 150;            % reactor length [m]


%% Calculations

S = pi/4*D^2;           % Crosse section Area [m2]
Yin = [FAin, FBin, Tin, Tein]';
[z, Y] = ode45(@PFR, [0, L], Yin, []);


%% Post processing
Fa = Y(:, 1);
Fb = Y(:, 2);
T = Y(:, 3);
Te = Y(:, 4);
X = (FAin-Fa)/FAin;

figure
yyaxis left
plot(z, T-273.15, z, Te-273.15)
legend("T", "Text")
% yyaxis right
% plot(z, X)

%% Parametric analysis
% [TODO]


%% ODE system (defined locally, requires at least MATLAB-2016b)
function dY = PFR(~,Y)

    global A;       % [1/s]
    global Ea;      % [cal/mol]
    global P;       % [Pa]
    global S;       % [m2]
    global U;       % [cal/m2/h/K]
    global D;       % [m]
    global deltaHR; % [cal/mol]
    global Cp;      % [cal/mol/K]
    global Fe;      % [kmol/h]
    global Cpe;     % [cal/mol/K]
    global Se;      % [m2]

    Fa = Y(1);      % [kmol/h]
    Fb = Y(2);
    T = Y(3);
    Te = Y(4);
    
    % Time is in hours
    k = A*exp(-Ea/1.987/T)*3600; % [1/h]
    
    Ftot = Fa+Fb;       % [kmol/h]
    Ctot = P/8314/T;    % [kmol/m3]
    CA = Fa/Ftot*Ctot;  % [kmol/m3]
    r = k*CA;           % [kmol/m3/h]
    Qr = -deltaHR*r;        % [cal/m3/h] (?)
    
    a = 4/D;            % [1/m]
    
    dFa = -r*S; % [kmol/h/m]
    dFb = r*S;
    dT = (Qr + U*(Te-T)*a)*S/Ftot/Cp; % [K/m]
    dTe = -U*(Te-T)*a*Se/Fe/Cpe;
    
    dY = [dFa, dFb, dT, dTe]';

end
