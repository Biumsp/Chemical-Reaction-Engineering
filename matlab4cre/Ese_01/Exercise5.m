close all;  clear variables;


%% Input data
T = 750+273.15; % [K]
P = 3*1e5;      % [Pa]
MW  = 25;       % [kg/kmol]
FA0 = 20;       % [kmol/h]
FB0 = 0;        % [kmol/h]
FC0 = 0;        % [kmol/h]
L = 100;        % [m]
d = 0.08;       % [m]

A1 = 2.e8;      % [1/s]
A2 = 4.e8;      % [1/s]
E1 = 40000;     % [cal/mol]
E2 = 26000;     % [cal/mol]

% Kinetic constants
k1 = A1*exp(-E1/1.987/T);   % [1/s]
k2 = A2*exp(-E2/1.987/T);   % [1/s]

% Initial conditions
Ctot = P/8314/T;        % [kmol/m3]
Ftot = FA0+FB0+FC0;     % [kmol/h]
CA0  = Ctot*FA0/Ftot;   % [kmol/m3]
CB0  = Ctot*FB0/Ftot;   % [kmol/m3]
CC0  = Ctot*FC0/Ftot;   % [kmol/m3]
Atot = pi/4*d^2;        % [m2]
Qtot = Ftot/Ctot;       % [m3/h]
v    = Qtot/Atot/3600;  % [m/h]
tau  = L/v;             % [s]


%% ODE solution (without QSSA)
C0 = [CA0 CB0 CC0]';
[t, C] = ode15s(@isothermalPFR, [0 tau], C0, [], k1,k2);

CA = C(:,1);    % [kmol/m3]
CB = C(:,2);    % [kmol/m3]
CC = C(:,3);    % [kmol/m3]

%% ODE solution (with QSSA)
C0qssa = [CA0 CC0]';
[tqssa, Cqssa] = ode45(@isothermalPFRwithQSSA, [0 tau], C0qssa, [], k1,k2);

CAqssa = Cqssa(:,1);    % [kmol/m3]
CCqssa = Cqssa(:,2);    % [kmol/m3]
CBqssa = k1/k2*CAqssa;  % [kmol/m3]


%% Post processing

% Comparison: withouth and with QSSA (linear scale)
figure;
plot(t, CB,'-o', tqssa,CBqssa,'-o');
legend('without QSSA', 'with QSSA');
xlabel('time [s]');
ylabel('C [kmol/m3]');
title('Concentration of species B');

% Comparison: withouth and with QSSA (log scale)
figure;
semilogx(t, CB,'-o', tqssa,CBqssa,'-o');
legend('without QSSA', 'with QSSA');
xlabel('time [s]');
ylabel('C [kmol/m3]');
title('Concentration of species B');


%% ODE system (defined locally, requires at least MATLAB-2016b)
function dC = isothermalPFR(t,C, k1,k2)

    CA = C(1);      % [kmol/m3]
    CB = C(2);      % [kmol/m3]
    CC = C(3);      % [kmol/m3]

    r1 = k1*CA;     % [kmol/m3/s]
    r2 = k2*CB;     % [kmol/m3/s]

    dCA = -r1;      % [kmol/m3/s]
    dCB = r1-r2;    % [kmol/m3/s]
    dCC = r2;       % [kmol/m3/s]

    dC = [dCA dCB dCC]';

end

%% ODE system (defined locally, requires at least MATLAB-2016b)
function dC = isothermalPFRwithQSSA(t,C, k1,k2)

    % Be careful, no differential equation for B is written
    CA = C(1);      % [kmol/m3]
    CC = C(2);      % [kmol/m3]
    
    % CB is evaluated on the basis of QSSA
    CB = k1/k2*CA;  % [kmol/m3]

    r1 = k1*CA;     % [kmol/m3/s]
    r2 = k2*CB;     % [kmol/m3/s]

    dCA = -r1;      % [kmol/m3/s]
    dCC = r2;       % [kmol/m3/s]

    dC = [dCA dCC]';

end