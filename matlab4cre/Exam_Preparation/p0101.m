% Practical 1
% Exercise  1

% Reactions in Series in a PFR

close all
clear variables
clc

% -------------------------------------------------------------------------
% Data  
% -------------------------------------------------------------------------

% Stoichiometry
SC = [-1 0; 1 -1; 0 1];

d = 8e-2;       % diameter [m]
L = 100;        % length [m]

T = 750+273.15; % Temperature [K]
P = 3e5;        % Pressure [Pa]

MWa = 25;       % Molecular Weight [kg/kmol]

Fin = 20;       % Inlet Molar Flow Rate [kmol/h]
FAin = 20;      
FBin = 0;
FCin = 0;

Ctot = P/8314/T;        % Total Concentration [kmol/m^3]
CAin = FAin/Fin*Ctot;
CBin = FBin/Fin*Ctot;
CCin = FCin/Fin*Ctot;

Q = Fin/Ctot/3600;      % Volumetric Flow Rate [m^3/s]
v = Q/(pi*d^2/4);       % Velocity [m/s]
tau = L/v;              % Residence Time [s]

k1 = 1.2e8*exp(-37e3/1.987/T);   % Kinetic Constant 1 [1/s]
k2 = 4.0e8*exp(-39e3/1.987/T);   % Kinetic Constant 2 [1/s]

% -------------------------------------------------------------------------
% Solution
% -------------------------------------------------------------------------

% Inlet 
Cin = [CAin, CBin, CCin];

[t, C] = ode23s(@function_1, [0 tau], Cin, [], k1, k2, SC);

CA = C(:,1);
CB = C(:,2);
CC = C(:,3);

Yb = CB/CAin;
Yc = CC/CAin;

Sb_c = (k1*CA - k2*CB)./(k2*CB);

% -------------------------------------------------------------------------
% Graphical Post Processing
% -------------------------------------------------------------------------

figure
plot(t, C)
title('Concentrations')
legend(['CA'; 'CB'; 'CC'])

figure
plot(t, Sb_c)
title('Selectivity B/C')

figure
plot(t, Yb, t, Yc)
title('Yields')
legend('Yb', 'Yc')

% Comparison with analytical result
CA = CAin*exp(-k1*t);
CB = k1/(k2 - k1)*(exp(-k1*t) - exp(-k2*t))*CAin;
CC = CAin - CA - CB;

Yb_an = CB/CAin;
Yc_an = CC/CAin;
Sb_c_an = (k1*CA - k2*CB)./(k2*CB);

figure
subplot(121)
plot(t, Yb, '-o', t, Yb_an, t, Yc, '-o', t, Yc_an)
title('Analytical vs. Numerical Yields')
legend('Yb', 'Yb_an', 'Yc', 'Yc_an')

subplot(122)
plot(t, Sb_c, '-o', t, Sb_c_an)
title('Analytical vs. Numerical Selectivity')
legend('Sb_c', 'Sb_c_an')

% -------------------------------------------------------------------------
% Functions
% -------------------------------------------------------------------------

function yy = function_1(t, C, k1, k2, SC)

    r  = [k1*C(1) k2*C(2)]';
    yy = SC*r;

end