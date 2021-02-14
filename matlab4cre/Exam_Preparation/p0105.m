% Practical 1
% Exercise  4

% Reactions in Series in a PFR and QSSA

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

k1 = 2.0e8*exp(-40e3/1.987/T);   % Kinetic Constant 1 [1/s]
k2 = 4.0e8*exp(-26e3/1.987/T);   % Kinetic Constant 2 [1/s]

% -------------------------------------------------------------------------
% Solution
% -------------------------------------------------------------------------

% Inlet 
Cin = [CAin, CBin, CCin];

% We ask matlab to interpolate the results to evaluate the maximum
% difference at lines 91 & 92
[~, C] = ode23s(@function_1, 0:0.1:tau, Cin, [], k1, k2, SC);

CA = C(:,1);
CB = C(:,2);
CC = C(:,3);

Yb = CB/CAin;
Yc = CC/CAin;

% Selectivity should always be approximately 0, because Cb ~ k1*Ca/k2
Sb_c = (k1*CA - k2*CB)./(k2*CB);

% -------------------------------------------------------------------------
% Solution with QSSA
% -------------------------------------------------------------------------

% inlet 
Cin = [CAin CCin];
[t, C_qssa] = ode23s(@function_2, 0:0.1:tau, Cin, [], k1, k2);

CA_qssa = C_qssa(:,1);
CB_qssa = CA*k1/k2;
CC_qssa = C_qssa(:,2);

% -------------------------------------------------------------------------
% Analytical Solution 
% -------------------------------------------------------------------------

CA_an = CAin*exp(-k1*t);
CB_an = k1/(k2 - k1)*(exp(-k1*t) - exp(-k2*t))*CAin;
CC_an = CAin - CA - CB;

Yb_an = CB/CAin;
Yc_an = CC/CAin;
Sb_c_an = (k1*CA - k2*CB)./(k2*CB);

% -------------------------------------------------------------------------
% Results
% -------------------------------------------------------------------------

maxError_num = max(abs(CB - CB_qssa));
maxError_an  = max(abs(CB_an - CB_qssa));

fprintf('Max Error (numerical) = %f [kmol/m^3]; max Error (analytical) = %f [kmol/m^3] \n',...
        maxError_num, maxError_an)
    
% -------------------------------------------------------------------------
% Graphical Post Processing
% -------------------------------------------------------------------------

figure
plot(t, C)
title('Concentrations')
legend(['CA'; 'CB'; 'CC'])

figure
plot(t, Sb_c)
title('Selectivity B/C (USELESS PLOT: line 60)')

figure
plot(t, Yb, t, Yc)
title('Yields')
legend('Yb', 'Yc')

figure
plot(t, Yb, '-o', t, Yb_an, t, Yc, '-o', t, Yc_an)
title('Analytical vs. Numerical Yields')
legend('Yb', 'Yb_an', 'Yc', 'Yc_an')

% -------------------------------------------------------------------------
% Functions
% -------------------------------------------------------------------------

function yy = function_1(t, C, k1, k2, SC)

    r  = [k1*C(1) k2*C(2)]';
    yy = SC*r;
    
end

function yy = function_2(t, C, k1, k2)

    r  = [k1*C(1) k2*C(2)]';
    yy = [-r(1) r(2)]';

end