% Practical 1
% Exercise  1

% Reactions in Series in a PFR

close all
clear variables
clc

% -------------------------------------------------------------------------
% Data  
% -------------------------------------------------------------------------

Ftot = 50/3.6;
Cp = 30*4.186;
P = 3*101325;

Dh0r = -5000*4.186;
U = 50*4.186/3.6;
Te = 573.15;
Tin = 573.15;
T_Limit = 380;

d = 8e-2;
L = 150;

V = pi*d^2/4*L;

Qin = Ftot*8.314*Tin/P;
CAin = Ftot/Qin;

k = @(T) 2e8*exp(-24e3/1.987/T);

% -------------------------------------------------------------------------
% Solution with governing equations
% -------------------------------------------------------------------------

inlet = [Ftot, 0, Tin];

[V_gov, var] = ode23s(@function_1, [0 V], inlet, [], Ftot, Cp, Qin, Dh0r, U, Te, d, k, Tin);

Fa = var(:, 1);
Fb = var(:, 2);
T_gov = var(:, 3);

Q  = Qin*T_gov/Tin;

CA = Fa./Q;
CB = Fb./Q;

X = (CAin - CA)/CAin;

% -------------------------------------------------------------------------
% Results 1
% -------------------------------------------------------------------------

disp('Governing Equations')
fprintf('Max T = %f [°C] \n', max(T_gov)-273.15);
fprintf('X = %f [-] \n', X(end))

% -------------------------------------------------------------------------
% Graphical Post Processing
% -------------------------------------------------------------------------

figure
plot(V_gov, T_gov)
title('Original U: T vs. V')

figure
plot(V_gov, CA, V_gov, CB)
title('Original U: Ca and Cb vs. V')

% -------------------------------------------------------------------------
% Find minimum U
% -------------------------------------------------------------------------


for U = U:0.1:2*U
    
    [V_gov, var] = ode23s(@function_1, [0 V], inlet, [], Ftot, Cp, Qin, Dh0r, U, Te, d, k, Tin);
    T_gov = var(:, 3);
    
    maxT = max(T_gov)-273.15;
    if maxT <= T_Limit
        break
    end
end
Fa = var(:, 1);
Q  = Qin*T_gov/Tin;
CA = Fa./Q;
X = (CAin - CA)/CAin;

% -------------------------------------------------------------------------
% Results 2
% -------------------------------------------------------------------------

disp('Minimum U')
fprintf('Min U = %f [W/m^2/K], Max T = %f [°C] \n', U, maxT);
fprintf('X = %f [-] \n', X(end))

% -------------------------------------------------------------------------
% Functions
% -------------------------------------------------------------------------

function yy = function_1(V, var, Ftot, Cp, Qin, Dh0r, U, Te, d, k, Tin)

    Fa = var(1);
    T  = var(3);
    
    Q  = Qin*T/Tin;
    
    CA  = Fa/Q;
    R   = k(T)*CA;
    Qr  = -Dh0r*R;
    Qex = U*4/d*(T - Te);
    
    dFa = -R;
    dFb =  R;
    dT  =  (Qr - Qex)/(Ftot*Cp);

    yy = [dFa dFb dT]';

end
