% Practical 2
% Exercise  1

% Design of an adiabatic PFR

close all
clear variables
clc

% -------------------------------------------------------------------------
% Data  
% -------------------------------------------------------------------------

kf_exp = 31.1/3600; % Forward Kinetic Constant  [1/s]    @ T = 360 [K]
Ea = 6.57e4;        % Activation Energy         [J/mol]
Keq_exp = 3.03;     % Equilibrium Constant      [-]      @ T = 333.15 [K]
Dh0r_exp = -6.9;    % Reaction Heat             [J/mol]  @ T = 300 [K]

Tin = 330;          % Temperature [K]
Fin = 163/3.6;      % Inlet Molar Flow Rate [mol/s]

CAin = 9.30e-3;     % [mol/m^3]
CBin = 0;           % [mol/m^3]
CIin = CAin/9;      % [mol/m^3]
Ctot_in = CAin + CBin + CIin;

Qin = Fin/Ctot_in;  % [m^3/s]

% [J/mol/K]
CpA = 131; 
CpB = 171;
CpI = 161; 

X_target = 0.6;     % [-]

% -------------------------------------------------------------------------
% Data Pre-Processing
% -------------------------------------------------------------------------

FAin = CAin/Ctot_in*Fin;    % [mol/s]
FBin = CBin/Ctot_in*Fin;    % [mol/s]
FIin = CIin/Ctot_in*Fin;    % [mol/s]

% Definitions
ThetaA = FAin/FAin; % [-]
ThetaB = FBin/FAin; % [-]
ThetaI = FIin/FAin; % [-]

CpIn = ThetaA*CpA + ThetaB*CpB + ThetaI*CpI;
DCp = -CpA + CpB; 

% Functions
Dh0r = @(T) Dh0r_exp + DCp*(T - 300);
Keq  = @(T) Keq_exp*exp(-(1./T - 1/333.15)*(Dh0r_exp - DCp*300)/8.314...
            + DCp/8.314*log(T/333.15));
kf   = @(T) kf_exp*exp(-Ea/8.314*(1./T - 1/360));
kb   = @(T) kf(T)./Keq(T);       

% Reaction Heat
Qr = @(T, CA, CB) -Dh0r(T).*(kf(T)*CA - kb(T)*CB);  

% -------------------------------------------------------------------------
% Solution: design equations
% -------------------------------------------------------------------------

V = FAin*integral(@(X)function_1(X, kf, kb, Dh0r, Tin, CpIn, DCp, CAin), 0, X_target);

% -------------------------------------------------------------------------
% Solution: verify using he governing equations
% -------------------------------------------------------------------------

% Inlet
inlet = [FAin FBin FIin Tin];

[V, var] = ode23s(@function_2, [0 V], inlet, [],...
                   Qin, Tin, Qr, CpIn, DCp, kf, kb, CAin);
               
Fa = var(:, 1);
Fb = var(:, 2);
T  = var(:, 3);

Q  = Qin*T/Tin;
CA = Fa./Q;
CB = Fb./Q;

X = (CAin - CA)/CAin;

% -------------------------------------------------------------------------
% Results
% -------------------------------------------------------------------------

fprintf('Target Conversion: %f [-] \n', X_target)
fprintf('Volume = %f [m^3]; Conversion = %f [-]\n', V(end), X(end));

% -------------------------------------------------------------------------
% Graphical Post Processing
% -------------------------------------------------------------------------

figure
plot(V, CA, V, CB)
title('Concentrations')
legend('CA', 'CB')

figure
plot(V, T)
title('Temperature')

figure
plot(V, X)
title('Conversion')

% -------------------------------------------------------------------------
% Functions
% -------------------------------------------------------------------------

function yy = function_1(X, kf, kb, Dh0r, Tin, CpIn, DCp, CAin)

    T = Tin + (-Dh0r(Tin)*X)./(CpIn + DCp*X);
    CA = CAin*(1 - X).*Tin/T;
    CB = CAin*X;

    RA = -kf(T).*CA + kb(T).*CB;
    yy = -1./RA;

end

function yy = function_2(V, var, Qin, Tin, Qr, CpIn, DCp, kf, kb, CAin)

    Fa = var(1);
    Fb = var(2);
    Fi = var(3);
    T  = var(4);
    
    Q  = Qin*T/Tin;
    CA = Fa/Q;
    CB = Fb/Q;
    
    Ftot = Fa + Fb + Fi;
    
    X = (CAin - CA)/CAin;    
    r = kf(T)*CA - kb(T)*CB;
    
    dFa = -r;
    dFb =  r;
    dFi =  0;
    dT  =  Qr(T, CA, CB)/Ftot/(CpIn + DCp*X);
    
    yy = [dFa dFb dFi dT]';

end