% Practical 1
% Exercise  2

% Reactions in Parallel in a Batch Reactor

close all
clear variables
clc

% -------------------------------------------------------------------------
% Data
% -------------------------------------------------------------------------

SC = [-1 -1; 1 0; 0 1];

MWa = 35;       % Molecular Weight [kg/kmol]
MWb = 35;
MWc = 35;

T = 373.15;     % [K]
rho = 800;      % [kg/m^3]

k1 = 8e4*exp(-8e3/1.987/T); % [1/h]
k2 = 1e5*exp(-1e4/1.987/T); % [1/h]

V = 1;              % [m^3]
X_target = 0.98;    % [-]

Ctot = rho/MWa;     % [kmol/m^3]
CAin = Ctot;

% Costs 
tD = 1;             % Downtime [h]
Cf = @(t) 100*tD + 25*t + 8000; % [$]
Gain_B = 15;        % [$/kg]

% -------------------------------------------------------------------------
% Solution
% -------------------------------------------------------------------------

% Inlet
tau = 3;        % [h]
Cin = [CAin, 0, 0];

[t, C] = ode23s(@function_1, [0 tau], Cin, [], k1, k2, SC);

CA = C(:,1);
CB = C(:,2);
CC = C(:,3);

Yb = CB/CAin;
Yc = CC/CAin;

X = (CAin - CA)/CAin;

Sb_c = (k1*CA - k2*CB)./(k2*CB);

% -------------------------------------------------------------------------
% Results
% -------------------------------------------------------------------------

ii = find(X >= X_target, 1);
CB_end = CB(ii);
ProdB = 24*CB_end*V/(t(ii) + 1)*MWb;

fprintf('X = %f [-] @ t = %f [h] : Daily Productivity of B = %f [kg/d] = %f [$/d]\n',...
        X(ii), t(ii), ProdB, Gain_B*ProdB);
    
% Maximize Gain
ProdB = 24*CB*V./(t + 1)*MWb;
Gain = ProdB*Gain_B;
Pain = Cf(t);
Edge = Gain - Pain;

[MaxEdge, ii_max] = max(Edge);
fprintf('Maximum Edge = %f [$]; Prod B = %f [kg/$]; X = %f [-] @ t = %f [h] \n',...
       MaxEdge, ProdB(ii_max), X(ii_max), t(ii_max));

% -------------------------------------------------------------------------
% Graphical Post Processing
% -------------------------------------------------------------------------

figure
yyaxis left
plot(t, X)
yyaxis right
plot(t, Edge)
title('Conversion and Edge')

% -------------------------------------------------------------------------
% Functions
% -------------------------------------------------------------------------

function yy = function_1(t, C, k1, k2, SC)

    r  = [k1*C(1) k2*C(1)]';
    yy = SC*r;

end