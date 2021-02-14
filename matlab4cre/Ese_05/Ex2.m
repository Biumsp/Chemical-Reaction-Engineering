close all; clear variables;

global k;            % [1/min]
global CAin;         % [kmol/m3]
global tauEff;       % [min]
global alpha;        % [-]

%% Experimental data                                                       
nExp = 10;
tExp = [0 0.5 1 1.5 2 2.5 3 3.5 4 5]';
FExp = [0.0421434 0.4344916 0.6686519 0.8122302 0.8811733 ...
        0.9367312 0.9622991 0.9786636 0.9827176 0.995]';

%% User data                                                               
V    = 0.7;    % total volume [l]   
Q    = 0.7;    % total volumetric flow rate [l/min]
CAin = 10;     % inlet concentration of A [kmol/m3]
k    = 5;      % kinetic constant [m3/kmol/min]
tau  = V/Q;    % nominal residence time [min]

%% Ideal CSTR solution                                                     %
CAcstr = (-1+sqrt(1+4*k*tau*CAin))/(2*k*tau);

%% Analysis of experimental data                                           %

% Linear regression
% [TODO]

% Statistics
SSres = (Y-X*a)'*(Y-X*a);
SStot = (Y-mean(Y))'*(Y-mean(Y));
R2 = 1 - SSres/SStot;

% Calculation of compartment model parameters
q = a(1);
m = a(2);
alpha = 1-1/exp(q);
beta = 1-(1-alpha)/m/tau;
Veff = (1-beta)*V;
Qeff = (1-alpha)*Q;
tauEff = Veff/Qeff;

% Print on the screen 
fprintf('alpha: %f \n', alpha);
fprintf('beta:  %f \n', beta);
fprintf('R2 coefficient: %f \n', R2);

%% NLS solution                                                         
% [TODO]

% Comparison
fprintf('A) Compartment model: %f - CSTR: %f [kmol/m3]\n', CAout, CAcstr);
fprintf('X) Compartment model: %f - CSTR: %f [kmol/m3]\n', 1-CAout/CAin, 1-CAcstr/CAin);


%% Non linear System (NLS)
function F = CSTR(C)

    % [TODO]

end