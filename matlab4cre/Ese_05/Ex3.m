close all; clear variables;

global k;            % [1/min]
global CAin;         % [kmol/m3]
global C2;
global C3;

%% Experimental data                                                       %
nExp = 10;
tExp = [0 0.5 1 1.5 2 2.5 3 3.5 4 5]';
FExp = [0.0421434 0.4344916 0.6686519 0.8122302 0.8811733 ...
        0.9367312 0.9622991 0.9786636 0.9827176 0.995]';

%% User data                                                               
V    = 0.7;    % total volume [l]   
Q    = 0.7;    % total volumetric flow rate [l/min]
CAin = 10;     % [kmol/m3]
k    = 5;      % kinetic constant [m3/kmol/min]
tau  = V/Q;    % nominal residence time [min]


%% Ideal CSTR and PFR solutions                                            
CAcstr = (-1+sqrt(1+4*k*tau*CAin))/(2*k*tau);
CApfr = CAin/(1+k*tau*CAin);


%% Analysis of experimental data                                           

figure; hold on;
plot(tExp,FExp, '-o');
plot(tExp,1-exp(-tExp/tau)/tau, '-x');
xlabel('time [min]'); ylabel('CDF');
legend('Real reactor', 'Ideal CSTR');

% Linear regression
Y = log(1 - FExp);
X = [ones(nExp, 1) tExp];
a = (X'*X)\(X'*Y);
C2 = exp(a(1));
C3 = -1/a(2);

% Statistics
SSres = (Y-X*a)'*(Y-X*a);
SStot = (Y-mean(Y))'*(Y-mean(Y));
R2 = 1 - SSres/SStot;

% Residence Time Distribution function
EExp = C2/C3*exp(-tExp/C3);

% Print on the screen 
fprintf('C2: %f (CSTR=1)\n',  C2);
fprintf('C3: %f (CSTR=%f)\n', C3, tau);
fprintf('R^2 coefficient: %f \n', R2);

%% Integration of SM 
% REMEMBER: in the integral we have to multiply E by the C_batch_react
% Which is the same of the PFR
% As the max value of t we should use the max t we have in the data, but
% here we can put any value, because we have the analytical solution. So we
% use 100*residenceTime, to be sure.
%CAbatch_times_E = @(t) (CAin/(1+k*tau*CAin)).*C2/C3*exp(-t/C3);
CAout = integral(@CAbatch_Times_E, 0, 100*tau);

%% Comparison
fprintf('A) Segregated model: %f - CSTR: %f - PFR: %f [kmol/m3]\n', ...
         CAout, CAcstr, CApfr);
fprintf('X) Segregated model: %f - CSTR: %f - PFR: %f \n', ...
         1-CAout/CAin, 1-CAcstr/CAin, 1-CApfr/CAin);

     
%% Function to be integrated
function F = CAbatch_Times_E(t)

    global k
    global CAin C2 C3 tau

    F = (CAin./(1+k.*tau.*CAin)).*C2./C3.*exp(-t./C3);

end
