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

%% User data                                                               %
V    = 0.7;    % total volume [l]   
Q    = 0.7;    % total volumetric flow rate [l/min]
CAin = 10;     % [kmol/m3]
k    = 5;      % kinetic constant [m3/kmol/min]
tau  = V/Q;    % nominal residence time [min]


%% Ideal CSTR and PFR solutions                                            %
CAcstr = (-1+sqrt(1+4*k*tau*CAin))/(2*k*tau);
CApfr = CAin/(1+k*tau*CAin);

% ------------------------------------------------------------------------%
% Fitting of F using non-linear regression analysis                       %
% This was already done in Ex. 3 using linear regression analysis         %
% Here we simply try an alternative method, but results are the same!     %
%-------------------------------------------------------------------------%

firstGuess = [1 tau];
nlm = fitnlm(tExp,FExp,@nonLinearModel,firstGuess);
C2   = nlm.Coefficients.Estimate(1);
C3   = nlm.Coefficients.Estimate(2);

% Quality evaluation
FModel = 1-C2*exp(-tExp/C3);
SSres = sum( (FExp-FModel).^2  );
SStot = sum( (FExp-mean(FExp)).^2 );
R2 = 1 - SSres/SStot;

% Residence Time Distribution function
EExp = C2/C3*exp(-tExp/C3);

% Print on the screen 
fprintf('C2: %f \n', C2);
fprintf('C3: %f \n', C3);
fprintf('R2 coefficient: %f \n', R2);


%% ODE solution
lambdaMax = 100*tau;
[lambda, CA] = ode23s(@ODEEquation, [lambdaMax 0], CAin);

CAout = CA(end);
plot(lambda, CA)

% Comparison
fprintf('A) Max. Mixedness: %f - CSTR: %f - PFR: %f [kmol/m3]\n', ...
         CAout, CAcstr, CApfr);
fprintf('X) Max. Mixedness: %f - CSTR: %f - PFR: %f \n', ...
         1-CAout/CAin, 1-CAcstr/CAin, 1-CApfr/CAin);
     
     
%% ODE system corresponding to the MM model
function dCA = ODEEquation(~,CA)

    global C3 CAin k
    
    RA = -k*CA^2;
    
    % E/(1-F) = 1/C3
    dCA = (1/C3)*(CA - CAin) - RA;

end


%% Function fitting the experimental CDF
function F = nonLinearModel(a,t)

    F = 1-a(1)*exp(-t/a(2));

end
