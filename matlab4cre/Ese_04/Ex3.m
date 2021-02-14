close all; clear variables;

%% Experimental data
Nexp = 9;
T =     [ 300        350        350        400        450        450        500        500        550        ]; % [K]
kappa = [ 2.5197e-04 3.4984e-02 3.8545e-02 1.5634e+00 2.9993e+01 2.5964e+01 3.3637e+02 2.9365e+02 2.1903e+03 ]; % [1/min]
      

%% Hypothesis: kappa=A*T^n*exp(-E/RT)
% Linearization: y = a0+a1*x1+a2*x2
%                y = ln(k), a0=ln(A), a1=n, a2=-E/R, x1=ln(T), x2=1/T
% [TODO]

% Recover original parameters
% [TODO]

% Statistics
Ymod = X*a;
SSres = (Yexp-Ymod)'*(Yexp-Ymod);
SStot = (Yexp-mean(Yexp))'*(Yexp-mean(Yexp));
R2 = 1 - SSres/SStot;

% Plot
% [TODO]

% Print on the screen 
fprintf('Frequency factor (1/min):    %e \n', A);
fprintf('Temperature exponent:        %f \n', n);
fprintf('Activation energy (cal/mol): %f \n', E);
fprintf('R2 coefficient:              %f \n', R2);
