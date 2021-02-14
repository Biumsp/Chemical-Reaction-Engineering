close all; clear variables;
global CA0;

%% Experimental data
Nexp = 9;
CA =  [ 10 8.0262 6.5575 5.4393 4.5708 3.8847 3.3346 2.8877 2.5203 ];   % mol/l
t  =  [ 0  5      10     15     20     25     30     35     40     ];   % min


%% Approach 1: linear regression analysis (squential procedure)
% Reaction orders to be tested
n = [1.05 1.10 1.15 1.20 1.25 1.30 1.35 1.40 1.45 1.50 1.55 1.60];

% Hypothesis: rA = k*CA^n
% Integration:   CA^(1-n) = CA0^(1-n) -k*(1-n)*t
% Linearization: y = q + m*x
%                y = CA^(1-n)-CA0^(1-n), q = 0, m = -k*(1-n), x = t
%                y = a1*x

for j=1:length(n)
    
    CA0 = CA(1);
    Yexp =  CA.^(1-n(j))-CA0.^(1-n(j));
    Yexp = Yexp';
    Xexp = [t' ];
    
    a = (Xexp'*Xexp)\(Xexp'*Yexp);
    k = a(1)/(1 - n(j));
    Ymod = Xexp*a;
    
    SSres = (Yexp - Ymod)'*(Yexp - Ymod);
    YexpAvg = mean(Yexp);
    SStot = (Yexp - YexpAvg)'*(Yexp - YexpAvg);
    R2 = 1- SSres/SStot;
    fprintf('n = %f, k = %f, R2 = %f\n', n(j), k, R2);
    
end


%% Approach 2: non-linear regression analysis

% Non linear regression
% [TODO]

% Results
% [TODO]

%% Function to be regressed
function y = nonLinearModel(a,x)

    global CA0;
    
    % [TODO]

end
