close all; clear variables;

%% Experimental data
Nexp = 9;
CA =  [ 10 8.0262 6.5575 5.4393 4.5708 3.8847 3.3346 2.8877 2.5203 ];   % mol/l
t  =  [ 0  5      10     15     20     25     30     35     40     ];   % min

plot(t, CA)
%% Calculation of rA=-dCA/dt
dCAdt(1) = (CA(2) - CA(1))/(t(2) - t(1));
for i = 2:Nexp-1
   dCAdt(i) = (CA(i + 1) - CA(i-1))/(t(i + 1) - t(i-1)); 
end
dCAdt(Nexp) = (CA(Nexp) - CA(Nexp - 1))/(t(Nexp) - t(Nexp-1));
rA = -dCAdt;

% %% Plot of raw data
% figure; title('Experimental data');
% scatter(CA,rA, 'o'); 
% xlabel('C_{A} [mol/l]'); ylabel('r_A [mol/l/min]');


% %% Hypothesis: rA = k*CA^n
% % Linearization: y = q + m*x
% %                y = ln(r), q = ln(k), m = n, x = ln(CA)
% y = log(rA);
% x = log(CA);

%% Linear regression analysis
% Yexp = a0 + a1*Xexp
% Yexp = ln(r), a0 = ln(k), a1 = n, Xexp = [1, ln(CA0)]
Yexp = log(rA)';
Xexp = [ones(Nexp, 1) log(CA)'];
A = Xexp'*Xexp;
b = Xexp'*Yexp;

% A*a = b
a = A\b;

% Recover original parameters
n = a(2); % Matlab counts from 1
k = exp(a(1));
fprintf('n: %f, k: %f\n', n, k)

% Statistics
Ymod = Xexp*a;
SSres = (Yexp - Ymod)'*(Yexp - Ymod);
YexpAvg = mean(Yexp);
SStot = (Yexp - YexpAvg)'*(Yexp - YexpAvg);
R2 = 1- SSres/SStot;
fprintf('R2 = %f\n', R2);

% Plot
CAmod = 2.52:0.1:10;
xmod = log(CAmod);
ymod = log(k) + n*log(CAmod);
figure;
plot(x, y, 'o'); hold on;
plot(xmod, ymod, '-')


% Print on the screen 
% [TODO]
