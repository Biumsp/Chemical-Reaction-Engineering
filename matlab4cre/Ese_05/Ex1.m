close all; clear variables;

%% Experimental data
N = 21;
t = [0 1 2 2.2 2.4 2.6 2.8 3 3.2 3.4 3.6 3.8 4 4.2 4.4 4.6 4.8 5 6 7 8]';
E = [0 0 0 0.010640797 0.042563186 0.085126373 0.127689559 0.255379118 ...
     0.681010982 0.851263728 0.851263728 0.766137355 0.427055382 ...
     0.284703588 0.170252746 0.127689559 0.085126373 0.042563186 ...
     0.021281593 0	0 ];

%% User data
k = 0.1;        % kinetic constant [m3/kmol/s]
L = 30;         % reactor length [m]
D = 0.04;       % reactor diameter [m]
v = 8.1;        % velocity [m/s]
CAin = 10;      % inlet concentration [kmol/m3]

%% Preliminary calculations
A = pi*D^2/4;   % cross section area [m2]
Q = v*A;        % volumetric flow rate [m3/s]
V = A*L;        % volume [m3]

%% Ideal PFR
tau = V/Q;
CApfr = CAin/(1+k*tau*CAin);

%% Analysis of RTD     


%% Normalization of E (to be sure that the area below is 1)
% We compute the moment of order zero 
% (just the measure (should be 1, but it won't because of errors)) and divide by it
% So we evaluate a numerical integral

m0 = n_moments(E, t, 0);
E = E/m0;

% Plot the RTD
figure;
plot(t,E, '-o');
xlabel('time [s]');
ylabel('E [1/s]');

%% Mean residence time: tm = int(t*E*dt)

% Mean Residence Time (tm): moment of ordr 1
tm = n_moments(E, t, 1);

%% Variance: sigma2 = m2 - int(t2*E*dt)
m2 = n_moments(E, t, 2);
sigma2 = m2 - tm^2;

% Adimensionalized sigma^2
sigmaTeta2 = sigma2/tm^2;

fprintf('From RTD: tm=%f [s] - sigma2=%f [s2] - sigmaTeta2=%f \n', ...
         tm, sigma2, sigmaTeta2);

%% Dispersion model: dispersion coefficient
% We use the approx for Pe big (We expect it to be big because the RTD
% looks like a Dirac's delta)
Pe = 2/sigmaTeta2;
GammaEff = L*v/Pe;

fprintf('DM: Pe=%f - GammaEff=%f [m2/s]\n', Pe, GammaEff);


%% Tanks in Series model
n = 1/sigmaTeta2;

% Plus
n_plus = ceil(n);
tau_plus = tm/n_plus;
CA_plus(1) = CAin;
for i = 2:n_plus+1
    
    CA_plus(i) = (-1 + sqrt(1 + 4*k*tau_plus*CA_plus(i-1)))/(2*k*tau_plus);
    
end

% Minus
n_minus = floor(n);
tau_minus = tm/n_minus;
CA_minus(1) = CAin;
for i = 2:n_minus+1
    
    CA_minus(i) = (-1 + sqrt(1 + 4*k*tau_minus*CA_minus(i-1)))/(2*k*tau_minus);
    
end

% Interpolation
CA_tis = CA_minus(end) + (CA_plus(end) - CA_minus(end))/(n_plus-n_minus)*(n-n_minus);

%Print on the screen
fprintf('A) TIS: %f - PFR: %f [kmol/m3]\n', CA_tis(end), CApfr);
fprintf('X) TIS: %f - PFR: %f [kmol/m3]\n', 1-CA_tis/CAin, 1-CApfr/CAin);
