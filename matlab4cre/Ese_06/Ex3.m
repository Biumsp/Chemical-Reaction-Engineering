%% Data struct
data.b=2;                        % stoichiometric coefficient species B
data.c=1;                        % stoichiometric coefficient species C
data.GammaA=1e-6/3600;           % diffusion coefficient [m2/s]
data.GammaB=1e-6/3600;           % diffusion coefficient [m2/s]
data.fL=0.98;                    % liquid fraction [-]
data.HA=1e5*1e3;                 % Henry's constant [Pa.m3/kmol]
data.a=20;                       % interfacial area per unti ov volume [m2/m3]
data.KLa=20/3600.;               % mass transfer coefficient (liquid side) [1/s]
data.KGa=0.01*(1e-3/3600);       % mass transfer coefficient (gas side) [kmol/m3/s/Pa]
data.kl=100e6/1e-6/3600.;        % kinetic constant (liquid volume basis) [m6/kmol^2/s]


%% Input data
T = 303;        % temperature [K]
pTot = 101325;  % total pressure (gas) [Pa]
pAin = 5000;    % inlet partial pressure of species A [Pa]
CBin = 0.00100; % inlet concentration of species B [kmol/m3]
CCin = 0.0;     % inlet concentration of species C [kmol/m3]
QL = 10/1000;   % liquid volumetric flow rate [m3/s]
Ctot = 0.301;   % total concentration liquid phase [kmol/m3]
QG = 0.4/1000;  % gas volumetric flow rate [m3/s]
L = 5;          % total length [m]
D = 0.40;       % diameter [m]


%% Counter-current Tower Reactor
FtotG = (pTot/8314/T)*QG;   % total molar flow rate (gas) [kmol/s]
FtotL = Ctot*QL;            % total molar flow rate (liquid) [kmol/s]
A = pi/4*D^2;               % cross section area [m2]
Vtot = A*L;                 % total volume [m3]


%% Shooting algorithm
maxiterations = 100;    % maximum number of iterations
alpha =0.10;            % underrelaxation factor
deltaPthreshold = 0.1;  % [Pa]

pAout = (1-0.80)*pAin;  % first guess solution (from co-current)
for i=1:maxiterations
    
    % ODE solution
    Yin = [pAout CBin CCin]';
    [V, Y] = ode45(@TowerCounterCurrent, [0 Vtot], Yin, [], pTot, Ctot, FtotG, FtotL, data);
    pA = Y(:,1);
    
    % Check solution
    error = pA(end) - pAin;
    
    % Print current status
    fprintf('Iteration: %d - pAGuess: %f Pa - pAin: %f Pa - Error: %f Pa \n', ...
             i, pA(end), pAin, error);
    
    % Exit the loop
    if (abs(error) < deltaPthreshold)
        break;
    end
         
    % New guess
    pAout = pAout - alpha*error;
    
end

z = V/A;        % axial coordinate [m]
pA = Y(:,1);    % partial pressure of A [Pa]
CB = Y(:,2);    % concentration of B [kmol/m3]
CC = Y(:,3);    % concentration of C [kmol/m3]


%% Post-processing

fprintf('A conversion: %f \n', 1-pA(1)/pAin);

figure;
xlabel('axial coordinate [m]');
hold on;
yyaxis left; ylabel('partial pressure of A [Pa]'); 
plot(z,pA);
yyaxis right; ylabel('concentration of B [mol/m3]'); 
plot(z,CB*1000);

figure;
xlabel('axial coordinate [m]');
hold on;
yyaxis left; ylabel('concentration of B [mol/m3]'); 
plot(z,CB*1000);
yyaxis right; ylabel('concentration of C [mol/m3]'); 
plot(z,CC*1000);


%% Reconstructing the overall rate of change and the enhancing factor
for i=1:size(pA,1)
    [rIIII(i),E(i)] = OverallRateOfChange(pA(i), CB(i), data);
end

figure;
xlabel('axial coordinate [m]');
hold on;
yyaxis left; ylabel('overall rate of change [mol/m3/hr]'); 
plot(z,rIIII*3600*1000);
yyaxis right; ylabel('enhancing factor'); 
plot(z,E);


%% ODE system (defined locally, requires at least MATLAB-2016b)
function dY = TowerCounterCurrent(V,Y, pTot, Ctot, FtotG, FtotL, data)

    pA = Y(1);      % partial pressure of A in gaseous phase [Pa]
    CB = Y(2);      % concentration of B in liquid phase [kmol/m3]
    
    [rIIII,E] = OverallRateOfChange(pA, CB, data);    % [kmol/m3/s]
        
    dpA =  pTot/FtotG*rIIII;                    % [Pa/m3]
    dCB = -Ctot/FtotL*data.b*rIIII;             % [kmol/m3/m3]
    dCC =  Ctot/FtotL*data.c*rIIII;             % [kmol/m3/m3]

    dY = [dpA dCB dCC]';

end
