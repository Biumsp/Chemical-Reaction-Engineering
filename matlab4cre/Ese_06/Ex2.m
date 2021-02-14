close all; clear variables;


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
CBin = 0.001;   % inlet concentration of species B [kmol/m3]
CCin = 0.0;     % inlet concentration of species C [kmol/m3]
QL = 10/1000;   % liquid volumetric flow rate [m3/s]
Ctot = 0.301;   % total concentration liquid phase [kmol/m3]
QG = 0.4/1000;  % gas volumetric flow rate [m3/s]
L = 5;          % total length [m]
D = 0.40;       % diameter [m]


%% Co-current Tower Reactor
FtotG = (pTot/8324/T)*QG;
FtotL = Ctot*QL;
A = pi/4*D^2;
Vtot = A*L;

%% ODE solution

ICs = [pAin, CBin, CCin]';
[V, Y] = ode45(@TowerCoCurrent, [0 Vtot], ICs, [], pTot, Ctot, FtotG, FtotL, data);

pA= Y(:, 1);
CB= Y(:, 2);
CC= Y(:, 3);
z = V/A;

%% Post-processing

fprintf('A conversion: %f \n', 1-pA(end)/pAin);

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

function dY = TowerCoCurrent(V,Y, pTot, Ctot, FtotG, FtotL, data)

    pA = Y(1);
    CB = Y(2);
    rIIII = OverallRateOfChange(pA, CB, data);
    
    dpA_over_dV = -pTot/FtotG*rIIII;
    dCB_over_dV = -Ctot/FtotL*data.b*rIIII;
    dCC_over_dV =  Ctot/FtotL*rIIII;
    
    dY = [dpA_over_dV, dCB_over_dV, dCC_over_dV]';

end
