clear variables

% -------------------------------------------------------------------------
% Data
% -------------------------------------------------------------------------

data = struct;


% You can modify from here on
leg = {'FA'; 'FB'; 'FC'; 'T'; 'P'; 'Tex'};

data.isothermal      = 1; % Boolean (1 or 0)
data.adiabatic       = 0; % Boolean (1 or 0)
data.isoperibolic    = 0; % Boolean (1 or 0)
data.isobaric        = 1; % Boolean (1 or 0)
data.countercurrent  = 0; % Boolean (1 or 0)
data.packed_bed      = 1; % Boolean (1 or 0)

data.phase = 'G'; % Phase ('L' of 'G')

data.rho_cat = 2000; % Catalyst density [kg_cat/m^3_cat]
data.void    = 0.4; % Void fraction of bed [-] (1 - data.void) is the fraction occupied by the catalyst

data.Tin  = 250+273.15;  % Inlet temperature [K]
data.Pin  = 500e3;  % Inlet pressure [Pa]
data.Fin  = [1 0 0]'*data.Pin/8.314/data.Tin*3;  % Inlet molar flow rates [mol/s]

data.SC    = [-1 1 2]';   % Matrix of stoichiometric coefficients
data.RO    = [ 2 0 0]';
data.Ea_R  = -log(0.05/2)/(1/data.Tin - 1/673.15);   % Activation Energies over R [K]
data.A     = 0.05/exp(-data.Ea_R/data.Tin);   % Frequency factors [1/s]

data.Cpmix = 15;  % If constant [J/mol/K]
data.Dh0r  = -800; % Reaction enthalpies [J/mol]

% -------------------------------------------------------------------------
% Solution
% -------------------------------------------------------------------------

inlet = [data.Fin' data.Tin data.Pin 0]';

V = 0.1;
[V, var] = ode23s(@PFR, [0 V], inlet, [], data);

Tex = var(:, end);
P   = var(:, end-1);
T   = var(:, end-2);
F   = var(:, 1:length(inlet)-3);

for ii = 1:length(inlet)-3
    C(:,ii) = F(:,ii)./sum(F,2).*P./T/8.314;
end

% -------------------------------------------------------------------------
% Results
% -------------------------------------------------------------------------

try 
P = data.pressure_fcn(data, T, F, V);
catch
end
 
X = (data.Fin(1) - F(:,1)) /data.Fin(1);
figure
plot(V, X)
disp(X(end))

% -------------------------------------------------------------------------
% Graphical Post Processing
% -------------------------------------------------------------------------

figure
yyaxis left
plot(V, C)
yyaxis right
% plot(V, T)
% legend(leg{1:end-3})
% legend('Location','best')