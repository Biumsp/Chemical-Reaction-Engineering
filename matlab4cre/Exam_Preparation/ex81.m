close all
clear variables

% -------------------------------------------------------------------------
% Data
% -------------------------------------------------------------------------

data = struct;

% You can modify from here on
leg = {'CA'; 'CB'; 'CC'; 'T'; 'P'; 'Tex'};

data.isothermal      = 1; % Boolean (1 or 0)
data.adiabatic       = 0; % Boolean (1 or 0)
data.isoperibolic    = 0; % Boolean (1 or 0)
data.isobaric        = 1; % Boolean (1 or 0)
data.countercurrent  = 0; % Boolean (1 or 0)
data.packed_bed      = 0; % Boolean (1 or 0)

data.phase = 'G'; % Phase ('L' of 'G')

data.Fin  = [6e4/3600 0 0 0]';  % Inlet molar flow rates [mol/s]
data.Tin  = 1033; % [K]
data.Pin  = 101325;

data.SC    = [-2 -1; 1 -1; 1 1; 0 1];   % Matrix of stoichiometric coefficients

data.rate_fcn = @rates; % handle to a function returning a row vector of
% reaction rates like r = func(T, C) (returns [mol/m^3/s] or [mol/kg_cat/s])

% -------------------------------------------------------------------------
% Solution
% -------------------------------------------------------------------------

inlet = [data.Fin' data.Tin data.Pin 0]';

V = 1.6;
[V, var] = ode23s(@PFR, [0 V], inlet, [], data);

Tex = var(:, end);
P   = var(:, end-1);
T   = var(:, end-2);
F   = var(:, 1:length(inlet)-3);

% -------------------------------------------------------------------------
% Results
% -------------------------------------------------------------------------

 try 
    P = data.pressure_fcn(data, T, F, V);
 catch
 end
 
X = (data.Fin(1) - F(:, 1))/data.Fin(1);
index = find(X >= 0.5, 1);
fprintf('V(X = 0.5) = %f [m^3]\n', V(index))
 
% -------------------------------------------------------------------------
% Graphical Post Processing
% -------------------------------------------------------------------------

yyaxis left
plot(V, F)
yyaxis right
plot(V, X)
legend(leg{1:end-3})
legend('Location','best')

function r = rates(~, C, ~)

    Keq1 = 0.31;
    Keq2 = 0.48;
    
    k1f = 7e5/36e5;
    k1b = k1f/Keq1;
    
    k2f = 4e5/36e5;
    k2b = k2f/Keq2;
    
    r1 = k1f*C(1)^2 - k1b*C(2)*C(3);
    r2 = k2f*C(1)*C(2) - k2b*C(3)*C(4);
    r = [r1 r2];

end