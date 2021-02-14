

close all
clear variables
clc

% -------------------------------------------------------------------------
% Data
% -------------------------------------------------------------------------

leg = {'CA'; 'CB'; 'CC'; 'T'};

data = struct;

data.Cin = [55e3 0 0]';
data.Tin = 80+273.15;
data.Pin = 1;
data.P   = 1;
data.SC  = [-1 -1; 1 0; 0 1];
data.RO  = 0;
data.Ea_R= [2e4 1.8e4]/1.987;
data.A   = [3e14, 2e13];
data.Cp  = 0;
data.isothermal = 1;
data.adiabatic  = 0;
data.dCp  = 0; 
data.beta = 0;
data.MexCpex = 0; 
data.Tex  = 0;
data.Qin  = 0;
data.tau  = 60;
data.MW   = 0;
data.V    = 0;
data.h0   = 0;
data.T0   = 0;
data.Lex  = 0;
data.phase = 'G';

% -------------------------------------------------------------------------
% Solution
% -------------------------------------------------------------------------

inlet = [data.Cin' data.Tin]'; % MUST be a column vector

[t, var] = ode23s(@CSTR, [0 0.001*data.tau], inlet, [], data);

T = var(:, end);
C = var(:, 1:length(inlet)-1);


% -------------------------------------------------------------------------
% Graphical Post Processing
% -------------------------------------------------------------------------

for ii = 1:length(leg)
    
    if ii == length(leg)
        fprintf('T = %f [°C] = %f [K] \n', T(end) - 273.15, T(end));
    else
        fprintf('%s = %f [mol/m^3] = %f [kmol/m^3] =  \n', leg{ii}, C(end, ii), C(end, ii)/1e3);
    end
    
end

yyaxis left
plot(t, C)
yyaxis right
plot(t, T)
legend(leg)
legend('Location','best')