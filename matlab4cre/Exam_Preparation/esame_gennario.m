% -------------------------------------------------------------------------
% Data
% -------------------------------------------------------------------------

% Species = [O2 C2H4 C2H4O CO2 H2O] = [A B C D E];
SC = [-1 -1; 
      -2 -1/3; 
       2 0; 
       0 2/3; 
       0 2/3];
   
MWc = 44e-3;
Prod = 100/3600/MWc; % [mol/s]

Fin  = [0.06 0.94 0 0 0]';

Tin  = 549;
Pin  = 9.65*101325;

k    = [0.255 0.361]; % [boh]

rho_cat = 1800;
Gamma = 4.64e-6;
radius = 75e-6;

phi = radius/3*sqrt(k/Gamma);
eta = tanh(phi)./phi;

k = k.*eta/rho_cat;

opts = odeset('RelTol', 1e-7, 'AbsTol', 1e-10);
[W, F] = ode23s(@funk, [0 100], Fin, opts, Pin, Tin, SC, k);

FA = F(:, 1);
X = (Fin(1)-FA)/Fin(1);
index = find(X >= 0.4, 1);

F = F(1:index, :);
Fin = Fin*Prod/F(end,3);

opts = odeset('RelTol', 1e-7, 'AbsTol', 1e-10);
[W, F] = ode23s(@funk, [0 1000], Fin, opts, Pin, Tin, SC, k);

FA = F(:, 1);
X = (Fin(1)-FA)/Fin(1);
index = find(X >= 0.4, 1);

F = F(1:index, :);
W = W(1:index);

figure
plot(W, F)

fprintf('Wcat = %f [kg]\n', W(end))



function dF = funk(W, F, Pin, Tin, SC, k)

    Ctot = Pin/8.314/Tin;
    
    CA = F(1)/sum(F)*Ctot;
    r = k*CA;
    
    dF = SC*r';    

end
