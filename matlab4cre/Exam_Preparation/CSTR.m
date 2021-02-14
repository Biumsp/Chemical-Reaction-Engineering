% REMEMBER
% Species are rows while reactions are columns:  
% es. Cin = [1 2 3]'; A = [1e8 2e9];

% -------------------------------------------------------------------------
% Data
% -------------------------------------------------------------------------

% data = struct;

% % You can modify from here on
% leg = {'CA'; 'CB'; 'CC'; 'T'};

% data.isothermal = 0; % Boolean (1 or 0)
% data.adiabatic  = 0; % Boolean (1 or 0)

% data.phase = ; % Phase ('L' of 'G')

% data.MW   = ;  % Absolute values are not relevant: only relative values are
% data.Cin = ;   % Inlet concentrations [mol/m^3]
% data.Tin = ;   % Inlet temperature [K]
% data.Qin  = ;  % Inlet volumetric flow rate [m^3/s]
% data.tau  = ;  % Residence time [s]
% data.V    = ;  % Volume [m^3]

% data.SC  = ;   % Matrix of stoichiometric coefficients
% data.RO  = 0;  % Matrix of reaction orders (if non-elementary)
% data.Ea_R  = ; % Activation Energies over R [K]
% data.A   = ;   % Frequency factors [1/s]
% data.k   = ;   % Kinetic constants (if constants) [proper SI units]
% data.reaction_rate = ; % handle to a function returning a row vector of
% % reaction rates like r = func(T, C) (if needed of course)

% data.Cp  = ;   % Vector of specific heat [J/mol/K]
% data.dCp = 0;  % Derivative of Cps with respect to T [J/mol/K^2]
% data.h0   = ;  % Reference enthalpies [J/mol]
% data.T0   = ;  % Reference T [K]

% data.MCpex = ; % Flow rate times Cp of external fluid [W/K]
% data.Tex  = ;  % Temperature of external fluid [K]
% data.beta = ;  % Beta factor (1 - exp(-L/l)) for heat exchange [-]

% data.Lex  = ;  % Work done on the system [W]
% data.Qex = ;   % External heat (if constant) (positive if removed) [W]

% % Beta Version ------------------------------------------------------------
% % % % data.reaction_rate = struct;
% % % % data.reaction_rate.b  = ;      % stoichiometric coefficient species B
% % % % data.reaction_rate.GammaA = ;  % diffusion coefficient [m2/s]
% % % % data.reaction_rate.GammaB = ;  % diffusion coefficient [m2/s]
% % % % data.reaction_rate.fL  = ;     % liquid fraction [-]
% % % % data.reaction_rate.HA  = ;     % Henry's constant [Pa.m3/kmol]
% % % % data.reaction_rate.a   = ;     % interfacial area per unti ov volume [m2/m3]
% % % % data.reaction_rate.KLa = ;     % mass transfer coefficient (liquid side) [1/s]
% % % % data.reaction_rate.KGa = ;     % mass transfer coefficient (gas side) [kmol/m3/s/Pa]
% % % % data.reaction_rate.kl  = ;     % kinetic constant (liquid volume basis) [m6/kmol^2/s]
% % Beta Version ------------------------------------------------------------

% % -------------------------------------------------------------------------
% % Solution
% % -------------------------------------------------------------------------
% 
% inlet = [data.Cin' data.Tin]'; % MUST be a column vector
% 
% [t, var] = ode23s(@CSTR, [0 0.001*data.tau], inlet, [], data);
% 
% T = var(:, end);
% C = var(:, 1:length(inlet)-1);
% 
% 
% % -------------------------------------------------------------------------
% % Graphical Post Processing
% % -------------------------------------------------------------------------
% 
% for ii = 1:length(leg)
%     
%     if ii == length(leg)
%         fprintf('T = %f [°C] = %f [K] \n', T(end) - 273.15, T(end));
%     else
%         fprintf('%s = %f [mol/m^3] = %f [kmol/m^3] =  \n', leg{ii}, C(end, ii), C(end, ii)/1e3);
%     end
%     
% end
% 
% yyaxis left
% plot(t, C)
% yyaxis right
% plot(t, T)
% legend(leg)
% legend('Location','best')


function yy = CSTR(t, var, data)

    % ---------------------------------------------------------------------
    % Recover the variables
    % ---------------------------------------------------------------------

    T = var(end);
    C = var(1:length(var)-1);
    
    Cin = data.Cin;
    Tin = data.Tin;
    SC = data.SC;
    RO = data.RO;
    Ea_R = data.Ea_R;
    A  = data.A;
    k = data.k;
    Cp = data.Cp;
    isothermal = data.isothermal;
    adiabatic = data.adiabatic;
    dCp = data.dCp;
    beta = data.beta;
    MCpex = data.MexCpex;
    Tex = data.Tex;
    Qin = data.Qin;
    tau = data.tau;
    MW  = data.MW;
    V = data.V;
    h0 = data.h0;
    T0 = data.T0;
    Lex = data.Lex;
    phase = data.phase;
    Qex = data.Qex;
    
    if and(Qin == 0, V == 0)
        V = tau;
        Qin = 1;
    elseif V == 0
        V = tau*Qin;
    elseif Qin == 0
        Qin = V/tau;
    end
    
    if not(reaction_rate)
        if data.RO == 0
            [rows, cols] = size(SC);
            RO = zeros(rows, cols);
            for ii = 1:rows
                for jj = 1:cols                
                    RO(ii, jj) = abs(min(0, SC(ii, jj)));
                end
            end
        end
    end
    
    if data.dCp == 0
        dCp = zeros(length(C), 1);
    end
    
    % ---------------------------------------------------------------------
    % Check if the matrices are correctly shaped
    % ---------------------------------------------------------------------
    
    size_C = size(C);
    if size_C(2) > 1
        disp('Error (1)')
    end
    size_SC = size(SC);
    size_RO = size(RO);
    if or(size_RO(1) ~= size_SC(1), size_RO(2) ~= size_SC(2))
        disp('Error (2)')
    end
    if not(isothermal)
        size_Cp = size(Cp);
        if or(size_Cp(2) > 1, size_Cp(1) ~= size_C(1))
            disp('Error (3)')
        end
        size_dCp = size(dCp);
        if or(size_dCp(2) > 1, size_dCp(1) ~= size_C(1))
            disp('Error (4)')
        end
    end
    size_Cin = size(Cin);
    if or(size_Cin(2) > 1, size_Cin(1) ~= size_C(1))
        disp('Error (5)')
    end
    if not(reaction_rate)
        if k == 0
            size_Ea = size(Ea_R);
            if size_Ea(1) > 1
                disp('Error (6)')
            end
            size_A = size(A);
            if or(size_A(1) > 1, size_A(2) ~= size_Ea(2))
                disp('Error (7)')
            end
        else
            Ea_R = zeros(1, size_SC(2));
            A    = zeros(1, size_SC(2));
        end
    else
        Ea_R = zeros(1, size_SC(2));
        A    = zeros(1, size_SC(2));
    end
    size_MW = size(MW);
    if size_MW(2) > 1
        disp('Error (8)')
    end
    if size_MW(1) ~= size_C(1)
        equimolar = 1;
        for ii = 1:length(A)
            if sum(SC(:, ii)) ~= 0
                disp('Error (9)')
                equimolar = 0;
            end
        end
        if equimolar
            MW = ones(length(C), 1);
        end
    end
    if not(isothermal)
        size_h0 = size(h0);
        if or(size_h0(2) > 1, size_h0(1) ~= size_C(1))
            disp('Error (10)')
        end
    end
    
    % ---------------------------------------------------------------------
    % Balances for species
    % ---------------------------------------------------------------------
    
    Ctot = sum(C);
    x    = C/Ctot;
    xin  = Cin/sum(Cin);
    
    if phase == 'L'
        gamma   = 1;
        MWin_MW = 1;
    elseif phase == 'G'
        gamma   = Tin/T; 
        MWin_MW = (xin'*MW)/(x'*MW);
    else
        disp('Error: data.phase should be either L or G')
    end
    
    if not(reaction_rate)
        k  = A.*exp(-Ea_R/T);

        r = zeros(1, length(A));
        for ii = 1:length(A)
            r(ii) = k(ii)*prod(C.^RO(:, ii));
        end
    else
        r = reation_rate(T, C);
    end
    
    dC = Qin*(Cin - C/gamma*MWin_MW)/V + SC*r';
    
    % ---------------------------------------------------------------------
    % Balance for temperature
    % ---------------------------------------------------------------------
    
    if isothermal
        dT = 0;
    else
        
        Fin = Qin*Cin;
        F   = Qin/gamma*MWin_MW*C;
        hin = h0 + Cp*(Tin - T0);
        h   = h0 + Cp*(T   - T0);
        
        Cp_mix  = x'*Cp;
        dCp_mix = x'*dCp;
        dn      = V*sum(SC*r');
        
        if Qex == 0
            if not(adiabatic)  
                Qex = MCpex*(T - Tex)*beta;
            end
        end

        dT = (Fin'*hin - h'*(F + x*dn) + 8.314*T*dn - Qex + Lex);
        dT =  dT/V/Ctot/(Cp_mix + T*dCp_mix);
    
    end
    
    % ---------------------------------------------------------------------
    % Output vector
    % ---------------------------------------------------------------------
    
    yy = [dC' dT]';

end