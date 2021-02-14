% REMEMBER
% Species are rows while reactions are columns:  
% es. Cin = [1 2 3]'; A = [1e8 2e9];

% % -------------------------------------------------------------------------
% % Data
% % -------------------------------------------------------------------------

% data = struct;


% % You can modify from here on
% leg = {'FA'; 'FB'; 'T'; 'P'; 'Tex'};

% data.isothermal      = 0; % Boolean (1 or 0)
% data.adiabatic       = 0; % Boolean (1 or 0)
% data.isoperibolic    = 0; % Boolean (1 or 0)
% data.isobaric        = 0; % Boolean (1 or 0)
% data.countercurrent  = 0; % Boolean (1 or 0)
% data.packed_bed      = 0; % Boolean (1 or 0)

% data.phase = ; % Phase ('L' of 'G')
% data.rho   = ; % Only if 'L' [kg/m^3]

% data.rho_bed = ; % Catalytic bed density [kg_catalyst/m^3_reactor]
% data.rho_cat = ; % Catalyst density [kg_cat/m^3_cat]
% data.void    = ; % Void fraction of bed [-] (1 - data.void) is the fraction occupied by the catalyst

% data.MW   = ;  % Absolute values are not relevant: only relative values are
% data.Fin  = ;  % Inlet molar flow rates [mol/s]
% data.Tin  = ;  % Inlet temperature [K]
% data.Pin  = ;  % Inlet pressure [Pa]

% data.SC    = ;   % Matrix of stoichiometric coefficients
% data.RO    = 0;  % Matrix of reaction orders (if non-elementary)
% data.Ea_R  = ;   % Activation Energies over R [K]
% data.A     = ;   % Frequency factors [1/s]
% data.k     = ;   % Kinetic constants (if constants) [proper SI units]
% data.rate_fcn = ; % handle to a function returning a row vector of
% % reaction rates like r = func(T, C, data)(returns [mol/m^3/s] or [mol/kg_cat/s])
% data.pressure_fcn  = ; % handle to a function that returns the P (if any) in Pa
% % and take as arguments (data, T, F, V): must handle vector arguments

% data.Cpmix = ;  % If constant [J/mol/K]
% data.Cp    = ;  % Vector of specific heat [J/mol/K]
% data.Dh0r  = ; % Reaction enthalpies [J/mol]
% data.h0    = ;  % Reference enthalpies [J/mol]
% data.T0    = ;  % Reference T [K]

% data.U     = ; % Global heat exchange coefficients [W/m^2/K]
% data.MCpex = ; % Flow rate times Cp of external fluid [W/K]
% data.Tex   = ; % Temperature of external fluid (if isoperibolic)[K]
% data.a     = ; % Surface per unit volume [1/m]
% data.Cross_S_ext = ;% External cross section [m^2]
% data.Cross_S = ; % Internal Cross Section [m^2]

% data.d = ; % Internal EQUIVALENT diameter [m]
% data.mi = ; % Dinamic viscosity [Pa*s]
% data.ni = ; % Kinematic viscosity [m^2/s]

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
% %
% % If ISOThermal but not adiabatic: if you need the Tex profile you get
% like this (maybe it's better to solve it analytically
% % Tex = (Qr - U*a*T)/U/a;


% % -------------------------------------------------------------------------
% % Solution
% % -------------------------------------------------------------------------
% 
% inlet = [data.Fin' data.Tin data.Pin data.Tex]';

% V = ;
% [V, var] = ode23s(@PFR, [0 V], inlet, [], data);
% 
% Tex = var(:, end);
% P   = var(:, end-1);
% T   = var(:, end-2);
% F   = var(:, 1:length(inlet)-3);

% % -------------------------------------------------------------------------
% % Results
% % -------------------------------------------------------------------------

%  try 
%     P = data.pressure_fcn(data, T, F, V);
%  catch
%  end

% % -------------------------------------------------------------------------
% % Graphical Post Processing
% % -------------------------------------------------------------------------
% 
% yyaxis left
% plot(V, F)
% yyaxis right
% plot(V, T)
% legend(leg{1:end-3})
% legend('Location','best')


function yy = PFR(V, var, data)

    % ---------------------------------------------------------------------
    % Recover the variables
    % ---------------------------------------------------------------------
    
    isobaric       = data.isobaric;
    isothermal     = data.isothermal;
    adiabatic      = data.adiabatic;
    isoperibolic   = data.isoperibolic;
    countercurrent = data.countercurrent;
    packed_bed     = data.packed_bed;
    
    
    Tex = var(end);
    P   = var(end-1);
    T   = var(end-2);
    F   = var(1:length(var)-3);
    
    try
        P = data.pressure_fcn(data, T, F, V);
    catch
    end
    if isoperibolic 
        Tex = data.Tex;
    end
           
    Fin = data.Fin;
    phase = data.phase;
    
    SC = data.SC;
    [Ns, Nr] = size(SC);
    try 
        rate = data.rate_fcn;
    catch
        try
            k = data.k;
            Nr = length(k);
        catch
            Ea_R = data.Ea_R;
            A  = data.A;
            
            try
                RO = data.RO;
            catch
                [Ns, Nr] = size(SC);
                RO = zeros(Ns, Nr);
                for ii = 1:Ns
                    for jj = 1:Nr                
                        RO(ii, jj) = abs(min(0, SC(ii, jj)));
                    end
                end
            end
            size_RO = size(RO);
            if or(size_RO(1) ~= Ns, size_RO(2) ~= Nr)
                disp('Error (1)')
            end
            
            size_Ea = size(Ea_R);
            if size_Ea(1) > 1
                disp('Error (2)')
            end
            size_A = size(A);
            if or(size_A(1) > 1, size_A(2) ~= size_Ea(2))
                disp('Error (3)')
            end
        end
    end
    
    [Ns, cols] = size(F);
    if cols > 1
        disp('Error (4)')
    end
    
    if not(adiabatic || isothermal)
        try 
            data.Qex;
        catch
            U     = data.U;
            a     = data.a;
            if not(isoperibolic) 
                Cross_S = data.Cross_S;
                Cross_S_ext = data.Cross_S_ext;
                MCpex = data.MCpex;
            end
        end
    end
    if not(isobaric)
        Cross_S = data.Cross_S;
    end
    try
        MW  = data.MW;
        size_MW = size(MW);
        if size_MW(2) > 1
            disp('Error (7)')
        end
        if size_MW(1) ~= Ns
            equimolar = 1;
            for ii = 1:Nr
                if sum(SC(:, ii)) ~= 0
                    disp('Error (8)')
                    equimolar = 0;
                end
            end
            if equimolar
                MW = ones(Ns, 1);
            end
        end
    catch
        if not(isobaric)
            try 
                data.pressure_fcn
            catch
                disp('Error: missing MWs')
            end
        end
    end
    if not(isothermal)
        try 
            data.Dh0r;
        catch
            h0  = data.h0;
            T0  = data.T0;
            size_h0 = size(h0);
            if or(size_h0(2) > 1, size_h0(1) ~= Ns)
                disp('Error (9)')
            end
            if exists(data.Cp, 'var')
                Cp = data.Cp;
                size_Cp = size(Cp);
                if or(size_Cp(2) > 1, size_Cp(1) ~= Ns)
                    disp('Error (5)')
                end
            else
                Cp = zeros(Ns, 1);
            end
        end
    end
       
    % ---------------------------------------------------------------------
    % Check if the matrices are correctly shaped
    % ---------------------------------------------------------------------
    
    
    size_Fin = size(Fin);
    if or(size_Fin(2) > 1, size_Fin(1) ~= Ns)
        disp('Error (6)')
    end
    
    % ---------------------------------------------------------------------
    % Balances for species
    % ---------------------------------------------------------------------
    
    Ftot = sum(F);    
    Ctot = P/8.314/T;
    C = F/Ftot*Ctot;
    x = F/Ftot;

    try
        r = rate(T, C, data);
    catch
        k  = A.*exp(-Ea_R/T);

        r = zeros(1, Nr);
        for ii = 1:Nr
            r(ii) = k(ii)*prod(C.^RO(:, ii));
        end        
    end
        
    dF = SC*r';
    
    % ---------------------------------------------------------------------
    % Balance for temperature
    % ---------------------------------------------------------------------
    
    if isothermal
        dT = 0;
    else
        try
            Dh0r = data.Dh0r;
        catch
            h   = h0 + Cp*(T   - T0);
            Dh0r = h'*SC; 
        end
        
        Qr = - Dh0r*r';
        
        try
            Cpmix = data.Cpmix;
        catch
            Cpmix  = x'*Cp;
        end
        try 
            Qex = data.Qex;
        catch
            if adiabatic
                Qex = 0;
            else
                Qex = U*a*(T - Tex);
            end
        end
        
        dT = (Qr - Qex)/Ftot/Cpmix;
    
    end
    
    % ---------------------------------------------------------------------
    % Balance for pressure
    % ---------------------------------------------------------------------
    
    if isobaric
        dP = 0;
    else
        try 
            data.pressure_fcn;
            dP = 0;
        catch
            v   = Ftot/Ctot/Cross_S;
            if phase == 'G'
                rho = x'*MW*P/8.314/T;
            elseif phase == 'L'
                rho = data.rho;
            end
            d  = data.d;
            try
                mi = data.mi;
            catch
                mi = data.ni*rho;
            end                
            Re = rho*v*d/mi;
            f  = 0.079*Re^-0.25;
            dP = - rho*v^2*(dT/T + a/2*f)/(1 - rho*v^2/P);
        end
    end
    
    % ---------------------------------------------------------------------
    % Balance for external temperature
    % ---------------------------------------------------------------------
    
    if not(isothermal)
        try
            try 
                Qex = data.Qex;
            catch
                Qex = U*a*(T - Tex);
            end
            dTex = (-1)^(countercurrent)*Cross_S_ext/Cross_S*Qex/MCpex;
        catch
            dTex = 0;
        end
    else
        dTex = 0;
    end
    
    % ---------------------------------------------------------------------
    % Output vector
    % ---------------------------------------------------------------------
   
    yy = [dF' dT dP dTex]';
    
    if packed_bed
        try 
            rho_bed = data.rho_bed;
        catch
            rho_cat = data.rho_cat;
            void    = data.void;
            rho_bed = (1 - void)*rho_cat;
        end
        
        yy = yy*rho_bed;
    end
    
end