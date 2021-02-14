function [rIIII,E,MH,resL,resG,resK] = OverallRateOfChange(pA, CB, data)
 
    % resL,resG,resK are the resistences in L, G and kinetics

    %% Extract data from data struct
    b = data.b;             % stoichiometric coefficient species B
    GammaA = data.GammaA;   % diffusion coefficient [m2/s]
    GammaB = data.GammaB;   % diffusion coefficient [m2/s]
    fL = data.fL;           % liquid fraction [-]
    HA = data.HA;           % Henry's constant [Pa.m3/kmol]
    a = data.a;             % interfacial area per unti ov volume [m2/m3]
    KLa = data.KLa;         % mass transfer coefficient (liquid side) [1/s]
    KGa = data.KGa;         % mass transfer coefficient (gas side) [kmol/m3/s/Pa]
    kl = data.kl;           % kinetic constant (liquid volume basis) [m6/kmol^2/s]
    
    
    %% Preliminary calculations
    klStar = kl*CB;
    klPrime = klStar*CB;
    KL = KLa/a;
    MH = sqrt(klStar*GammaA*CB)/KL;
    resG = 1/KGa;
    resK = HA/(klPrime*fL);
    % for resL we need E

    %% Enhancement factor
    if (MH<1)               % Regions I and II
        
        % We approximate
        if (MH < 0.3)
            E = 1;
        else
            E = 1+MH^(2/3);
        end
        
        % We don't need any iterative procedure
        resL = HA/KLa/E;
        rIIII = pA/(resG+resl+resK);

    else                    % Region III
        
        % FGS
        pAI = pA;
        
        max_iter = 100;
        for i = 1:max_iter
            
            Ei = 1 + GammaB/GammaA*CB*HA/b/pAI;
            
            if (MH < 5*Ei)
                E = MH;
            else
                E = Ei;
            end
            
            resL = HA/KLa/E;
            rIIII = pA/(resG+resL+resK);
            rIIIIGas = KGa*(pA - pAI);
            
            if (abs(rIIIIGas - rIIII)/rIIII < 0.001) 
                break
            end
            
            pAI = pA - rIIII/KGa;
         
        end
        
    end
    
end
