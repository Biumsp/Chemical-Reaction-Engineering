from .reaction_mechanism import ReactionMechanism
from .species import Species

class Mixture:

    name = None              
    mixtures = 0         
    R = 8.314462            # Gas constant [J/mol/K]
    RM = None               # Reaction mechanism
    comp = None             # Composition: list of species
    C = None                # List of concentrations [kmol/m^3]
    x = None                # Molar fractions [-]
    w = None                # Mass fractions [-]
    MW = None               # List of Molecular Weights [kmol/kg]
    MWmix = None            # Average Molecular Weight [kmol/kg]
    T = None                # Temperature [K]
    P = None                # Pressure [Pa]
    rho = None              # Density [kg/m^3]
    Ctot = None             # Total concentration [kmol/m^3]
    VFR = None              # Volumetric Flow Rate [m^3/s]
    MFR = None              # Mass Flow Rate [kg/s]
    phase = "G"             # Phase ('L' or 'G')
    Ns = 0                  # Number of species [-]
    Nr = 0                  # Number of reactions [-]

    variables =  {'name' : (int, str),
                  'R' : (float,),
                  'RM' : (ReactionMechanism,),
                  'comp' : (list, Species),
                  'C' : (list,),
                  'x' : (list,),
                  'w' : (list,),
                  'MW' : (list, float, int),
                  'MWmix' : (float, int),
                  'T' : (float, int),
                  'P' : (float, int),
                  'rho' : (float, int),
                  'Ctot' : (float, int),
                  'VFR' : (float, int),
                  'MFR' : (float, int),
                  'phase' : (str,),
                  'Ns' : (int,),
                  'Nr' : (int,)}

    def __init__(self, name = None):
        if name:
            self.name = name
        else:
            self.name = Mixture.mixtures
            Mixture.mixtures += 1

    def __str__(self):
        return f"Mixture '{self.name}'"
        
    def __repr__(self):
        return "<class 'Mixture'>"
    
    def __setattr__(self, name, value):
        self._check_validity(name, value)
        self._check_coherence(name, value)
        
        super.__setattr__(self, name, value)
        
        self._update(name, value)
        
    def _check_validity(self, name, value):
        if name not in self.variables:
            raise ValueError(f"{type(self).__name__} object has no attribute '{name}'")  
        elif type(value) not in self.variables[name]:
            raise TypeError(f"Wrong type {type(value)} for attribute '{name}'")
            
        return True
            
    def _check_coherence(self, name, value):
        
        if self.Ns:
            if name == "x" or name == "w" or name == "comp":
                    if self.Ns != len(value):
                        raise ValueError(f"The length of the variable {name} doesn't match the number of species already evaluated")

            if name == "RM":
                if value.Ns != self.Ns:
                        raise ValueError(f"The length of the variable {name} doesn't match the number of species already evaluated")

        if name == "MFR" and self.MFR:
            if abs((value - self.MFR)/self.MFR) > 0.0001:
                raise ValueError("The Mass Flow Rate provided is different from the value calculated from the density and the Volumetric Flow Rate")
        elif name == "rho" and self.rho:
            if abs((value - self.rho)/self.rho) > 0.0001:
                raise ValueError("The density provided is different from the value calculated from the Mass Flow Rate and the Volumetric Flow Rate")
        elif name == "VFR" and self.VFR:
            if abs((value - self.VFR)/self.VFR) > 0.0001:
                raise ValueError("The Volumetric Flow Rate provided is different from the value calculated from the Mass Flow Rate and the density")
        elif name == "Ctot" and self.Ctot:
            if abs((value - self.Ctot)/self.Ctot) > 0.0001:
                raise ValueError("The Volumetric Flow Rate provided is different from the value calculated from the Mass Flow Rate and the density")
        elif name == "x" and self.x:
            if sum(value) != 1:
                raise ValueError("The molar fractions don't add up to 1")
            for i in range(len(value)):
                if abs((value[i] - self.x[i])/self.x[i]) > 0.0001:
                    raise ValueError("The molar fractions provided are different from those calculated from the concentrations")
        elif name == "w" and self.w:
            if sum(value) != 1:
                raise ValueError("The mass fractions don't add up to 1")
            for i in range(len(value)):
                if abs((value[i] - self.w[i])/self.w[i]) > 0.0001:
                    raise ValueError("The mass fractions provided are different from those calculated from the molar fractions and the Molecular Weights")
                
    def _update(self, name, value):
        
        if name == "VFR" and self.rho and not self.MFR:
            self.MFR = self.VFR*self.rho
        elif name == "VFR" and self.MFR and not self.rho:
            self.rho = self.MFR/self.VFR
        elif name == "rho" and self.VFR and not self.MFR:
            self.MFR = self.VFR*self.rho
        elif name == "MFR" and self.VFR and not self.rho:
            self.rho = self.MFR/self.VFR
        elif name == "MFR" and self.rho and not self.VFR:
            self.VFR = self.MFR/self.rho
        elif name == "rho" and self.MFR and not self.VFR:
            self.VFR = self.MFR/self.rho
            
        elif name == "C" and not self.Ctot:
            self.Ctot = sum(self.C)
            self.x = [Ci/self.Ctot for Ci in self.C]
            
        if name == "comp" and len(self.comp) == 1:
            self.comp == [self.comp]
        elif name == "comp" and not self.MW:
            MW = [s.MW for s in self.comp]
            if None not in MW:
                self.MW = MW 
                
        elif name == "MW" and type(self.MW) != list:
            self.MW = [self.MW]
            
        elif name == "MW" and self.x and not self.MWmix:
            self.MWmix = sum([self.x[i]*self.MW[i] for i in range(self.Ns)])
        elif name == "MW" and self.w and not self.MWmix:
            self.MWmix = 1/sum([self.w[i]/self.MW[i] for i in range(self.Ns)])
                
        elif name == "Ctot" and self.x and not self.C:
            self.C = [xi*self.Ctot for xi in self.x]
            
        elif name == "x" and self.MW and self.MWmix and not self.w:
            self.w = [self.x[i]*self.MW[i]/self.MWmix for i in range(self.Ns)]
        elif name == "w" and self.MW and self.MWmix and not self.x:
            self.x = [self.w[i]/self.MW[i]*self.MWmix for i in range(self.Ns)]
                                       
        elif name == "x" and self.MW and not self.MWmix:
            self.MWmix = sum([self.x[i]*self.MW[i] for i in range(self.Ns)])
        elif name == "w" and self.MW and not self.MWmix:
            self.MWmix = 1/sum([self.w[i]/self.MW[i] for i in range(self.Ns)])
            
        elif name == "MWmix" and self.x and not self.w:
            self.w = [self.x[i]*self.MW[i]/self.MWmix for i in range(self.Ns)]
        elif name == "MWmix" and self.w and not self.x:
            self.x = [self.w[i]/self.MW[i]*self.MWmix for i in range(self.Ns)]
                               
        if not self.Ns:
            if name == "x" or name == "w" or name == "comp":
                self.Ns = len(value)
            if name == "RM":
                self.Ns = self.RM.Ns
                               
        if not self.Nr:
            if name == "RM":
                self.Nr = self.RM.Nr