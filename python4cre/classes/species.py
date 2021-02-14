from .formation_enthalpy import FormationEnthalpy
from .formation_entropy import FormationEntropy
from .formation_gibbs import FormationGibbs
from .specific_heat import SpecificHeat

class Species:

    name   = None
    Dh0f   = None     # Formation enthalpy 
    Ds0f   = None     # Formation entropy
    Dg0f   = None     # Gibbs' energy of formation
    Cp     = None     # Specific heat 
    Pitzer = None     # Pitzer's acentric factor 
    Diff   = None     # Diffusivity in the mixture
    
    variables = {"name" : (str, int),
                 "Dh0f" : (FormationEnthalpy,),
                 "Ds0f" : (FormationEntropy,),
                 "Dg0f" : (FormationGibbs,),
                 "Cp" : (SpecificHeat, int, float),
                 "Pitzer" : (float, int),
                 "Diff" : (float, int)}

    def __init__(self, name):
        self.name = name
        
    def __str__(self):
        return f"{type(self)}, name = {self.name}"
    
    def __repr__(self):
        return "<class Species>"
        
    def __setattr__(self, name, value):
        self._check_validity(name, value)
        self._check_coherence(name, value)
        
        super.__setattr__(self, name, value)
        
        self._update(name)
        
    def _check_validity(self, name, value):
        if name not in self.variables:
            raise ValueError(f"{type(self).__name__} object has no attribute '{name}'")  
        elif type(value) not in self.variables[name]:
            raise TypeError(f"Wrong type {type(value)} for attribute '{name}'")
            
        return True
            
    def _check_coherence(self, name, value):
        pass
    
    def _update(self, name):
        if name == "Cp":
            if self.Dh0f and self.Ds0f and (self.Dh0f.Tref != self.Ds0f.Tref):
                self._update_Dh0f()                
        elif name == "Dh0f" and self.Cp:
            if self.Ds0f and (self.Dh0f.Tref != self.Ds0f.Tref):
                self._update_Dh0f() 
        elif name == "Ds0f" and self.Cp:
            if self.Dh0f and (self.Dh0f.Tref != self.Ds0f.Tref):
                self._update_Dh0f() 
                
        if self.Dh0f and self.Ds0f and not self.Dg0f:
            if self.Dh0f.Tref == self.Ds0f.Tref:
                Dg0f = self.Dh0f - self.Ds0f*self.Ds0f.Tref
                self.Dg0f = FormationGibbs(Dg0f, self.Ds0f.Tref)
            
    def _update_Dh0f(self):
        print(f"""formation enthalpy and formation entropy of species {self.name} have different reference temperatures""")  
        print("Updating Dh0f...")
        
        value = self.Cp*self.Ds0f.Tref - self.Cp*self.Dh0f.Tref
        new_Dh0f = FormationEnthalpy(value, self.Ds0f.Tref)
        
        print(f"Dh0f updated from '{self.Dh0f}' \nto '{new_Dh0f}'")
        
        self.Dh0f = new_Dh0f