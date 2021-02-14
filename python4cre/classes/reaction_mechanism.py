class ReactionMechanism(object):
    
    mechanisms = 0
    name = None
    A    = None
    n    = None
    Ea   = None
    SC   = None
    RO   = None
    Ns   = None
    Nr   = None
    
    variables = {"name": (str, int),
                 "A" : (list, float, int),
                 "n" : (list, float, int),
                 "Ea" : (list, float, int),
                 "SC" : (list, float, int),
                 "RO" : (list, float, int),
                 "Ns" : (int,),
                 "Nr" : (int,)}
                 

    def __init__(self, name = None):

        if name:
            self.name = name
        else:
            self.name = ReactionMechanism.mechanisms
            ReactionMechanism.mechanisms += 1
        
    def __str__(self):
        return f"{type(self)}, name = {self.name}"
    
    def __repr__(self):
        return "<class ReactionMechanism>"
        
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
        
        def check_length(N, v):
            if (N != len(v)):
                raise ValueError(f"The length of the list '{name}' doesn't match the length of the lists previously defined")
        
        if name == "A":
            if type(value) != list:
                value = [value]
            if self.Nr:
                check_length(self.Nr, value) 
            else:
                self.Nr = len(value)
        
        elif name == "n":
            if type(value) != list:
                value = [value]
            if self.Nr:
                check_length(self.Nr, value) 
            else:
                self.Nr = len(value)
        
        elif name == "Ea":
            if type(value) != list:
                value = [value]
            if self.Nr:
                check_length(self.Nr, value) 
            else:
                self.Nr = len(value)
                
        elif name == "SC":
            if type(value[0]) != list:
                value = [[c] for c in value]
            if self.Ns:
                check_length(self.Ns, value) 
            else:
                self.Ns = len(value)
                
            if self.Nr:
                for vv in value:
                    check_length(self.Nr, vv) 
            else:
                self.Nr = len(value[0])
                self.SC = value
        
        elif name == "RO":
            if type(value[0]) != list:
                value = [[c] for c in value]
            if self.Ns:
                check_length(self.Ns, value) 
            else:
                self.Ns = len(value)
                
            if self.Nr:
                for vv in value:
                    check_length(self.Nr, vv) 
            else:
                self.Nr = len(value[0])
                self.RO = value
    
    def _update(self, name):
        if name == "SC" and not self.RO:
            self.RO = self.SC