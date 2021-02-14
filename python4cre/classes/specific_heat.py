import numpy as np

class SpecificHeat(object):

    def __init__(self, values, functions = 0):

        if type(values) != list:
            values = [values]
        self.values = values

        if type(functions) != list:
            functions = [functions]
        self.functions = functions
        
        self.int_functions = []
        self._get_integral()

    def __str__(self):
        function = [str(round(c, 2)) + "T^" + str(e) for c, e in zip(self.values, self.functions)]
        return "Cp = " +  " + ".join(function) + " [J/mol/K]"
    
    def __repr__(self):
        return "<class 'SpecificHeat'>"

    def __mul__(self, other):
        """Integrates the specific heat:\n
           + If Cp*value returns the integral from 0 to value of Cp\n
           + If Cp*(upper, lower) or Cp*[upper, lower] returns the integral from lower to upper of Cp"""
        
        if type(other) == float or type(other) == int:
            return sum([fun(other) for fun in self.int_functions])
        
        elif type(other) == tuple or type(other) == list:
            if len(other) == 2:
                return self*other[0] - self*other[1]
            else:
                raise ValueError("the extremes of integrations must be 2 numbers")

    def __rmul__(self, other):
        return SpecificHeat([other*v for v in self.values], self.functions)
    
    def __add__(self, other):
        if type(other) == SpecificHeat:
            functions = {}
            for c, f in zip(self.values + other.values, self.functions + other.functions):
                if f in functions:
                    functions[f] += c
                else:
                    functions.update({f : c})
                    
            return SpecificHeat(list(functions.values()), list(functions.keys()))
        
        else:
            raise TypeError(f"unsupported operand type(s) for +: '{type(other).__name__}' and 'SpecificHeat'")

    def _get_integral(self):

        for e, c in zip(self.functions, self.values):
            if e == -1:
                self.int_functions.append(lambda x, e = e, c = c: c*np.log(x))
            else:
                self.int_functions.append(lambda x, e = e, c = c: (c*x**(e + 1))/(e + 1))