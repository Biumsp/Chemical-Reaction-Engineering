class FormationGibbs(float):
    def __new__(cls, value, Tref = 298.15, Pref = 1e5):
        return float.__new__(cls, value)

    def __init__(self, value, Tref = 298.15, Pref = 1e5):
        self.value = value
        self.Tref = Tref
        self.Pref = Pref
        
    def __str__(self):
        return f"Dg0f = {round(self.value, 2)} [J/mol] @(T = {self.Tref} [K], P = {self.Pref*1e-5} [bar])"
    
    def __repr__(self):
        return "<class 'FormationGibbs'>"