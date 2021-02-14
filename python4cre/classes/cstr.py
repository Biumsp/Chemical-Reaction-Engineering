import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt 
from .constants import BIGFILLER, FILLER
import time

class CSTR:

    CSTRs = 0                # Number of active CSTRs
    inlet = None             # Inlet mixture
    outlet = None            # Inlet mixture
    ext = None               # External (cooling7heating) mixture

    EoS = "IG"               # Equation Of State (Default: Ideal Gas)
    method = "BDF"           # Numerical method for the integration (Default: BDF)
    display_opt = None       # Display Options

    US = None
    Lex = 0
    Qadd = 0

    adiabatic = False       
    isothermal = False      
    isoperibolic = False
    constant_heat_exchange = False

    def __init__(self, name = None):
        if name:
            self.name = name
        else:
            name = CSTR.CSTRs
            CSTR.CSTRs += 1

    def __repr__(self):
        return f"<CSTR {self.name}>"


        """Set display options: opt = Display option (default: None)\n
                                opt Values: - 'off' (turns off the display of stuff),\n
                                            - 'detailed' (shows detail about the numerical integration and the internal assumptions)
                                            -  None (shows the internal assumptions but doesn't show details about the numerical integration)  """
        """Set general parameters: tau  = Residence time (default: V/W0)\n
                                   V    = Recator Volume (default: None)\n"""
        """Set the equation of state and the numerical integration method:\n
           EoS = Equation of State;\n
           EoS Options: "IG", "RK", "RKS", "PR" (default: Ideal Gas ("IG"))\n
           Tc = Critical temperatures of all species (default: None)\n
           Pc = Critical pressures of all species (default: None)\n
           omega = Pitzer acentric factors of all species (default: None)\n
           method = Numerical integration method;\n
           method Options: "RK45", "RK23", "LSODA", "BDF", "Radau", "DOP853" (default: "RK45")\n\n
           For info on the numerical integration methods see:\n
           'https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html#scipy.integrate.solve_ivp'"""
        

    def geometry(self, tau = None, V = None):
        
        self.tau = tau
        self.V = V
        if not self.tau:
            if (not self.W0) | (not self.V):
                self._error("The user must provide the residence time or the inlet flow rate and the volume")
            else:
                self.tau = self.V/self.W0

        self.US = US
        self.Tex0 = Tex0
        self.Qadd = Qadd
        self.Lex = Lex

        if not self.US:
            self.adiabatic = True
        if (self.adiabatic and Qadd):
            self.constant_heat_exchange = True
        

    def _display_assumptions(self):
        if self.display_opt != "off":
            print(FILLER)
            print("CSTR", self.name)
            print(FILLER)
            print()
            if self.isothermal:
                print("Isothermal CSTR assumed")
            else:
                print("Non-Isothermal CSTR assumed")
            if self.adiabatic:
                print("Adiabatic CSTR assumed")
            else:
                print("Non-Adiabatic CSTR assumed")
            if self.constant_heat_exchange:
                print("Only constant heat exchange assumed")
            if self.phase == "L":
                print("Liquid phase assumed")
            if self.phase == "G":
                print("Gas phase assumed")
                print("Using Ideal Gas EoS")
            print()
            print(BIGFILLER, "\n")


    def solve(self, alpha = 50):
        """Solves the CSTR using a false-transient method""" 

        if alpha == 50:
            self._display_assumptions()

        t = (0, alpha*self.tau)
        ICs = [self.T0] + self.C0
        sol = solve_ivp(ODE_CSTR, t, ICs, args = (self,), method = self.method)
        y = sol.y

        if sum([abs(y[i, -1] - y[i, -2]) for i in range(self.Ns)]) > 1e-6:
            if alpha > 1000:
                self._error("False-transient method couldn't reach steady-state")
            else:
                alpha *= 2
                self.solve(alpha = alpha) 
        else:
            self.C = y[1:]
            if len(self.C) == 1:
                self.C = [self.C]
            self.T = y[0] 

        self.X = [(self.C0[i] - self.C[i])/self.C0[i] for i in range(self.Ns)]
        self.X = [i for i in self.X if i > 0]

        print("CSTR {}, Solution Found!".format(self.name))


def ODE_CSTR(t, y, r):

    if r.EoS != "IG":
        print("Non ideal EoS still unavailable: assuming Ideal Gas")

    T = y[0]
    C = y[1:]

    # Molar fractions
    xin = [r.C0[i]/sum(r.C0) for i in range(r.Ns)]
    x   = [C[i] /sum(C)  for i in range(r.Ns)] 

    # Kinetic constants and reaction rates
    k = [r.A[i]*(T**(r.n[i]))*np.exp(-r.Ea[i]/r.R/T) for i in range(r.Nr)]
    if not r.RO:
        Ri = [sum([r.SC[i][j] * k[j] * np.prod([C[m]**abs(min(r.SC[m][j], 0)) for m in range(r.Ns)]) for j in range(r.Nr)]) 
                  for i in range(r.Ns)]
    else:
        Ri = [sum([r.SC[i][j] * k[j] * np.prod([C[m]**r.RO[m][j] for m in range(r.Ns)]) for j in range(r.Nr)]) 
                  for i in range(r.Ns)]

    if r.phase == "L":
        Ai = [(r.C0[i] - C[i])/r.tau + Ri[i] for i in range(r.Ns)]
    elif r.phase == "G":
        gammaMW = r.P0*T/r.P/r.T0 * sum([xin[i]*r.MW[i] for i in range(r.Ns)])/sum([x[i]*r.MW[i] for i in range(r.Ns)])
        Ai = [(r.C0[i] - gammaMW*C[i])/r.tau + Ri[i] for i in range(r.Ns)]
    
    if r.isothermal:
        dTdt = 0
        return [dTdt] + Ai
    elif r.adiabatic | r.constant_heat_exchange:
        Qex = r.Qadd
    else:
        beta = 1 - np.exp(- r.US/r.Wex/r.Cpex)
        Qex = -r.Wex*r.Cpex*beta*(T - r.Tex0)

    Cpmix = sum([x[i]*r.Cp[i] for i in range(r.Ns)])

    L1 = [r.SC[i][0]*x[i]*(r.Dh0f[i] - 2*int(r.phase == "G")*r.R*T + integral_polinomial(r.Cp[i], r.Cp_fun[i], r.Tr[i], T)) for i in range(r.Ns)]
    L2 = [xin[i]*(r.Dh0f[i] + integral_polinomial(r.Cp[i], r.Cp_fun[i], r.Tr[i], r.T0)) for i in range(r.Ns)]
    L3 = [x[i]*(r.Dh0f[i] + integral_polinomial(r.Cp[i], r.Cp_fun[i], r.Tr[i], T)) for i in range(r.Ns)]
    if r.phase == "L":
        if Qex + r.Lex == 0:
            dTdt = - sum(L1)*sum(Ai)/(Cpmix - int(r.phase == "G")*r.R)/sum(C) + (sum(L2)*sum(r.C0) - sum(L3)*sum(C))/r.tau/(Cpmix - int(r.phase == "G")*r.R)/sum(C)
        else:
            dTdt = (Qex + r.Lex - sum(L1)*sum(Ai)*r.V)/r.V/(Cpmix - int(r.phase == "G")*r.R)/sum(C) + (sum(L2)*sum(r.C0) - sum(L3)*sum(C))/r.tau/(Cpmix - int(r.phase == "G")*r.R)/sum(C)

    else:
        if Qex + r.Lex == 0:
            dTdt = - sum(L1)*sum(Ai)/(Cpmix - int(r.phase == "G")*r.R)/sum(C) + (sum(L2)*sum(r.C0) - gammaMW*sum(L3)*sum(C))/r.tau/(Cpmix - int(r.phase == "G")*r.R)/sum(C)
        else:
            dTdt = (Qex + r.Lex - sum(L1)*sum(Ai)*r.V)/r.V/(Cpmix - int(r.phase == "G")*r.R)/sum(C) + (sum(L2)*sum(r.C0) - gammaMW*sum(L3)*sum(C))/r.tau/(Cpmix - int(r.phase == "G")*r.R)/sum(C)

    return [dTdt] + Ai
