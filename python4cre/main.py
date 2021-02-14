from models.cstr import CSTR
import numpy as np


names = ["Methane", "Ethane", "Propane"]

C0   = [1.345, 14.24, 243.54]
MW   = [16, 30, 44]

Cp   = [SpecificHeat([1.2, 3, 0.02], [0, -1, 2]),
        SpecificHeat([0.7, 5, 0.02], [0, -1, 2]),
        SpecificHeat([1, 3.5, 0.02], [0, -1, 2])]
Dh0f = [FormationEnthalpy(-1324),
        FormationEnthalpy(-1004, 300),
        FormationEnthalpy(-3456, 257)]
Ds0f = [FormationEntropy(-14),
        FormationEntropy(-4, 320),
        FormationEntropy(-3.6, 267)]

species = [Species(name) for name in names]
for s, i in zip(species, range(len(species))):
    s.Cp = Cp[i]
    s.Dh0f = Dh0f[i]
    s.Ds0f = Ds0f[i]
    s.MW = MW[i]

A    = [1.2e13, 1.1e14, 2e11]
n    = [0.1, 0.3, -0.8]
Ea   = [10000, 20000, 40000]

SC   = [[-1, -1,  2],
        [ 1, -1,  1],
        [ 0,  1, -2]]
RO   = [[ 1,  0,  0],
        [ 1,  0,  1],
        [ 0,  1,  0]]

RM = ReactionMechanism()
RM.A = A
RM.Ea = Ea
RM.n = n
RM.SC = SC
RM.RO = RO

T = 456
P = 2e5
VFR = 12
MFR = 14

mix = Mixture("Inlet Mixture")
mix.RM = RM
mix.comp = species
mix.C = C0
mix.T = T
mix.P = P

ext = ExternalFluid()

cstr1 = CSTR(1)

SC = [[-1], [-1], [1], [0]]
RO = [[1], [0], [0], [0]]

Win = (2.640+6.6)/3600
V = 1.135
T0 = 24+273.15

N0 = [19.5/3.6, 364/3.6, 0, 32.6/3.6]
C0 = [i/Win for i in N0]

A = [16.96e12/3600]
Ea = [72000]
Cp = [146, 75, 192, 82]
Dh0f = [-148918, -275000, -505360, 0]
Tr = [293.15, 293.15, 293.15, 0]

cstr1.inlet(C0, T0, W0 = Win, phase = "L")
cstr1.parameters(V = V, isoT = True)
cstr1.reactions(SC, A, RO = RO, Ea = Ea, Cp = Cp, Dh0f = Dh0f, Tr = Tr)
cstr1.solve()

