import numpy 
from matplotlib import pyplot as plt
import copy

class CSTR:
    def __init__(self, tau, dt = 1/50):
        self.tau = tau
        self.Cout = {}
        self.dt = dt

    def _convolute(self, Cin):
        if len(Cin[0]) == 1:
            t_s = [i*self.dt for i in range(200)]
            Cout_s = [self._E(t)*Cin[1][0] for t in t_s]
            return Flow([t_s, Cout_s])

        t_s = Cin[0] + [i*self.dt + Cin[0][-1] for i in range(300)]
        Cin[1] = Cin[1] + [1e-7 for _ in range(300)]

        Cout_s = []
        for t in t_s:
            Cout = 0
            
            for ii in range(int(t//self.dt)):
                if ii >= len(Cin[1]):
                    break
                Cout += dt*Cin[1][ii]*self._E(t - ii*self.dt)

            Cout_s.append(Cout)

        self.Cout = [t_s, Cout_s]
        self.Cout = self._clean(self.Cout)
        return Flow(copy.deepcopy(self.Cout))

    def __rmul__(self, other):
        if type(other) != list:
            Cin = other.Cout
        else:
            Cin = other
        return self._convolute(Cin)

    def _E(self, t):
        return numpy.exp(-t/self.tau)/self.tau

    def _clean(self, Cout):
        remove = []
        for i in range(len(Cout[0])):
            if Cout[1][i] <= 1e-7:
                remove.append(i)
        
        for index in sorted(remove, reverse=True):
            del Cout[0][index]
            del Cout[1][index]

        return Cout

class PFR:
    def __init__(self, tau):
        self.tau = tau
        self.Cout = {}

    def _convolute(self, Cin):
        t = []
        for i in Cin[0]:
            t.append(i + self.tau)

        self.Cout = [t, Cin[1]]
        self.Cout = self._clean(self.Cout)

        return Flow(copy.deepcopy(self.Cout))

    def __rmul__(self, other):
        if type(other) != list:
            Cin = other.Cout
        else:
            Cin = other
        return self._convolute(Cin)

    def _clean(self, Cout):
        remove = []
        for i in range(len(Cout[0])):
            if Cout[1][i] <= 1e-7:
                remove.append(i)
        
        for index in sorted(remove, reverse=True):
            del Cout[0][index]
            del Cout[1][index]

        return Cout


class LFR:
    def __init__(self, tau, dt = 1/50):
        self.tau = tau
        self.Cout = {}
        self.dt = dt

    def _convolute(self, Cin):
        if len(Cin[0]) == 1:
            t_s = [i*self.dt for i in range(200)]
            Cout_s = [self._E(t)*Cin[1][0] for t in t_s]
            return Flow([t_s, Cout_s])

        t_s = Cin[0] + [i*self.dt + Cin[0][-1] for i in range(300)]
        Cin[1] = Cin[1] + [1e-7 for _ in range(300)]

        Cout_s = []
        for t in t_s:
            Cout = 0
            
            for ii in range(int(t//self.dt)):
                if ii >= len(Cin[1]):
                    break
                Cout += self.dt*Cin[1][ii]*self._E(t - ii*self.dt)

            Cout_s.append(Cout)

        self.Cout = [t_s, Cout_s]
        self.Cout = self._clean(self.Cout)
        return Flow(copy.deepcopy(self.Cout))

    def __rmul__(self, other):
        if type(other) != list:
            Cin = other.Cout
        else:
            Cin = other
        return self._convolute(Cin)

    def _E(self, t):
        if t < self.tau/2:
            return 0
            
        E = self.tau**2/2/t**3

        return max(0, E)

    def _clean(self, Cout):
        remove = []
        for i in range(len(Cout[0])):
            if Cout[1][i] <= 1e-7:
                remove.append(i)
        
        for index in sorted(remove, reverse=True):
            del Cout[0][index]
            del Cout[1][index]

        return Cout


def pulse():
    t = [0]
    C = [1]
    
    Cout = [t, C]
    return Flow(Cout)


class Flow:
    def __init__(self, C):
        self.Cout = C
    
    def __rmul__(self, other):
        if type(other) != float:
            raise TypeError(f"Cannot multiply Flow and {type(other)}")
        else:
            C = [other*i for i in self.Cout[1]]
            return Flow([copy.deepcopy(self.Cout[0]), C])

    def __add__(self, other):
        t_s = other.Cout[0]
        new_t = []
        C = []
        for t, i in zip(t_s, list(range(len(t_s)))):
            new_t.append(t)
            try:
                index = self.Cout[0].index(t)
                C.append(self.Cout[1][index] + other.Cout[1][i])
            except:
                C.append(other.Cout[1][i])

        for t, i in zip(self.Cout[0], list(range(len(self.Cout[0])))):
            if t in new_t:
                continue

            new_t.append(t)
            try:
                index = other.Cout[0].index(t)
                C.append(other.Cout[1][index] + self.Cout[1][i])
            except:
                C.append(self.Cout[1][i])

        dd = dict(zip(C, new_t))
        C.sort(key = lambda x: dd[x])
        new_t = sorted(new_t)

        return Flow([new_t, C])


dt = 1/50 

p = pulse()
print(type(p))
cstr = CSTR(1, dt)
pfr1 = PFR(1)
pfr2 = PFR(5)
lfr = LFR(1, dt)

RTD = p*lfr*cstr
#print("Accuracy:", sum(RTD.Cout[1])*dt)

plt.plot(RTD.Cout[0], RTD.Cout[1])
plt.show()