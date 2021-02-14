Ctot = 101325/8.314/298;
CH0 = Ctot/2;

FHin = 10;
Dh0r = -22.06*2*1000;
Cpin = 8.8 + 11.2;
DCp = 7.6*2 - 8.8 - 11.2;
X = 0.9;

T = 298 + (-Dh0r*X)/(DCp*X + Cpin)
Rh = -1e12*exp(-20000/T)*1e-6*CH0^2*(1-X)^2*298^2/T^2;
V = 10*X/-Rh