A = 250/60/1e12;
E_R = 6000/1.987;
Dhr = -3000*4.186;
Cain = 0.05*1e6;
FAin = 0.25*25/3600;
ThetaB = 3;
Tin = 25+273.15;
Cp_mix = 10*4.186;
X_Target = 0.6;

T = @(X) Tin + (-Dhr*X/Cp_mix);
r = @(X) A*exp(-E_R./T(X)).*Cain^3.*(1-X).*(ThetaB - 2*X).^2;

funk = @(X) 1./r(X);

V = FAin*integral(funk, 0, X_Target)