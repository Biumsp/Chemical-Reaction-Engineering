CA0 = 2e3;
ThetaB = 1.5e3/CA0;

tau1 = 0.8*3600;
X1 = 0.5;
tau2 = 5*3600;
X2 = 0.6;

alpha = CA0*X1^2/((1-X1)*(ThetaB-X1));
beta = X1/tau1/((1-X1)*(ThetaB-X1));

gamma = CA0*(1-X2)*(ThetaB-X2);
kb = (X2/tau2 - gamma*beta)/(alpha*gamma - CA0*X2^2)

kf = kb*alpha + beta

Keq = kf/kb