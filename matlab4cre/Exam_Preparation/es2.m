X = 0.9999;
Af = 1;
Ab = 9.89;
Ef = 5555;
Eb = 19444;

r = @(T) Af*exp(-Ef./T)*(1 - X) - Ab*exp(-Eb./T)*X;

T = 873.15:1173.15;

plot(T, r(T))