X = 0:0.01:1;
ThetaB = 10/90;
Ca0 = 90;
Ca = Ca0*(1-X);
Cb = -Ca0*(X) + Ca0*(ThetaB + 2*X);
R = Ca.*Cb;
plot(X, 1./R)
X(R == max(R))