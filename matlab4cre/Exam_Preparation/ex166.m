T = [40 40 40 60]'+273.15;
CAin = [2.5 2 1.5 2.5]';
CA = [1 0.873 0.727 0.750]';

Qin = 1/60;
V = 0.5e-3;

LHS = Qin./CA.^2.*(CAin - CA.*(0.9 + 0.1*CA)./(0.9 + 0.1*CAin));

y = log(LHS)
X = [ones(4,1) 1./T]

A = X'*X;
b = X'*y;

a = A\b;

A = exp(a(1))
EA = -a(2)*8.314

% Statistics
Ymod = X*a;
SSres = (y - Ymod)'*(y - Ymod);
YexpAvg = mean(y);
SStot = (y - YexpAvg)'*(y - YexpAvg);
R2 = 1- SSres/SStot;
fprintf('R2 = %f\n', R2);
