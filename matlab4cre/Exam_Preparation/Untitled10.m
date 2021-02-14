close
C0 = [1000 0];

[t, C] = ode23s(@funk, [0, 70000], C0);
Ca = C(:, 1);
Cb = C(:,2);

X = (C0(1) - Ca)/C0(1);
figure 
plot(t, X)

A = [0.5 25000];
Ea_R = [20e3 60e3]/8.314;
T = 300;

k   =  A.*exp(-Ea_R/T);
U = 500;
S_V = 4/0.42;
r = [k(1)*Ca k(2)*Cb];
Dh0r = [-150e3 -40e3];

Tex = (r*Dh0r' + U*S_V*T)/U/S_V;
figure
plot(X, Tex)


function yy = funk(t, C)

    Ca = C(1);
    Cb = C(2);
    
    A = [0.5 25000];
    Ea_R = [20e3 60e3]/8.314;
    T = 300;
    
    k   =  A.*exp(-Ea_R/T);
    dCa = -k(1)*Ca;
    dCb =  k(1)*Ca - k(2)*Cb;
    yy = [dCa dCb]';

end