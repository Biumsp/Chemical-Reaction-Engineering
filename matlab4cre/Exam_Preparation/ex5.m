clear
close


A = 8e8;
E_R = 2e4;
Dh0r = -80e3;
P = 3*101325;
Qin = 0.25;

L = 12;
rho_cp = 20000;
Tin = 720;
U = 80;
a = 4/sqrt(4/pi);


FNin = P/8.314/Tin*Qin/4;
FHin = 3*FNin;
SC = [-1; -3; 2];


inlet = [FNin FHin 0 Tin Tin];

[x, var] = ode23s(@funk,[0 L], inlet,[], Qin, A, E_R, P, Dh0r, rho_cp, SC,  U, a, Tin);

F = var(:, 1:3);
T = var(:, 4);
Tex = var(:, 5);

figure
plot(x, F)
figure
plot(x, T)
figure
plot(x, Tex)


function yy = funk(x, var, Qin, A, E_R, P, Dh0r, rho_cp, SC, U, a, Tin)

    F = var(1:3);
    T = var(4);
    Tex = var(5);
    
    C = F/sum(F)*P/8.314/T;
    
    r = A*exp(-E_R/T)*C(1)*C(2);
    dF = r*SC;
    
    
    Qex = (T - Tex)*U*a;
    Qr = -Dh0r*r;
    
    dT = (Qr - Qex)/rho_cp/Qin;
    
    Q = Qin*Tex/Tin;
    dTex = -Qex/rho_cp/Q;

    yy = [dF' dT dTex]';

end