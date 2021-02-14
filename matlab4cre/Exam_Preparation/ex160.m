U = 1360;
S = 3.3;
A = 4e6;
Ea_R = 7900;
V = 5;
rhoCp = 4.2e6;
Dh0r = -1.67e6*100e-3;
CA0 = 1e3;
T0 = 293.15;

inlet = [CA0, T0]';
[t1, var1] = ode23(@funk_1, [0 100000], inlet, [], U, S, A, Ea_R, V, rhoCp, Dh0r);
CA = var1(:,1);
T = var1(:,2);

index = find(T >= 90+273.15, 1);

T1 = T(1:index);
CA1 = CA(1:index);
t1 = t1(1:index);
CAsteam = CA1(end);
Msteam = rhoCp*V*(35)/2257e3;

fprintf('M steam = %f [kg]\n', Msteam');

inlet = [CA1(end), T1(end)]';
[t2, var2] = ode23(@funk_2, [t1(index) t1(index)+10000], inlet, [], U, S, A, Ea_R, V, rhoCp, Dh0r);
CA = var2(:,1);
T = var2(:,2);

index = find(T <= 45+273.15, 1);

T2 = T(1:index);
CA2 = CA(1:index);
t2 = t2(1:index);

CA = [CA1' CA2'];
T  = [T1' T2'];
t = [t1' t2'];


figure;
yyaxis left
plot(t, CA)
yyaxis right
plot(t, T-273)

fprintf('total t = %f [s]\n', t(end));


function yy = funk_1(t, var, U, S, A, Ea_R, V, rhoCp, Dh0r)

    CA = var(1);
    T = var(2);
    
    if T < 55+273.15
        Tex = 120+273.15;
    elseif T > 55+273.15
        Tex = T;
    end
    
    k   = A*exp(-Ea_R/T);
    Qr  = -Dh0r*k*CA;
    
    
    dCA = -k*CA;
    dT  = (Qr - U*S*(T - Tex))/rhoCp;
    
    yy = [dCA dT]';

end

function yy = funk_2(t, var, U, S, A, Ea_R, V, rhoCp, Dh0r)

    CA = var(1);
    T = var(2);
    
    Tex = 15 + 273.15;
    
    k   = A*exp(-Ea_R/T);
    Qr  = -Dh0r*k*CA;
    
    
    dCA = -k*CA;
    dT  = (Qr - U*S*(T - Tex))/rhoCp;
    
    yy = [dCA dT]';

end