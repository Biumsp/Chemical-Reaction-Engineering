

[x, X] = ode23s(@ex, [0 0.1], 0, [], data);

plot(x, X)







function yy = ex(x, X, data)

    rho_cat = data.rho_cat;
    void    = data.void;
    k       = data.A*exp(-data.Ea_R/(250+273.15));
    CA0     = data.Pin/8.314/data.Tin;
    u0      = 3;

    yy = rho_cat*(1-void)*k*(1-X)^2*CA0^2/(1+2*X)^2/u0/CA0;
end