clear variables

Da = 1;

for ii = 1:length(Da)
    
    tau = Da(ii);
    
    y(ii) = integral(@(t)funz(tau, t), 0, 10*tau);
    y(ii) = integral(@(t)funz_2(y(ii), tau, t), 0, 10*Da(ii));
    
end



function yy = funz(tau, t)

    E = 1/tau.*exp(-t./tau);
    C = 1./(1 + t);

    yy = C.*E;
end

function yy = funz_2(Cin, tau, t)

    E = 1/tau.*exp(-t./tau);
    C = Cin./(1 + t*Cin);

    yy = C.*E;

end
    
    
    