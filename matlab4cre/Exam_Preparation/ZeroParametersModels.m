function [sigmaTheta2, Pe] = ZeroParametersModels(t, E)

    function m = moments(x, y, n)
        m = 0;
        for ii = 1:length(x)-1
            m = m + 1/2*(x(ii)^n*y(ii) + x(ii+1)^n*y(ii+1))*(x(ii+1) - x(ii));
        end
    end

    m0 = moments(t, E, 0);
    E  = E/m0;
    
    m1 = moments(t, E, 1);
    m2 = moments(t, E, 2);
    
    sigmaTheta2 = (m2 - m1^2)/m1^2; 
    Pe = (2+sqrt(4 + 32*sigmaTheta2))/2/sigmaTheta2;
    
end