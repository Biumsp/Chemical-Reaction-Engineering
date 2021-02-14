function [a, R2] = LinearRegression(y, X)

    A = X'*X;
    b = X'*y;
    
    a = A\b;
    
    % Statistics
    Ymod = X*a;
    SSres = (y - Ymod)'*(y - Ymod);
    YexpAvg = mean(y);
    SStot = (y - YexpAvg)'*(y - YexpAvg);
    R2 = 1- SSres/SStot;
    fprintf('R2 = %f\n', R2);

end