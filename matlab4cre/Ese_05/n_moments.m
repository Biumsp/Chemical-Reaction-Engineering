function yy = n_moments(F, x, n)
% Computes the moment of order n of the function F(x)

    yy = 0;
    for i = 1:length(F)-1
        yy = yy + ((x(i+1)^n)*F(i + 1) + (x(i)^n)*F(i))/2 * (x(i+1) - x(i));
    end

end