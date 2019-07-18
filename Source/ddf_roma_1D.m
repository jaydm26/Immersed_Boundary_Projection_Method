function ddf = ddf_roma_1D(r,params)
    %
    % Calculates the value of the discrete delta function for the input r
    % using the Roma et.al. (1999) version for the discrete delta 
    % function.

    r = abs(r/params.dx);
    if r < 0.5
        ddf = 1/3 * (1 + sqrt(1-3*r^2));
    elseif r >= 0.5 && r < 1.5
        ddf = 1/6 * (5 - 3*r - sqrt(1-3*(1-r)^2));
    else
        ddf = 0;
    end
end