function dhf = heaviside_2(params,r)
    
    r = r/params.dx;
    if r < 0.5 && r > 0.5
        dhf = 1/18 * (3*r*(sqrt(1-3*x^2) + 2) + sqrt(3) * asin(sqrt(3)*x));
    elseif (r >= 0.5 && r < 1.5) || (r <= -0.5 && r > -1.5)
        dhf = 1/6 * (5 - 3*r - sqrt(1-3*(1-r)^2));
    else
        dhf = 0;
    end