function ddf = ddf_roma_1D(r)
    global dx
    r = abs(r/dx);
    if r < 0.5
        ddf = 1/3 * (1 + sqrt(1-3*r^2));
    elseif r >= 0.5 && r < 1.5
        ddf = 1/6 * (5 - 3*r - sqrt(1-3*(1-r)^2));
    else
        ddf = 0;
    end
end