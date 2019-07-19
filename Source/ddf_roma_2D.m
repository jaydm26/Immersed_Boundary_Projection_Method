function ddf = ddf_roma_2D(params,x,y)
    %DDF_ROMA_2D Execute the discrete delta function in 2 dimensions.
    %
    % ddf = ddf_roma_2D(params,x,y)
    %
    % Created by Jay Mehta (18 July 2019)
    
    ddf = ddf_roma_1D(params,x) * ddf_roma_1D(params,y);
end