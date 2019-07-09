function y = fast_idst(x,dim)
    
    % Perform a IDST for x in the required dimension
    
    [Nx,Ny] = size(x);
    y = zeros(Nx,Ny);
    
    switch dim
        case 1
            y = 2/(Ny+1) .* fast_dst(x,1);
        case 2
            y = 2/(Nx+1) .* fast_dst(x,2);
    end
    