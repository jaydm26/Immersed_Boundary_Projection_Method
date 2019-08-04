function u = L_inv_operation(f,g_hat)
    %L_INV_OPERATION Solves the Poisson equation: 
    %      d^2 u / dx^2 + d^2 u / dy^2 = f 
    % using the Lattice Green's Function operator.
    %
    % u = L_inv_operation(f,g_hat)
    %
    % Variable lookup:
    % 
    % u: solution for the Poisson's Equation.
    %
    % f: input for the Poisson's Equation.
    %
    % g_hat: FFT2 of the L_inv operator. Obtained from the function L_inv.
    %
    % Created by Jay Mehta (18 July 2019)
    
    Nx = f.size(1);
    Ny = f.size(2);
    
    switch f.data
        case "node"
            u = NodeData(Nx,Ny);
            rhs_hat = fft2(f.x,2*Nx+1,2*Ny+1);
            out_hat = rhs_hat .* g_hat;
            out_hat = ifft2(out_hat);
            u.x = out_hat(Nx+1:end,Ny+1:end);
        case "cell"
            u = CellData(Nx,Ny);
            rhs_hat =  fft2(f.x,2*Nx+2,2*Ny+2);
            out_hat = rhs_hat .* g_hat;
            out_hat = ifft2(out_hat);
            u.x = out_hat(Nx+1:end,Ny+1:end);
    end
end