function out = L_inv_operation(rhs,DataType)
    
    global Nx Ny g_hat_n g_hat_c
    switch DataType
        case "node"
            out = NodeData(Nx,Ny);
            rhs_hat = fft2(rhs,2*Nx+1,2*Ny+1);
            out_hat = rhs_hat .* g_hat_n;
            out_hat = ifft2(out_hat);
            out.x = out_hat(Nx+1:end,Ny+1:end);
        case "cell"
            out = CellData(Nx,Ny);
            rhs_hat =  fft2(rhs,2*Nx+2,2*Ny+2);
            out_hat = rhs_hat .* g_hat_c;
            out_hat = ifft2(out_hat);
            out.x = out_hat(Nx+1:end,Ny+1:end);
    end
end