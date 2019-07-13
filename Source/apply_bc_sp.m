function gamma = apply_bc_sp(gamma,gamma0,velocity)
    
    global Nx Ny
    
    gamma.x(:,1) = 0;
    gamma.x(1,2:Ny) = gamma.x(2,2:Ny);
    gamma.x(Nx+1,2:Ny) = gamma.x(Nx,2:Ny);
%     gamma.x(:,Ny+1) = 0;
    U_inf = interpol(NodeData(Nx,Ny),velocity,1);
    U_inf = U_inf.x';
    for i = 1:Nx+1
        for j = Ny+1
            gamma.x(i,j) = gamma0.x(i,j) - U_inf(i,j-1) * (gamma0.x(i,j)-gamma0.x(i,j-1));
        end
    end
end