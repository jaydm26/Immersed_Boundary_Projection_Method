function T = apply_bc_temp(T,T0,velocity)
    
    global Nx Ny
    
    T.x(:,1) = 0;
    T.x(1,2:Ny) = T.x(2,2:Ny);
    T.x(Nx+1,2:Ny) = T.x(Nx,2:Ny);
    U_inf = interpol(CellData(Nx,Ny),velocity,1);
    U_inf = U_inf.x';
    for i = 1:Nx+1
        for j = Ny+1
            T.x(i,j) = T0.x(i,j) - U_inf(i,j-1) * (T0.x(i,j)-T0.x(i,j-1));
        end
    end
end
    