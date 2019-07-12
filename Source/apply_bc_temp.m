function T = apply_bc_temp(T,T0)
    
    global Nx Ny Co
    
    T.x(:,1) = 0;
    T.x(1,2:Ny) = T.x(2,2:Ny);
    T.x(Nx+1,2:Ny) = T.x(Nx,2:Ny);
    for i = 1:Nx+1
        for j = Ny+1
            T.x(i,j) = T0.x(i,j) - Co * (T0.x(i,j)-T0.x(i,j-1));
        end
    end
end
    