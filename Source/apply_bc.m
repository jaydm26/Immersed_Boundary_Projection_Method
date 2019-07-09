function U = apply_bc(U,t)
    global Nx Ny uL uR uB uT vL vR vB vT
    
    switch U.data
        case 'edge'
            for i = 2:Ny+1
                U.x(1,i)    = uL(i,t);
                U.x(Nx+1,i) = uR(i,t);
            end
            for i = 2:Nx
                U.x(i,1)    = -U.x(i,2)     + 2*uB(i,t);
                U.x(i,Ny+2) = -U.x(i,Ny+1)  + 2*uT(i,t);
            end
            for i = 2:Ny
                U.y(1,i)    = -U.y(2,i)     + 2*vL(i,t);
                U.y(Nx+2,i) = -U.y(Nx+1,i)  + 2*vR(i,t);
            end
            for i = 2:Nx+1
                U.y(i,1)    = vB(i,t);
                U.y(i,Ny+1) = vT(i,t);
            end
    end
end