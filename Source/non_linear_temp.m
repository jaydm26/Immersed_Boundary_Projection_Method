function nlt = non_linear_temp(U,T)
    %% Four terms to evaluate for d(uu) and d(vu) for X-direction
    
    global Nx Ny dx dy
    
    nlt = CellData(Nx,Ny);
    
    T1 = interpol(EdgeData(Nx,Ny),T,1);
    uT = EdgeData(Nx,Ny);
    uT.x = U.x .* T1.x;
    duTdx = div(CellData(Nx,Ny),uT,1);
    duTdx.x = duTdx.x/dx;
    
    T2 = interpol(EdgeData(Nx,Ny),T,2);
    vT = EdgeData(Nx,Ny);
    vT.y = U.y .* T2.y;
    dvTdy = div(CellData(Nx,Ny),vT,2);
    dvTdy.x = dvTdy.x/dy;
    
    nlt.x = duTdx.x + dvTdy.x;
    nlt.x(:,1) = 0;
    nlt.x(:,end) = 0;
    nlt.x(1,:) = 0;
    nlt.x(end,:) = 0;
end