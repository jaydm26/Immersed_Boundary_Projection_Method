function nltt = non_linear_temp(domain,velocity,T)
    %NON_LINEAR_TEMP Computes the non-linear term duT/dx + dvT/dy present in the 
    % Energy equation.
    %
    % nltt = non_linear_temp(velocity,T)
    %
    % Variable lookup:
    %
    % velocity: Velocity field (EdgeData).
    %
    % T: Temperature field (CellData).
    %
    % domain: domain parameters.
    %
    % nltt: computed non-linear terms.
    %
    % Created by Jay Mehta (18 July 2019)
    
    %% Four terms to evaluate for d(uu) and d(vu) for X-direction
    
    Nx = domain.Nx;
    Ny = domain.Ny;
    dx = domain.dx;
    dy = domain.dy;
    
    nltt = CellData(Nx,Ny);
    
    T1 = interpol(EdgeData(Nx,Ny),T,1);
    uT = EdgeData(Nx,Ny);
    uT.x = velocity.x .* T1.x;
    duTdx = div(CellData(Nx,Ny),uT,1);
    duTdx.x = duTdx.x/dx;
    
    T2 = interpol(EdgeData(Nx,Ny),T,2);
    vT = EdgeData(Nx,Ny);
    vT.y = velocity.y .* T2.y;
    dvTdy = div(CellData(Nx,Ny),vT,2);
    dvTdy.x = dvTdy.x/dy;
    
    nltt.x = duTdx.x + dvTdy.x;
    nltt.x(:,1) = 0;
    nltt.x(:,end) = 0;
    nltt.x(1,:) = 0;
    nltt.x(end,:) = 0;
end