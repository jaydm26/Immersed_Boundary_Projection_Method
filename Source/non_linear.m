function nlt = non_linear(domain,velocity)
    %NON_LINEAR Computes the non-linear term u.du/dx + v.du/dy and 
    % u.dv/dx + v.dv/dy present in the Navier-Stokes equation.
    %
    % nlt = non_linear(domain,velocity)
    %
    % Variable lookup:
    %
    % velocity: velocity field (EdgeData).
    %
    % domain: domain parameters.
    %
    % nlt: computed non-linear terms.
    
    %% Four terms to evaluate for d(uu) and d(vu) for X-direction
    
    Nx = domain.Nx;
    Ny = domain.Ny;
    dx = domain.dx;
    dy = domain.dy;
    
    nlt = EdgeData(Nx,Ny);
    
    u = interpol(CellData(Nx,Ny),velocity,1);
    u.x = u.x .* u.x;
    du2dx = div(EdgeData(Nx,Ny),u,1);
    du2dx.x = du2dx.x/dx;
    
    un = interpol(NodeData(Nx,Ny),velocity,1);
    vn = interpol(NodeData(Nx,Ny),velocity,2);
    uv = NodeData(Nx,Ny);
    uv.x = un.x .* vn.x;
    duvdy = div(EdgeData(Nx,Ny),uv,2);
    duvdy.x = duvdy.x/dy;
    
    duvdx = div(EdgeData(Nx,Ny),uv,1);
    duvdx.y = duvdx.y/dx;
    
    v = interpol(CellData(Nx,Ny),velocity,2);
    v.x = v.x .* v.x;
    dv2dy = div(EdgeData(Nx,Ny),v,2);
    dv2dy.y = dv2dy.y/dy;
    
    nlt.x = du2dx.x + duvdy.x;
    nlt.x(:,1) = 0;
    nlt.x(:,end) = 0;
    nlt.y = duvdx.y + dv2dy.y;
    nlt.y(1,:) = 0;
    nlt.y(end,:) = 0;
end