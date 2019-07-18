function nlt = non_linear(params,domain,velocity)
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
    dx = params.dx;
    dy = params.dx;
    
    nlt = EdgeData(Nx,Ny);
    
    u = interpol(velocity,CellData(Nx,Ny),1);
    u.x = u.x .* u.x;
    du2dx = div(u,EdgeData(Nx,Ny),1);
    du2dx.x = du2dx.x/dx;
    
    un = interpol(velocity,NodeData(Nx,Ny),1);
    vn = interpol(velocity,NodeData(Nx,Ny),2);
    uv = NodeData(Nx,Ny);
    uv.x = un.x .* vn.x;
    duvdy = div(uv,EdgeData(Nx,Ny),2);
    duvdy.x = duvdy.x/dy;
    
    duvdx = div(uv,EdgeData(Nx,Ny),1);
    duvdx.y = duvdx.y/dx;
    
    v = interpol(velocity,CellData(Nx,Ny),2);
    v.x = v.x .* v.x;
    dv2dy = div(v,EdgeData(Nx,Ny),2);
    dv2dy.y = dv2dy.y/dy;
    
    nlt.x = du2dx.x + duvdy.x;
    nlt.x(:,1) = 0;
    nlt.x(:,end) = 0;
    nlt.y = duvdx.y + dv2dy.y;
    nlt.y(1,:) = 0;
    nlt.y(end,:) = 0;
end