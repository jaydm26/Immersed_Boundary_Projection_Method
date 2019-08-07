function nlt = non_linear_alt(params,velocity)
    %NON_LINEAR_ALT Computes the non-linear term u.du/dx + v.du/dy and 
    % u.dv/dx + v.dv/dy present in the Navier-Stokes equation.
    %
    % nlt = non_linear_alt(params,velocity)
    %
    % Variable lookup:
    %
    % velocity: velocity field (EdgeData).
    %
    % params: flow parameters
    %
    % nlt: computed non-linear terms.
    %
    % Created by Jay Mehta (18 July 2019)
    
    %% Four terms to evaluate for d(uu) and d(vu) for X-direction
    
    Nx = velocity.size(1);
    Ny = velocity.size(2);
    dx = params.dx;
    dy = params.dx;
    
    nlt = EdgeData(Nx,Ny);
    
%     w = min(1.2 * params.dt * max(max(max(abs(velocity.x))),...
%         max(max(abs(velocity.y)))),1);
    w = 0;

    u = interpol(velocity,CellData(Nx,Ny),1);
    u_upwind = upwinding(velocity,CellData(Nx,Ny),1);
    u.x = u.x .* (u.x - w * u_upwind.x);
    du2dx = div(u,EdgeData(Nx,Ny),1);
    du2dx.x = du2dx.x/dx;
    
    un = interpol(velocity,NodeData(Nx,Ny),1);
    vn = interpol(velocity,NodeData(Nx,Ny),2);
    u_upwind = upwinding(velocity,NodeData(Nx,Ny),1);
    uv = NodeData(Nx,Ny);
    uv.x = (un.x - w * u_upwind.x) .* vn.x;
    duvdy = div(uv,EdgeData(Nx,Ny),2);
    duvdy.x = duvdy.x/dy;
    
    v_upwind = upwinding(velocity,NodeData(Nx,Ny),2);
    uv.x = un.x .* (vn.x - w * v_upwind.x);
    duvdx = div(uv,EdgeData(Nx,Ny),1);
    duvdx.y = duvdx.y/dx;
    
    v = interpol(velocity,CellData(Nx,Ny),2);
    v_upwind = upwinding(velocity,CellData(Nx,Ny),2);
    v.x = v.x .* (v.x - w * v_upwind.x);
    dv2dy = div(v,EdgeData(Nx,Ny),2);
    dv2dy.y = dv2dy.y/dy;
    
    nlt.x = du2dx.x + duvdy.x;
    nlt.x(:,1) = 0;
    nlt.x(:,end) = 0;
    nlt.y = duvdx.y + dv2dy.y;
    nlt.y(1,:) = 0;
    nlt.y(end,:) = 0;
end