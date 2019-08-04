function nlt = non_linear(params,velocity)
    %NON_LINEAR Computes the non-linear term u.du/dx + v.du/dy and 
    % u.dv/dx + v.dv/dy present in the Navier-Stokes equation.
    %
    % nlt = non_linear(params,velocity)
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
    
    % d(vu) = v.du + u.dv
    % d(uu) = u.du + u.du. They are calculated differently to maintain
    % conservativeness of the equations in discrete space.
    
    t1_1 = interpol(velocity,NodeData(Nx,Ny),2);
    t1_2 = div(velocity,NodeData(Nx,Ny),1);
    t1_2.x = t1_2.x/dy;
    t1 = NodeData(Nx,Ny);
    t1.x = t1_1.x .* t1_2.x;
    t1 = interpol(t1,EdgeData(Nx,Ny),1); % v.du
    
    t2_1 = div(velocity,CellData(Nx,Ny),2);
    t2_1.x = t2_1.x/dy;
    t2_1 = interpol(t2_1,EdgeData(Nx,Ny),1);
    t2_2 = velocity.x;
    t2 = EdgeData(Nx,Ny);
    t2.x = t2_1.x .* t2_2; % u.dv
    
    t3_1 = interpol(velocity,CellData(Nx,Ny),1);
    t3_2 = div(velocity,CellData(Nx,Ny),1);
    t3_2.x = t3_2.x/dx;
    t3 = CellData(Nx,Ny);
    t3.x = t3_1.x .* t3_2.x;
    t3 = interpol(t3,EdgeData(Nx,Ny),1); % u.du

    t4_1 = div(velocity,CellData(Nx,Ny),1);
    t4_1.x = t4_1.x/dx;
    t4_1 = interpol(t4_1,EdgeData(Nx,Ny),1);
    t4_2 = velocity.x;
    t4 = EdgeData(Nx,Ny);
    t4.x = t4_1.x .* t4_2; % u.du
    
    nlt = EdgeData(Nx,Ny);
    nlt.x = t1.x + t2.x + t3.x + t4.x;
    
    %% Four terms to evaluate for d(uv) and d(vv) for Y-direction
    
    % d(uv) = u.dv + v.du
    % d(vv) = v.dv + v.dv
    
    t1_1 = interpol(velocity,NodeData(Nx,Ny),1);
    t1_2 = div(velocity,NodeData(Nx,Ny),2);
    t1_2.x = t1_2.x/dx;
    t1 = NodeData(Nx,Ny);
    t1.x = t1_1.x .* t1_2.x;
    t1 = interpol(t1,EdgeData(Nx,Ny),2); %u.dv
    
    t2_1 = div(velocity,CellData(Nx,Ny),1);
    t2_1.x = t2_1.x/dx;
    t2_1 = interpol(t2_1,EdgeData(Nx,Ny),2);
    t2_2 = velocity.y;
    t2 = EdgeData(Nx,Ny);
    t2.y = t2_1.y .* t2_2; % v.du
    
    t3_1 = interpol(velocity,CellData(Nx,Ny),2);
    t3_2 = div(velocity,CellData(Nx,Ny),2);
    t3_2.x = t3_2.x/dy;
    t3 = CellData(Nx,Ny);
    t3.x = t3_1.x .* t3_2.x;
    t3 = interpol(t3,EdgeData(Nx,Ny),2); % v.dv
    
    t4_1 = div(velocity,CellData(Nx,Ny),2);
    t4_1.x = t4_1.x/dy;
    t4_1 = interpol(t4_1,EdgeData(Nx,Ny),2);
    t4_2 = velocity.y;
    t4 = EdgeData(Nx,Ny);
    t4.y = t4_1.y .* t4_2; % v.dv
    
    nlt.y = t1.y + t2.y + t3.y + t4.y;
end