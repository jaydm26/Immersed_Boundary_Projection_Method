function nltt = non_linear_temp(params,velocity,T)
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
    
    %% Four terms to evaluate for d(uT) and d(vT) for X-direction
    
    Nx = T.size(1);
    Ny = T.size(2);
    dx = params.dx;
    dy = params.dy;
    
   % d(vT) = v.dT + T.dv
    % d(uT) = u.dT + T.du. They are calculated differently to maintain
    % conservativeness of the equations in discrete space.
    
    t1_1 = interpol(velocity,NodeData(Nx,Ny),2);
    t1_2 = div(T,NodeData(Nx,Ny),1);
    t1_2.x = t1_2.x/dy;
    t1 = NodeData(Nx,Ny);
    t1.x = t1_1.x .* t1_2.x;
    t1 = interpol(t1,EdgeData(Nx,Ny),1); % v.dT
    
    t2_1 = div(velocity,CellData(Nx,Ny),2);
    t2_1.x = t2_1.x/dy;
    t2_1 = interpol(t2_1,EdgeData(Nx,Ny),1);
    t2_2 = T.x;
    t2 = EdgeData(Nx,Ny);
    t2.x = t2_1.x .* t2_2; % T.dv
    
    t3_1 = interpol(velocity,CellData(Nx,Ny),1);
    t3_2 = div(T,CellData(Nx,Ny),1);
    t3_2.x = t3_2.x/dx;
    t3 = CellData(Nx,Ny);
    t3.x = t3_1.x .* t3_2.x;
    t3 = interpol(t3,EdgeData(Nx,Ny),1); % u.dT

    t4_1 = div(velocity,CellData(Nx,Ny),1);
    t4_1.x = t4_1.x/dx;
    t4_1 = interpol(t4_1,EdgeData(Nx,Ny),1);
    t4_2 = T.x;
    t4 = EdgeData(Nx,Ny);
    t4.x = t4_1.x .* t4_2; % T.du
    
    nltt = EdgeData(Nx,Ny);
    nltt.x = t1.x + t2.x + t3.x + t4.x;
end