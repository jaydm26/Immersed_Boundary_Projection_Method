function T = apply_bc_temp(params,T,T0,velocity)
    %APPLY_BC_TEMP  Applies Boundary Condtions of Inlet Temperature to the 
    % left wall, no temperature gradient on the top and bottom walls in the 
    % vertical direction, and outlet to the right wall of a data field 
    % stored in the Cell Space. This function is specifically created for 
    % temperature related problems in the Immersed Boundary Projection 
    % Method. Refer to the reference for further explanation.
    % 
    % This works for three cases:
    %
    % If two argument is provided, the convective outlet is replaced
    % by the solid wall.
    % 
    % If three arguments are provided, the convective outlet uses the
    % reference velocity U (through the Courant Number) as the convecting
    % velocity.
    % 
    % If four arguments are provided, the convective outlet uses the
    % X-velocity field as the convecting velocity. Note that this is
    % unstable. *Do not use*
    %
    % T = APPLY_BC_TEMP(params,T,T0,velocity)
    %
    % Variable lookup:
    % 
    % T: Temperature field (CellData) on which the boundary conditions
    % are being applied.
    %
    % params: flow parameters.
    %
    % T0: Temperature field (CellData) from the previous time step. Used
    % to calculate the convective outlet boundary condition.
    %
    % velocity: Current Velocity field (EdgeData). Used to calculate the
    % convective outlet boundary condition.
    %
    % Nx: Number of divisions in the X-direction.
    %
    % Ny: Number of divisions in the Y-direction.
    %
    % dt: Size of each time step moved during one iteration.
    % 
    % dx: Size of each cell. *Ensure that dx = dy. The code will give
    % erroneous results otherwise.*
    %
    % Co: Courant Number defined by Co = U * dt / dx. Refer to reference
    % for further explanation.
    %
    % Created by Jay Mehta (18 July 2019)
    
    if T.data ~= "cell"
        error("T is not a Cell field.")
    end
    if T0.data ~= "cell"
        error("T0 is not a Cell field.")
    end
    if velocity.data ~= "edge"
        error("velocity is not an Edge field.")
    end
    
    Nx = T.size(1);
    Ny = T.size(2);
    
    switch nargin
        case 2
            T.x(:,1) = 0;
            T.x(1,2:Ny+1) = T.x(2,2:Ny+1);
            T.x(Nx+2,2:Ny+1) = T.x(Nx+1,2:Ny+1);
            T.x(:,Ny+2) = 0;
        case 3
            T.x(:,1) = 0;
            T.x(1,2:Ny+1) = T.x(2,2:Ny+1);
            T.x(Nx+2,2:Ny+1) = T.x(Nx+1,2:Ny+1);
            for i = 1:Nx+2
                for j = Ny+2
                    T.x(i,j) = T0.x(i,j) - params.Co * (T0.x(i,j)-T0.x(i,j-1));
                end
            end
        case 4
            T.x(:,1) = 0;
            T.x(1,2:Ny+1) = T.x(2,2:Ny+1);
            T.x(Nx+2,2:Ny+1) = T.x(Nx+1,2:Ny+1);
            U_inf = interpol(NodeData(Nx,Ny),velocity,1);
            U_inf = U_inf.x';
            for i = 1:Nx+2
                for j = Ny+2
                    T.x(i,j) = T0.x(i,j) - U_inf(i,j-1) * params.dt/params.dx * (T0.x(i,j)-T0.x(i,j-1));
                end
            end
    end
end
    