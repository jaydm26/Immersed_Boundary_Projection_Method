function gamma = apply_bc_sp(params,gamma,gamma0,velocity)
    %APPLY_BC_SP Applies Boundary Condtions of Inlet to the left wall, no 
    % traction to the top and bottom wall, and convective outlet to the 
    % right wall of a data field stored in Node Space. This function is 
    % specifically created for the Nullspace Method. Refer to the reference
    % for further explanation.
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
    % gamma = APPLY_BC_SP(params,gamma,gamma0,velocity)
    %
    % Variable lookup:
    % 
    % gamma: Vorticity field (NodeData) on which the boundary conditions
    % are being applied.
    %
    % params: flow parameters.
    %
    % gamma0: Vorticity field (NodeData) from the previous time step. Used
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
    
    if nargin == 2 && gamma.data ~= "node"
        error("gamma is not a Node field.")
    end
    if nargin == 3 && gamma0.data ~= "node"
        error("gamma0 is not a Node field.")
    end
    if nargin == 4 && velocity.data ~= "edge"
        error("velocity is not an Edge field.")
    end
    Nx = gamma.size(1);
    Ny = gamma.size(2);
    
    switch nargin
        case 2
            gamma.x(:,1) = 0;
            gamma.x(1,2:Ny) = gamma.x(2,2:Ny);
            gamma.x(Nx+1,2:Ny) = gamma.x(Nx,2:Ny);
            gamma.x(:,Ny+1) = 0;
        case 3
            gamma.x(:,1) = 0;
            gamma.x(1,2:Ny) = gamma.x(2,2:Ny);
            gamma.x(Nx+1,2:Ny) = gamma.x(Nx,2:Ny);
            for i = 1:Nx+1
                for j = Ny+1
                    gamma.x(i,j) = gamma0.x(i,j) - params.Co * (gamma0.x(i,j)-gamma0.x(i,j-1));
                end
            end
        case 4
            gamma.x(:,1) = 0;
            gamma.x(1,2:Ny) = gamma.x(2,2:Ny);
            gamma.x(Nx+1,2:Ny) = gamma.x(Nx,2:Ny);
            U_inf = interpol(NodeData(Nx,Ny),velocity,1);
            U_inf = U_inf.x';
            for i = 1:Nx+1
                for j = Ny+1
                    gamma.x(i,j) = gamma0.x(i,j) - U_inf(i,j-1) * params.dt/params.dx * (gamma0.x(i,j)-gamma0.x(i,j-1));
                end
            end
    end
end