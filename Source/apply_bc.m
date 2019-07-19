function velocity = apply_bc(bc,velocity,t)
    %APPLY_BC Applies Dirichlet or Neumann Boundary conditions to the 
    % velocity fields. Note that the boundary conditions must be predefined
    % and declared as global variables.
    %
    % It is important to note here that these conditions have to be defined
    % by the user. Examples available in the reference.
    %
    % velocity = APPLY_BC(velocity,t,bc)
    %
    % Variable lookup:
    %
    % velocity: Velocity Field (EdgeData) on which the boundary conditions are
    % being applied.
    %
    % t: time at which the boundary conditons are being applied. Used when
    % boundary condition is time dependent.
    %
    % Nx: Number of divisions in the X-direction.
    %
    % Ny: Number of divisions in the Y-direction.
    %
    % uL: Boundary Condition on the left wall for velocity in the 
    % X-direction.
    %
    % uR: Boundary Condition on the right wall for velocity in the 
    % X-direction.
    %
    % uB: Boundary Condition on the bottom wall for velocity in the 
    % X-direction.
    %
    % uT: Boundary Condition on the top wall for velocity in the 
    % X-direction.
    %
    % vL: Boundary Condition on the left wall for velocity in the 
    % Y-direction.
    %
    % vR: Boundary Condition on the right wall for velocity in the 
    % Y-direction.
    %
    % vB: Boundary Condition on the bottom wall for velocity in the 
    % Y-direction.
    %
    % vT: Boundary Condition on the top wall for velocity in the 
    % Y-direction.
    %
    % Created by Jay Mehta (18 July 2019)
    
    Nx = velocity.size(1);
    Ny = velocity.size(2);
    
    switch velocity.data
        case "edge"
            for i = 2:Ny+1
                velocity.x(1,i)    = bc.uL(i,t);
                velocity.x(Nx+1,i) = bc.uR(i,t);
            end
            for i = 2:Nx
                velocity.x(i,1)    = -velocity.x(i,2)     + 2*bc.uB(i,t);
                velocity.x(i,Ny+2) = -velocity.x(i,Ny+1)  + 2*bc.uT(i,t);
            end
            for i = 2:Ny
                velocity.y(1,i)    = -velocity.y(2,i)     + 2*bc.vL(i,t);
                velocity.y(Nx+2,i) = -velocity.y(Nx+1,i)  + 2*bc.vR(i,t);
            end
            for i = 2:Nx+1
                velocity.y(i,1)    = bc.vB(i,t);
                velocity.y(i,Ny+1) = bc.vT(i,t);
            end
        otherwise
            error("U is not a Edge field.")
    end
end