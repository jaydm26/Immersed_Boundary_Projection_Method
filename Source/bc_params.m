function bc = bc_params(uL,uR,uB,uT,vL,vR,vB,vT)
    %
    % Creates a data structure to pass to the apply_bc function. This data
    % structure contains the data for the boundary conditions.
    %
    % Variable lookup:
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
    
    bc = struct;
    
    bc.uL = uL;
    bc.uR = uR;
    bc.uB = uB;
    bc.uT = uT;
    
    bc.vL = vL;
    bc.vR = vR;
    bc.vB = vB;
    bc.vT = vT;
end