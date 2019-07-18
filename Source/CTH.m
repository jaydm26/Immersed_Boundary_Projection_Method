function gamma = CTH(Fx,Fy)
    %
    % Executes the C^T H operation on the Lagrangian body forces.
    %
    % Variable lookup:
    %
    % Fx: Lagrangian Forces in the X-direction.
    %
    % Fy: Lagrangian Forces in the Y-direction.
    %
    % gamma: Result of the C^T H operation. H places data on the Edge
    % field. Curl of Edge field yields a Node field.
    
    q = H_operation("edge",Fx,Fy);
    
    gamma = curl_2(q);
end