function op = CTH(params,domain,xi,eta,Fx,Fy)
    %CTH Executes the C^T H operation on the Lagrangian body forces.
    %
    % op = CTH(params,domain,xi,eta,Fx,Fy)
    %
    % Variable lookup:
    %
    % params: flow parameters.
    %
    % domain: domain parameters.
    %
    % xi: X-coordinate of the Lagrangian points.
    %
    % eta: Y-corrdinate of the Lagrangian points.
    %
    % Fx: Lagrangian Forces in the X-direction.
    %
    % Fy: Lagrangian Forces in the Y-direction.
    %
    % op: Result of the C^T H operation. H places data on the Edge
    % field. Curl of Edge field yields a Node field.
    %
    % Created by Jay Mehta (18 July 2019)
    
    q = H_operation(params,domain,"edge",xi,eta,Fx,Fy);
    
    op = curl_2(q);
end