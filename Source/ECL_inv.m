function [x,y] = ECL_inv(params,domain,g_hat,xi,eta,gamma)
    %ECL_INV Executes the ECL^-1 operation on Node Data. Refer to reference for
    % further explanation.
    %
    % [x,y] = ECL_inv(params,domain,xi,eta,gamma)
    %
    % Variable lookup:
    %
    % x, y : output on the Lagrangian points.
    %
    % params: flow parameters.
    %
    % domain: data structure containing all domains.
    %
    % xi: X-coordinate of the Lagrangian points.
    %
    % eta: Y-corrdinate of the Lagrangian points.
    %
    % gamma: input vorticity.
    
    sf = L_inv_operation(gamma.x,g_hat);
     
    q = curl_2(sf);
    
    [x,y] = E_operation(params,domain,xi,eta,q);
end