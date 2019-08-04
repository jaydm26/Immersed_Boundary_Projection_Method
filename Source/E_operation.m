function [op1,op2] = E_operation(params,domain,xi,eta,ip)
    %E_OPERATION Executes the interpolation operation to take data from the
    % field on to the Lagrangian points. Refer to reference for further 
    % explanation.
    %
    % [op1,op2] = E_operation(params,domain,xi,eta,ip)
    %
    % Variable lookup:
    %
    % op1, op2 : output on the Lagrangian points.
    %
    % params: flow parameters.
    %
    % domain: data structure containing all domains.
    %
    % xi: X-coordinate of the Lagrangian points.
    %
    % eta: Y-corrdinate of the Lagrangian points.
    %
    % ip: input
    %
    % Created by Jay Mehta (18 July 2019)
    
    switch ip.data
        case "node"
            op1 = zeros(length(xi),1);
            for k = 1:length(xi)
                op1(k,1) = E_op(params,domain,xi(k),eta(k),ip);
            end
        case "cell"
            op1 = zeros(length(xi),1);
            for k = 1:length(xi)
                op1(k,1) = E_op(params,domain,xi(k),eta(k),ip);
            end
        case "edge"
            op1 = zeros(length(xi),1);
            op2 = zeros(length(xi),1);
            for k = 1:length(xi)
                op1(k,1) = E_op(params,domain,xi(k),eta(k),ip,1);
                op2(k,1) = E_op(params,domain,xi(k),eta(k),ip,2);
            end
    end
end
            
    