function op = H_op(params,xi,eta,x,y,ip)
    %H_OP Regularizes the data from the Lagrangian points to the field.
    % Refer to reference for further explanation.
    %
    % op = H_op(params,xi,eta,x,y,ip)
    %
    % Variable lookup:
    %
    % op: output on the flow field.
    %
    % params: flow parameters.
    %
    % xi: X-coordinate of the Lagrangian points.
    %
    % eta: Y-corrdinate of the Lagrangian points.
    %
    % x: X-coordinate of the flow field
    %
    % y: Y-coordinate of the flow field
    %
    % ip: input on the Lagrangian points.
    %
    % Created by Jay Mehta (18 July 2019)
    
    op = zeros(length(ip),1);
    
    for k = 1:length(ip)
        if ip(k) == 0 || (x-xi(k))/params.dx > 1.5 || (y-eta(k))/params.dx > 1.5
            op(k) = 0;
        else
            op(k) = ip(k) * ddf_roma_2D(params,x-xi(k),y-eta(k));
        end
    end
    op = sum(op);
end