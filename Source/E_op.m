function op = E_op(params,domain,xi,eta,ip,dim)
    %E_OP Interpolates the field variable to the Lagrangian Points. Refer 
    % to the reference for further explanation.
    %
    % out = E_op(xi,eta,s,dim)
    %
    % Variable lookup:
    %
    % params: flow parameters.
    %
    % op: output of the flow field on the Lagrangian points.
    %
    % domain: data structure containing all domains.
    %
    % xi: X-coordinate of the Lagrangian points.
    %
    % eta: Y-corrdinate of the Lagrangian points.
    %
    % ip: input
    %
    % dim: dimension along which the Edge field is operated.
    %   1: X field of the Edge Field
    %   2: Y field of the Edge Field
    %
    % Created by Jay Mehta (18 July 2019)
    
    Nx = ip.size(1);
    Ny = ip.size(2);
    
    X_n = domain.X_n;
    Y_n = domain.Y_n;
    X_e_x = domain.X_e_x;
    Y_e_x = domain.Y_e_x;
    X_e_y = domain.X_e_y;
    Y_e_y = domain.Y_e_y;
    X_c = domain.X_c;
    Y_c = domain.Y_c;
    
    if ip.data == "node"
        if nargin == 5
            X_n = X_n';
            Y_n = Y_n';
            for i = 1:Nx+1
                for j = 1:Ny+1
                    if ip.x(i,j) == 0 || (X_n(i,j)-xi)/params.dx > 1.5 || (Y_n(i,j)-eta)/params.dx > 1.5
                        op(i,j) = 0;
                    else
                        op(i,j) = ip.x(i,j) * ddf_roma_2D(params,X_n(i,j)-xi,Y_n(i,j)-eta);
                    end
                end
            end
            op = sum(sum(op));
%             X_n = X_n';
%             Y_n = Y_n';
        else
            error("Too many inputs.")
        end
    elseif ip.data == "cell"
        if nargin == 5
            X_c = X_c';
            Y_c = Y_c';
            for i = 1:Nx+2
                for j = 1:Ny+2
                    if ip.x(i,j) == 0 || (X_c(i,j)-xi)/params.dx > 1.5 || (Y_c(i,j)-eta)/params.dx > 1.5
                        op(i,j) = 0;
                    else
                        op(i,j) = ip.x(i,j) * ddf_roma_2D(params,X_c(i,j)-xi,Y_c(i,j)-eta);
                    end
                end
            end
            op = sum(sum(op));
%             X_c = X_c';
%             Y_c = Y_c';
        else
            error("Too many inputs.")
        end
    elseif ip.data == "edge"
        if nargin == 6
            if dim == 1
                X_e_x = X_e_x';
                Y_e_x = Y_e_x';
                for i = 1:Nx+1
                    for j = 1:Ny+2
                        if ip.x(i,j) == 0 || (X_e_x(i,j)-xi)/params.dx > 1.5 || (Y_e_x(i,j)-eta)/params.dx > 1.5
                            op(i,j) = 0;
                        else
                            op(i,j) = ip.x(i,j) * ddf_roma_2D(params,X_e_x(i,j)-xi,Y_e_x(i,j)-eta);
                        end
                    end
                end
                op = sum(sum(op));
%                 X_e_x = X_e_x';
%                 Y_e_x = Y_e_x';
            elseif dim == 2
                X_e_y = X_e_y';
                Y_e_y = Y_e_y';                
                for i = 1:Nx+2
                    for j = 1:Ny+1
                        if ip.y(i,j) == 0 || (X_e_y(i,j)-xi)/params.dx > 1.5 || (Y_e_y(i,j)-eta)/params.dx > 1.5
                            op(i,j) = 0;
                        else
                            op(i,j) = ip.y(i,j) * ddf_roma_2D(params,X_e_y(i,j)-xi,Y_e_y(i,j)-eta);
                        end
                    end
                end
                op = sum(sum(op));
%                 X_e_y = X_e_y';
%                 Y_e_y = Y_e_y';
            else
                error("Dimension exceeds 2.")
            end
        else
            error("Incorrect inputs.")
        end
    else
        error("Data Type Not Supported.")
    end
end