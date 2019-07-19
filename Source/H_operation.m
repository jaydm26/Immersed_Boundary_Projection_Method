function op = H_operation(params,domain,DataType,xi,eta,ip1,ip2)
    %H_OPERATION Executes the regularization operation to take data from 
    % the Lagrangian points to the flow field. Refer to reference for 
    % further explanation.
    %
    % op = H_operation(params,domain,DataType,xi,eta,ip1,ip2)
    %
    % Variable lookup:
    %
    % op: output on the flow field
    %
    % params: flow parameters.
    %
    % domain: data structure containing all domains.
    %
    % DataType: Data Type on which the data is to be regularized.
    %   "cell" to regularize onto Cell Space.
    %   "edge" to regularize onto Edge Space (requires 2 inputs).
    %   "node" to regularize onto Node Space.
    %
    % xi: X-coordinate of the Lagrangian points.
    %
    % eta: Y-corrdinate of the Lagrangian points.
    %
    % ip1, ip2: inputs on the Lagrangian points.
    %
    % Created by Jay Mehta (18 July 2019)
    
    Nx = domain.Nx;
    Ny = domain.Ny;
    X_n = domain.X_n;
    Y_n = domain.Y_n;
    X_e_x = domain.X_e_x;
    Y_e_x = domain.Y_e_x;
    X_e_y = domain.X_e_y;
    Y_e_y = domain.Y_e_y;
    X_c = domain.X_c;
    Y_c = domain.Y_c;
    
    switch DataType
        case "node"
            op = NodeData(Nx,Ny);
            X_n = X_n';
            Y_n = Y_n';
            for i = 1:Nx+1
                for j = 1:Ny+1
                    op.x(i,j) = H_op(params,xi,eta,X_n(i,j),Y_n(i,j),ip1);
                end
            end
%             X_n = X_n';
%             Y_n = Y_n';
        case "cell"
            op = CellData(Nx,Ny);
            X_c = X_c';
            Y_c = Y_c';
            for i = 1:Nx+2
                for j = 1:Ny+2
                    op.x(i,j) = H_op(params,xi,eta,X_c(i,j),Y_c(i,j),ip1);
                end
            end
%             X_c = X_c';
%             Y_c = Y_c';
        case "edge"
            if nargin == 7
                op = EdgeData(Nx,Ny);

                X_e_x = X_e_x';
                Y_e_x = Y_e_x';

                for i = 1:Nx+1
                    for j = 1:Ny+2
                        op.x(i,j) = H_op(params,xi,eta,X_e_x(i,j),Y_e_x(i,j),ip1);
                    end
                end

                X_e_y = X_e_y';
                Y_e_y = Y_e_y';

                for i = 1:Nx+2
                    for j = 1:Ny+1
                        op.y(i,j) = H_op(params,xi,eta,X_e_y(i,j),Y_e_y(i,j),ip2);
                    end
                end
%                 X_e_x = X_e_x';
%                 Y_e_x = Y_e_x';
%                 X_e_y = X_e_y';
%                 Y_e_y = Y_e_y';
            else
                error("Too many inputs.")
            end
    end
end