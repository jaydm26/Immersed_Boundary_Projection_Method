function out = H_operation(DataType,Fx,Fy)
    
    global Nx Ny X_n Y_n X_e_x Y_e_x X_e_y Y_e_y X_c Y_c
    
    switch DataType
        case "node"
            out = NodeData(Nx,Ny);
            for i = 1:Nx+1
                for j = 1:Ny+1
                    out.x(i,j) = H_op(X_n(i,j),Y_n(i,j),Fx);
                end
            end
            
        case "cell"
            out = CellData(Nx,Ny);
            for i = 1:Nx+2
                for j = 1:Ny+2
                    out.x(i,j) = H_op(X_c(i,j),Y_c(i,j),Fx);
                end
            end
            
        case "edge"
            if nargin == 3
                out = EdgeData(Nx,Ny);

                X_e_x = X_e_x';
                Y_e_x = Y_e_x';

                for i = 1:Nx+1
                    for j = 1:Ny+2
                        out.x(i,j) = H_op(X_e_x(i,j),Y_e_x(i,j),Fx);
                    end
                end

                X_e_y = X_e_y';
                Y_e_y = Y_e_y';

                for i = 1:Nx+2
                    for j = 1:Ny+1
                        out.y(i,j) = H_op(X_e_y(i,j),Y_e_y(i,j),Fy);
                    end
                end

                X_e_x = X_e_x';
                Y_e_x = Y_e_x';
                X_e_y = X_e_y';
                Y_e_y = Y_e_y';
            end
    end
end