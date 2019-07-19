function op = laplacian_2(ip,dir)
    %LAPLACIAN_2 Computes the laplacian of the input field. If two 
    % arguments are provided, it will compute a directional laplacian.
    %
    % op = laplacian_2(ip,dir)
    %
    % Variable lookup:
    %
    % op: laplacian of the input field.
    %
    % ip: input field.
    %
    % dir: direction of the laplacian- 1: X-direction
    %                                  2: Y-direction
    %
    % Created by Jay Mehta (18 July 2019)
    
    if nargin == 1
        switch ip.data
            case 'cell'
                Nx = ip.size(1);
                Ny = ip.size(2);
                op = CellData(Nx,Ny);
                for i = 2:Nx+1
                    for j = 2:Ny+1
                        op.x(i,j) = ip.x(i-1,j) + ip.x(i+1,j) + ip.x(i,j+1) + ip.x(i,j-1) - 4 * ip.x(i,j);
                    end
                end
            case 'node'
                Nx = ip.size(1);
                Ny = ip.size(2);
                op = NodeData(Nx,Ny);
                for i = 2:Nx
                    for j = 2:Ny
                        op.x(i,j) = -4*ip.x(i,j) + ip.x(i-1,j) + ip.x(i+1,j) + ip.x(i,j+1) + ip.x(i,j-1);
                    end
                end
            case 'edge'
                Nx = ip.size(1);
                Ny = ip.size(2);
                op = EdgeData(Nx,Ny);
                for i = 2:Nx
                    for j = 2:Ny+1
                        op.x(i,j) = -4*ip.x(i,j) + ip.x(i-1,j) + ip.x(i+1,j) + ip.x(i,j+1) + ip.x(i,j-1);
                    end
                end
                for i = 2:Nx+1
                    for j = 2:Ny
                        op.y(i,j) = -4*ip.y(i,j) + ip.y(i-1,j) + ip.y(i+1,j) + ip.y(i,j+1) + ip.y(i,j-1);
                    end
                end
        end
    elseif nargin == 2
        switch ip.data
            case 'cell'
                Nx = ip.size(1);
                Ny = ip.size(2);
                op = CellData(Nx,Ny);
                if dir == 1
                    for i = 2:Nx+1
                        for j = 2:Ny+1
                            op.x(i,j) = ip.x(i-1,j) + ip.x(i+1,j) - 2 * ip.x(i,j);
                        end
                    end
                elseif dir == 2
                    for i = 2:Nx+1
                        for j = 2:Ny+1
                            op.x(i,j) = ip.x(i,j-1) + ip.x(i,j+1) - 2 * ip.x(i,j);
                        end
                    end
                else
                    error("Cannot operate in a dimension greater than 2")
                end
            case 'node'
                Nx = ip.size(1);
                Ny = ip.size(2);
                op = NodeData(Nx,Ny);
                if dir == 1
                    for i = 2:Nx
                        for j = 2:Ny
                            op.x(i,j) = -2*ip.x(i,j) + ip.x(i-1,j) + ip.x(i+1,j);
                        end
                    end
                elseif dir == 2
                    for i = 2:Nx
                        for j = 2:Ny
                            op.x(i,j) = -2*ip.x(i,j) + ip.x(i,j-1) + ip.x(i,j+1);
                        end
                    end
                else
                    error("Cannot operate in a dimension greater than 2")
                end
            case 'edge'
                Nx = ip.size(1);
                Ny = ip.size(2);
                op = EdgeData(Nx,Ny);
                if dir == 1
                    for i = 2:Nx
                        for j = 2:Ny+1
                            op.x(i,j) = -2*ip.x(i,j) + ip.x(i-1,j) + ip.x(i+1,j);
                        end
                    end
                    for i = 2:Nx+1
                        for j = 2:Ny
                            op.y(i,j) = -2*ip.y(i,j) + ip.y(i-1,j) + ip.y(i+1,j);
                        end
                    end
                elseif dir == 2
                    for i = 2:Nx
                        for j = 2:Ny+1
                            op.x(i,j) = -2*ip.x(i,j) + ip.x(i,j-1) + ip.x(i,j+1);
                        end
                    end
                    for i = 2:Nx+1
                        for j = 2:Ny
                            op.y(i,j) = -2*ip.y(i,j) + ip.y(i,j-1) + ip.y(i,j+1);
                        end
                    end
                else
                    error("Cannot operate in a dimension greater than 2")
                end
        end
    end
end