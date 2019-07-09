function A = laplacian(U,dir)
    if nargin == 1
        switch U.data
            case 'cell'
                Nx = U.size(1);
                Ny = U.size(2);
                A = CellData(Nx,Ny);
                for i = 2:Nx+1
                    for j = 2:Ny+1
                        A.x(i,j) = U.x(i-1,j) + U.x(i+1,j) + U.x(i,j+1) + U.x(i,j-1) - 4 * U.x(i,j);
                    end
                end
            case 'node'
                Nx = U.size(1);
                Ny = U.size(2);
                A = NodeData(Nx,Ny);
                for i = 2:Nx
                    for j = 2:Ny
                        A.x(i,j) = -4*U.x(i,j) + U.x(i-1,j) + U.x(i+1,j) + U.x(i,j+1) + U.x(i,j-1);
                    end
                end
            case 'edge'
                Nx = U.size(1);
                Ny = U.size(2);
                A = EdgeData(Nx,Ny);
                for i = 2:Nx
                    for j = 2:Ny+1
                        A.x(i,j) = -4*U.x(i,j) + U.x(i-1,j) + U.x(i+1,j) + U.x(i,j+1) + U.x(i,j-1);
                    end
                end
                for i = 2:Nx+1
                    for j = 2:Ny
                        A.y(i,j) = -4*U.y(i,j) + U.y(i-1,j) + U.y(i+1,j) + U.y(i,j+1) + U.y(i,j-1);
                    end
                end
        end
    elseif nargin == 2
        switch U.data
            case 'cell'
                Nx = U.size(1);
                Ny = U.size(2);
                A = CellData(Nx,Ny);
                if dir == 1
                    for i = 2:Nx+1
                        for j = 2:Ny+1
                            A.x(i,j) = U.x(i-1,j) + U.x(i+1,j) - 2 * U.x(i,j);
                        end
                    end
                elseif dir == 2
                    for i = 2:Nx+1
                        for j = 2:Ny+1
                            A.x(i,j) = U.x(i,j-1) + U.x(i,j+1) - 2 * U.x(i,j);
                        end
                    end
                else
                    error("Cannot operate in a dimension greater than 2")
                end
            case 'node'
                Nx = U.size(1);
                Ny = U.size(2);
                A = NodeData(Nx,Ny);
                if dir == 1
                    for i = 2:Nx
                        for j = 2:Ny
                            A.x(i,j) = -2*U.x(i,j) + U.x(i-1,j) + U.x(i+1,j);
                        end
                    end
                elseif dir == 2
                    for i = 2:Nx
                        for j = 2:Ny
                            A.x(i,j) = -2*U.x(i,j) + U.x(i,j-1) + U.x(i,j+1);
                        end
                    end
                else
                    error("Cannot operate in a dimension greater than 2")
                end
            case 'edge'
                Nx = U.size(1);
                Ny = U.size(2);
                A = EdgeData(Nx,Ny);
                if dir == 1
                    for i = 2:Nx
                        for j = 2:Ny+1
                            A.x(i,j) = -2*U.x(i,j) + U.x(i-1,j) + U.x(i+1,j);
                        end
                    end
                    for i = 2:Nx+1
                        for j = 2:Ny
                            A.y(i,j) = -2*U.y(i,j) + U.y(i-1,j) + U.y(i+1,j);
                        end
                    end
                elseif dir == 2
                    for i = 2:Nx
                        for j = 2:Ny+1
                            A.x(i,j) = -2*U.x(i,j) + U.x(i,j-1) + U.x(i,j+1);
                        end
                    end
                    for i = 2:Nx+1
                        for j = 2:Ny
                            A.y(i,j) = -2*U.y(i,j) + U.y(i,j-1) + U.y(i,j+1);
                        end
                    end
                else
                    error("Cannot operate in a dimension greater than 2")
                end
        end
    end
end