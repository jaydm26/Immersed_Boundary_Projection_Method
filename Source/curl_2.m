function A = curl(U,dim)
    switch U.data 
        case 'edge' 
            Nx = U.size(1);
            Ny = U.size(2);
            A = NodeData(Nx,Ny);
            for i = 1:Nx+1
                for j = 1:Ny+1
                    A.x(i,j) = U.x(i,j) - U.x(i,j+1) + U.y(i+1,j) - U.y(i,j);
                end
            end
        case 'node'
            if nargin == 1
                Nx = U.size(1);
                Ny = U.size(2);
                A = EdgeData(Nx,Ny);
                for i = 1:Nx+1
                    for j = 2:Ny+1
                        A.x(i,j) = U.x(i,j) - U.x(i,j-1);
                    end
                end
                for i = 2:Nx+1
                    for j = 1:Ny+1
                        A.y(i,j) = U.x(i-1,j) - U.x(i,j);
                    end
                end
            elseif nargin == 2
                Nx = U.size(1);
                Ny = U.size(2);
                A = EdgeData(Nx,Ny);
                if dim == 1
                    for i = 1:Nx+1
                        for j = 2:Ny+1
                            A.x(i,j) = U.x(i,j) - U.x(i,j-1);
                        end
                    end
                elseif dim == 2
                    for i = 2:Nx+1
                        for j = 1:Ny+1
                            A.y(i,j) = U.x(i,j) - U.x(i-1,j);
                        end
                    end
                end
            end
    end    
end