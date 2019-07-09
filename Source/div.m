function A = div(A,U,dim)
    
    global Nx Ny
    
    if nargin == 2
        if U.data ~= "edge"
            error("Check Data Type. Data Type can only be edge")
        else
            Nx = U.size(1);
            Ny = U.size(2);
            A = CellData(Nx,Ny);
            for i = 2:Nx+1
                for j = 2:Ny+1
                    A.x(i,j) = U.x(i,j) - U.x(i-1,j) + U.y(i,j) - U.y(i,j-1);
                end
            end
        end
    elseif nargin == 3
        if U.data == "edge" 
            switch lower(A.data)
                case "node"
                    Nx = U.size(1);
                    Ny = U.size(2);
                    if dim == 1
                        for i = 1:Nx+1
                            for j = 1:Ny+1
                                A.x(i,j) = U.x(i,j+1) - U.x(i,j);
                            end
                        end
                    elseif dim == 2
                        for i = 1:Nx+1
                            for j = 1:Ny+1
                                A.x(i,j) = U.y(i+1,j) - U.y(i,j);
                            end
                        end
                    end
                case "cell"
                    Nx = U.size(1);
                    Ny = U.size(2);
                    if dim == 1
                        for i = 2:Nx+1
                            for j = 1:Ny+2
                                A.x(i,j) = U.x(i,j) - U.x(i-1,j);
                            end
                        end
                    elseif dim == 2
                        for i = 1:Nx+2
                            for j = 2:Ny+1
                                A.x(i,j) = U.y(i,j) - U.y(i,j-1);
                            end
                        end
                    end
            end
        elseif U.data == "cell"
            switch lower(A.data)
                case "edge"
                    if dim == 1
                        for i = 1:Nx+1
                            for j = 1:Ny+2
                                A.x(i,j) = U.x(i+1,j) - U.x(i,j);
                            end
                        end
                    elseif dim == 2
                        for i = 1:Nx+2
                            for j = 1:Ny+1
                                A.y(i,j) = U.x(i,j+1) - U.x(i,j);
                            end
                        end
                    end
            end
        elseif U.data == "node"
            switch lower(A.data)
                case "edge"
                    if dim == 1
                        for i = 2:Nx+1
                            for j = 1:Ny+1
                                A.y(i,j) = U.x(i,j) - U.x(i-1,j);
                            end
                        end
                    elseif dim == 2
                        for i = 1:Nx+1
                            for j = 2:Ny+1
                                A.x(i,j) = U.x(i,j) - U.x(i,j-1);
                            end
                        end
                    end
            end
        end
    end
end