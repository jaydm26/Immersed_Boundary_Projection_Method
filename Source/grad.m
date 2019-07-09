function A = grad(P,dim)
    
    switch nargin
        case 1
             if P.data ~= "cell"
                error("Check Data Type. Data Type can only be cell")
            else
                Nx = P.size(1);
                Ny = P.size(2);
                A = EdgeData(Nx,Ny);
                for i = 1:Nx+1
                    for j = 1:Ny+2
                        A.x(i,j) = P.x(i+1,j) - P.x(i,j);
                    end
                end
                for i = 1:Nx+2
                    for j = 1:Ny+1
                        A.y(i,j) = P.x(i,j+1) - P.x(i,j);
                    end
                end
             end
             
        case 2
            switch dim
                case 1
                    Nx = P.size(1);
                    Ny = P.size(2);
                    A = EdgeData(Nx,Ny);
                    for i = 1:Nx+1
                        for j = 1:Ny+2
                            A.x(i,j) = P.x(i+1,j) - P.x(i,j);
                        end
                    end
                case 2
                    Nx = P.size(1);
                    Ny = P.size(2);
                    A = EdgeData(Nx,Ny);
                    for i = 1:Nx+2
                        for j = 1:Ny+1
                            A.y(i,j) = P.x(i,j+1) - P.x(i,j);
                        end
                    end
            end
    end
end