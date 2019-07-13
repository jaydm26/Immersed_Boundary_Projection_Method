function U = interpol(U,V,dim)
    
    % Interpolates data in V to U.
    if U.size ~= V.size
        error("Data sizes do not match")
    else
        if U.data == "edge"
            switch lower(V.data)
                case "cell"
                    Nx = U.size(1);
                    Ny = U.size(2);
                    if dim == 1
                        for i = 1:Nx+1
                            for j = 1:Ny+2
                                U.x(i,j) = 0.5*(V.x(i+1,j) + V.x(i,j));
                            end
                        end
                    elseif dim == 2
                        for i = 1:Nx+2
                            for j = 1:Ny+1
                                U.y(i,j) = 0.5*(V.x(i,j+1) + V.x(i,j));
                            end
                        end
                    end
                case "node"
                    Nx = U.size(1);
                    Ny = U.size(2);
                    if dim == 1
                        for i = 1:Nx+1
                            for j = 2:Ny+1
                                U.x(i,j) = 0.5*(V.x(i,j) + V.x(i,j-1));
                            end
                        end
                    elseif dim == 2
                        for i = 2:Nx+1
                            for j = 1:Ny+1
                                U.y(i,j) = 0.5*(V.x(i,j) + V.x(i-1,j));
                            end
                        end
                    end
                case "edge"
                    Nx = U.size(1);
                    Ny = U.size(2);
                    for i = 1:Nx+1
                        for j = 2:Ny+1
                            U.x(i,j) = 0.25 * (V.y(i,j-1) + V.y(i+1,j-1) + V.y(i,j) + V.y(i+1,j));
                        end
                    end
                    for i = 2:Nx+1
                        for j = 1:Ny+1
                            U.y(i,j) = 0.25 * (V.x(i-1,j) + V.x(i-1,j+1) + V.x(i,j) + V.x(i,j+1));
                        end
                    end
            end
        elseif U.data == "cell"
            switch lower(V.data)
                case "cell"
                case "edge"
                    Nx = U.size(1);
                    Ny = U.size(2);
                    if dim == 1
                        for i = 2:Nx+1
                            for j = 1:Ny+2
                                U.x(i,j) = 0.5 * (V.x(i,j) + V.x(i-1,j));
                            end
                        end
                    elseif dim == 2
                        for i = 1:Nx+2
                            for j = 2:Ny+1
                                U.x(i,j) = 0.5 * (V.y(i,j) + V.y(i,j-1));
                            end
                        end
                    end
                case "node"
            end
        elseif U.data == "node"
            switch lower(V.data)
                case "cell"
                case "edge"
                    if dim == 1
                        Nx = U.size(1);
                        Ny = U.size(2);
                        for i = 1:Nx+1
                            for j = 1:Ny+1
                                U.x(i,j) = 0.5 * (V.x(i,j+1) + V.x(i,j));
                            end
                        end
                    elseif dim == 2
                        Nx = U.size(1);
                        Ny = U.size(2);
                        for i = 1:Nx+1
                            for j = 1:Ny+1
                                U.x(i,j) = 0.5 * (V.y(i+1,j) + V.y(i,j));
                            end
                        end
                    end
                case "node"
            end
        end
    end
end 