function op = interpol(ip,op,dim)
    %INTERPOL Interpolates data in ip to op. Refer to reference for further
    % explanation.
    %
    % op = interpol(ip,op,dim)
    %
    % Variable lookup:
    %
    % op: interpolated output.
    %
    % ip: input that has to be interpolated
    %
    % dim: dimension to interpolate over.
    %
    % Created by Jay Mehta (18 July 2019)
    
    if op.size ~= ip.size
        error("Data sizes do not match")
    else
        if op.data == "edge"
            switch lower(ip.data)
                case "cell"
                    Nx = op.size(1);
                    Ny = op.size(2);
                    if dim == 1
                        for i = 1:Nx+1
                            for j = 1:Ny+2
                                op.x(i,j) = 0.5*(ip.x(i+1,j) + ip.x(i,j));
                            end
                        end
                    elseif dim == 2
                        for i = 1:Nx+2
                            for j = 1:Ny+1
                                op.y(i,j) = 0.5*(ip.x(i,j+1) + ip.x(i,j));
                            end
                        end
                    end
                case "node"
                    Nx = op.size(1);
                    Ny = op.size(2);
                    if dim == 1
                        for i = 1:Nx+1
                            for j = 2:Ny+1
                                op.x(i,j) = 0.5*(ip.x(i,j) + ip.x(i,j-1));
                            end
                        end
                    elseif dim == 2
                        for i = 2:Nx+1
                            for j = 1:Ny+1
                                op.y(i,j) = 0.5*(ip.x(i,j) + ip.x(i-1,j));
                            end
                        end
                    end
                case "edge"
                    Nx = op.size(1);
                    Ny = op.size(2);
                    for i = 1:Nx+1
                        for j = 2:Ny+1
                            op.x(i,j) = 0.25 * (ip.y(i,j-1) + ip.y(i+1,j-1) + ip.y(i,j) + ip.y(i+1,j));
                        end
                    end
                    for i = 2:Nx+1
                        for j = 1:Ny+1
                            op.y(i,j) = 0.25 * (ip.x(i-1,j) + ip.x(i-1,j+1) + ip.x(i,j) + ip.x(i,j+1));
                        end
                    end
            end
        elseif op.data == "cell"
            switch lower(ip.data)
                case "cell"
                case "edge"
                    Nx = op.size(1);
                    Ny = op.size(2);
                    if dim == 1
                        for i = 2:Nx+1
                            for j = 1:Ny+2
                                op.x(i,j) = 0.5 * (ip.x(i,j) + ip.x(i-1,j));
                            end
                        end
                    elseif dim == 2
                        for i = 1:Nx+2
                            for j = 2:Ny+1
                                op.x(i,j) = 0.5 * (ip.y(i,j) + ip.y(i,j-1));
                            end
                        end
                    end
                case "node"
            end
        elseif op.data == "node"
            switch lower(ip.data)
                case "cell"
                case "edge"
                    if dim == 1
                        Nx = op.size(1);
                        Ny = op.size(2);
                        for i = 1:Nx+1
                            for j = 1:Ny+1
                                op.x(i,j) = 0.5 * (ip.x(i,j+1) + ip.x(i,j));
                            end
                        end
                    elseif dim == 2
                        Nx = op.size(1);
                        Ny = op.size(2);
                        for i = 1:Nx+1
                            for j = 1:Ny+1
                                op.x(i,j) = 0.5 * (ip.y(i+1,j) + ip.y(i,j));
                            end
                        end
                    end
                case "node"
            end
        end
    end
end 