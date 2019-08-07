function op = upwinding(ip,op,dim)
    %UPWINDING Evaluates the value of the terms by using the upwinding
    %scheme.
    %
    % op = upwinding(ip,op,dim)
    %
    % Variable lookup:
    %
    % ip: input field.
    %
    % op: output field that contains the upwinded values.
    %
    % dim: dimension to operate on.
    %
    % Created by Jay Mehta (18 July 2019)
    
    Nx = ip.size(1);
    Ny = ip.size(2);    
    
    switch ip.data
        case "edge"
            switch op.data
                case "cell"
                    if dim == 1
                        for i = 2:Nx+1
                            for j = 1:Ny+2
                                op.x(i,j) = 0.5 * (ip.x(i,j) - ip.x(i-1,j));
                            end
                        end
                    elseif dim == 2
                        for i = 1:Nx+2
                            for j = 2:Ny+1
                                op.x(i,j) = 0.5 * (ip.y(i,j) - ip.y(i,j-1));
                            end
                        end
                    end
                case "node"
                    if dim == 1
                        for i = 1:Nx+1
                            for j = 1:Ny+1
                                op.x(i,j) = 0.5 * (ip.x(i,j+1) - ip.x(i,j));
                            end
                        end
                    elseif dim == 2
                        for i = 1:Nx+1
                            for j = 1:Ny+1
                                op.x(i,j) = 0.5 * (ip.y(i+1,j) - ip.y(i,j));
                            end
                        end
                    end
            end
        case "cell"
            switch op.data
                case "edge"
                    if dim == 1
                        for i = 1:Nx+1
                            for j = 1:Ny+2
                                op.x(i,j) = 0.5*(ip.x(i+1,j) - ip.x(i,j));
                            end
                        end
                    elseif dim == 2
                        for i = 1:Nx+2
                            for j = 1:Ny+1
                                op.y(i,j) = 0.5*(ip.x(i,j+1) - ip.x(i,j));
                            end
                        end
                    end
            end
    end
end