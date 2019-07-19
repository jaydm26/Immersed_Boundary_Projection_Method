function op = div(ip,op,dim)
    %DIV Executes the divergence operation on the data in Edge Space if one
    % argument is provided.
    %
    % op = DIV(ip,op,dim)
    %
    % Variable lookup:
    %
    % ip: Input.
    %
    % op: Output.
    %
    % Nx: Number of divisions in the X-direction.
    %
    % Ny: Number of divisions in the Y-direction.
    %
    % dim: dimension to operate on.
    %
    % Created by Jay Mehta (18 July 2019)
    
    Nx = ip.size(1);
    Ny = ip.size(2);
    
    % This is the divergence operation.
    if nargin == 1
        if ip.data ~= "edge"
            error("Check Data Type. Data Type can only be edge")
        else
            op = CellData(Nx,Ny);
            for i = 2:Nx+1
                for j = 2:Ny+1
                    op.x(i,j) = ip.x(i,j) - ip.x(i-1,j) + ip.y(i,j) - ip.y(i,j-1);
                end
            end
        end
        
    elseif nargin == 2
        error("No direction provided.")
        
    elseif nargin == 3
        % This is the differencing operations.
        if ip.data == "edge" 
            switch lower(op.data)
                case "node"
                    if dim == 1
                        for i = 1:Nx+1
                            for j = 1:Ny+1
                                op.x(i,j) = ip.x(i,j+1) - ip.x(i,j);
                            end
                        end
                    elseif dim == 2
                        for i = 1:Nx+1
                            for j = 1:Ny+1
                                op.x(i,j) = ip.y(i+1,j) - ip.y(i,j);
                            end
                        end
                    end
                case "cell"
                    if dim == 1
                        for i = 2:Nx+1
                            for j = 1:Ny+2
                                op.x(i,j) = ip.x(i,j) - ip.x(i-1,j);
                            end
                        end
                    elseif dim == 2
                        for i = 1:Nx+2
                            for j = 2:Ny+1
                                op.x(i,j) = ip.y(i,j) - ip.y(i,j-1);
                            end
                        end
                    end
            end
        elseif ip.data == "cell"
            switch lower(op.data)
                case "edge"
                    if dim == 1
                        for i = 1:Nx+1
                            for j = 1:Ny+2
                                op.x(i,j) = ip.x(i+1,j) - ip.x(i,j);
                            end
                        end
                    elseif dim == 2
                        for i = 1:Nx+2
                            for j = 1:Ny+1
                                op.y(i,j) = ip.x(i,j+1) - ip.x(i,j);
                            end
                        end
                    end
            end
        elseif ip.data == "node"
            switch lower(op.data)
                case "edge"
                    if dim == 1
                        for i = 2:Nx+1
                            for j = 1:Ny+1
                                op.y(i,j) = ip.x(i,j) - ip.x(i-1,j);
                            end
                        end
                    elseif dim == 2
                        for i = 1:Nx+1
                            for j = 2:Ny+1
                                op.x(i,j) = ip.x(i,j) - ip.x(i,j-1);
                            end
                        end
                    end
            end
        end
    end
end