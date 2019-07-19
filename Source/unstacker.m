function op = unstacker(domain,ip,type)
    %UNSTACKER Unstacks the ip data structure in a column for POD analysis.
    %
    % op = UNSTACKER(domain,ip,type)
    % 
    % Variable lookup:
    %
    % domain: domain parameters.
    %
    % ip: input stacked column.
    %
    % op: output field variable.
    %
    % type: location of field variable- "cell" for Cell Space
    %                                   "edge" for Edge Space
    %                                   "node" for Node Space
    %
    % Created by Jay Mehta (18 July 2019)
    
    Nx = domain.Nx;
    Ny = domain.Ny;
    
    switch type
        case "cell"
            op = CellData(Nx,Ny);
            i = 1;
            for j = 1:Ny+2
                op.x(:,j) = ip(i:i+Nx+2-1,1);
                i = i + Nx+2;
            end
        case "edge"
            op = EdgeData(Nx,Ny);
            i = 1;
            for j = 1:Ny+2
                op.x(:,j) = ip(i:i+Nx+1-1,1);
                i = i + Nx+1;
            end
            for j = 1:Ny+1
                op.y(:,j) = ip(i:i+Nx+2-1,1);
                i = i + Nx+2;
            end
        case "node"
            op = NodeData(Nx,Ny);
            i = 1;
            for j = 1:Ny+1
                op.x(:,j) = ip(i:i+Nx+1-1,1);
                i = i + Nx+1;
            end
    end
end
    