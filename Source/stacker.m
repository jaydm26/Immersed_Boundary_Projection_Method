function op = stacker(ip)
    %STACKER Stacker stacks the ip data structure in a column for POD analysis.
    %
    % Variable lookup:
    %
    % ip: input field variable.
    %
    % op: output stacked column.
    %
    % Created by Jay Mehta (18 July 2019)
    
    Nx = ip.size(1);
    Ny = ip.size(2);
    
    switch ip.data
        case "edge"
            op = zeros(((Nx+1)*(Ny+2) + (Nx+2)*(Ny+1)),1);
            i = 1;
            for j = 1:Ny+2
                op(i:i+Nx+1-1,1) = ip.x(:,j);
                i = i + Nx+1;
            end
            for j = 1:Ny+1
                op(i:i+Nx+2-1,1) = ip.y(:,j);
                i = i + Nx+2;
            end
        case "cell"
            op = zeros(((Nx+2)*(Ny+2)),1);
            i = 1;
            for j = 1:Ny+2
                op(i:i+Nx+2-1,1) = ip.x(:,j);
                i = i + Nx+2;
            end
        case "node"
            op = zeros(((Nx+1)*(Ny+1)),1);
            i = 1;
            for j = 1:Ny+1
                op(i:i+Nx+1-1,1) = ip.x(:,j);
                i = i + Nx+1;
            end
    end
end