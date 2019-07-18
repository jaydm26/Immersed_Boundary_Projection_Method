function Y = stacker(X)
    
    % X is a data structure. Stacker stacks the variables in columns for
    % POD analysis.
    global Nx Ny
    
    switch X.data
        case "edge"
            Y = zeros(((Nx+1)*(Ny+2) + (Nx+2)*(Ny+1)),1);
            i = 1;
            for j = 1:Ny+2
                Y(i:i+Nx+1-1,1) = X.x(:,j);
                i = i + Nx+1;
            end
            for j = 1:Ny+1
                Y(i:i+Nx+2-1,1) = X.y(:,j);
                i = i + Nx+2;
            end
        case "cell"
            Y = zeros(((Nx+2)*(Ny+2)),1);
            i = 1;
            for j = 1:Ny+2
                Y(i:i+Nx+2-1,1) = X.x(:,j);
                i = i + Nx+2;
            end
        case "node"
            Y = zeros(((Nx+1)*(Ny+1)),1);
            i = 1;
            for j = 1:Ny+1
                Y(i:i+Nx+1-1,1) = X.x(:,j);
                i = i + Nx+1;
            end
    end
end