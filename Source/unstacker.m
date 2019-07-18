function Y = unstacker(X,type)
    
    global Nx Ny
    
    switch type
        case "cell"
            Y = CellData(Nx,Ny);
            i = 1;
            for j = 1:Ny+2
                Y.x(:,j) = X(i:i+Nx+2-1,1);
                i = i + Nx+2;
            end
        case "edge"
            Y = EdgeData(Nx,Ny);
            i = 1;
            for j = 1:Ny+2
                Y.x(:,j) = X(i:i+Nx+1-1,1);
                i = i + Nx+1;
            end
            for j = 1:Ny+1
                Y.y(:,j) = X(i:i+Nx+2-1,1);
                i = i + Nx+2;
            end
        case "node"
            Y = NodeData(Nx,Ny);
            i = 1;
            for j = 1:Ny+1
                Y.x(:,j) = X(i:i+Nx+1-1,1);
                i = i + Nx+1;
            end
    end
end
    