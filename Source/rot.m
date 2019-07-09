function A = rot(U)
    if U.data ~= "edge"
        error("Check Data Type. Data Type can only be edge")
    else
        Nx = U.size(1);
        Ny = U.size(2);
        A = NodeData(Nx,Ny);
        for i = 1:Nx+1
            for j = 1:Ny+1
                A.x(i,j) = U.x(i,j) - U.x(i,j+1) + U.y(i+1,j) - U.y(i,j);
            end
        end
    end    
end