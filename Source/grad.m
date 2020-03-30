function op = grad(ip)
    %GRAD Takes the gradient of the field contained in ip. Here, field must
    % be cell only. 
    %
    % op = grad(ip)
    %
    % Variable lookup:
    %
    % ip: input field.
    %
    % op: gradient of the input field.
    %
    % Created by Jay Mehta (18 July 2019)
    
    Nx = ip.size(1);
    Ny = ip.size(2);
    
    switch ip.data
        case "cell"
            op = EdgeData(Nx,Ny);
            for i = 1:Nx+1
                for j = 1:Ny+2
                    op.x(i,j) = ip.x(i+1,j) - ip.x(i,j);
                end
            end
            for i = 1:Nx+2
                for j = 1:Ny+1
                    op.y(i,j) = ip.x(i,j+1) - ip.x(i,j);
                end
            end
        otherwise
            error("Check Data Type. Data Type can only be cell")
    end

end