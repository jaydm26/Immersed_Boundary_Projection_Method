function A = integrate(f)
    %
    % Created by Jay Mehta (18 July 2019)
    
    switch f.data
        case 'cell'
            Nx = f.size(1);
            Ny = f.size(2);
            f2 = CellData(Nx,Ny);
            f2.x = ones(Nx+2,Ny+2);
            A = dot_prod(f,f2);
        case 'node'
            Nx = f.size(1);
            Ny = f.size(2);
            f2 = CellData(Nx,Ny);
            f2.x = ones(Nx+1,Ny+1);
            A = dot_prod(f,f2);
        otherwise
            error("Check data type. Only cell and node are accepted")
    end
end