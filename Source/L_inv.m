function g_hat = L_inv(domain,DataType)
    %L_INV Create the L_inv operator (including FFT). This is coupled with 
    % the L_inv operation. Used for Node Data and Cell Data to solve the
    % Poisson's equation: d^2 u / dx^2 + d^2 u / dy^2 = f. Uses the Lattice
    % Green's Function.
    %
    % g_hat = L_inv(domain,x_range,y_range,DataType)
    %
    % Variable lookup:
    %
    % g_hat: fft2 of the Lattice Green's Function.
    %
    % domain: data structure containing all domains.
    %
    % x_range: minimum and maximum values of the X-corrdinate in the
    % domain.
    %
    % y_range: minimum and maximum values of the Y-corrdinate in the
    % domain.
    %
    % DataType: Type of L_inv operator- "cell" for Cell Space
    %                                   "node" for Node Space
    %
    % Created by Jay Mehta (18 July 2019)
    
    Nx = domain.Nx;
    Ny = domain.Ny;
    x_range = domain.x_range;
    y_range = domain.y_range;
    params = flow_parameters_init;
    domain2 = domain_parameters_init;
    domain2.x_range = x_range;
    domain2.y_range = y_range;
    
    switch DataType
        case "node"
            Nx = 2*Nx; % Had to update due to global declaration
            Ny = 2*Ny; % Reset after this setup
            domain2.Nx = Nx;
            domain2.Ny = Ny;
            params.dx = (x_range(2) - x_range(1))/Nx;

            [X_n2, Y_n2] = DomainSetup(params,domain2,"node");
            domain2.X_n = X_n2;
            domain2.Y_n = Y_n2;

            i_center = Nx/2 + 1;
            j_center = Ny/2 + 1;

            g_0 = NodeData(Nx,Ny);
            X_n2 = X_n2';
            Y_n2 = Y_n2';
            % Define boundaries of g_0 using the Green's Function for Laplace in 2D
            G = @(x,y) 1/(4*pi) * log(x^2+y^2);

            for i = [1,Nx+1]
                for j = 1:Ny+1
                    g_0.x(i,j) = G(X_n2(i,j),Y_n2(i,j));
                end
            end

            for i = 1:Nx+1
                for j = [1,Ny+1]
                    g_0.x(i,j) = G(X_n2(i,j),Y_n2(i,j));
                end
            end

            f = NodeData(Nx,Ny);
            f.x(i_center,j_center) = 1;

            rhs = laplacian_2(g_0);
            rhs.x = -rhs.x + f.x;

            g_0 = smoothing(g_0,rhs,"...","dst");
            g_hat = fft2(g_0.x); % In preparation of the FFT operator
            
        case "cell"
            Nx = 2*Nx; % Had to update due to global declaration
            Ny = 2*Ny; % Reset after this setup
            domain2.Nx = Nx;
            domain2.Ny = Ny;
            params.dx = (x_range(2) - x_range(1))/Nx;

            [X_c2, Y_c2] = DomainSetup(params,domain2,"cell");

            i_center = Nx/2 + 1;
            j_center = Ny/2 + 1;

            g_0 = CellData(Nx,Ny);
            X_c2 = X_c2';
            Y_c2 = Y_c2';
            % Define boundaries of g_0 using the Green's Function for Laplace in 2D
            G = @(x,y) 1/(4*pi) * log(x^2+y^2);

            for i = [1,Nx+2]
                for j = 1:Ny+2
                    g_0.x(i,j) = G(X_c2(i,j),Y_c2(i,j));
                end
            end

            for i = 1:Nx+2
                for j = [1,Ny+2]
                    g_0.x(i,j) = G(X_c2(i,j),Y_c2(i,j));
                end
            end

            f = CellData(Nx,Ny);
            f.x(i_center,j_center) = 1;

            rhs = laplacian_2(g_0);
            rhs.x = -rhs.x + f.x;

            g_0 = smoothing(g_0,rhs,"...","dst");
            g_hat = fft2(g_0.x); % In preparation of the FFT operator
    end
end
