function g_hat = L_inv(DataType)
    % Create the L_inv operator (including FFT). Coupled with the L_inv
    % operation.
    
    global Nx Ny x_range y_range
    
    switch DataType
        case "node"
            Nx = 2*Nx; % Had to update due to global declaration
            Ny = 2*Ny; % Reset after this setup

            [X_n2, Y_n2] = DomainSetup(x_range,y_range,Nx,Ny,"node");

            i_center = Nx/2 + 1;
            j_center = Ny/2 + 1;

            g_0 = NodeData(Nx,Ny);
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

            Nx = Nx/2; % Resetting everything
            Ny = Ny/2;
            
        case "cell"
            Nx = 2*Nx; % Had to update due to global declaration
            Ny = 2*Ny; % Reset after this setup

            [X_c2, Y_c2] = DomainSetup(x_range,y_range,Nx,Ny,"cell");

            i_center = Nx/2 + 1;
            j_center = Ny/2 + 1;

            g_0 = CellData(Nx,Ny);
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

            Nx = Nx/2; % Resetting everything
            Ny = Ny/2;
    end
end
