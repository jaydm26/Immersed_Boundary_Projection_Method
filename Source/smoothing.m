function u = smoothing(u,f,method,type,TypeVal)
    %SMOOTHING Conventional methods to the solve the Poisson's Equation using the
    % Jacobi Method, Gauss-Siedel, Successive Over-Relaxation, Discrete
    % Cosine Transforms and Discrete Sine Transforms.
    %
    % u = smoothing(u,f,method,type,TypeVal)
    %
    % Varibale lookup:
    %
    % u: solution for the Poisson's Equation.
    %
    % f: input for the Poisson's Equation.
    %
    % method: selects the method for solving-
    %   "jacobi" for Jacobi Method
    %   "gs" for Gauss-Siedel Method
    %   "sor" for Successive Over-Relaxation
    %   "dct" for Discrete Cosine Transform
    %   "dst" for Discrete Sine Transform
    %
    % type: selects the type for solving-
    %   "tol" for tolerance based solver
    %   "iter" for iteration based solver
    %   "..." for Fourier Transform Methods.
    %
    % TypeVal: Value of the type variable. Controls the tolerance or number
    % of iterations.
    %
    % Created by Jay Mehta (18 July 2019)
     
    Nx = u.size(1);
    Ny = u.size(2);
    
    switch lower(type)
        case {"tol","tolerance"}
            switch lower(method)
                case "jacobi"
                    %% Jacobi Method
                    switch lower(u.data)
                        case "edge"
                            U_temp = EdgeData(Nx,Ny);
                            conv = 1;
                            tol = TypeVal;
                            while conv > tol
                                for i = 2:Nx
                                    for j = 2:Ny+1
                                        U_temp.x(i,j) = 0.25 * (u.x(i+1,j) + u.x(i-1,j) + u.x(i,j+1) + u.x(i,j-1) - f.x(j,i));
                                    end
                                end
                                U_temp = apply_bc(U_temp,1);
                                conv = max(max(abs(u.x-U_temp.x)));
                            end
                        case "cell"
                            U_temp = CellData(Nx,Ny);
                            conv = 1;
                            tol = TypeVal;
                            while conv > tol
                                for i = 2:Nx+1
                                    for j = 2:Ny+1
                                        U_temp.x(i,j) = 0.25 * (u.x(i+1,j) + u.x(i-1,j) + u.x(i,j+1) + u.x(i,j-1) - f.x(j,i));
                                    end
                                end
                                conv = max(max(abs(u.x-U_temp.x)));
                            end
                    end

                case {"gauss-seidel","gs"}
                    %% Gauss-Siedel Method
                    switch lower(u.data)
                        case "edge"
                            U_temp = EdgeData(Nx,Ny);
                            conv = 1;
                            tol = TypeVal;
                            while conv > tol
                                U_temp.x = u.x;
                                for i = 2:Nx
                                    for j = 2:Ny+1
                                        u.x(i,j) = 0.25 * (u.x(i+1,j) + u.x(i-1,j) + u.x(i,j+1) + u.x(i,j-1) - f.x(i,j));
                                    end
                                end
                                conv = max(max(abs(u.x-U_temp.x)));
                            end
                        case "cell"
                            U_temp = CellData(Nx,Ny);
                            conv = 1;
                            tol = TypeVal;
                            n_iter = 0;
                            while conv > tol
                                U_temp.x = u.x;
                                for i = 2:Nx+1
                                    for j = 2:Ny+1
                                        u.x(i,j) = 0.25 * (u.x(i+1,j) + u.x(i-1,j) + u.x(i,j+1) + u.x(i,j-1) - f.x(i,j));
                                    end
                                end
                                conv = max(max(abs(u.x-U_temp.x)));
                                n_iter = n_iter + 1;
                                if n_iter > 500
                                    break
                                end
                            end
                        case "node"
                            U_temp = NodeData(Nx,Ny);
                            conv = 1;
                            tol = TypeVal;
                            n_iter = 0;
                            while conv > tol
                                U_temp.x = u.x;
                                for i = 2:Nx
                                    for j = 2:Ny
                                        u.x(i,j) = 0.25 * (u.x(i+1,j) + u.x(i-1,j) + u.x(i,j+1) + u.x(i,j-1) - f.x(i,j));
                                    end
                                end
                                conv = max(max(abs(u.x-U_temp.x)));
                                n_iter = n_iter + 1;
                                if n_iter > 500
                                    break
                                end
                            end
                    end

                case {"sor","successive over-relaxation"}
                    %% Successive Over-Relaxation Method
                    switch lower(u.data)
                        case "edge"
                            U_temp = EdgeData(Nx,Ny);
                            conv = 1;
                            tol = TypeVal;
                            l11 = 1 - 2 * ((pi/(Nx+1))^2 + (pi/(Ny+1))^2);
                            w= 2/(1+sqrt(1-l11^2));
                            while conv > tol
                                U_temp.x = u.x;
                                for i = 2:Nx
                                    for j = 2:Ny+1
                                        u.x(i,j) = 0.25 * w * (u.x(i+1,j) + u.x(i-1,j) + u.x(i,j+1) + u.x(i,j-1) - f.x(i,j)) + (1-w) * U_temp.x(i,j);
                                    end
                                end
                                conv = max(max(abs(u.x-U_temp.x)));
                            end
                        case "cell"
                            U_temp = CellData(Nx,Ny);
                            conv = 1;
                            tol = TypeVal;
                            l11 = 1 - 2 * ((pi/(Nx+1))^2 + (pi/(Ny+1))^2);
                            w= 2/(1+sqrt(1-l11^2));
                            while conv > tol
                                U_temp.x = u.x;
                                for i = 2:Nx+1
                                    for j = 2:Ny+1
                                        u.x(i,j) = 0.25 * w * (u.x(i+1,j) + u.x(i-1,j) + u.x(i,j+1) + u.x(i,j-1) - f.x(i,j)) + (1-w) * U_temp.x(i,j);
                                    end
                                end
                                conv = max(max(abs(u.x-U_temp.x)));
                            end
                        case "node"
                            U_temp = NodeData(Nx,Ny);
                            conv = 1;
                            tol = TypeVal;
                            l11 = 1 - 2 * ((pi/(Nx+1))^2 + (pi/(Ny+1))^2);
                            w= 2/(1+sqrt(1-l11^2));
                            while conv > tol
                                U_temp.x = u.x;
                                for i = 2:Nx
                                    for j = 2:Ny
                                        u.x(i,j) = 0.25 * w * (u.x(i+1,j) + u.x(i-1,j) + u.x(i,j+1) + u.x(i,j-1) - f.x(i,j)) + (1-w) * U_temp.x(i,j);
                                    end
                                end
                                conv = max(max(abs(u.x-U_temp.x)));
                            end
                    end
                            
            end
        case {"iter","iteration"}
            switch lower(method)
                case "jacobi"
                    %% Jacobi Method
                    U_temp = EdgeData(Nx,Ny);
                    n_iter = TypeVal;
                    for n = 1:n_iter
                        for i = 2:Nx
                            for j = 2:Ny+1
                                U_temp.x(i,j) = 0.25 * (u.x(i+1,j) + u.x(i-1,j) + u.x(i,j+1) + u.x(i,j-1) - f.x(i,j));
                            end
                        end                        
                    end
                case {"gauss-seidel","gs"}
                    %% Gauss-Siedel Method
                    U_temp = EdgeData(Nx,Ny);
                    n_iter = TypeVal;
                    for n = 1:n_iter
                        U_temp.x = u.x;
                        for i = 2:Nx
                            for j = 2:Ny+1
                                u.x(i,j) = 0.25 * (u.x(i+1,j) + u.x(i-1,j) + u.x(i,j+1) + u.x(i,j-1) - f.x(i,j));
                            end
                        end
                    end
                case {"sor","successive over-relaxation"}
                    %% Successive Over-Relaxation Method
                    U_temp = EdgeData(Nx,Ny);
                    l11 = 1 - 2 * ((pi/(Nx+1))^2 + (pi/(Ny+1))^2);
                    w= 2/(1+sqrt(1-l11^2));
                    n_iter = TypeVal;
                    for n = 1:n_iter
                        U_temp.x = u.x;
                        for i = 2:Nx
                            for j = 2:Ny+1
                                u.x(i,j) = 0.25 * w *  (u.x(i+1,j) + u.x(i-1,j) + u.x(i,j+1) + u.x(i,j-1) - f.x(i,j)) + (1-w) * U_temp.x(i,j);
                            end
                        end
                    end
            end
        case {"multigrid","mg"}
            %% Add Multigrid
        case "dct"
            %% Discrete Cosine Transforms
            u = f.x(2:Nx+1,2:Ny+1);
            u = dct(u,[],1);
            u = dct(u,[],2);

            lam = cos(pi*(0:Nx-1)'/Nx) + cos(pi*(0:Ny-1)/Ny) - 2 * ones(Ny,Ny);
            lam(1,1) = 1;
            u(1,1) = 0;
            u = 0.5 * u ./ lam;

            u = idct(u,[],1);
            u = idct(u,[],2);
        case "dst"
            %% Discrete Sine Transforms
            switch lower(u.data)
                case "cell"
                    f = f.x(2:Nx+1,2:Ny+1);
                    f = dst(f);
                    f = dst(f');
                    f = f';

                    lam = -4 * sin(pi*(1:Nx)'/(2*(Nx+1))); % Reduced domain leads to a drop
                    for i = 1:Nx
                        for j = 1:Ny
                            f(i,j) = f(i,j)/(lam(i)+lam(j));
                        end
                    end
                    f = idst(f);
                    f = idst(f');
                    f = f';
                    u.x(2:Nx+1,2:Ny+1) = f;
                case "node"
                    f = f.x(2:Nx,2:Ny);
                    f = dst(f);
                    f = dst(f');
                    f = f';
                    for i = 1:Nx-1
                        cosn = cos(pi*i/Nx);
                        for j = 1:Ny-1
                            cosm = cos(pi*j/Ny);
                            f(i,j) = 0.5 * f(i,j)/(cosm+cosn-2);
                        end
                    end
                    f = idst(f);
                    f = idst(f');
                    f = f';
                    u.x(2:Nx,2:Ny) = f;
            end
    end
end