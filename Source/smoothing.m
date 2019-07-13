function U = smoothing(U,f,method,type,TypeVal)
     
    global Nx Ny
    
    switch lower(type)
        case {"tol","tolerance"}
            switch lower(method)
                case "jacobi"
                    %% Jacobi Method
                    switch lower(U.data)
                        case "edge"
                            U_temp = EdgeData(Nx,Ny);
                            conv = 1;
                            tol = TypeVal;
                            while conv > tol
                                for i = 2:Nx
                                    for j = 2:Ny+1
                                        U_temp.x(i,j) = 0.25 * (U.x(i+1,j) + U.x(i-1,j) + U.x(i,j+1) + U.x(i,j-1) - f.x(j,i));
                                    end
                                end
                                U_temp = apply_bc(U_temp,1);
                                conv = max(max(abs(U.x-U_temp.x)));
                            end
%                             U.x = -U_temp.x;
                        case "cell"
                            U_temp = CellData(Nx,Ny);
                            conv = 1;
                            tol = TypeVal;
                            while conv > tol
                                for i = 2:Nx+1
                                    for j = 2:Ny+1
                                        U_temp.x(i,j) = 0.25 * (U.x(i+1,j) + U.x(i-1,j) + U.x(i,j+1) + U.x(i,j-1) - f.x(j,i));
                                    end
                                end
                                conv = max(max(abs(U.x-U_temp.x)));
                            end
%                             U.x = -U.x;
                    end

                case {"gauss-seidel","gs"}
                    %% Gauss-Siedel Method
                    switch lower(U.data)
                        case "edge"
                            U_temp = EdgeData(Nx,Ny);
                            conv = 1;
                            tol = TypeVal;
                            while conv > tol
                                U_temp.x = U.x;
                                for i = 2:Nx
                                    for j = 2:Ny+1
                                        U.x(i,j) = 0.25 * (U.x(i+1,j) + U.x(i-1,j) + U.x(i,j+1) + U.x(i,j-1) - f.x(i,j));
                                    end
                                end
                                conv = max(max(abs(U.x-U_temp.x)));
                            end
%                             U.x = -U.x;
                        case "cell"
                            U_temp = CellData(Nx,Ny);
                            conv = 1;
                            tol = TypeVal;
                            n_iter = 0;
                            while conv > tol
                                U_temp.x = U.x;
                                for i = 2:Nx+1
                                    for j = 2:Ny+1
                                        U.x(i,j) = 0.25 * (U.x(i+1,j) + U.x(i-1,j) + U.x(i,j+1) + U.x(i,j-1) - f.x(i,j));
                                    end
                                end
                                conv = max(max(abs(U.x-U_temp.x)));
                                n_iter = n_iter + 1;
                                if n_iter > 500
                                    break
                                end
                            end
%                             U.x = -U.x;
                        case "node"
                            U_temp = NodeData(Nx,Ny);
                            conv = 1;
                            tol = TypeVal;
                            n_iter = 0;
                            while conv > tol
                                U_temp.x = U.x;
                                for i = 2:Nx
                                    for j = 2:Ny
                                        U.x(i,j) = 0.25 * (U.x(i+1,j) + U.x(i-1,j) + U.x(i,j+1) + U.x(i,j-1) - f.x(i,j));
                                    end
                                end
                                conv = max(max(abs(U.x-U_temp.x)));
                                n_iter = n_iter + 1;
                                if n_iter > 500
                                    break
                                end
                            end
%                             U.x = -U.x;
                    end

                case {"sor","successive over-relaxation"}
                    %% Successive Over-Relaxation Method
                    switch lower(U.data)
                        case "edge"
                            U_temp = EdgeData(Nx,Ny);
                            conv = 1;
                            tol = TypeVal;
                            l11 = 1 - 2 * ((pi/(Nx+1))^2 + (pi/(Ny+1))^2);
                            w= 2/(1+sqrt(1-l11^2));
                            while conv > tol
                                U_temp.x = U.x;
                                for i = 2:Nx
                                    for j = 2:Ny+1
                                        U.x(i,j) = 0.25 * w * (U.x(i+1,j) + U.x(i-1,j) + U.x(i,j+1) + U.x(i,j-1) - f.x(i,j)) + (1-w) * U_temp.x(i,j);
                                    end
                                end
                                conv = max(max(abs(U.x-U_temp.x)));
                            end
%                             U.x = -U.x;
                        case "cell"
                            U_temp = CellData(Nx,Ny);
                            conv = 1;
                            tol = TypeVal;
                            l11 = 1 - 2 * ((pi/(Nx+1))^2 + (pi/(Ny+1))^2);
                            w= 2/(1+sqrt(1-l11^2));
%                             n_iter = 0;
                            while conv > tol
                                U_temp.x = U.x;
                                for i = 2:Nx+1
                                    for j = 2:Ny+1
                                        U.x(i,j) = 0.25 * w * (U.x(i+1,j) + U.x(i-1,j) + U.x(i,j+1) + U.x(i,j-1) - f.x(i,j)) + (1-w) * U_temp.x(i,j);
                                    end
                                end
                                conv = max(max(abs(U.x-U_temp.x)));
%                                 n_iter = n_iter + 1;
%                                 if n_iter > 500
%                                     break
%                                 end
                            end
%                             U.x = -U.x;
                        case "node"
                            U_temp = NodeData(Nx,Ny);
                            conv = 1;
                            tol = TypeVal;
                            l11 = 1 - 2 * ((pi/(Nx+1))^2 + (pi/(Ny+1))^2);
                            w= 2/(1+sqrt(1-l11^2));
%                             n_iter = 0;
                            while conv > tol
                                U_temp.x = U.x;
                                for i = 2:Nx
                                    for j = 2:Ny
                                        U.x(i,j) = 0.25 * w * (U.x(i+1,j) + U.x(i-1,j) + U.x(i,j+1) + U.x(i,j-1) - f.x(i,j)) + (1-w) * U_temp.x(i,j);
                                    end
                                end
                                conv = max(max(abs(U.x-U_temp.x)));
                            end
%                             U.x = -U.x;
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
                                U_temp.x(i,j) = 0.25 * (U.x(i+1,j) + U.x(i-1,j) + U.x(i,j+1) + U.x(i,j-1) - f.x(i,j));
                            end
                        end                        
                    end
%                     U.x = -U_temp.x;
                case {"gauss-seidel","gs"}
                    %% Gauss-Siedel Method
                    U_temp = EdgeData(Nx,Ny);
                    n_iter = TypeVal;
                    for n = 1:n_iter
                        U_temp.x = U.x;
                        for i = 2:Nx
                            for j = 2:Ny+1
                                U.x(i,j) = 0.25 * (U.x(i+1,j) + U.x(i-1,j) + U.x(i,j+1) + U.x(i,j-1) - f.x(i,j));
                            end
                        end
                    end
%                     U.x = -U.x;
                case {"sor","successive over-relaxation"}
                    %% Successive Over-Relaxation Method
                    U_temp = EdgeData(Nx,Ny);
                    l11 = 1 - 2 * ((pi/(Nx+1))^2 + (pi/(Ny+1))^2);
                    w= 2/(1+sqrt(1-l11^2));
                    n_iter = TypeVal;
                    for n = 1:n_iter
                        U_temp.x = U.x;
                        for i = 2:Nx
                            for j = 2:Ny+1
                                U.x(i,j) = 0.25 * w *  (U.x(i+1,j) + U.x(i-1,j) + U.x(i,j+1) + U.x(i,j-1) - f.x(i,j)) + (1-w) * U_temp.x(i,j);
                            end
                        end
                    end
%                     U.x = -U.x;
            end
        case {"multigrid","mg"}
            %% Add Multigrid when available
            
        case "dct"
            %% Discrete Cosine Transforms
            U = f.x(2:Nx+1,2:Ny+1);
            U = dct(U,[],1);
            U = dct(U,[],2);

            lam = cos(pi*(0:Nx-1)'/Nx) + cos(pi*(0:Ny-1)/Ny) - 2 * ones(Ny,Ny);
            lam(1,1) = 1;
            U(1,1) = 0;
            U = 0.5 * U ./ lam;

            U = idct(U,[],1);
            U = idct(U,[],2);
            
        case "dst"
            %% Discrete Sine Transforms
            switch lower(U.data)
                case "cell"
                    f = f.x(2:Nx+1,2:Ny+1);
                    f = fast_dst(f,2);
                    f = fast_idst(f,1);

                    lam = -4 * sin(pi*(1:Nx)'/(2*(Nx+1))); % Reduced domain leads to a drop
                    for i = 1:Nx
                        for j = 1:Ny
                            f(i,j) = f(i,j)/(lam(i)+lam(j));
                        end
                    end

                    f = fast_idst(f,2);
                    f = fast_dst(f,1);
                    U.x(2:Nx+1,2:Ny+1) = f;
                case "node"
%                     f = f.x(2:Nx,2:Ny);
%                     f = fast_dst(f,2);
%                     f = fast_idst(f,1);
% 
%                     lam = -4 * sin(pi*(1:Nx-1)'/(2*(Nx))).^2; % Reduced domain leads to a drop
%                     for i = 1:Nx-1
%                         for j = 1:Ny-1
%                             f(i,j) = f(i,j)/(lam(i)+lam(j));
%                         end
%                     end
% 
%                     f = fast_idst(f,2);
%                     f = fast_dst(f,1);
%                     U.x(2:Nx,2:Ny) = f;
%% EDIT
                    f = f.x(2:Nx,2:Ny);
                    f = dst(f);
                    f = dst(f');
                    f = f';

%                     lam = -4 * sin(pi*(1:Nx-1)'/(2*(Nx))).^2; % Reduced domain leads to a drop
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
                    U.x(2:Nx,2:Ny) = f;
%                     U.x = -U.x;

            end

    end
end