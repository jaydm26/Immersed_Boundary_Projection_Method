function out = E_op(xi,eta,s,dim)
    global Nx Ny X_n Y_n X_e_x Y_e_x X_e_y Y_e_y X_c Y_c dx dy
    if s.data == "node"
        if nargin == 3
            X_n = X_n';
            Y_n = Y_n';
            for i = 1:Nx+1
                for j = 1:Ny+1
                    if s.x(i,j) == 0 || (X_n(i,j)-xi)/dx > 1.5 || (Y_n(i,j)-eta)/dx > 1.5
                        out(i,j) = 0;
                    else
                        out(i,j) = s.x(i,j) * ddf_roma_2D(X_n(i,j)-xi,Y_n(i,j)-eta);
                    end
                end
            end
            out = sum(sum(out));
            X_n = X_n';
            Y_n = Y_n';
        end
    elseif s.data == "cell"
        if nargin == 3
            X_c = X_c';
            Y_c = Y_c';
            for i = 1:Nx+2
                for j = 1:Ny+2
                    if s.x(i,j) == 0 || (X_c(i,j)-xi)/dx > 1.5 || (Y_c(i,j)-eta)/dx > 1.5
                        out(i,j) = 0;
                    else
                        out(i,j) = s.x(i,j) * ddf_roma_2D(X_c(i,j)-xi,Y_c(i,j)-eta);
                    end
                end
            end
            out = sum(sum(out));
            X_c = X_c';
            Y_c = Y_c';
        end
    elseif s.data == "edge"
        if nargin == 4
            if dim == 1
                X_e_x = X_e_x';
                Y_e_x = Y_e_x';
                for i = 1:Nx+1
                    for j = 1:Ny+2
                        if s.x(i,j) == 0 || (X_e_x(i,j)-xi)/dx > 1.5 || (Y_e_x(i,j)-eta)/dx > 1.5
                            out(i,j) = 0;
                        else
                            out(i,j) = s.x(i,j) * ddf_roma_2D(X_e_x(i,j)-xi,Y_e_x(i,j)-eta);
                        end
                    end
                end
                out = sum(sum(out));
                X_e_x = X_e_x';
                Y_e_x = Y_e_x';
            elseif dim == 2
                X_e_y = X_e_y';
                Y_e_y = Y_e_y';                
                for i = 1:Nx+2
                    for j = 1:Ny+1
                        if s.y(i,j) == 0 || (X_e_y(i,j)-xi)/dx > 1.5 || (Y_e_y(i,j)-eta)/dx > 1.5
                            out(i,j) = 0;
                        else
                            out(i,j) = s.y(i,j) * ddf_roma_2D(X_e_y(i,j)-xi,Y_e_y(i,j)-eta);
                        end
                    end
                end
                out = sum(sum(out));
                X_e_y = X_e_y';
                Y_e_y = Y_e_y';
            else
                error("Dimension exceeds 2.")
            end
        end
    else
        error("Data Type Not Supported.")
    end
end