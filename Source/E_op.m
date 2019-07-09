function out = E_op(xi,eta,s,dim)
    global Nx Ny X_n Y_n x_range y_range dx
    if s.data == "node"
        if nargin == 3
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
        end
    elseif s.data == "cell"
    elseif s.data == "edge"
        if nargin == 4
            if dim == 1
                [X_e_x, Y_e_x] = DomainSetup(x_range,y_range,Nx,Ny,"xe");
                X_e_x = transpose(X_e_x);
                Y_e_x = transpose(Y_e_x);
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
            elseif dim == 2
                [X_e_y, Y_e_y] = DomainSetup(x_range,y_range,Nx,Ny,"ye");
                X_e_y = transpose(X_e_y);
                Y_e_y = transpose(Y_e_y);                
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
            else
                error("Dimension exceeds 2.")
            end
        end
    else
        error("Data Type Not Supported.")
    end
end