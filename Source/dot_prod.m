function A = dot_prod(f,g)
    %DOT_PROD Takes the dot product between data field f and g if they are 
    % the same type.
    %
    % A = dot_prod(f,g)
    %
    % Created by Jay Mehta (18 July 2019)
    
    % Check what data type:
    data1 = f.data;
    data2 = g.data;
    if data1 ~= data2
        error("Data Types do not match. Check inputs.")
    else
        switch data1
            case 'cell'
                if f.size ~= g.size
                    error('Vector dimensions do not agree')
                else
                    Nx = f.size(1);
                    Ny = f.size(2);
                    % interiors
                    A = sum(dot(f.x(2:Nx+1,2:Ny+1),g.x(2:Nx+1,2:Ny+1)));
                    % boundaries
                    A = A + 0.5 * (dot(f.x(1,2:Ny+1),g.x(1,2:Ny+1)) + ...
                        dot(f.x(Nx+2,2:Ny+1),g.x(Nx+2,2:Ny+1)) + ...
                        dot(f.x(2:Nx+1,1),g.x(2:Nx+1,1)) + ...
                        dot(f.x(2:Nx+1,Ny+2),g.x(2:Nx+1,Ny+2)));
                    % corners
                    A = A + 0.25 * (f.x(1,1)*g.x(1,1) + ...
                        f.x(Nx+2,1)*g.x(Nx+2,1) + ...
                        f.x(1,Ny+2)*g.x(1,Ny+2) + ...
                        f.x(Nx+2,Ny+2)*g.x(Nx+2,Ny+2));
                    A = A /(Nx*Ny);
                end
            case 'node'
                if f.size ~= g.size
                    error('Vector dimensions do not agree')
                else
                    Nx = f.size(1);
                    Ny = f.size(2);
                    % interiors
                    A = sum(dot(f.x(2:Nx,2:Ny),g.x(2:Nx,2:Ny)));
                    % boundaries
                    A = A + 0.5 * (dot(f.x(1,2:Ny),g.x(1,2:Ny)) + ...
                        dot(f.x(Nx+1,2:Ny),g.x(Nx+1,2:Ny)) + ...
                        dot(f.x(2:Nx,1),g.x(2:Nx,1)) + ...
                        dot(f.x(2:Nx,Ny+1),g.x(2:Nx,Ny+1)));
                    % corners
                    A = A + 0.25 * (f.x(1,1)*g.x(1,1) + ...
                        f.x(Nx+1,1)*g.x(Nx+1,1) + ...
                        f.x(1,Ny+1)*g.x(1,Ny+1) + ...
                        f.x(Nx+1,Ny+1)*g.x(Nx+1,Ny+1));
                    A = A /(Nx*Ny);
                end
                
            case 'edge'
                if f.size ~= g.size
                    error('Vector dimensions do not agree')
                else
                    Nx = f.size(1);
                    Ny = f.size(2);
                    % interior
                    A = sum(dot(f.x(2:Nx,2:Ny+1),g.x(2:Nx,2:Ny+1))) + ...
                        sum(dot(f.y(2:Nx+1,2:Ny),g.y(2:Nx+1,2:Ny)));
                    % boundaries
                    A = A + 0.5 * (dot(f.x(1,2:Ny+1),g.x(1,2:Ny+1)) + ...
                        dot(f.x(Nx+1,2:Ny+1),g.x(Nx+1,2:Ny+1)) + ...
                        dot(f.y(2:Nx+1,1),g.y(2:Nx+1,1)) + ...
                        dot(f.y(2:Nx+1,Ny+1),g.y(2:Nx+1,Ny+1)));
                    % corners
                    A = A + 0.25 * (f.x(1,1)*g.x(1,1) + ...
                        f.x(Nx+1,1)*g.x(Nx+1,1) + ...
                        f.x(1,Ny+2)*g.x(1,Ny+2) + ...
                        f.x(Nx+1,Ny+2)*g.x(Nx+1,Ny+2) +...
                        f.y(1,1)*g.y(1,1) + ...
                        f.y(Nx+2,1)*g.y(Nx+2,1) + ...
                        f.y(1,Ny+1)*g.y(1,Ny+1) + ...
                        f.y(Nx+2,Ny+1)*g.y(Nx+2,Ny+1));
                    A = A /(Nx*Ny);
                end
        end
    end
end