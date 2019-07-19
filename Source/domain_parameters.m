function domain = domain_parameters(X_c,Y_c,X_e_x,Y_e_x,X_e_y,Y_e_y,...
        X_n,Y_n,Nx,Ny,x_range,y_range)
    %DOMAIN_PARAMETERS Creates a data structure that contains all the 
    % domains required for the flow.
    %
    % domain = domain_parameters(X_c,Y_c,X_e_x,Y_e_x,X_e_y,Y_e_y,...
    %   X_n,Y_n,Nx,Ny,x_range,y_range)
    %
    % Variable lookup:
    %
    % X_c,Y_c = X and Y coordinates for Cell Space data.
    %
    % X_e_x,Y_e_x = X and Y coordinates for X direction Edge Space data.
    %
    % X_e_y,Y_e_y = X and Y coordinates for Y direction Edge Space data.
    %
    % X_n,Y_n = X and Y coordinates for Node Space data.
    %
    % Nx: Number of divisions in the X-direction.
    %
    % Ny: Number of divisions in the Y-direction.
    %
    % x_range: minimum and maximum values of the X-corrdinate in the
    % domain.
    %
    % y_range: minimum and maximum values of the Y-corrdinate in the
    % domain.
    %
    % Created by Jay Mehta (18 July 2019)
    
    domain = struct;
    
    domain.X_c = X_c;
    domain.Y_c = Y_c;
    domain.X_e_x = X_e_x;
    domain.Y_e_x = Y_e_x;
    domain.X_e_y = X_e_y;
    domain.Y_e_y = Y_e_y;
    domain.X_n = X_n;
    domain.Y_n = Y_n;
    domain.Nx = Nx;
    domain.Ny = Ny;
    domain.x_range = x_range;
    domain.y_range = y_range;
end