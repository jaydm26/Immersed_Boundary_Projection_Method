function u = EdgeData(Nx,Ny)
    %EDGEDATA Create a data structure that stores data that would be stored in the
    % Edge Space of size (Nx+1,Ny+2) for the X-direction and (Nx+2,Ny+1) 
    % for the Y-direction. Types of field variables that could be stored in
    % the cell space is velocity.
    %
    % u = EdgeData(Nx,Ny)
    %
    % Variable lookup:
    %
    % Nx: Number of divisions in the X-direction.
    %
    % Ny: Number of divisions in the Y-direction.
    % 
    % Example:
    % 
    % velocity = EdgeData(64,64) creates a data structure with the field of:
    %
    % x: Stores the values for the X-direction.
    %
    % y: Stores the values for the Y-direction.
    %
    % data: Refers the the data type for the structure.
    %
    % size: Refers to the size used to create the data field.
    %
    % Created by Jay Mehta (18 July 2019)
    
    %% Edge Data Builder
    u = struct;
    u.x = zeros(Nx+1,Ny+2);
    u.y = zeros(Nx+2,Ny+1);
    u.data = 'edge';
    u.size = [Nx Ny];
end