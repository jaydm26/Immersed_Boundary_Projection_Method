function w = NodeData(Nx,Ny)
    %NODEDATA Create a data structure that stores data that would be stored 
    % in the Node Space of size (Nx+1,Ny+1). Types of field variables that
    % could be stored in the cell space are temperature, pressure, kinetic 
    % energy, etc.
    %
    % w = NodeData(Nx,Ny)
    %
    % Variable lookup:
    %
    % Nx: Number of divisions in the X-direction.
    %
    % Ny: Number of divisions in the Y-direction.
    % 
    % Example:
    % 
    % T = NodeData(64,64) creates a data structure with the field of:
    %
    % x: Stores the values.
    %
    % data: Refers the the data type for the structure.
    %
    % size: Refers to the size used to create the data field.
    %
    % Created by Jay Mehta (18 July 2019)
    
    %% Node Data Builder
    w = struct;
    w.x = zeros(Nx+1,Ny+1);
    w.data = 'node';
    w.size = [Nx Ny];
end