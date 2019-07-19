function p = CellData(Nx,Ny)
    %CELLDATA Create a data structure that stores data that would be stored
    % in the Cell Space of size (Nx+2,Ny+2). Types of field variables that 
    % could be stored in the cell space are temperature, pressure, kinetic 
    % energy, etc.
    %
    % p = CELLDATA(Nx,Ny)
    %
    % Variable lookup:
    %
    % Nx: Number of divisions in the X-direction.
    %
    % Ny: Number of divisions in the Y-direction.
    % 
    % Example:
    % 
    % T = CELLDATA(64,64) creates a data structure with the field of:
    %
    % x: Stores the values.
    %
    % data: Refers the the data type for the structure.
    %
    % size: Refers to the size used to create the data field.
    %
    % Created by Jay Mehta (18 July 2019)
    
    %% Cell Data Builder
    p = struct;
    p.x = zeros(Nx+2,Ny+2);
    p.data = 'cell';
    p.size = [Nx Ny];
end