function p = CellData(Nx,Ny)
    %
    % Create a data structure that stores data that would be stored in the
    % Cell Space. Types of field variables that could be stored in the cell
    % space are temperature, pressure, kinetic energy, etc.
    % 
    % Example:
    % 
    % T = CellData(64,64) creates a data structure with the field of:
    %
    % x: Stores the values.
    %
    % data: Refers the the data type for the structure.
    %
    % size: Refers to the size used to create the data field.
    
    %% Cell Data Builder
    p = struct;
    p.x = zeros(Nx+2,Ny+2);
    p.data = 'cell';
    p.size = [Nx Ny];
end