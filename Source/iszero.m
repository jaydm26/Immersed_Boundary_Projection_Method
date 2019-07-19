function [flag,n] = iszero(A)
    %ISZERO Check if the input matrix A is all zeros or not. Useful for checking
    % the nature of very large but sparse matrices.
    % flag = 1 corresponds to at least one non-zero value inside the
    % matrix.
    % n gives the number of non-zero values.
    %
    % [flag,n] = iszero(A)
    %
    % Example:
    % A = zeros(10,10);
    % [flag,n] = iszero(A);
    %
    % This gives flag = 0 and n = 0
    %
    % A = ones(10,10);
    % [flag,n] = iszero(A);
    % 
    % This gives flag = 1 and n = 100
    %
    % Created by Jay Mehta (18 July 2019)
    
    [x,y] = size(A);
    n = 0;
    flag = 0;
    for i = 1:x
        for j = 1:y
            if A(i,j) ~= 0
                flag = 1;
                n = n + 1;
            else
            end
        end
    end
end