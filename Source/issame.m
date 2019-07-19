function [flag] = issame(A,B)
    %ISSAME Check whether two matrices contain the same values.
    % flag = 0 says that the A = B completely.
    % flag = 1 means that A ~= B.
    %
    % [flag] = issame(A,B)
    %
    % Example:
    %
    % A = zeros(3,3);
    % B = 1;
    % flag = issame(A,B)
    %
    % flag = 1 for this case.
    %
    % A = zeros(3,3);
    % B = ones(3,3);
    % flag = issame(A,B)
    % 
    % flag = 1 for this case.
    %
    % A = ones(3,3);
    % B = ones(3,3);
    % flag = issame(A,B);
    % 
    % flag = 0 for this case.
    %
    % Created by Jay Mehta (18 July 2019)
    
    [X1,Y1] = size(A);
    [X2,Y2] = size(B);
    flag = 0;
    
    if X1 ~= X2 || Y1 ~= Y2
        flag = 1;
    else
    end
    
    [flag,~] = iszero(A-B);
end