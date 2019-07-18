function A = MatrixA_Generator(xi,eta,Type)
    %MATRIXA_GENERATOR Creates the matrix A which is used to replace 
    % ECL^-1(EC)^T or EE^T matrix used to speed up the conjugate gradient 
    % solution.
    %
    % A = MatrixA_Generator(xi,eta,Type)
    %
    % Variable lookup:
    %
    % xi: X-coordinate of the Lagrangian points.
    %
    % eta: Y-corrdinate of the Lagrangian points.
    %
    % A: ECL^-1(EC)^T or EE^T operation expressed as a matrix.
    %
    % Type: Flow variable to be solved- "vel" for velocity/vorticity field.
    %                                   "temp" for temperature field.
    
    switch Type
        case "vel"

            k = length(xi);
            A = zeros(2*k,2*k);

            for i = 1:2*k
                x = zeros(2*k,1);
                x(i) = 1;
                X = afun(x);
                A(:,i) = X;
            end
        case "temp"
            k = length(xi);
            A = zeros(k,k);
            
            for i = 1:k
                x = zeros(k,1);
                x(i) = 1;
                X = afun_temp(x);
                A(:,i) = X;
            end
    end
end