function AX = afun(X,xi,eta)
    %AFUN Function to create Matrix A such that A = ECL^-1(EC)^T. This 
    % matrix is used to solve the no-slip condition on the body by 
    % calculating the force delta_f. Refer to the reference for further 
    % explanation.
    %
    % This matrix is only used for solving the fluid flow using the
    % nullspace method. For temperature, refer to afun_temp.
    %
    % AX = AFUN(X,body_map)
    %
    % Variable lookup:
    %
    % xi,eta: X and Y coordinate of the Langrangian Points.
    
    k = length(xi);
    Fx = X(1:k);
    Fy = X(k+1:end);
    
    gamma = CTH(Fx,Fy);
    
    [Ax,Ay] = ECL_inv(gamma);
    
    AX = [-Ax;-Ay];
end