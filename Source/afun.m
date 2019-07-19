function AX = afun(params,domain,g_hat,xi,eta,X)
    %AFUN Function to create Matrix A such that A = ECL^-1(EC)^T. This 
    % matrix is used to solve the no-slip condition on the body by 
    % calculating the force delta_f. Refer to the reference for further 
    % explanation.
    %
    % This matrix is only used for solving the fluid flow using the
    % nullspace method. For temperature, refer to afun_temp.
    %
    % AX = afun(params,domain,g_hat,xi,eta,X)
    %
    % Variable lookup:
    %
    % xi,eta: X and Y coordinate of the Langrangian Points.
    %
    % params: flow parameters
    %
    % domain: domain parameters
    %
    % g_hat: FFT2 of Lattice Green's Function.
    %
    % Created by Jay Mehta (18 July 2019)
    
    k = length(xi);
    Fx = X(1:k);
    Fy = X(k+1:end);
    
    gamma = CTH(params,domain,xi,eta,Fx,Fy);
    
    [Ax,Ay] = ECL_inv(params,domain,g_hat,xi,eta,gamma);
    
    AX = [-Ax;-Ay];
end