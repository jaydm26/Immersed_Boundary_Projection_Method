function AX = afun_temp(params,domain,xi,eta,X)
    %AFUN_TEMP Function to create Matrix A such that A = E(E)^T. This 
    % matrix is used to solve the T_w = T_f at the body-fluid interface. 
    % This is done by calculating a forcing field on the body such that T_f
    % is always equal to T_w. Refer to the reference for further
    % explanation.
    %
    % This matrix is only used for solving the temperature field. For
    % temperature, refer to afun.
    %
    % AX = afun_temp(params,domain,xi,eta,X)
    %
    % Variable lookup:
    %
    % xi,eta: X and Y coordinate of the Langrangian Points.
    %
    % params: flow parameters
    %
    % domain: domain parameters
    %
    % Created by Jay Mehta (18 July 2019)
    
    HX = H_operation(params,domain,"cell",xi,eta,X);
    
    AX = E_operation(params,domain,xi,eta,HX);
end