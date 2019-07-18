function AX = afun_temp(X)
    %AFUN_TEMP Function to create Matrix A such that A = E(E)^T. This 
    % matrix is used to solve the T_w = T_f at the body-fluid interface. 
    % This is done by calculating a forcing field on the body such that T_f
    % is always equal to T_w. Refer to the reference for further
    % explanation.
    %
    % This matrix is only used for solving the temperature field. For
    % temperature, refer to afun.
    %
    % AX = AFUN_TEMP(X)
    
    HX = H_operation("cell",X);
    
    AX = E_operation("cell",HX);
end