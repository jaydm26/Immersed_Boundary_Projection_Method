function AX = afun(X)
    
    global body_map
    
    k = length(body_map(:,1));
    Fx = X(1:k);
    Fy = X(k+1:end);
    
    gamma = CTH(Fx,Fy);
    
    [Ax,Ay] = ECL_inv(gamma);
    
    AX = [-Ax;-Ay];
end