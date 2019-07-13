function AX = afun_temp(X)
    
    HX = H_operation("cell",X);
    
    AX = E_operation("cell",HX);
end