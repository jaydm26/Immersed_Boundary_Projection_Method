function out = E_op_simple(xi,eta,s)
    global Nx Ny X_n Y_n
    for i = 1:Nx+1
        for j = 1:Ny+1
            out(i,j) = s(i,j) * ddf_roma_2D(X_n(i,j)-xi,Y_n(i,j)-eta);
        end
    end
    out = sum(sum(out));
end