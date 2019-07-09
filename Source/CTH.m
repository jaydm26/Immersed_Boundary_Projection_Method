function gamma = CTH(Fx,Fy)
    
    global X_e_x Y_e_x X_e_y Y_e_y Nx Ny
    
    q = EdgeData(Nx,Ny);
    
    X_e_x = X_e_x';
    Y_e_x = Y_e_x';
    
    for i = 1:Nx+1
        for j = 1:Ny+2
            q.x(i,j) = H_op(X_e_x(i,j),Y_e_x(i,j),Fx);
        end
    end
    
    X_e_y = X_e_y';
    Y_e_y = Y_e_y';
    
    for i = 1:Nx+2
        for j = 1:Ny+1
            q.y(i,j) = H_op(X_e_y(i,j),Y_e_y(i,j),Fy);
        end
    end
    
    X_e_x = X_e_x';
    Y_e_x = Y_e_x';
    X_e_y = X_e_y';
    Y_e_y = Y_e_y';
    
    gamma = curl_2(q);
end