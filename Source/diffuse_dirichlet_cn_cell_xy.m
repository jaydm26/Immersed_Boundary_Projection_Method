function [t,T] = diffuse_dirichlet_cn_cell_xy(t,T,T0,rhs,dt)

    global Nx Ny Fo
    T_bc = CellData(Nx,Ny);
    T_bc = apply_bc_temp(T_bc,T0); % At time t
    rhs.x = rhs.x + T_bc.x * 0.5 * Fo;
    T_bc = apply_bc_temp(T_bc,T0); % At time t+dt
    rhs.x = rhs.x + T_bc.x * 0.5 * Fo;
    
    %% Round 1
    
    LC = zeros(Nx,Nx);
    for i = 1:Nx
        LC(i,i) = -2;
        if i == 1
            LC(i+1,i) = 1;
        elseif i == Nx
            LC(i-1,i) = 1;
        else
            LC(i+1,i) = 1;
            LC(i-1,i) = 1;
        end
    end
    
    % Now construct the AX = B problem
    
    A = eye(Nx) - 0.5 * Fo * LC;
    for j = 2:Ny+1
        B = rhs.x(2:Nx+1,j);
        a = zeros(length(A)-1,1);
        b = zeros(length(A),1);
        c = zeros(length(A)-1,1); 
        for i = 1:length(A)
            if i == 1
                b(i,1) = A(i,i);
                c(i,1) = A(i,i+1);
            elseif i == length(A)
                a(i-1,1) = A(i,i-1);
                b(i,1) = A(i,i);
            else
                a(i-1,1) = A(i,i-1);
                b(i,1) = A(i,i);
                c(i,1) = A(i,i+1);
            end
        end
        T.x(2:Nx+1,j) = trisolve(a,b,c,B,'reg');
    end
    
    
    %% Round 2
    
    rhs.x = T.x';

    LC = zeros(Ny,Ny);
    for i = 1:Ny
        LC(i,i) = -2;
        if i == 1
            LC(i+1,i) = 1;
        elseif i == Ny
            LC(i-1,i) = 1;
        else
            LC(i+1,i) = 1;
            LC(i-1,i) = 1;
        end
    end
    
    A = eye(Ny) - 0.5*Fo*LC;
    for j = 2:Nx+1
        B = rhs.x(2:Ny+1,j);
        a = zeros(length(A)-1,1);
        b = zeros(length(A),1);
        c = zeros(length(A)-1,1);
        for i = 1:length(A)
            if i == 1
                b(i,1) = A(i,i);
                c(i,1) = A(i,i+1);
            elseif i == length(A)
                a(i-1,1) = A(i,i-1);
                b(i,1) = A(i,i);
            else
                a(i-1,1) = A(i,i-1);
                b(i,1) = A(i,i);
                c(i,1) = A(i,i+1);
            end
        end
        T.x(2:Ny+1,j) = trisolve(a,b,c,B,'reg');
    end

    %% Undoing the interchange and transpose
    
    T.x = T.x';
    
    T = apply_bc_temp(T,T0);
    
    t = t + dt;
    
end
