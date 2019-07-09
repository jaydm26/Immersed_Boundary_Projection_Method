function [t,U] = diffuse_dirichlet_cn_edge_xy(t,U,rhs,dt)

    global Nx Ny Fo
    diff_u = laplacian(U);
    rhs.x = rhs.x + diff_u.x * Fo;
    rhs.y = rhs.y + diff_u.y * Fo;
    U_bc = EdgeData(Nx,Ny);
    U_bc = apply_bc(U_bc,1);
    rhs.x = rhs.x + U_bc.x * 0.5 * Fo;
    rhs.y = rhs.y + U_bc.y * 0.5 * Fo;
    U_bc = apply_bc(U_bc,1);
    rhs.x = rhs.x + U_bc.x * 0.5 * Fo;
    rhs.y = rhs.y + U_bc.y * 0.5 * Fo;
    
    
    %% Round 1 (X & Y are in their place)
    
    %% For X-direction
    LE = zeros(Nx-1,Nx-1);
    for i = 1:Nx-1
        LE(i,i) = -2;
        if i == 1
            LE(i+1,i) = 1;
        elseif i == Nx-1
            LE(i-1,i) = 1;
        else
            LE(i+1,i) = 1;
            LE(i-1,i) = 1;
        end
    end
    
    % Now constructing the AX=B problem.

    A = eye(Nx-1) - 0.5* Fo*LE;
    for j = 2:Ny+1
        B = rhs.x(2:Nx,j);% + 0.5* Fo * B2(:,j);
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
        U.x(2:Nx,j) = trisolve(a,b,c,B,'reg');
    end

%% For Y-direction
    
    LE = zeros(Nx,Nx);
    for i = 1:Nx
        LE(i,i) = -2;
        if i == 1
            LE(i+1,i) = 1;
            LE(i,i) = -3;
        elseif i == Nx
            LE(i-1,i) = 1;
            LE(i,i) = -3;
        else
            LE(i+1,i) = 1;
            LE(i-1,i) = 1;
        end
    end
    
    % Now constructing the AX=B problem.

    A = eye(Nx) - 0.5*Fo*LE;
    for j = 2:Ny
        B = rhs.y(2:Nx+1,j);% + 2*Fo*B2(:,j);
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
        U.y(2:Nx+1,j) = trisolve(a,b,c,B,'reg');
    end

    %% Round 2 (X & Y are interchanged and transposed)
    
    rhs.x = U.y';
    rhs.y = U.x';
    
    %% For X-direction
    LE = zeros(Ny-1,Ny-1);
    for i = 1:Ny-1
        LE(i,i) = -2;
        if i == 1
            LE(i+1,i) = 1;
        elseif i == Ny-1
            LE(i-1,i) = 1;
        else
            LE(i+1,i) = 1;
            LE(i-1,i) = 1;
        end
    end

    A = eye(Ny-1) - 0.5*Fo*LE;
    for j = 2:Nx+1
        B = rhs.x(2:Ny,j);
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
        U.x(2:Ny,j) = trisolve(a,b,c,B,'reg');
    end

%% For Y-direction
    
    LE = zeros(Ny,Ny);
    for i = 1:Ny
        LE(i,i) = -2;
        if i == 1
            LE(i+1,i) = 1;
            LE(i,i) = -3;
        elseif i == Ny
            LE(i-1,i) = 1;
            LE(i,i) = -3;
        else
            LE(i+1,i) = 1;
            LE(i-1,i) = 1;
        end
    end

    A = eye(Ny) - 0.5*Fo*LE;
    for j = 2:Nx
        B = rhs.y(2:Ny+1,j);
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
            U.y(2:Ny+1,j) = trisolve(a,b,c,B,'reg');
    end
    
    %% Undoing the interchange and transpose
    
    U_temp = U;
    U.x = U_temp.y';
    U.y = U_temp.x';
    
    U.x(1,2:Ny+1) = 0; % Bottom
    U.x(end,2:Ny+1) = 0; % Top
    U.x(2:Nx,1) = -U.x(2:Nx,2);
    U.x(2:Nx,end) = -U.x(2:Nx,end-1);
    
    U.y(1,2:Ny) = -U.y(2,2:Ny); % Bottom
    U.y(end,2:Ny) = -U.y(end-1,2:Ny); % Top
    U.y(2:Nx+1,1) = 0;
    U.y(2:Nx+1,end) = 0;
    
    t = t + dt;
    
end
