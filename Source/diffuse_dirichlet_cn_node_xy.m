function [t,gamma] = diffuse_dirichlet_cn_node_xy(t,rhs,gamma,gamma0,velocity)

    global Nx Ny Fo dt
    
    if nargin == 5
        gamma_bc = NodeData(Nx,Ny);
        gamma_bc = apply_bc_sp(gamma_bc,gamma0,velocity); % At time t
        rhs.x = rhs.x + gamma_bc.x * 0.5 * Fo;
        gamma_bc = apply_bc_sp(gamma_bc,gamma0,velocity); % At time t+dt
        rhs.x = rhs.x + gamma_bc.x * 0.5 * Fo;
    elseif nargin == 4
        gamma_bc = NodeData(Nx,Ny);
        gamma_bc = apply_bc_sp(gamma_bc,gamma0); % At time t
        rhs.x = rhs.x + gamma_bc.x * 0.5 * Fo;
        gamma_bc = apply_bc_sp(gamma_bc,gamma0); % At time t+dt
        rhs.x = rhs.x + gamma_bc.x * 0.5 * Fo;
    elseif nargin == 3
        gamma_bc = NodeData(Nx,Ny);
        gamma_bc = apply_bc_sp(gamma_bc); % At time t
        rhs.x = rhs.x + gamma_bc.x * 0.5 * Fo;
        gamma_bc = apply_bc_sp(gamma_bc); % At time t+dt
        rhs.x = rhs.x + gamma_bc.x * 0.5 * Fo;
    end
    %% Round 1
    
    LN = zeros(Nx-1,Nx-1);
    for i = 1:Nx-1
        LN(i,i) = -2;
        if i == 1
            LN(i+1,i) = 1;
        elseif i == Nx-1
            LN(i-1,i) = 1;
        else
            LN(i+1,i) = 1;
            LN(i-1,i) = 1;
        end
    end
    
    % Now construct the AX = B problem
    
    A = eye(Nx-1) - 0.5 * Fo * LN;
    for j = 2:Ny
        B = rhs.x(2:Nx,j);
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
        gamma.x(2:Nx,j) = trisolve(a,b,c,B,'reg');
    end
    
    
    %% Round 2
    
    rhs.x = gamma.x';

    LN = zeros(Ny-1,Ny-1);
    for i = 1:Ny-1
        LN(i,i) = -2;
        if i == 1
            LN(i+1,i) = 1;
        elseif i == Ny-1
            LN(i-1,i) = 1;
        else
            LN(i+1,i) = 1;
            LN(i-1,i) = 1;
        end
    end
    
    A = eye(Ny-1) - 0.5*Fo*LN;
    for j = 2:Nx
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
        gamma.x(2:Ny,j) = trisolve(a,b,c,B,'reg');
    end

    %% Undoing the interchange and transpose
    
    gamma.x = gamma.x';
    if nargin == 5
        gamma = apply_bc_sp(gamma,gamma0,velocity);
    elseif nargin == 4
        gamma = apply_bc_sp(gamma,gamma0);
    elseif nargin == 3
        gamma = apply_bc_sp(gamma);
    end
    
    t = t + dt;
    
end
