function [t,T] = diffuse_dirichlet_cn_cell_xy(params,t,rhs,T,T0,velocity)
    %DIFFUSE_DIRICHLET_CN_CELL_XY Solves the diffusion problem 
    % A * delta_T^* = rhs1 on the Cell Space. Refer to references for 
    % further explanations.
    %
    % [t,T] = diffuse_dirichlet_cn_cell_xy(params,t,rhs,T,T0,velocity)
    %
    % Variable lookup:
    %
    % Nx: Number of divisions in the X-direction.
    %
    % Ny: Number of divisions in the Y-direction.
    %
    % t: Current time.
    %
    % params: flow parameters.
    %
    % rhs: Right Hand Side (CellData) for the diffusion problem.
    %
    % T: Temperature field (CellData) on which diffusion takes place.
    %
    % T0: Temperature field (CellData) from the previous time step.
    %
    % velocity: Current Velocity field (EdgeData). Used to calculate the
    % convective outlet boundary condition.
    %
    % Created by Jay Mehta (18 July 2019)
    
    Nx = T.size(1);
    Ny = T.size(2);
    
    if nargin == 6
        T_bc = CellData(Nx,Ny);
        T_bc = apply_bc_temp(params,T_bc,T0,velocity); % At time t
        diff_T_bc = laplacian_2(T_bc);
        rhs.x = rhs.x + diff_T_bc.x * 0.5 * params.Fo_t;
        T_bc = apply_bc_temp(params,T_bc,T0,velocity); % At time t+dt
        diff_T_bc = laplacian(T_bc);
        rhs.x = rhs.x + diff_T_bc.x * 0.5 * params.Fo_t;
    elseif nargin == 5
        T_bc = CellData(Nx,Ny);
        T_bc = apply_bc_temp(params,T_bc,T0); % At time t
        diff_T_bc = laplacian_2(T_bc);
        rhs.x = rhs.x + diff_T_bc.x * 0.5 * params.Fo_t;
        T_bc = apply_bc_temp(params,T_bc,T0); % At time t+dt
        diff_T_bc = laplacian_2(T_bc);
        rhs.x = rhs.x + diff_T_bc.x * 0.5 * params.Fo_t;
    elseif nargin == 4
        T_bc = CellData(Nx,Ny);
        T_bc = apply_bc_temp(params,T_bc); % At time t
        diff_T_bc = laplacian_2(T_bc);
        rhs.x = rhs.x + diff_T_bc.x * 0.5 * params.Fo_t;
        T_bc = apply_bc_temp(params,T_bc); % At time t+dt
        diff_T_bc = laplacian_2(T_bc);
        rhs.x = rhs.x + diff_T_bc.x * 0.5 * params.Fo_t;
    else
        error("Too many inputs.")
    end
    
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
    
    A = eye(Nx) - 0.5 * params.Fo_t * LC;
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
    T.x = T.x';

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
    
    A = eye(Ny) - 0.5 * params.Fo_t * LC;
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
    
    if nargin == 6
        T = apply_bc_temp(params,T,T0,velocity);
    elseif nargin == 5
        T = apply_bc_temp(params,T,T0);
    elseif nargin == 4
        T = apply_bc_temp(params,T);
    else
        error("Too many inputs.")
    end
    
    t = t + params.dt;
    
end
