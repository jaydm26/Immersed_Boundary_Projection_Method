function [t,gamma] = diffuse_dirichlet_cn_node_xy(params,t,rhs,gamma,gamma0,velocity)
    %DIFFUSE_DIRICHLET_CN_NODE_XY Solves the diffusion problem 
    % A * delta_gamma^* = rhs1 on the Node Space. Refer to references for 
    % further explanations.
    %
    % [t,gamma] = diffuse_dirichlet_cn_node_xy(params,t,rhs,gamma,gamma0,velocity)
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
    % rhs: Right Hand Side (NodeData) for the diffusion problem.
    %
    % T: Temperature field (NodeData) on which diffusion takes place.
    %
    % T0: Temperature field (NodeData) from the previous time step.
    %
    % velocity: Current Velocity field (EdgeData). Used to calculate the
    % convective outlet boundary condition.
    %
    % Created by Jay Mehta (18 July 2019)
    
    Nx = gamma.size(1);
    Ny = gamma.size(2);
    
    if nargin == 6
        gamma_bc = NodeData(Nx,Ny);
        gamma_bc = apply_bc_sp(params,gamma_bc,gamma0,velocity); % At time t
        diff_gamma_bc = laplacian_2(gamma_bc);
        rhs.x = rhs.x + diff_gamma_bc.x * 0.5 * params.Fo;
        gamma_bc = apply_bc_sp(params,gamma_bc,gamma0,velocity); % At time t+dt
        diff_gamma_bc = laplacian_2(gamma_bc);
        rhs.x = rhs.x + diff_gamma_bc.x * 0.5 * params.Fo;
    elseif nargin == 5
        gamma_bc = NodeData(Nx,Ny);
        gamma_bc = apply_bc_sp(params,gamma_bc,gamma0); % At time t
        diff_gamma_bc = laplacian_2(gamma_bc);
        rhs.x = rhs.x + diff_gamma_bc.x * 0.5 * params.Fo;
        gamma_bc = apply_bc_sp(params,gamma_bc,gamma0); % At time t+dt
        diff_gamma_bc = laplacian_2(gamma_bc);
        rhs.x = rhs.x + diff_gamma_bc.x * 0.5 * params.Fo;
    elseif nargin == 4
        gamma_bc = NodeData(Nx,Ny);
        gamma_bc = apply_bc_sp(params,gamma_bc); % At time t
        diff_gamma_bc = laplacian_2(gamma_bc);
        rhs.x = rhs.x + diff_gamma_bc.x * 0.5 * params.Fo;
        gamma_bc = apply_bc_sp(params,gamma_bc); % At time t+dt
        diff_gamma_bc = laplacian_2(gamma_bc);
        rhs.x = rhs.x + diff_gamma_bc.x * 0.5 * params.Fo;
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
    
    A = eye(Nx-1) - 0.5 * params.Fo * LN;
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
    gamma.x = gamma.x';

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
    
    A = eye(Ny-1) - 0.5 * params.Fo * LN;
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
    
    if nargin == 6
        gamma = apply_bc_sp(params,gamma,gamma0,velocity);
    elseif nargin == 5
        gamma = apply_bc_sp(params,gamma,gamma0);
    elseif nargin == 4
        gamma = apply_bc_sp(params,gamma);
    end
    
    t = t + params.dt;
    
end
