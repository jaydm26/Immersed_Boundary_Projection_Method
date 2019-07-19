function [t,velocity] = diffuse_dirichlet_cn_edge_xy(params,bc,t,rhs,velocity)
    %DIFFUSE_DIRICHLET_CN_EDGE_XY Solves the diffusion problem 
    % A * delta_q^* = rhs1 on the Edge Space. Refer to references for 
    % further explanations.
    %
    % [t,velocity] = diffuse_dirichlet_cn_edge_xy(params,bc,t,rhs,velocity)
    %
    % Variable lookup:
    %
    % Nx: Number of divisions in the X-direction.
    %
    % Ny: Number of divisions in the Y-direction.
    %
    % t: Current time.
    %
    % velocity: Current Velocity field (EdgeData).
    %
    % params: flow parameters.
    % 
    % bc: Boundary conditions for the Edge Field.
    %
    % rhs: Right Hand Side (EdgeData) for the diffusion problem.
    %
    % Created by Jay Mehta (18 July 2019)
    
    Nx = velocity.size(1);
    Ny = velocity.size(2);
    
    velocity_bc = EdgeData(Nx,Ny);
    velocity_bc = apply_bc(bc,velocity_bc,t);
    rhs.x = rhs.x + velocity_bc.x * 0.5 * params.Fo;
    rhs.y = rhs.y + velocity_bc.y * 0.5 * params.Fo;
    velocity_bc = apply_bc(bc,velocity_bc,t);
    rhs.x = rhs.x + velocity_bc.x * 0.5 * params.Fo;
    rhs.y = rhs.y + velocity_bc.y * 0.5 * params.Fo;
    
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

    A = eye(Nx-1) - 0.5* params.Fo * LE;
    for j = 2:Ny+1
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
        velocity.x(2:Nx,j) = trisolve(a,b,c,B,'reg');
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

    A = eye(Nx) - 0.5 * params.Fo * LE;
    for j = 2:Ny
        B = rhs.y(2:Nx+1,j);
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
        velocity.y(2:Nx+1,j) = trisolve(a,b,c,B,'reg');
    end

    %% Round 2 (X & Y are interchanged and transposed)
    
    rhs.x = velocity.y';
    rhs.y = velocity.x';
    
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

    A = eye(Ny-1) - 0.5 * params.Fo * LE;
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
        velocity.x(2:Ny,j) = trisolve(a,b,c,B,'reg');
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

    A = eye(Ny) - 0.5 * params.Fo * LE;
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
            velocity.y(2:Ny+1,j) = trisolve(a,b,c,B,'reg');
    end
    
    %% Undoing the interchange and transpose
    
    velocity_temp = velocity;
    velocity.x = velocity_temp.y';
    velocity.y = velocity_temp.x';
    
    velocity = apply_bc(bc,velocity,t);
    
    t = t + dt;
    
end
