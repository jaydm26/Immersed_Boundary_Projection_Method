function [ u ] = poisson2d( f, griddata )
%POISSON2D Solve a Poisson problem on a 2-d staggered grid
%   u = POISSON2D(f,gridparams) solves the 2d Poisson problem
%    
%            (d^2/dx^2 + d^2/dy^2) u = f
%   
%   for the scalar field u, given the
%   force field f. The grid information is provided in griddata. It
%   is assumed that grid spacing is uniform and the same in both
%   directions.
%
%   If the size of f is (M+1)x(N+1), where M and N are
%   described in griddata, then it is assumed that f and u are node-based
%   variables. In this form, it is solved as a Dirichlet problem with
%   homogeneous boundary conditions (and any non-zero boundary values are
%   assumed to already be contained in the forcing f.
%
%   If f is M x N, then f and u are assumed cell-centered variables and the
%   problem is treated as a Neumann problem with homogeneous boundary
%   values. It is assumed that non-zero Neumann data is included in f.
%
%   If f is neither of these dimensions, then the routine returns with an
%   error.
%
%   MAE 259A
%   J. D. Eldredge
%   3/11/2014

%% Unpack some information from gridparams
M = griddata.M;  % Number of cells in x direction
N = griddata.N;  % Number of cells in y direction

dxsq = griddata.dx*griddata.dy;

%% Check the size of f and determine problem type.
size1 = size(f,1);
size2 = size(f,2);
padx = 0; pady = 0;
if (size1==M+1),
    % Node-based problem
    f = f(2:M,:);
    padx = 1;
    solvetype = 1; % Discrete sine transform
elseif (size1==M),
    % Cell-based problem
    solvetype = 2; % Discrete cosine transform
else
    disp('Error in poisson2d. Inconsistent size.');
    return
end
if (size2==N+1),
    % Node-based problem
    f = f(:,2:N);
    pady = 1;
    if (solvetype ~= 1),
        disp('Error in poisson2d. Inconsistent data organization.');
        return
    end
elseif (size2==N),
    % Cell-based problem
    if (solvetype ~= 2),
        disp('Error in poisson2d. Inconsistent data organization.');
        return
    end
else
    disp('Error in poisson2d. Inconsistent size.');
    return
end
   
if (solvetype == 1),
    %% Dirichlet type
    
    
    %% Carry out the sine transform of the force field
    % In x direction
    f = dst(f);
    % In y direction. Need to transpose to ensure columns are transformed
    f = dst(f');
    f = f';
    
    %% Divide by eigenvalues
    u = zeros(size(f));
    for n = 1:N-1
        cosn = cos(pi*n/N);
        for m = 1:M-1
            cosm = cos(pi*m/M);
            u(m,n) = 0.5*dxsq*f(m,n)/(cosm+cosn-2);
        end
    end
    
    %% Transform back
    % In x direction
    u = idst(u);
    % In y direction. Need to transpose to ensure columns are transformed
    u = idst(u');
    u = u';
    
    %% Pad with boundary values
    if (padx),
        u = [zeros(size(u,1),1) u zeros(size(u,1),1)];
    end
    if (pady),
        u = [zeros(1,size(u,2)); u; zeros(1,size(u,2))];
    end

elseif (solvetype == 2),
    %% Neumann type
    
    %% Carry out the cosine transform of the force field
    % In x direction
    f = dctstag(f);
    % In y direction. Need to transpose to ensure columns are transformed
    f = dctstag(f');
    f = f';
    
    %% Divide by eigenvalues. Skip the (0,0) eigenvalue.
%     u = zeros(size(f));
%     for n = 0:N-1
%         cosn = cos(pi*n/N);
%         for m = 0:M-1
%             cosm = cos(pi*m/M);
%             if (m==0 && n==0), continue; end;
%             u(m+1,n+1) = 0.5*dxsq*f(m+1,n+1)/(cosm+cosn-2);
%         end
%     end
    lam = repmat(cos(pi*(0:M-1)'/M),1,N) ...
        + repmat(cos(pi*(0:N-1)/N),M,1) ...
        - 2*ones(M,N);
    lam(1,1) = 1; f(1,1) = 0;
    u = 0.5*dxsq*f./lam;

    %% Transform back
    % In x direction
    u = idctstag(u);
    % In y direction. Need to transpose to ensure columns are transformed
    u = idctstag(u');
    u = u';
    
end
