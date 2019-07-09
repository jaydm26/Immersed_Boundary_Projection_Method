profile on
clearvars
clc

% Global Variable Declaration

global Nx Ny uB uT uL uR vB vT vL vR dx dy dt Fo

% Problem Setup
Nx = 128;
Ny = 128;

x_range = [0 1];
y_range = [0 1];

[X_xe,Y_xe] = DomainSetup(x_range,y_range,Nx,Ny,"xe");
X_topo = X_xe';
Y_topo = Y_xe';
U0 = 100;
uB = zeros(Nx+1,1); % Bottom Edge is Nx+1 long
uT = U0*ones(Nx+1,1); % Top Edge is Nx+1 long
uL = zeros(Ny+2,1); % Left Edge is Ny+2 long
uR = zeros(Ny+2,1); % Right Edge is Ny+2 long

vB = zeros(Ny+2,1);
vT = zeros(Ny+2,1);
vL = zeros(Nx+1,1);
vR = zeros(Nx+1,1);

nu = 1;
Re = U0* 1 /nu;
% Fo = 1e-4;
% dt = Fo * dx^2/nu;
dt = 0.1*dx/U0;
Fo = nu * dt/dx^2;
tf = 100*dt;
t = 0:dt:tf;

% Initialize the variables

p = CellData(Nx,Ny);
U = EdgeData(Nx,Ny);

U = apply_bc(U,1);
U_array = EdgeData(Nx,Ny);
U_array.x = zeros(Nx+1,Ny+2,length(t));
U_array.y = zeros(Nx+2,Ny+1,length(t));
U_array.x(:,:,1) = U.x;
U_array.y(:,:,1) = U.y;
dU = EdgeData(Nx,Ny);


for i = 2:length(t)
    
            % Update the Right Hand Side
            diff_u = laplacian(U);
            diff_u.x = diff_u.x * Fo;
            diff_u.y = diff_u.y * Fo;
            non_lin1 = non_linear_alt(U);
            non_lin1.x = dt * non_lin1.x;
            non_lin1.y = dt * non_lin1.y;

            rhs1 = EdgeData(Nx,Ny);
            rhs1.x = U.x + diff_u.x - non_lin1.x;
            rhs1.y = U.y + diff_u.y - non_lin1.y;
            
            % Diffusion Problem (Step 1): Using ADI
            U_star = rhs1;

            % Solve Pressure Poisson (Step 2)
            p = CellData(Nx,Ny);
%           pass dp (initial guess), rhs
            rhs2 = div(CellData(Nx,Ny),U_star);
            rhs2.x = dx/dt * rhs2.x;
            p = smoothing(p,rhs2,"sor","tol",1e-3);
            
            % Pressure Obtained above is the correct pressure (Step 3)
            
            % Correct the Velocity (Step 4)
            rhs3 = grad(p);
            rhs3.x = rhs3.x * dt/dx;
            rhs3.y = rhs3.y * dt/dx;
            U.x = U_star.x - rhs3.x;
            U.y = U_star.y - rhs3.y;
            
            % Add the correction to the true velocities (Step 5)
            
%             U.x = U_array.x(:,:,i-1) + dU.x; 
%             U.y = U_array.y(:,:,i-1) + dU.y;
            U = apply_bc(U,1);
            U_array.x(:,:,i) = U.x;
            U_array.y(:,:,i) = U.y;
            
            Un = interpol(NodeData(Nx,Ny),U,1);
            U_1 = EdgeData(Nx,Ny);
            U_1.x = U_array.x(:,:,i-1);
            Un_1 = interpol(NodeData(Nx,Ny),U_1,1);
            conv = max(max(abs(Un.x - Un_1.x)));
            if conv < 1e-6
                break
            end
end
% profile viewer 
figure
surf(X_topo,Y_topo,U_array.x(:,:,end));
% figure
% contour(X_topo,Y_topo,U_array.x(:,:,end));