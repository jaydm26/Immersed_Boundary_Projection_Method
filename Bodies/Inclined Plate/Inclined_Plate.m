%% Clear Everything
clear all
clc
rmpath('/Users/jaymehta/Desktop/Research UCLA/IBPM/Source')
rmpath('/Users/jaymehta/Desktop/Research UCLA/IBPM/Source/dst_idst')
%% Add Paths
addpath('/Users/jaymehta/Desktop/Research UCLA/IBPM/Source')
addpath('/Users/jaymehta/Desktop/Research UCLA/IBPM/Source/dst_idst')
%% Set up the problem domain and the problem object

global Nx Ny dx dy dt X_n Y_n body_map x_range y_range g_hat_n U...
    X_e_x Y_e_x X_e_y Y_e_y Fo Co
 
% Domain

Nx = 64;
Ny = 64;

x_range = [-5 5];
y_range = [-5 5];

%% Create the L^-1 operator using Lattice Green's function

g_hat_n = L_inv("node");

%% Setting up the Domain of the Problem

[X_n, Y_n] = DomainSetup(x_range,y_range,Nx,Ny,"node");
[X_e_x, Y_e_x] = DomainSetup(x_range,y_range,Nx,Ny,"xe");
[X_e_y, Y_e_y] = DomainSetup(x_range,y_range,Nx,Ny,"ye");
[X_c,Y_c] = DomainSetup(x_range,y_range,Nx,Ny,"cell");

%% Setting up the Object

xL = -0.5 * sind(30);
yL = 0.5 * cosd(30);
xR = 0.5 * sind(30);
yR = -0.5 * cosd(30);

[body_map,L,theta] = Line_Builder(xL,xR,yL,yR);

%% Forming the A matrix for PCG. Load it from the side bar

A = MatrixA_Generator("vel");

%% Initializing the variables pre-solving

% Employ the Streamfunction - Vorticity Method to avoid the coupled forcing
% functions which show up in the Pressure-Velocity Method.

% Of course, in this method, we will employ the usual CN/AB2
% discretization, and then use the Block LU Decomposition. This will allow
% us to use the delta formulation.
U = 1;
V = 0;

Re = 40;
nu = U * L / Re;

Co = 0.1;
Fo = 5;
dt = min([Fo * dx^2/nu,Co*dx]);

Fo = nu * dt/dx^2;
Co = dt/dx;
t_steady = 1/nu;
tf = t_steady;
t = 0:dt:tf;

%% Pre-Setup

velocity = EdgeData(Nx,Ny); % Velocity Field
gamma = NodeData(Nx,Ny); % Vorticity
Fx = zeros(length(body_map(:,1)),1);
Fy = zeros(length(body_map(:,2)),1);
ub = zeros(length(body_map(:,1)),1); % X-component of Velocity on the body
vb = zeros(length(body_map(:,1)),1); % Y-component of Velocity on the body
sf = NodeData(Nx,Ny);
tol = 1e-1;
Drag = 0;
Lift = 0;

%% Boundary Conditions on Gamma

gamma0 = gamma;
gamma = apply_bc_sp(gamma,gamma0,velocity);

%% Starter (Euler)
for a = 2
    % Set up R1
    gamma0 = gamma;
    diff_gamma = laplacian_2(gamma);
    nl = non_linear_alt(velocity);
    nl.x = -nl.x;
    nl.y = -nl.y;
    nl = curl_2(nl);
    
    ff = CTH(Fx,Fy);
    
    rhs1 = NodeData(Nx,Ny);
    rhs1.x = Fo * diff_gamma.x + dt * nl.x + dt * ff.x;
    
    % Solve the diffusion problem
    
    delta_gamma_star = rhs1;
    
    % Now set up R2. Note here that R2 will have two components.
    
    temp = NodeData(Nx,Ny);
    temp.x = gamma.x + delta_gamma_star.x;
    [ub2,vb2] = ECL_inv(temp);
    
    % Forming the two rhs2's
    
    rhs2_x = (ub - U * ones(length(body_map(:,1)),1)) + ub2;
    rhs2_y = (vb - V * ones(length(body_map(:,2)),1)) + vb2;
    
    rhs2 = [rhs2_x;rhs2_y];
    % Now obtain the delta_f in the two directions by pcg
    delta_f = pcg(A, rhs2, 1e-1);
    delta_f = delta_f./dt;
    
    % Now delta_f_star in the X and Y direction are the correct delta_f in the X and Y direction 
    
    delta_f_x = delta_f(1:length(body_map(:,1)));
    delta_f_y = delta_f(length(body_map(:,1))+1:end);
    
    % Now we correct for the vorticity to ensure we have no slip on the body
    
    gamma_c = CTH(delta_f_x,delta_f_y);
    
    delta_gamma = delta_gamma_star.x + dt * gamma_c.x;
    
    % Update gamma and the forcing function % Adding the unstable component
    
    gamma.x = gamma.x + delta_gamma;% + gamma_inst.x;
    gamma = apply_bc_sp(gamma,gamma0,velocity);
    Fx = Fx + delta_f_x;
    Fy = Fy + delta_f_y;
    
    % Obtain the streamfunction and velocity
    
    sf = L_inv_operation(-gamma.x,"node"); 
    
    velocity = curl_2(sf);
    velocity.x = velocity.x + U * ones(Nx+1,Ny+2);
    velocity.y = velocity.y;
    
    % Calculate Lift and Drag
    
    Hq = H_operation("edge",Fx,Fy);
    
    Drag(a) = -sum(sum(Hq.x))*dx^2;
    Lift(a) = sum(sum(Hq.y))*dx^2;
    
    % Checking the residual for each step and breaking the loop if convergence has reached
    
    conv = 1/Nx * norm(delta_gamma(2:Nx,2:Ny))/dt / norm(gamma.x(2:Nx,2:Ny))/dt;
    if conv <= tol
        break
    end
end

%% CN-AB2

for a = 3:length(t)
    
    % Set up R1
    gamma0 = gamma;
    diff_gamma = laplacian_2(gamma);
    rhs1 = NodeData(Nx,Ny);
    rhs1.x = -dt * 0.5 * nl.x;
    nl = non_linear_alt(velocity);
    nl.x = -nl.x;
    nl.y = -nl.y;
    nl = curl_2(nl);
    
    ff = CTH(Fx,Fy);
    
    rhs1.x = Fo * diff_gamma.x + 1.5 * dt * nl.x + dt * ff.x;
    
    % Solve the diffusion problem
    
    [~,delta_gamma_star] = diffuse_dirichlet_cn_node_xy(t(a),gamma,gamma0,rhs1,dt,velocity);
    
    % Now set up R2. Note here that R2 will have two components.
    
    temp = NodeData(Nx,Ny);
    temp.x = gamma.x + delta_gamma_star.x;
    [ub2,vb2] = ECL_inv(temp);
    
    rhs2_x = (ub - U * ones(length(body_map(:,1)),1)) + ub2;
    rhs2_y = (vb - V * ones(length(body_map(:,2)),1)) + vb2;
    
    rhs2 = [rhs2_x;rhs2_y];
    % Now obtain the delta_f in the two directions by pcg
    delta_f = pcg(A, rhs2, 1e-1,40);
    delta_f = delta_f./dt;
    % Now delta_f_star in the X and Y direction are the correct delta_f in the X and Y direction
    
    % Assigning
    
    delta_f_x = delta_f(1:length(body_map(:,1)));
    delta_f_y = delta_f(length(body_map(:,1))+1:end);
    
    % Now we correct for the vorticity to ensure we have no slip on the body
    
    gamma_c = CTH(delta_f_x,delta_f_y);
    
    delta_gamma = delta_gamma_star.x + dt * gamma_c.x;
    
    % Update gamma and the forcing function
    
    gamma.x = gamma.x + delta_gamma;% + gamma_inst.x;
    gamma = apply_bc_sp(gamma,gamma0,velocity);
    Fx = Fx + delta_f_x;
    Fy = Fy + delta_f_y;
    
    % Obtain the streamfunction and velocity
    
    sf = L_inv_operation(-gamma.x,"node");
    
    velocity = curl_2(sf);
    velocity.x = velocity.x + U * ones(Nx+1,Ny+2);
    velocity.y = velocity.y;
    
    % Calculate Lift and Drag
    
    Hq = H_operation("edge",Fx,Fy);
    
    Drag(a) = -sum(sum(Hq.x))*dx^2;
    Lift(a) = sum(sum(Hq.y))*dx^2;
    
    % Checking the residual for each step and breaking the loop if convergence has reached
    
    conv = 1/Nx * norm(delta_gamma(2:Nx,2:Ny))/dt / norm(gamma.x(2:Nx,2:Ny))/dt;
    if conv <= tol
        break
    end
end
%% Streamlines

xi = body_map(:,1);
eta = body_map(:,2);
f1 = figure;
contour(X_n./(L),Y_n./(L),(sf.x'),-3:0.4:3)
hold on
plot(xi./(L),eta./(L),"LineWidth",2);
hold off
pbaspect([1 1 1])
title({'Streamlines for a Flat Plate of L = ',num2str(L),' flow is of uniform velocity of U = ',num2str(U)})
xlabel("X/L")
ylabel("Y/L")
f1.WindowState = 'fullscreen';

%% Quiver Plot

xi = body_map(:,1);
eta = body_map(:,2);
f4 = figure;
qx = NodeData(Nx,Ny);
qy = NodeData(Nx,Ny);
qx = interpol(qx,velocity,1);
qy = interpol(qy,velocity,2);
quiver(X_n./(L),Y_n./(L),qx.x',qy.x');
hold on
plot(xi./(L),eta./(L),"LineWidth",2);
hold off
pbaspect([1 1 1])
title({'Velocity for a Flat Plate of L = ',num2str(L),' flow is of uniform velocity of U = ',num2str(U)})
xlabel("X/L")
ylabel("Y/L")
axis tight
f4.WindowState = 'fullscreen'