%% Clear Everything
clear all
clc
rmpath('/Users/jaymehta/Desktop/Research UCLA/IBPM/Source')
rmpath('/Users/jaymehta/Desktop/Research UCLA/IBPM/Source/dst_idst')

%% Add Paths
addpath('/Users/jaymehta/Desktop/Research UCLA/IBPM/Source')
addpath('/Users/jaymehta/Desktop/Research UCLA/IBPM/Source/dst_idst')
%% Set up the problem domain and the problem object

global Nx Ny dx dy dt X_n Y_n body_map x_range y_range g_hat_n g_hat_c U...
    X_e_x Y_e_x X_e_y Y_e_y X_c Y_c Fo Co Fo_t
 
% Domain

Nx = 64;
Ny = 64;

x_range = [-10 10];
y_range = [-10 10];

%% Create the L^-1 operator using Lattice Green's function

g_hat_n = L_inv("node");
g_hat_c = L_inv("cell");

%% Setting up the Domain of the Problem

[X_n, Y_n] = DomainSetup(x_range,y_range,Nx,Ny,"node");
[X_e_x, Y_e_x] = DomainSetup(x_range,y_range,Nx,Ny,"xe");
[X_e_y, Y_e_y] = DomainSetup(x_range,y_range,Nx,Ny,"ye");
[X_c,Y_c] = DomainSetup(x_range,y_range,Nx,Ny,"cell");

%% Setting up the Object

R = 1;
xc = 0;
yc = 0;

N_theta = floor(2*pi*R/dx);
d_theta = 2*pi/N_theta;

theta_range = 0:d_theta:2*pi-d_theta;

body_map = zeros(length(theta_range),2);

for i = 1:length(theta_range)
    body_map(i,1) = xc + R*cos(theta_range(i));
    body_map(i,2) = yc + R*sin(theta_range(i));
end

%% Forming the A matrix for PCG. Load it from the side bar

A = MatrixA_Generator("vel");
A_temp = MatrixA_Generator("temp");

%% Initializing the variables pre-solving

% Employ the Streamfunction - Vorticity Method to avoid the coupled forcing
% functions which show up in the Pressure-Velocity Method.

% Of course, in this method, we will employ the usual CN/AB2
% discretization, and then use the Block LU Decomposition. This will allow
% us to use the delta formulation.
U = 1;
V = 0;

Re = 40;
Pr = 1;
nu = U * R / Re;
alpha = nu/Pr;

Co = 0.1;
Fo = 5;
Fo_t = 5;
dt = min([Fo * dx^2/nu,Co*dx,Fo_t * dx^2/alpha]);

Fo = nu * dt/dx^2;
Fo_t = alpha * dt/dx^2;
Co = U * dt/dx;
t_steady = max([(2*R)^2/nu,(2*R)^2/alpha]);
tf = t_steady;
time_range = 0:dt:tf;

%% Pre-Setup

velocity = EdgeData(Nx,Ny); % Velocity Field
gamma = NodeData(Nx,Ny); % Vorticity
Fx = zeros(length(body_map(:,1)),1);
Fy = zeros(length(body_map(:,2)),1);
ub = zeros(length(body_map(:,1)),1); % X-component of Velocity on the body
vb = zeros(length(body_map(:,1)),1); % Y-component of Velocity on the body
sf = NodeData(Nx,Ny);

T = CellData(Nx,Ny); % Temperature Field
FT = zeros(length(body_map(:,1)),1);
Tb = 1;

tol = 1e-2;
Drag = 0;
Lift = 0;

%% Boundary Conditions on Gamma

gamma0 = gamma;
gamma = apply_bc_sp(gamma,gamma0,velocity);

T0 = T;
T = apply_bc_temp(T,T0,velocity);

%% Starter (Euler)
for t = 2
    % Set up R1
    gamma0 = gamma;
    diff_gamma = laplacian_2(gamma);
    nl = non_linear_alt(velocity);
    nl = curl_2(nl);
    
    Hf = CTH(Fx,Fy);
    
    rhs1 = NodeData(Nx,Ny);
    rhs1.x = Fo * diff_gamma.x - dt * nl.x + dt * Hf.x;
    
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
    
    Drag(t) = -sum(sum(Hq.x))*dx^2;
    Lift(t) = sum(sum(Hq.y))*dx^2;
    
    % After solving for the flow field, we can solve for the temperature
    % field.
    diff_T = laplacian_2(T);
    nlt = non_linear_temp(velocity,T);
    HfT = H_operation("cell",FT);
    
    rhs3 = CellData(Nx,Ny);
    rhs3.x = Fo_t * diff_T.x - dt * nlt.x + dt * HfT.x;
    
    % Solve the diffusion problem
    
    delta_T_star = rhs3;
    
    temp = CellData(Nx,Ny);
    temp.x = T.x + delta_T_star.x;
    Tb2 = E_operation("cell",temp);
    
    rhs4 = Tb - Tb2;
    delta_FT = pcg(A_temp, rhs4, 1e-1);
    delta_FT = delta_FT./dt;
    
    HDFT = H_operation("cell",delta_FT);
    delta_T = delta_T_star.x + dt * HDFT.x;
    
    % Assigning
    T.x = T.x + delta_T;
    T = apply_bc_temp(T,T0,velocity);
    FT = FT + delta_FT;
    
    % Checking the residual for each step and breaking the loop if convergence has reached
    
    conv_flow = 1/Nx * norm(delta_gamma(2:Nx,2:Ny))/dt / norm(gamma.x(2:Nx,2:Ny))/dt;
    conv_therm = 1/Nx * norm(delta_T(2:Nx+1,2:Ny+1))/dt / norm(T.x(2:Nx+1,2:Ny+1))/dt;
    if conv_flow <= tol && conv_therm <= tol
        break
    end
end

%% CN-AB2

for t = 96:length(time_range)
    
    % Set up R1
    gamma0 = gamma;
    diff_gamma = laplacian_2(gamma);
    rhs1 = NodeData(Nx,Ny);
    rhs1.x = dt * 0.5 * nl.x;
    nl = non_linear_alt(velocity);
    nl = curl_2(nl);
    
    Hf = CTH(Fx,Fy);
    
    rhs1.x = rhs1.x + Fo * diff_gamma.x - 1.5 * dt * nl.x + dt * Hf.x;
    
    % Solve the diffusion problem
    
    [~,delta_gamma_star] = diffuse_dirichlet_cn_node_xy(time_range(t),gamma,gamma0,rhs1,dt,velocity);
    
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
    
    Drag(t) = -sum(sum(Hq.x))*dx^2;
    Lift(t) = sum(sum(Hq.y))*dx^2;
    
    % After solving for the flow field, we can solve for the temperature
    % field.
    diff_T = laplacian_2(T);
    rhs3.x = 0.5 * dt * nlt.x;
    nlt = non_linear_temp(velocity,T);
    HfT = H_operation("cell",FT);
    
    rhs3 = CellData(Nx,Ny);
    rhs3.x = rhs3.x + Fo_t * diff_T.x - 1.5 * dt * nlt.x + dt * HfT.x;
    
    % Solve the diffusion problem
    
    [~,delta_T_star] = diffuse_dirichlet_cn_cell_xy(time_range(t),T,T0,rhs3,dt,velocity);
    
    % Now setup R4
    temp = CellData(Nx,Ny);
    temp.x = T.x + delta_T_star.x;
    Tb2 = E_operation("cell",temp);
    
    rhs4 = Tb - Tb2;
    delta_FT = pcg(A_temp, rhs4, 1e-1);
    delta_FT = delta_FT./dt;
    
    HDFT = H_operation("cell",delta_FT);
    delta_T = delta_T_star.x + dt * HDFT.x;
    
    % Assigning
    T.x = T.x + delta_T;
    T = apply_bc_temp(T,T0,velocity);
    FT = FT + delta_FT;
    
    % Checking the residual for each step and breaking the loop if convergence has reached
    
    conv_flow = 1/Nx * norm(delta_gamma(2:Nx,2:Ny))/dt / norm(gamma.x(2:Nx,2:Ny))/dt;
    conv_therm = 1/Nx * norm(delta_T(2:Nx+1,2:Ny+1))/dt / norm(T.x(2:Nx+1,2:Ny+1))/dt;
    if conv_flow <= tol && conv_therm <= tol
        break
    end
end
%% Streamlines

xi = body_map(:,1);
xi = [xi;xi(1)];
eta = body_map(:,2);
eta = [eta;eta(1)];
f1 = figure;
contour(X_n./(2*R),Y_n./(2*R),(sf.x'))
hold on
plot(xi./(2*R),eta./(2*R),"LineWidth",2);
hold off
pbaspect([1 1 1])
title({'Streamlines for a Flat Plate of D = ',num2str(2*R),' flow is of uniform velocity of U = ',num2str(U)})
xlabel("X/D")
ylabel("Y/D")
f1.WindowState = 'fullscreen';

%% Vorticity

xi = body_map(:,1);
xi = [xi;xi(1)];
eta = body_map(:,2);
eta = [eta;eta(1)];
f2 = figure;
contour(X_n./(2*R),Y_n./(2*R),(gamma.x'))
hold on
plot(xi./(2*R),eta./(2*R),"LineWidth",2);
hold off
pbaspect([1 1 1])
title({'Streamlines for a Flat Plate of D = ',num2str(2*R),' flow is of uniform velocity of U = ',num2str(U)})
xlabel("X/D")
ylabel("Y/D")
f2.WindowState = 'fullscreen';

%% Quiver Plot

xi = body_map(:,1);
xi = [xi;xi(1)];
eta = body_map(:,2);
eta = [eta;eta(1)];
f4 = figure;
qx = NodeData(Nx,Ny);
qy = NodeData(Nx,Ny);
qx = interpol(qx,velocity,1);
qy = interpol(qy,velocity,2);
quiver(X_n./(2*R),Y_n./(2*R),qx.x',qy.x');
hold on
plot(xi./(2*R),eta./(2*R),"LineWidth",2);
hold off
pbaspect([1 1 1])
title({'Velocity for a Flat Plate of D = ',num2str(2*R),' flow is of uniform velocity of U = ',num2str(U)})
xlabel("X/D")
ylabel("Y/D")
axis tight
f4.WindowState = 'fullscreen';

%% Temperature Field

xi = body_map(:,1);
xi = [xi;xi(1)];
eta = body_map(:,2);
eta = [eta;eta(1)];
f3 = figure;
contourf(X_c./(2*R),Y_c./(2*R),(T.x'))
hold on
plot(xi./(2*R),eta./(2*R),"LineWidth",2);
hold off
pbaspect([1 1 1])
title({'Temperature for a Flat Plate of D = ',num2str(2*R),' flow is of uniform velocity of U = ',num2str(U)})
xlabel("X/D")
ylabel("Y/D")
f3.WindowState = 'fullscreen';