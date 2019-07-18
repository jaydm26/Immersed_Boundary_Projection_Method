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
    X_e_x Y_e_x X_e_y Y_e_y Fo Co R
 
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

R = 0.5;
xc = 0;
yc = 0;
body_function = @(x,y) (x-xc)^2 + (y-yc)^2 - (R)^2;

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

%% Initializing the variables pre-solving

% Employ the Streamfunction - Vorticity Method to avoid the coupled forcing
% functions which show up in the Pressure-Velocity Method.

% Of course, in this method, we will employ the usual CN/AB2
% discretization, and then use the Block LU Decomposition. This will allow
% us to use the delta formulation.
U = 1;
V = 0;

Re = 100;
nu = U * R / Re;

Co = 0.1;
Fo = 5;
dt = min([Fo * dx^2/nu,Co*dx]);

Fo = nu * dt/dx^2;
Co = dt/dx;
t_steady = (2*R)^2/nu;
tf = t_steady;
time_range = 0:dt:tf;

%% Pre-Setup

velocity = EdgeData(Nx,Ny); % Velocity Field
velocity_x_log = zeros(length(time_range),Nx+1,Ny+2);
velocity_y_log = zeros(length(time_range),Nx+2,Ny+1);
sf_log = zeros(length(time_range),Nx+1,Ny+1);
gamma_log = zeros(length(time_range),Nx+1,Ny+1);
velocity_stack = zeros(((Nx+1)*(Ny+2) + (Nx+2)*(Ny+1)),length(time_range));
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
gamma = apply_bc_sp(gamma,gamma0);

%% Starter (Euler)
for t = 2
    % Set up R1
    gamma0 = gamma;
    diff_gamma = laplacian_2(gamma);
    nl = non_linear_alt(velocity);
    nl = curl_2(nl);
    
    ff = CTH(Fx,Fy);
    
    rhs1 = NodeData(Nx,Ny);
    rhs1.x = Fo * diff_gamma.x - dt * nl.x + dt * ff.x;
    
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
    gamma = apply_bc_sp(gamma,gamma0);
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
    
    sf_log(t,:,:) = sf.x;
    gamma_log(t,:,:) = gamma.x;
    velocity_x_log(t,:,:) = velocity.x;
    velocity_y_log(t,:,:) = velocity.y;
    velocity_stack(:,t) = stacker(velocity);
    
    % Checking the residual for each step and breaking the loop if convergence has reached
    
    conv = 1/Nx * norm(delta_gamma(2:Nx,2:Ny))/dt / norm(gamma.x(2:Nx,2:Ny))/dt;
    if conv <= tol
        break
    end
    
end

%% CN-AB2 Inital

for t = 3:25
    
    % Set up R1
    gamma0 = gamma;
    diff_gamma = laplacian_2(gamma);
    rhs1 = NodeData(Nx,Ny);
    rhs1.x = dt * 0.5 * nl.x;
    nl = non_linear_alt(velocity);
    nl = curl_2(nl);
    
    ff = CTH(Fx,Fy);
    
    rhs1.x = rhs1.x + Fo * diff_gamma.x - 1.5 * dt * nl.x + dt * ff.x;
    
    % Solve the diffusion problem
    
    [~,delta_gamma_star] = diffuse_dirichlet_cn_node_xy(time_range(t),rhs1,gamma,gamma0);
    
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
    
    gamma.x = gamma.x + delta_gamma;
    gamma = apply_bc_sp(gamma,gamma0);
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
    
    sf_log(t,:,:) = sf.x;
    gamma_log(t,:,:) = gamma.x;
    velocity_x_log(t,:,:) = velocity.x;
    velocity_y_log(t,:,:) = velocity.y;
    velocity_stack(:,t) = stacker(velocity);
    
    % Checking the residual for each step and breaking the loop if convergence has reached
    
    conv = 1/Nx * norm(delta_gamma(2:Nx,2:Ny))/dt / norm(gamma.x(2:Nx,2:Ny))/dt;
    if conv <= tol
        sf_log = sf_log(1:t,:,:);
        gamma_log = gamma_log(1:t,:,:);
        velocity_x_log = velocity_x_log(1:t,:,:);
        velocity_y_log = velocity_y_log(1:t,:,:);
        velocity_stack = velocity_stack(:,1:t);
        break
    end
end

%% CN-AB2 Body Forces

for t = 2049:4096
    
    % Set up R1
    gamma0 = gamma;
    diff_gamma = laplacian_2(gamma);
    rhs1 = NodeData(Nx,Ny);
    rhs1.x = dt * 0.5 * nl.x;
    nl = non_linear_alt(velocity);
    nl = curl_2(nl);
    
    ff = CTH(Fx,Fy);
    
    rhs1.x = rhs1.x + Fo * diff_gamma.x - 1.5 * dt * nl.x + dt * ff.x;
    
    % Solve the diffusion problem
    
    [~,delta_gamma_star] = diffuse_dirichlet_cn_node_xy(time_range(t),rhs1,gamma,gamma0);
    
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
    delta_f_x(3) = delta_f_x(3) + 10;
    delta_f_y = delta_f(length(body_map(:,1))+1:end);
    delta_f_y(3) = delta_f_y(3) + 10;
    
    % Now we correct for the vorticity to ensure we have no slip on the body
    
    gamma_c = CTH(delta_f_x,delta_f_y);
    
    delta_gamma = delta_gamma_star.x + dt * gamma_c.x;
    
    % Update gamma and the forcing function
    
    gamma.x = gamma.x + delta_gamma;
    gamma = apply_bc_sp(gamma,gamma0);
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
    
    sf_log(t,:,:) = sf.x;
    gamma_log(t,:,:) = gamma.x;
    velocity_x_log(t,:,:) = velocity.x;
    velocity_y_log(t,:,:) = velocity.y;
    velocity_stack(:,t) = stacker(velocity);
    
    % Checking the residual for each step and breaking the loop if convergence has reached
    
    conv = 1/Nx * norm(delta_gamma(2:Nx,2:Ny))/dt / norm(gamma.x(2:Nx,2:Ny))/dt;
    if conv <= tol
        sf_log = sf_log(1:t,:,:);
        gamma_log = gamma_log(1:t,:,:);
        velocity_x_log = velocity_x_log(1:t,:,:);
        velocity_y_log = velocity_y_log(1:t,:,:);
        velocity_stack = velocity_stack(:,1:t);
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
% title({'Streamlines for a Flat Plate of D = ',num2str(2*R),' flow is of uniform velocity of U = ',num2str(U)})
title(strcat('Streamlines for a Flat Plate of D = ',num2str(2*R),' flow is of uniform velocity of U = ',num2str(U)));
xlabel("X/D")
ylabel("Y/D")
f1.WindowState = 'fullscreen';

%% Vorticity

xi = body_map(:,1);
xi = [xi;xi(1)];
eta = body_map(:,2);
eta = [eta;eta(1)];
f1 = figure;
contour(X_n./(2*R),Y_n./(2*R),(gamma.x'))
hold on
plot(xi./(2*R),eta./(2*R),"LineWidth",2);
hold off
pbaspect([1 1 1])
title({'Streamlines for a Flat Plate of D = ',num2str(2*R),' flow is of uniform velocity of U = ',num2str(U)})
xlabel("X/D")
ylabel("Y/D")
f1.WindowState = 'fullscreen';


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

%% Video Generator

video = video_generator("streamlines2",sf_log);

%% Principal Orthogonal Decomposition

% Find the Largest time scale from the Lift plot. Use the Data between two
% peaks. Load that data as gamma_log_2.

gam = NodeData(Nx,Ny);
stacked_gamma = zeros(((Nx+1)*(Ny+1)),length(gamma_log_2));

for i = 1:length(gamma_log_2)
    gam.x(:,:) = gamma_log_2(i,:,:);
    stacked_gamma(:,i) = stacker(gam);
end

m_stacked_gamma = mean(stacked_gamma,2);

for i = 1:length(gamma_log_2)
    stacked_gamma(:,i) = stacked_gamma(:,i) - m_stacked_gamma;
end

C = stacked_gamma * stacked_gamma';
[ef,ev] = eig(C);
ev = max(ev);
ev = ev';

% Which modes are useful?

sum_ev = sort(cumsum(sort(ev,"descend"))/sum(ev),"descend");

for i = 1:length(sum_ev)
    if sum_ev(end-i+1) >= 0.99
        max_req_ev = i;
        break
    end
end

for i = 1:max_req_ev
    p_gam = unstacker(ef(:,end-i+1),"node");
    figure
    contourf(X_n./(2*R),Y_n./(2*R),(p_gam.x'))
    hold on
    plot(xi./(2*R),eta./(2*R),'r',"LineWidth",2);
    hold off
    pbaspect([1 1 1])
    title(strcat('Enstrophy (Mode ',num2str(i),') for a Cylinder of D = ',num2str(2*R),' in a flow of uniform velocity of U = ',num2str(U)));
    xlabel("X/D")
    ylabel("Y/D")
end
    