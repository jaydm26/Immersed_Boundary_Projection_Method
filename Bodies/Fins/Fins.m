%% Immersed Boundary Projection Method
%
% Code written by Jay Mehta (July 2019).
% 
% This code was initially developed as a part of a Numerical Methods for
% Incompressible Flows class (MECH&AE 250H) at UCLA.
% 
% The code is developed from the works of Prof. Tim Colonius (CalTech),
% Prof. Sam Taira (UCLA) and Prof. Jeff Eldredge (UCLA). Other authors have been
% referenced when used.
%
% This code solves the problem for a viscous flow past a 2D cylinder.
% Initally the code is about pre-processing to generate the domain and the
% body. The middle part solves the equations. The trailing portion of the
% code is for post-processing. Each part (up to solving the problem) has 
% been explained in Taira, Colonius (JCP 2007) and Kajishima and Taira
% (Springer, 2017).

%% Clear Everything
clear all
clc
rmpath('Source')
rmpath('Source/dst_idst')

%% Add Paths
addpath('Source')
addpath('Source/dst_idst')

%% Set up the problem domain and the problem object

% Set up the parameters that have to be passed

params = flow_parameters_init;
domain = domain_parameters_init;

% Domain
Nx = 48;
Ny = 24;
domain.Nx = Nx;
domain.Ny = Ny;

x_range = [-120e-3 120e-3];
y_range = [-60e-3 60e-3];
domain.x_range = x_range;
domain.y_range = y_range;

dx = (x_range(2)-x_range(1))/Nx;
params.dx = dx;
dy = (y_range(2)-y_range(1))/Ny;
params.dy = dy;

%% Setting up the Domain of the Problem

[X_n, Y_n] = DomainSetup(params,domain,"node");
[X_e_x, Y_e_x] = DomainSetup(params,domain,"xe");
[X_e_y, Y_e_y] = DomainSetup(params,domain,"ye");
[X_c,Y_c] = DomainSetup(params,domain,"cell");

domain.X_n = X_n;
domain.Y_n = Y_n;
domain.X_e_x = X_e_x;
domain.Y_e_x = Y_e_x;
domain.X_e_y = X_e_y;
domain.Y_e_y = Y_e_y;
domain.X_c = X_c;
domain.Y_c = Y_c;

%% Create the L^-1 operator using Lattice Green's function (Liska, Colonius, 2016)

g_hat = L_inv(domain,"node");

%% Setting up the Object

% Line 1
[xi1,eta1,~,~] = Line_Builder(params,-20e-3,-20e-3,-12.7e-3 + 1.12*dx,12.7e-3 - 1.12*dx);
xi1 = xi1(1:end);
eta1 = eta1(1:end);

% Line 2
[xi2,eta2,~,~] = Line_Builder(params,-20e-3,-5e-3*sqrt(5),12.7e-3,12.7e-3);
xi2 = xi2(1:end-1);
eta2 = eta2(1:end-1);

% % Line 3
% [xi3,eta3,~,~] = Line_Builder(params,20e-3,20e-3,12.7e-3 - 1.12*dx,-12.7e-3 + 1.12*dx);
% xi3 = xi3(1:end);
% eta3 = eta3(1:end);

% Dimple Line 3
R = 15e-3;
xc = 0;
yc = 10e-3+12.7e-3;

N_theta = floor(2*acos((yc-12.7e-3)/R)*R/dx);
d_theta = 2*acos((yc-12.7e-3)/R)/N_theta;

theta_range = 1.5*pi-acos(abs((yc-12.7e-3))/R):d_theta:1.5*pi + acos((yc-12.7e-3)/R);

xi3 = zeros(length(theta_range),1);
eta3 = zeros(length(theta_range),1);

for i = 1:length(theta_range)
    xi3(i) = xc + R*cos(theta_range(i));
    eta3(i) = yc + R*sin(theta_range(i));
end

% Line 4
[xi4,eta4,~,~] = Line_Builder(params,5e-3*sqrt(5),20e-3,12.7e-3,12.7e-3);
xi4 = xi4(2:end);
eta4 = eta4(2:end);

% Line 5
[xi5,eta5,~,~] = Line_Builder(params,20e-3,20e-3,12.7e-3 - 1.12 * dx,-12.7e-3 + 1.12 * dx);
xi5 = xi5(1:end);
eta5 = eta5(1:end);

% Line 6
[xi6,eta6,~,~] = Line_Builder(params,20e-3,5e-3*sqrt(5),-12.7e-3,-12.7e-3);
xi6 = xi6(1:end-1);
eta6 = eta6(1:end-1);

% Line 7
R = 15e-3;
xc = 0;
yc = -10e-3-12.7e-3;

N_theta = floor(2*acos(abs((yc+12.7e-3))/R)*R/dx);
d_theta = 2*acos(abs((yc+12.7e-3))/R)/N_theta;

theta_range = 0.5*pi-acos(abs((yc+12.7e-3))/R):d_theta:0.5*pi + acos(abs((yc+12.7e-3))/R);

xi7 = zeros(length(theta_range),1);
eta7 = zeros(length(theta_range),1);

for i = 1:length(theta_range)
    xi7(i) = xc + R*cos(theta_range(i));
    eta7(i) = yc + R*sin(theta_range(i));
end

% Line 8
[xi8,eta8,~,~] = Line_Builder(params,-5e-3*sqrt(5),-20e-3,-12.7e-3,-12.7e-3);
xi8 = xi8(2:end);
eta8 = eta8(2:end);

char_L = 10*sqrt(5) * 1e-3;
params.char_L = char_L;

xi = [xi1;xi2;xi3;xi4;xi5;xi6;xi7;xi8];
eta = [eta1;eta2;eta3;eta4;eta5;eta6;eta7;eta8];

%% Forming the A matrix for PCG. Load it from the side bar

A = MatrixA_Generator(params,domain,xi,eta,"vel",g_hat);
AT = MatrixA_Generator(params,domain,xi,eta,"temp");

%% Initializing the variables pre-solving

% Employ the Streamfunction - Vorticity Method to avoid the coupled forcing
% functions which show up in the Pressure-Velocity Method.

% Of course, in this method, we will employ the usual CN/AB2
% discretization, and then use the Block LU Decomposition. This will allow
% us to use the delta formulation.

U = 1;
V = 0;
params.U = U;

Re = 40;
nu = U * char_L / Re;
params.nu = nu;

Pr = 1;
alpha = nu/Pr;
params.alpha = alpha;

Co = 1e-3;
Fo = 5;
Fo_t = Fo;

dt = min([Fo * (dx)^2/nu,Co*(dx)/U,Fo_t * (dx)^2/alpha]);
params.dt = dt;

Fo = nu * dt/(dx)^2;
Fo_t = alpha * dt/(dx)^2;
Co = U * dt/(dx);

params.Fo = Fo;
params.Fo_t = Fo_t;
params.Co = Co;

t_steady = max([(char_L)^2/nu,(char_L)^2/alpha]);
tf = t_steady;
time_range = 0:dt:t_steady;

%% Pre-Setup

velocity = EdgeData(Nx,Ny); % Velocity Field
gamma = NodeData(Nx,Ny); % Vorticity
Fx = zeros(length(xi),1);
Fy = zeros(length(eta),1);
ub = zeros(length(xi),1); % X-component of Velocity on the body
vb = zeros(length(eta),1); % Y-component of Velocity on the body
sf = NodeData(Nx,Ny);
T = EdgeData(Nx,Ny);
Tb = ones(length(xi),1);
FTX = zeros(length(xi),1);
FTY = zeros(length(eta),1);
tol = 1e-1;
Drag = 0;
Lift = 0;
rhs1 = NodeData(Nx,Ny);
rhs3 = EdgeData(Nx,Ny);

%% Boundary Conditions on Gamma

gamma0 = gamma;
gamma = apply_bc_sp(params,gamma,gamma0);
T0 = T;
T = apply_bc_temp(params,T,T0);

%% Starter (Euler)

for t = 2
    
    % Solving the for the flow velocities first
    
    % Set up R1
    gamma0 = gamma;
    velocity0 = velocity;
    diff_gamma = laplacian_2(gamma);
    nl = non_linear(params,velocity);
    nl = curl_2(nl);
    
    ff = CTH(params,domain,xi,eta,Fx,Fy);
    
    rhs1.x = Fo * diff_gamma.x - dt * nl.x + dt * ff.x;
    
    % Solve the diffusion problem
    
    delta_gamma_star = rhs1;
    
    % Now set up R2. Note here that R2 will have two components.
    
    temp = NodeData(Nx,Ny);
    temp.x = gamma.x + delta_gamma_star.x;
    [ub2,vb2] = ECL_inv(params,domain,g_hat,xi,eta,temp);
    
    % Forming the two rhs2's
    
    rhs2_x = ub - U * ones(length(xi),1) + ub2;
    rhs2_y = vb - V * ones(length(eta),1) + vb2;
    
    rhs2 = [rhs2_x;rhs2_y];
    % Now obtain the delta_f in the two directions by pcg
    delta_f = pcg(A, rhs2, 1e-1);
    delta_f = delta_f./dt;
    
    % Now delta_f_star in the X and Y direction are the correct delta_f in the X and Y direction 
    
    delta_f_x = delta_f(1:length(xi));
    delta_f_y = delta_f(length(xi)+1:end);
    
    % Now we correct for the vorticity to ensure we have no slip on the body
    
    gamma_c = CTH(params,domain,xi,eta,delta_f_x,delta_f_y);
    
    delta_gamma = delta_gamma_star.x + dt * gamma_c.x;
    
    % Update gamma and the forcing function % Adding the unstable component
    
    gamma.x = gamma.x + delta_gamma;
    gamma = apply_bc_sp(params,gamma,gamma0);
    Fx = Fx + delta_f_x;
    Fy = Fy + delta_f_y;
    
    % Obtain the streamfunction and velocity
    
    gamma.x = -gamma.x;
    sf = L_inv_operation(gamma,g_hat);
    gamma.x = -gamma.x;
    
    velocity = curl_2(sf);
    velocity.x = velocity.x + U * ones(Nx+1,Ny+2);
    velocity.y = velocity.y;
    
    % Calculate Lift and Drag
    
    Hq = H_operation(params,domain,"edge",xi,eta,Fx,Fy);
    
    Drag(t) = -sum(sum(Hq.x))*dx^2;
    Lift(t) = sum(sum(Hq.y))*dx^2;
    
    % Solving for the temperature field
    
    % Set R3
    T0 = T;
    diff_T = laplacian_2(T);
    nlt = non_linear_temp(params,velocity0,T);
    
    fT = H_operation(params,domain,"edge",xi,eta,FTX,FTY);
    
    rhs3 = CellData(Nx,Ny);
    rhs3.x = Fo_t * diff_T.x - dt .* nlt.x + dt * fT.x;
    rhs3.y = Fo_t * diff_T.y - dt .* nlt.y + dt * fT.y;
    
    % Solve the diffusion problem
    
    delta_T_star = rhs3;
    
    % Now set up R4.
    
    temp = EdgeData(Nx,Ny);
    temp.x = T.x + delta_T_star.x;
    temp.y = T.y + delta_T_star.y;
    [Tb_x2,Tb_y2] = E_operation(params,domain,xi,eta,temp);
    
    rhs4 = Tb - Tb_x2;
    
    delta_FT = pcg(AT,rhs4,1e-1);
    delta_FTX = delta_FT./dt;
    delta_FTY = zeros(length(eta),1);
    
    T_c = H_operation(params,domain,"edge",xi,eta,delta_FTX,delta_FTY);
    delta_T_x = delta_T_star.x + dt * T_c.x;
    delta_T_y = delta_T_star.y + dt * T_c.y;
    
    T.x = T.x + delta_T_x;
    T.y = T.y + delta_T_y;
    T = apply_bc_temp(params,T,T0);
    FTX = FTX + delta_FTX;
    
    % Checking the residual for each step and breaking the loop if convergence has reached
    
    conv_flow = 1/Nx * norm(delta_gamma(2:Nx,2:Ny))/dt / norm(gamma.x(2:Nx,2:Ny))/dt;
    conv_therm = 1/Nx * norm(delta_T_x(2:Nx+1,2:Ny+1))/dt / norm(T.x(2:Nx+1,2:Ny+1))/dt;
    if conv_flow <= tol && conv_therm <= tol
        break
    end
end

%% CN-AB2 Initial

for t = 3:200
    
    % Set up R1
    gamma0 = gamma;
    velocity0 = velocity;
    diff_gamma = laplacian_2(gamma);
    rhs1 = NodeData(Nx,Ny);
    rhs1.x = dt * 0.5 * nl.x;
    nl = non_linear(params,velocity);
    nl = curl_2(nl);
    
    ff = CTH(params,domain,xi,eta,Fx,Fy);
    
    rhs1.x = rhs1.x + Fo * diff_gamma.x - 1.5 * dt * nl.x + dt * ff.x;
    
    % Solve the diffusion problem
    
    [~,delta_gamma_star] = diffuse_dirichlet_cn_node_xy(params,time_range(t),rhs1,gamma,gamma0);
    
    % Now set up R2. Note here that R2 will have two components.
    
    temp = NodeData(Nx,Ny);
    temp.x = gamma.x + delta_gamma_star.x;
    [ub2,vb2] = ECL_inv(params,domain,g_hat,xi,eta,temp);
    
    rhs2_x = ub - U * ones(length(xi),1) + ub2;
    rhs2_y = vb - V * ones(length(eta),1) + vb2;
    
    rhs2 = [rhs2_x;rhs2_y];
    % Now obtain the delta_f in the two directions by pcg
    delta_f = pcg(A, rhs2, 1e-1,500);
    delta_f = delta_f./dt;
    % Now delta_f_star in the X and Y direction are the correct delta_f in the X and Y direction
    
    % Assigning
    
    delta_f_x = delta_f(1:length(xi));
    delta_f_y = delta_f(length(xi)+1:end);
    
    % Now we correct for the vorticity to ensure we have no slip on the body
    
    gamma_c = CTH(params,domain,xi,eta,delta_f_x,delta_f_y);
    delta_gamma = delta_gamma_star.x + dt * gamma_c.x;
    
    % Update gamma and the forcing function
    
    gamma.x = gamma.x + delta_gamma;
    gamma = apply_bc_sp(params,gamma,gamma0);
    Fx = Fx + delta_f_x;
    Fy = Fy + delta_f_y;
    
    % Obtain the streamfunction and velocity
    gamma.x = -gamma.x;
    sf = L_inv_operation(gamma,g_hat);
    gamma.x = -gamma.x;
    
    velocity = curl_2(sf);
    velocity.x = velocity.x + U * ones(Nx+1,Ny+2);
    velocity.y = velocity.y;
    
    % Calculate Lift and Drag
    
    Hq = H_operation(params,domain,"edge",xi,eta,Fx,Fy);
    
    Drag(t) = -sum(sum(Hq.x))*dx^2;
    Lift(t) = sum(sum(Hq.y))*dx^2;
    
    % Solving for the temperature field
    
    % Set R3
    rhs3.x = 0.5 * dt .* nlt.x;
    diff_T = laplacian_2(T);
    nlt = non_linear_temp(params,velocity0,T);
    
    fT = H_operation(params,domain,"edge",xi,eta,FTX,FTY);
    
    rhs3.x = rhs3.x + Fo_t * diff_T.x - 1.5 * dt .* nlt.x + dt * fT.x;
    rhs3.y = rhs3.y + Fo_t * diff_T.y - 1.5 * dt .* nlt.y + dt * fT.y;
    
    % Solve the diffusion problem
    
    [~,delta_T_star] = diffuse_dirichlet_cn_cell_xy(params,time_range(t),rhs3,T);
    
    % Now set up R4.
    
    temp = EdgeData(Nx,Ny);
    temp.x = T.x + delta_T_star.x;
    temp.y = T.y + delta_T_star.y;
    [Tb_x2,Tb_y2] = E_operation(params,domain,xi,eta,temp);
    
    rhs4 = Tb - Tb_x2;
    
    delta_FT = pcg(AT,rhs4,1e-1,100);
    delta_FTX = delta_FT./dt;
    delta_FTY = zeros(length(eta),1);
    
    T_c = H_operation(params,domain,"edge",xi,eta,delta_FTX,delta_FTY);
    delta_T_x = delta_T_star.x + dt * T_c.x;
    delta_T_y = delta_T_star.y + dt * T_c.y;
    
    T.x = T.x + delta_T_x;
    T.y = T.y + delta_T_y;
    T = apply_bc_temp(params,T);
    FTX = FTX + delta_FTX;
    FTY = FTY + delta_FTY;
    
    % Checking the residual for each step and breaking the loop if convergence has reached
    
    conv_flow = 1/Nx * norm(delta_gamma(2:Nx,2:Ny))/dt / norm(gamma.x(2:Nx,2:Ny))/dt;
    conv_therm = 1/Nx * norm(delta_T_x(2:Nx,2:Ny+1))/dt / norm(T.x(2:Nx,2:Ny+1))/dt;
    
    if conv_flow <= tol && conv_therm <= tol
        break
    end
end

%% Streamlines

xi_plot = [xi;xi(1)];
eta_plot = [eta;eta(1)];
f1 = figure;
contour(X_n./(char_L),Y_n./(char_L),(sf.x'))
hold on
plot(xi_plot./(char_L),eta_plot./(char_L),"r","LineWidth",2);
hold off
pbaspect([1 (y_range(2)-y_range(1))/(x_range(2)-x_range(1)) 1])
title(strcat('Streamlines for a Closed Body of characteristic length = '...
    ,num2str(char_L),' in a flow of uniform velocity of U = ',num2str(U)));
xlabel("X/L")
ylabel("Y/L")
f1.WindowState = 'fullscreen';

%% Vorticity

xi_plot = [xi;xi(1)];
eta_plot = [eta;eta(1)];
f2 = figure;
contour(X_n./(char_L),Y_n./(char_L),(gamma.x'))
hold on
plot(xi_plot./(char_L),eta_plot./(char_L),"r","LineWidth",2);
hold off
pbaspect([1 (y_range(2)-y_range(1))/(x_range(2)-x_range(1)) 1])
title(strcat('Vorticity for a Closed Body of characteristic length = '...
    ,num2str(char_L),' in a flow of uniform velocity of U = ',num2str(U)));
xlabel("X/L")
ylabel("Y/L")
f2.WindowState = 'fullscreen';

%% Quiver Plot

xi_plot = [xi;xi(1)];
eta_plot = [eta;eta(1)];
f4 = figure;
qx = interpol(velocity,NodeData(Nx,Ny),1);
qy = interpol(velocity,NodeData(Nx,Ny),2);
quiver(X_n./(char_L),Y_n./(char_L),qx.x',qy.x');
hold on
plot(xi_plot./(char_L),eta_plot./(char_L),"r","LineWidth",2);
hold off
pbaspect([1 (y_range(2)-y_range(1))/(x_range(2)-x_range(1)) 1])
title(strcat('Velocity for a Closed Body of characteristic length = '...
    ,num2str(char_L),' in a flow of uniform velocity of U = ',num2str(U)));
xlabel("X/L")
ylabel("Y/L")
axis tight
f4.WindowState = 'fullscreen';

%% Temperature Plot

xi_plot = [xi;xi(1)];
eta_plot = [eta;eta(1)];
T1 = interpol(T,CellData(Nx,Ny),1);
f3 = figure;
contourf(X_c./(char_L),Y_c./(char_L),(T1.x'),250,"LineColor","none")
hold on
plot(xi_plot./(char_L),eta_plot./(char_L),"r","LineWidth",2);
hold off
pbaspect([1 (y_range(2)-y_range(1))/(x_range(2)-x_range(1)) 1])
title(strcat('Temperature for a Closed Body of characteristic length = '...
    ,num2str(char_L),' in a flow of uniform velocity of U = ',num2str(U)));
xlabel("X/L")
ylabel("Y/L")
f3.WindowState = 'fullscreen';

%% Video Generator

video = video_generator(params,domain,"von-karman street",xi,eta,sf_log);

%% Shedding Frequency

fre = 1/dt * (2:t)/t;
fft_vel = mean(real(fft(transpose(velocity_stack(:,2:t)))),2);
plot(fre(2:22),fft_vel(2:22)) % Looking for a low frequency peak

%% Principal Orthogonal Decomposition

% Find the Largest time scale from the Lift plot. Use the Data between two
% peaks. Load that data as gamma_log_2.

gamma_log_2 = gamma_log(:,:,:);

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

xi1 = [xi;xi(1)];
eta1 = [eta;eta(1)];

for i = 1:max_req_ev
    p_gam = unstacker(domain,ef(:,end-i+1),"node");
    figure
    contourf(X_n./(2*R),Y_n./(2*R),(p_gam.x'))
    hold on
    plot(xi1./(2*R),eta1./(2*R),'r',"LineWidth",2);
    hold off
    pbaspect([1 1 1])
    title(strcat('Enstrophy (Mode ',num2str(i),') for a Cylinder of D = ',num2str(2*R),' in a flow of uniform velocity of U = ',num2str(U)));
    xlabel("X/D")
    ylabel("Y/D")
end

%% Body Function

function out = body_function(x,y)
    if (x>=-20e-3 && x<= 20e-3) && (y>=-12.7e-3 && y<=12.7e-3) % Blanket x and y -condition
        out = -1;
        if (x-0)^2 + (y-22.7e-3)^2 - (15e-3)^2 < 0 || (x-0)^2 + (y+22.7e-3)^2 - (15e-3)^2 < 0
            out = 1;
        end
    else
        out = 1;
    end
end