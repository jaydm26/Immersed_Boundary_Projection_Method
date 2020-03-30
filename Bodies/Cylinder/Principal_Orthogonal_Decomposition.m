%% Principal Orthogonal Decomposition (POD) for Reduced Order Modelling (ROM)
% Take the velocities from the simulation done before. Use them to find
% important modes over the spatio-temporal domain. Using these modes,
% project the full N-S onto these modes and solve the problem for
% remarkably less computational effort.

% Find the Largest time scale from the Lift plot. Use the Data between two
% peaks. Load that data as gamma_log_2.

% velocity_log_2 = velocity_log(:,:,:);

% It may be that POD may not be feasible for edge space data. It would be
% prudent to interpolate the data on to the cell or node space and then
% perform this operation.

vel = EdgeData(Nx,Ny);
size_log = size(velocity_x_log);
n_time_steps = size_log(1);
stacked_velocity = zeros(((Nx+2)*(Ny+1) + (Nx+1)*(Ny+2)),n_time_steps);

for i = 1:n_time_steps
    vel.x(:,:) = velocity_x_log(i,:,:);
    vel.y(:,:) = velocity_y_log(i,:,:);
    stacked_velocity(:,i) = stacker(vel);
end

m_stacked_velocity = mean(stacked_velocity,2);

for i = 1:n_time_steps
    stacked_velocity(:,i) = stacked_velocity(:,i) - m_stacked_velocity;
end

C = stacked_velocity * stacked_velocity';
[ef,ev] = eig(C);
ev = max(ev);
ev = ev';

% Which modes are useful?

sum_ev = sort(cumsum(sort(ev,"descend"))/sum(ev),"descend");

% Recovery Percentage = 90%

for i = 1:length(sum_ev)
    if sum_ev(end-i+1) >= 0.90
        max_req_ev = i;
        break
    end
end

xi1 = [xi;xi(1)];
eta1 = [eta;eta(1)];

%% Form the Eigenfunctions for the important modes and store them as tensor

for i = 1:max_req_ev
    if i == 1
        phi_vel = [unstacker(domain,ef(:,end-i+1),"edge")];
    else
        phi_vel = [phi_vel,unstacker(domain,ef(:,end-i+1),"edge")];
    end
end

%% Plot the important modes
for i = 1:max_req_ev
    p_vel = phi_vel(i);
    figure
    contourf(X_e_x./(2*R),Y_e_x./(2*R),(p_vel.x'))
    hold on
    plot(xi1./(2*R),eta1./(2*R),'r',"LineWidth",2);
    hold off
    pbaspect([1 1 1])
    title(strcat('Kinetic Energy (Mode ',num2str(i),') for a Cylinder of D = ',num2str(2*R),' in a flow of uniform velocity of U = ',num2str(U)));
    xlabel("X/D")
    ylabel("Y/D")
end

%% Using the important modes to recreate the flow
% First precalculate the tensors F and G.
% F_ijk = <phi_i,phi_j . grad(phi_k)>
% This will require me to obtain a method to get the gradient of a velocity
% field which will result in a tensor. Could and should use nested
% structures.