function params = flow_parameters(dx, dt, U, nu, alpha, Fo, Co, char_L, Fo_t)
    %FLOW_PARAMETERS Creates a params data structure which contains all 
    % important flow parameters.
    %
    % params = flow_parameters(dx, dt, U, nu, Fo, Co, char_L, Fo_t)
    %
    % Variable lookup:
    % 
    % dx: Size of each cell. *Ensure that dx = dy. The code will give
    % erroneous results otherwise.*
    %
    % dt: Size of each time step moved during one iteration.
    % 
    % U: Characteristic velocity for the flow.
    %
    % nu: Kinematic Viscosity of the flow.
    %
    % Fo: Fourier Number defined by Fo = nu * dt / dx^2. Refer to reference
    % for further explanation.
    %
    % Co: Courant Number defined by Co = U * dt / dx. Refer to reference
    % for further explanation.
    % 
    % char_L: Characteristic length for the flow.
    %
    % Created by Jay Mehta (18 July 2019)
    
    params = struct;
    
    params.dx = dx;
    params.dt = dt;
    params.U = U;
    params.nu = nu;
    params.alpha = alpha;
    params.Fo = Fo;
    params.Co = Co;
    params.char_L = char_L;
    params.Fo_t = Fo_t;
    
end