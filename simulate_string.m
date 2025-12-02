%% Simulate string and weights system for modal analysis
% simulates system for a given input frequency for modal analysis
% all initial displacements and velocities are zero

% INPUTS
%   omega_Uf: a frequency for the forcing function uf(t) = A*cos(w*t)
%   tspan: the time span to simulate over. tspan = [t_start, t_end]
%   save_vid: a bool controling whether the simulation is saved as a file

function simulate_string(omega_Uf, damping_coeff, tspan, save_vid)
    
    % assign system params (all arbitratily chosen)
    num_masses = 3; 
    total_mass = 4.5; 
    tension_force = 5; 
    string_length = 5;
    dx = string_length/(num_masses+1);
    amplitude_Uf = 0.1; 
    
    % construct the forcing function and its derivative:
    % uf = A*cos(w*t), uf' = -w*A*sin(w*t)
    Uf_func = @(t_in) amplitude_Uf*cos(omega_Uf*t_in);
    dUfdt_func = @(t_in) -omega_Uf*amplitude_Uf*sin(omega_Uf*t_in);
    
    % generate the struct
    string_params = struct();
    string_params.n = num_masses;
    string_params.M = total_mass;
    string_params.Uf_func = Uf_func;
    string_params.dUfdt_func = dUfdt_func;
    string_params.Tf = tension_force;
    string_params.L = string_length;
    string_params.c = damping_coeff;
    string_params.dx = dx;
    
    % load string_params into rate function
    my_rate_func = @(t_in, V_in) string_rate_func01(t_in, V_in, string_params);
    
    % initial conditions
    V0 = zeros(2*num_masses, 1);  
    
    % run the integration
    [tlist, Vlist] = ode45(my_rate_func, tspan, V0); 
    
    % animate the system
    animate_string(Vlist, tlist, string_params, save_vid)
end

