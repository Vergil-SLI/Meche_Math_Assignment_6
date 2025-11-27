%% simulate string and weights system

function string_simulation()
    
    % assign system params
    num_masses = 3; % per later analysis
    total_mass = 4.5; % arbitrary, may need to change
    tension_force = 5; % arbitrary, may need to change
    string_length = 5; % selected arbitrarily
    damping_coeff = 0.001; % selected per Orion's suggestion 0.1-0.01
    dx = string_length/(num_masses+1);
    amplitude_Uf = 0.1; % CHOSEN
    omega_Uf = 0.3; % CHOSEN
    
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
    my_rate_func = @(t_in, V_in) string_rate_func01_V2(t_in, V_in, string_params);
    
    % initial conditions
    U0 = zeros(num_masses, 1); % CHOSEN
    dUdt0 = zeros(num_masses, 1); % CHOSEN
    V0 = [U0; dUdt0];
    tspan = [0, 15]; % CHOSEN
    
    % run the integration
    [tlist, Vlist] = ode45(my_rate_func, tspan, V0); % may want to specify abstol?
    
    % animate the system
    save_vid = false;
    animate_string_sys_V2(Vlist, tlist, string_params, save_vid)
end

