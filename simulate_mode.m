%% Simulate the string system via differential equations for a specified
% environment / input. Here, the simulate code is refactored so that
% the user can specify a mode shape vector "U_r", the corresponding 
% resonant frequency "omega", a desired damping coefficient "damping" 
% and the function will simulate the predefined system (string_params)
% over a specified timespan "tspan"
% treats all ICs as zero (pos and vel)

% simulates movement of masses via U_r coefficients
function simulate_mode(omega, damping, tspan, string_params, save_vid)
    
    % assign system params
    num_masses = string_params.n;
    amplitude_Uf = 1; % ignore Uf coefficient
    omega_Uf = omega; % use the resonant frequency
    
    % construct the forcing function and its derivative:
    % uf = A*cos(w*t), uf' = -w*A*sin(w*t)
    Uf_func = @(t_in) amplitude_Uf*cos(omega_Uf*t_in);
    dUfdt_func = @(t_in) -omega_Uf*amplitude_Uf*sin(omega_Uf*t_in);
    
    % assign updated variable the struct
    string_params.Uf_func = Uf_func;
    string_params.dUfdt_func = dUfdt_func;
    string_params.c = damping;

    % create anonymized function
    my_rate_func = @(t_in, V_in) string_rate_func01(t_in, V_in, string_params);
    
    % initial conditions
    V0 = zeros(2*num_masses, 1);
    
    % run the integration
    [tlist, Vlist] = ode45(my_rate_func, tspan, V0); 
    
    % animate the system
    animate_string(Vlist, tlist, string_params, save_vid)

end